"""
Hard genotype: step 1 - precompute the per-chromosome graph structure used
to genotype edges.

This step only uses graph topology (the GFA and the Chrom/Edge columns of the
raw edge-depth matrix); no sample depths are read. It therefore runs once per
chromosome, and the output is shared by every sample in the genotyping step.

Steps:
  1. Read the GRCh38 walk (W-line with sample name "GRCh38") from the GFA to
     get the reference backbone node order.
  2. Build an undirected edge graph from the depth-matrix edges. Remove bridge
     edges that are not on GRCh38 (no ploidy information to genotype them).
  3. Find articulation points (plus the backbone start/end nodes) and split
     the graph into biconnected components:
       - "anchored" components : exactly 2 boundary articulation nodes, both
                                 on GRCh38 -> genotyped as a ploidy-2 flow DAG
       - "partially anchored"  : >=1 boundary node on GRCh38 but not a clean
                                 2-node GRCh38 boundary -> every edge genotyped
                                 with the copy-number model, using the local
                                 depth around the GRCh38 anchor node (the one
                                 with the lowest GRCh38 position)
       - no GRCh38 boundary    : left for step 6 (recovery)
  4. Split every haplotype walk in the GFA into segments between consecutive
     GRCh38 articulation nodes (ordered by GRCh38 position, not node ID).
     Partial walks that touch only one boundary, and walks that touch no
     boundary but overlap a component, are also recorded (over-inclusive;
     nodes outside the component are dropped in step 5).
  5. For each anchored component, build a directed multigraph from the walk
     segments, collapse directed cycles (strongly connected components) into
     single nodes to get a DAG, record the cyclic edges (genotyped with the
     copy-number model), and precompute the DAG traversal.
  6. Recovery: every depth-matrix edge that steps 3-5 cannot genotype (skipped
     components, removed off-GRCh38 bridge edges, DAG edges unreachable from
     the component's start node) is grouped by connected component of the
     missing edges. Each group gets one GRCh38 anchor node: its own backbone
     node with the lowest GRCh38 position if it has one, otherwise the
     nearest backbone node by breadth-first search through the full edge
     graph (bridges included). These edges are genotyped with the copy-number
     model around the anchor. Groups with no reachable backbone node are
     dropped (reported in the log).

Input:
    --chrom     : chromosome name (e.g. chr21)
    --gfa       : chromosome GFA file (W-lines are read; the GRCh38 walk must
                  be present exactly once)
    --depthfile : per-chromosome raw edge-depth matrix (Chrom, Edge, sample...),
                  plain or .gz, with or without a header line. Only the Chrom
                  and Edge columns are used.

Output:
    --out : pickle (dict) with keys:
              chr, path_edges, b38_node_order, start_node, end_node,
              articulation_nodes, outgoing_paths, incoming_paths,
              stem_positions, components, recovery_groups

Usage:
    python prep_component_structure.py \\
        --chrom     chr21 \\
        --gfa       chr21.gfa \\
        --depthfile chr21.depth_per_edge.txt.gz \\
        --out       chr21.structure.pkl
"""

import argparse
import pickle
import re
import sys
from collections import defaultdict, deque

import networkx as nx

from graph_utils import (
    open_maybe_gzip,
    parse_path,
    canonical_edge,
    reverse_complement_path,
    precompute_component_traversal,
)

NODE_SPLIT = re.compile(r">|<")
TOKEN_PATTERN = re.compile(r"[><]([^><]+)")


# ---------------------------------------------------------------------------
# Step 1: GRCh38 backbone from the GFA
# ---------------------------------------------------------------------------
def read_grch38_backbone(gfa):
    """Return the GRCh38 node ids in path order (exactly one GRCh38 W-line required)."""
    grch38_walks = []
    with open_maybe_gzip(gfa) as gfa_file:
        for line in gfa_file:
            if not line.startswith('W'):
                continue
            parts = line.rstrip('\n').split('\t')
            if parts[1] == "GRCh38":
                grch38_walks.append(parts[6])
    if len(grch38_walks) != 1:
        sys.exit(f"ERROR: expected exactly 1 GRCh38 W-line in {gfa}, found {len(grch38_walks)}")
    return NODE_SPLIT.split(grch38_walks[0])[1:]


def read_walks(gfa):
    """Return every W-line walk string in the GFA."""
    walks = []
    with open_maybe_gzip(gfa) as gfa_file:
        for line in gfa_file:
            if line.startswith('W'):
                walks.append(line.rstrip('\n').split('\t')[6])
    return walks


# ---------------------------------------------------------------------------
# Step 2: undirected edge graph from the depth matrix (topology only)
# ---------------------------------------------------------------------------
def read_edge_graph(depthfile):
    """
    Returns:
        G               : undirected graph of all depth-matrix edges
        raw_edge_lookup : (node1, node2) sorted -> original directed edge string
        all_edges       : [(node1, node2, canonical_edge_id), ...] for every edge
    """
    G = nx.Graph()
    raw_edge_lookup = {}
    all_edges = []
    with open_maybe_gzip(depthfile) as depth_file:
        for line in depth_file:
            edge = line.split('\t', 2)[1]
            if edge[:1] not in ('>', '<'):   # header line
                continue
            nodes = NODE_SPLIT.split(edge)[1:]
            if nodes[0] > nodes[1]:
                nodes[0], nodes[1] = nodes[1], nodes[0]
            G.add_edge(nodes[0], nodes[1])
            raw_edge_lookup[(nodes[0], nodes[1])] = edge
            handles = parse_path(edge)
            all_edges.append((nodes[0], nodes[1], canonical_edge(handles[0], handles[1])))
    return G, raw_edge_lookup, all_edges


def remove_off_backbone_bridges(G, grch38_node_set):
    """Remove bridge edges not on GRCh38 in place; return node -> [neighbors via removed bridges]."""
    removed_adj = defaultdict(list)
    for u, v in set(nx.bridges(G)):
        if u in grch38_node_set and v in grch38_node_set:
            continue
        G.remove_edge(u, v)
        removed_adj[u].append(v)
        removed_adj[v].append(u)
    return removed_adj


# ---------------------------------------------------------------------------
# Step 4: split haplotype walks into segments between GRCh38 articulation nodes
# ---------------------------------------------------------------------------
def orient_touchless_path(path, b38_node_order):
    """
    Orient a walk that touches no boundary node: use the GRCh38 order of its
    first vs last backbone node if they differ; otherwise use whichever handle
    direction (>/<) dominates. Reverse-complement if it runs backward.
    """
    first_b38 = last_b38 = None
    for match in TOKEN_PATTERN.finditer(path):
        order = b38_node_order.get(match.group(1))
        if order is None:
            continue
        if first_b38 is None:
            first_b38 = order
        last_b38 = order
    if first_b38 is not None and first_b38 != last_b38:
        return path if last_b38 > first_b38 else reverse_complement_path(path)
    return reverse_complement_path(path) if path.count("<") > path.count(">") else path


def split_paths(paths, split_nodes, b38_node_order, node_to_pairs):
    """
    For each pair of consecutive split nodes (sorted by GRCh38 position), collect
    the segments of every walk that lie between them, oriented start -> end.

      - Full match: walk contains both nodes. The largest range (first
        occurrence of start to last occurrence of end) is kept; walks that
        traverse the pair in reverse are reverse-complemented.
      - Partial walk: from the walk's first touched split node back to the walk
        start, and from its last touched split node to the walk end, recorded
        for BOTH pairs that node bounds (oriented to the node's role in each).
      - Touchless walk: touches no split node; recorded whole for every
        anchored component whose nodes it overlaps (node_to_pairs).

    Partial and touchless segments may include nodes outside the component;
    those are dropped when building each component's directed graph.

    Returns: {(start_node, end_node): set(segment strings)}
    """
    split_nodes = sorted(split_nodes, key=lambda x: b38_node_order[x])
    pairs = list(zip(split_nodes, split_nodes[1:]))
    split_node_set = set(split_nodes)

    # each split node bounds up to 2 adjacent pairs
    node_pairs = defaultdict(list)
    for pair in pairs:
        node_pairs[pair[0]].append(pair)
        node_pairs[pair[1]].append(pair)

    results = {pair: set() for pair in pairs}

    for path in paths:
        path = path.strip()
        first_pos = {}           # node -> (start, end) char span of first occurrence
        last_pos = {}            # node -> (start, end) char span of last occurrence
        earliest_touch = None    # first split node seen in the walk
        latest_touch = None      # last split node seen in the walk
        for match in TOKEN_PATTERN.finditer(path):
            node = match.group(1)
            if node not in split_node_set:
                continue
            if node not in first_pos:
                first_pos[node] = (match.start(), match.end())
            last_pos[node] = (match.start(), match.end())
            if earliest_touch is None:
                earliest_touch = node
            latest_touch = node

        # touchless walk
        if not first_pos:
            overlapping_pairs = set()
            token_count = 0
            for match in TOKEN_PATTERN.finditer(path):
                token_count += 1
                overlapping_pairs.update(node_to_pairs.get(match.group(1), ()))
            # a single-handle walk has no edge to contribute
            if not overlapping_pairs or token_count < 2:
                continue
            oriented = orient_touchless_path(path, b38_node_order)
            for pair in overlapping_pairs:
                results[pair].add(oriented)
            continue

        # full matches (forward and reverse). '>=' because adjacent tokens abut.
        relevant_pairs = set()
        for node in first_pos:
            relevant_pairs.update(node_pairs[node])

        for start, end in relevant_pairs:
            start_first = first_pos.get(start)
            end_last = last_pos.get(end)
            if start_first is not None and end_last is not None and end_last[0] >= start_first[1]:
                results[(start, end)].add(path[start_first[0]:end_last[1]])

            end_first = first_pos.get(end)
            start_last = last_pos.get(start)
            if end_first is not None and start_last is not None and start_last[0] >= end_first[1]:
                results[(start, end)].add(reverse_complement_path(path[end_first[0]:start_last[1]]))

        # partial pieces at the walk's outer ends; skipped when the touched
        # node is itself the first/last token (nothing on the open side)
        if first_pos[earliest_touch][0] > 0:
            before_piece = path[:last_pos[earliest_touch][1]]   # walk start -> node (ends with node)
            before_piece_rc = reverse_complement_path(before_piece)
            for pair in node_pairs[earliest_touch]:
                results[pair].add(before_piece_rc if pair[0] == earliest_touch else before_piece)

        if last_pos[latest_touch][1] < len(path):
            after_piece = path[first_pos[latest_touch][0]:]     # node -> walk end (starts with node)
            after_piece_rc = reverse_complement_path(after_piece)
            for pair in node_pairs[latest_touch]:
                results[pair].add(after_piece if pair[0] == latest_touch else after_piece_rc)

    return {pair: segments for pair, segments in results.items() if segments}


# ---------------------------------------------------------------------------
# Step 5: per anchored component - directed graph, cyclic edges, DAG traversal
# ---------------------------------------------------------------------------
def build_anchored_component(meta, segments):
    bicomponent = meta["bicomponent"]
    start_component_node = meta["start_component_node"]

    # distinct directed (left, right) handle pairs with both nodes in the component
    raw_edges = set()
    for segment in segments:
        handles = parse_path(segment)
        for left, right in zip(handles, handles[1:]):
            if left[1:] in bicomponent and right[1:] in bicomponent:
                raw_edges.add((left, right))

    D = nx.MultiDiGraph()
    observed_edges = []
    for left, right in raw_edges:
        source, target = left[1:], right[1:]
        D.add_edge(source, target, key=left + right)
        observed_edges.append((source, target, canonical_edge(left, right)))

    # directed cycles = SCCs with >1 node or a self-loop
    cyclic_components = [
        scc for scc in nx.strongly_connected_components(D)
        if len(scc) > 1 or any(D.has_edge(n, n) for n in scc)
    ]
    component_by_node = {}
    for scc_id, scc in enumerate(cyclic_components):
        for node in scc:
            component_by_node[node] = scc_id

    cyclic_edges = set()
    for source, target, edge_id in observed_edges:
        if (source in component_by_node and target in component_by_node
                and component_by_node[source] == component_by_node[target]):
            cyclic_edges.add(edge_id)

    if len(cyclic_components) == 1 and len(cyclic_components[0]) == len(bicomponent):
        # whole component is one cycle -> copy-number model only
        genotype_graph = None
        dag_seed_node = start_component_node
    elif cyclic_components:
        # collapse each cycle into a single node "C{index}" -> DAG
        dag = nx.MultiDiGraph()
        node_mapping = {node: node for node in D.nodes}
        for index, scc in enumerate(cyclic_components):
            for node in scc:
                node_mapping[node] = f"C{index}"
        for source, target, key in D.edges(keys=True):
            new_source, new_target = node_mapping[source], node_mapping[target]
            if new_source != new_target:
                dag.add_edge(new_source, new_target, key=key)
        assert nx.is_directed_acyclic_graph(dag)
        genotype_graph = dag
        # the start node may have been collapsed into a cycle node
        dag_seed_node = node_mapping.get(start_component_node, start_component_node)
    else:
        genotype_graph = D
        dag_seed_node = start_component_node

    if genotype_graph is not None:
        incoming, outgoing, edge_depth_keys = precompute_component_traversal(genotype_graph)
    else:
        incoming = outgoing = edge_depth_keys = None

    return {
        "component_index": meta["component_index"],
        "bicomponent": bicomponent,
        "start_component_node": start_component_node,
        "end_component_node": meta["end_component_node"],
        "start_component_node_ploidy": 2,
        "dag_seed_node": dag_seed_node,
        "incoming": incoming,
        "outgoing": outgoing,
        "edge_depth_keys": edge_depth_keys,
        "cyclic_edges": cyclic_edges,
    }


# ---------------------------------------------------------------------------
# Step 6: recovery of edges the component genotyping cannot reach
# ---------------------------------------------------------------------------
def reachable_dag_edges(incoming, outgoing, seed):
    """
    Edges the genotyping BFS (genotype_model.genotype_dag) will assign a dosage
    to. Which edges are reached depends only on the DAG, not on sample depths:
    every outgoing edge of a processed node is genotyped, and a node is
    processed once all of its incoming edges are genotyped.
    """
    typed = set()
    known = {seed}
    processed = set()
    queue = deque([seed])
    while queue:
        node = queue.popleft()
        if node in processed:
            continue
        node_outgoing = outgoing.get(node, [])
        typed.update(node_outgoing)
        processed.add(node)
        for edge in node_outgoing:
            target = edge[1]
            if target not in known and all(e in typed for e in incoming.get(target, [])):
                known.add(target)
                queue.append(target)
    return typed


def genotyped_edge_set(components):
    """Canonical ids of every edge the component genotyping can output."""
    genotyped = set()
    for entry in components:
        genotyped.update(entry["cyclic_edges"])
        if entry["incoming"] is not None:
            for edge in reachable_dag_edges(entry["incoming"], entry["outgoing"], entry["dag_seed_node"]):
                genotyped.add(entry["edge_depth_keys"][edge])
    return genotyped


def nearest_backbone_node(start_nodes, neighbors, b38_node_order):
    """Multi-source BFS from start_nodes; first GRCh38 node reached, or None."""
    visited = set(start_nodes)
    queue = deque(start_nodes)
    while queue:
        node = queue.popleft()
        if node in b38_node_order:
            return node
        for neighbor in neighbors(node):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    return None


def build_recovery_groups(missing_edges, neighbors, b38_node_order):
    """
    Group missing edges by connected component of the graph they induce and
    pick one GRCh38 anchor per group.

    Returns (groups, n_unresolved_edges, n_direct, n_bfs) where groups is
    [(anchor_node, [canonical_edge_id, ...]), ...].
    """
    missing_graph = nx.Graph()
    edges_by_node = defaultdict(list)
    for a, b, canon in missing_edges:
        missing_graph.add_edge(a, b)
        edges_by_node[a].append((a, b, canon))
        edges_by_node[b].append((a, b, canon))

    groups = []
    n_unresolved_edges = n_direct = n_bfs = 0
    for component_nodes in nx.connected_components(missing_graph):
        group_edges = sorted({e for node in component_nodes for e in edges_by_node[node]})
        backbone_nodes = [n for n in component_nodes if n in b38_node_order]
        if backbone_nodes:
            anchor = min(backbone_nodes, key=lambda n: b38_node_order[n])
            n_direct += 1
        else:
            anchor = nearest_backbone_node(sorted(component_nodes), neighbors, b38_node_order)
            if anchor is None:
                n_unresolved_edges += len(group_edges)
                continue
            n_bfs += 1
        groups.append((anchor, [canon for _a, _b, canon in group_edges]))
    return groups, n_unresolved_edges, n_direct, n_bfs


def add_recovery_groups(structure, all_edges, neighbors):
    """Compute and store structure["recovery_groups"]; returns a log message."""
    genotyped = genotyped_edge_set(structure["components"])
    missing_edges = [e for e in all_edges if e[2] not in genotyped]
    groups, n_unresolved, n_direct, n_bfs = build_recovery_groups(
        missing_edges, neighbors, structure["b38_node_order"]
    )
    structure["recovery_groups"] = groups
    return (f"{len(missing_edges)} of {len(all_edges)} edge(s) not reachable by component genotyping; "
            f"{len(groups)} recovery group(s) ({n_direct} with own GRCh38 node, {n_bfs} anchored by BFS); "
            f"{n_unresolved} edge(s) dropped (no GRCh38 node reachable)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Precompute per-chromosome graph structure (biconnected components, "
                    "cyclic edges, DAG traversal, recovery groups) for hard genotyping."
    )
    parser.add_argument("--chrom", required=True, help="Chromosome to process (e.g. chr21).")
    parser.add_argument("--gfa", required=True, help="Chromosome GFA file.")
    parser.add_argument(
        "--depthfile", required=True,
        help="Per-chromosome raw edge-depth matrix (Chrom, Edge, sample...); only Chrom/Edge are used."
    )
    parser.add_argument("--out", required=True, help="Output structure pickle.")
    args = parser.parse_args()

    # Step 1
    path_edges = read_grch38_backbone(args.gfa)   # GRCh38 node ids in path order
    grch38_node_set = set(path_edges)
    start_node = path_edges[0]
    end_node = path_edges[-1]
    b38_node_order = {node: index for index, node in enumerate(path_edges)}

    # Step 2
    G, raw_edge_lookup, all_edges = read_edge_graph(args.depthfile)
    removed_adj = remove_off_backbone_bridges(G, grch38_node_set)

    articulation_nodes = set(nx.articulation_points(G))
    articulation_nodes.update([start_node, end_node])

    # Step 3: classify biconnected components
    component_meta = []   # anchored components, finished in step 5
    components = []       # final per-component records
    n_skipped = n_skipped_edges = n_partial = 0

    for component_index, bicomponent in enumerate(nx.biconnected_components(G)):
        boundary_nodes = [n for n in bicomponent if n in articulation_nodes]
        backbone_boundary_nodes = [n for n in boundary_nodes if n in b38_node_order]

        # anchored: both boundary nodes on GRCh38 -> ploidy-2 flow DAG (step 5)
        if len(boundary_nodes) == 2 and len(backbone_boundary_nodes) == 2:
            i, j = sorted(b38_node_order[n] for n in boundary_nodes)
            component_meta.append({
                "component_index": component_index,
                "bicomponent": bicomponent,
                "start_component_node": path_edges[i],
                "end_component_node": path_edges[j],
            })
            continue

        # partially anchored: copy-number model for every edge, local depth
        # baseline taken around the GRCh38 boundary node with the lowest position
        if backbone_boundary_nodes:
            anchor_node = min(backbone_boundary_nodes, key=lambda n: b38_node_order[n])
            cyclic_edges = set()
            for u, v in G.subgraph(bicomponent).edges():
                handles = parse_path(raw_edge_lookup[(u, v) if u <= v else (v, u)])
                cyclic_edges.add(canonical_edge(handles[0], handles[1]))
            n_partial += 1
            components.append({
                "component_index": component_index,
                "bicomponent": bicomponent,
                "start_component_node": anchor_node,
                "end_component_node": anchor_node,
                "start_component_node_ploidy": 2,
                "dag_seed_node": anchor_node,
                "incoming": None,
                "outgoing": None,
                "edge_depth_keys": None,
                "cyclic_edges": cyclic_edges,
            })
            continue

        # no GRCh38 boundary node -> handled by recovery (step 6)
        n_skipped += 1
        n_skipped_edges += G.subgraph(bicomponent).number_of_edges()

    # node -> (start, end) pair(s) of the anchored component(s) containing it
    node_to_pairs = defaultdict(list)
    for meta in component_meta:
        pair = (meta["start_component_node"], meta["end_component_node"])
        for node in meta["bicomponent"]:
            node_to_pairs[node].append(pair)

    # Step 4
    backbone_articulation_nodes = [n for n in articulation_nodes if n in b38_node_order]
    path_split_by_articulation = split_paths(
        read_walks(args.gfa), backbone_articulation_nodes, b38_node_order, node_to_pairs
    )

    # first edge out of / last edge into each articulation node, used to
    # compute the per-sample depth around each node
    outgoing_paths = defaultdict(set)
    incoming_paths = defaultdict(set)
    for (pair_a, pair_b), segments in path_split_by_articulation.items():
        upstream, downstream = sorted((pair_a, pair_b), key=lambda n: b38_node_order[n])
        for segment in segments:
            handles = parse_path(segment)
            outgoing_paths[upstream].add(canonical_edge(handles[0], handles[1]))
            if downstream == end_node:
                incoming_paths[downstream].add(canonical_edge(handles[-2], handles[-1]))

    stem_positions = sorted(b38_node_order[n] for n in articulation_nodes if n in b38_node_order)

    # Step 5
    for meta in component_meta:
        segments = path_split_by_articulation.get(
            (meta["start_component_node"], meta["end_component_node"]), set()
        )
        components.append(build_anchored_component(meta, segments))

    structure = {
        "chr": args.chrom,
        "path_edges": path_edges,
        "b38_node_order": b38_node_order,
        "start_node": start_node,
        "end_node": end_node,
        "articulation_nodes": articulation_nodes,
        "outgoing_paths": outgoing_paths,
        "incoming_paths": incoming_paths,
        "stem_positions": stem_positions,
        "components": components,
    }

    # Step 6: BFS through the full edge graph (bridges removed in step 2 included)
    def neighbors(node):
        return list(G.neighbors(node)) + removed_adj.get(node, []) if node in G else removed_adj.get(node, [])

    recovery_msg = add_recovery_groups(structure, all_edges, neighbors)

    with open(args.out, "wb") as fh:
        pickle.dump(structure, fh, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"{args.chrom}: {len(components)} component(s) ({n_partial} partially anchored); "
          f"{n_skipped} component(s) / {n_skipped_edges} edge(s) without a GRCh38 boundary node")
    print(f"{args.chrom}: {recovery_msg}")
    print(f"{args.chrom}: wrote {args.out}")


if __name__ == "__main__":
    main()
