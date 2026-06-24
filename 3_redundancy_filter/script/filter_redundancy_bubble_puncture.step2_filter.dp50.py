"""
Step 2 of redundancy filtering via the bubble-puncture method.

** speed up filter process by only check subgraph within 50 steps **

Usage:
    python filter_redundancy_bubble_puncture.step2_filter.dp50.py \\
        --gml   chr21_subgraphs/subgraph_1.graph.gml \\ # subgraph
        --depth chr21.hprc-v2.0-mc-grch38.depth_per_edge.txt(.gz) \\ # to calculate avg depth per edge 
        --gfa   chr21.gfa \\ # to identify GRCh38 edges
        --out   chr21_results/subgraph_1.independent_edge_list.dp50.txt

Output:
    --out: one edge per line in original format (e.g. >12345<67890),
           listing all independent (kept) non-reference edges in this subgraph.

"""
import argparse
import gzip
import os
import re

import networkx as nx



# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Apply bubble-puncture filter to one biconnected-component subgraph."
)
parser.add_argument(
    "--gml", required=True,
    help="Subgraph GML file produced by Step 1."
)
parser.add_argument(
    "--depth", required=True,
    help="Edge depth file (plain text or .gz); columns: chr, edge_id, sample1, sample2, ..."
)
parser.add_argument(
    "--gfa", required=True,
    help="GFA file for the chromosome; GRCh38 reference walk is read from W lines."
)
parser.add_argument(
    "--out", required=True,
    help="Output file for independent (kept) non-reference edges."
)
args = parser.parse_args()

if os.path.dirname(args.out):
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

# prefix is used to detect the subgraph_remaining case
prefix = os.path.basename(args.gml).split('.')[0] # subgraph_100


# ---------------------------------------------------------------------------
# Load subgraph
# ---------------------------------------------------------------------------
G = nx.read_gml(args.gml)


# ---------------------------------------------------------------------------
# Read GRCh38 reference edges from GFA walk lines
# ---------------------------------------------------------------------------
grch38_edge_set = set()  # to store GRCh38 edges, canonical (node1, node2) tuples, smaller node first

with open(args.gfa, 'r') as gfa_file:
    for line in gfa_file:
        if not line.startswith('W'):
            continue
        parts = line.strip().split('\t')
        if parts[1] != "GRCh38":
            continue
        path = parts[6]
        # path is like: >102643231>102643232>102643233>102643234>102643235
        path_edges = re.split(">|<", path)[1:]
        # recording all edges in the path
        for i in range(len(path_edges) - 1):
            # smaller node first 
            node1 = path_edges[i]
            node2 = path_edges[i+1]
            if node1 > node2:
                node1, node2 = node2, node1
            grch38_edge_set.add((node1, node2)) # (102643231,102643232)


# ---------------------------------------------------------------------------
# Read edge depths for edges present in this subgraph
# ---------------------------------------------------------------------------
# edges that in subgraph 
graph_edge_list = set(G.edges())
edge_ave_depth = {}  # original edge string -> average depth across samples
node_edge_dict = {}  # (node1, node2) -> original edge string e.g. '>12345<67890'

opener = gzip.open if args.depth.endswith('.gz') else open
with opener(args.depth, 'rt') as edge_raw_file:
    # #all edges read depth raw count file 
    for line in edge_raw_file:
        parts = line.strip().split('\t')
        edge = parts[1]
        depth = parts[2:]
        nodes = re.split(r">|<", edge)[1:]
        # canonical order: smaller node first
        if nodes[0] > nodes[1]:
            nodes[0], nodes[1] = nodes[1], nodes[0]
        # if edge in graph_edge_list, then save the edge and its depth
        if (nodes[0], nodes[1]) in graph_edge_list or (nodes[1], nodes[0]) in graph_edge_list:
            edge_ave_depth[edge] = sum(int(d) for d in depth) / len(depth) # '>69345828<69345829': 21.89069
            node_edge_dict[(nodes[0], nodes[1])] = edge # ('69345828','69345829'): '>69345828<69345829'

# ---------------------------------------------------------------------------
# apply bubble-puncture filter to subgraph - within 50 steps
# ---------------------------------------------------------------------------

# find edges to consider (not in GRCh38 edges) 
edge_to_consider = set(node_edge_dict.keys()) - grch38_edge_set # ('69345828','69345829')
# order edges by average depth, descending
edge_to_consider = sorted(edge_to_consider, key=lambda x: edge_ave_depth[node_edge_dict[x]], reverse=True)
# loop through edges in edge_to_consider, and find independent edges
with open(args.out, 'w') as out_file:
    while edge_to_consider:
        # pop the first edge 
        edge = edge_to_consider.pop(0)  # ('69345828','69345829')
        # Get the local neighborhood (BFS) around both endpoints: depth < 50
        start_node = edge[0]
        local_nodes = set(nx.single_source_shortest_path_length(G, start_node, cutoff=50).keys())
        # get local subgraph 
        local_G = G.subgraph(local_nodes)
        local_G = nx.Graph(local_G) # unfreeze the graph 
        # remove bridge edges in the local subgraph
        local_bridge_edges = set(nx.bridges(local_G))
        # delete edge if edge is in local bridge edges
        local_bridge_edges.discard(edge)  # remove the edge itself if it is a bridge edge
        local_bridge_edges.discard((edge[1], edge[0]))  # remove the reverse edge if it exists
        local_G.remove_edges_from(local_bridge_edges) 
        # remove the edge from the local graph 
        local_G.remove_edge(*edge)  # remove the edge itself
        G.remove_edge(*edge)  
        # write edge to output file
        out_file.write(f"{node_edge_dict[edge]}\n")  # '>69345828<69345829'
        out_file.flush()
        # find bridge edges in the remaining subgraph (dependent edges)
        graph_bridge_edges = list(nx.bridges(local_G)) 
        # remove bridge edges from edge_to_consider and graph 
        for bri_edge in graph_bridge_edges:
            G.remove_edge(*bri_edge)  # remove from graph
            if bri_edge in edge_to_consider:
                edge_to_consider.remove(bri_edge)
            elif (bri_edge[1], bri_edge[0]) in edge_to_consider:
                edge_to_consider.remove((bri_edge[1], bri_edge[0]))
            else:
                continue


            


