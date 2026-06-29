"""
Replace biallelic edges by allele balance: step 4 - for snarls forming a
triangle (3-edge) or diamond (4-edge) structure on the GRCh38 path, find the
single non-reference edge that survived the bubble-puncture redundancy filter
and replace it with its allele balance against the (kept) GRCh38 edge(s).

A snarl qualifies if:
  - its edges form a 3-edge ("triangle") or 4-edge ("diamond") structure with
    exactly 2 distinct (direction-agnostic) traversals, and
  - enough of its nodes lie on the GRCh38 path (>=2 for triangle, ==3 for diamond), and
  - exactly one of its non-reference edge(s) was kept by the bubble-puncture filter.

For a qualifying snarl, the allele balance per sample is:
    AB = depth(kept non-reference edge) / (depth(kept non-reference edge) + depth(GRCh38 edge))
(for the diamond case, the GRCh38 edge with higher average depth is used)

Input:
    --edge_raw     : raw per-sample edge depth count file (after redundancy filtering)
    --gfa          : chromosome GFA file; GRCh38 path nodes are read from the W line
                      where column 2 == "GRCh38" (path string in column 7)
    --snarl_vcf    : step 3 output for snarls with a VCF-derived traversal
    --snarl_notvcf : step 3 output for snarls without a VCF-derived traversal
    --kept_edges   : edge list kept after the bubble-puncture redundancy filter

Output:
    --out : edge_id, allele_balance_per_sample (tab-separated)

Usage:
    python compute_edge_nonredundant_AB.py \\
        --edge_raw     chr1.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.txt.gz \\ 
        --gfa          chr1.gfa \\
        --snarl_vcf    chr1.snarl_ps_lv_at_edge.high_lv_edge.txt.gz \\
        --snarl_notvcf chr1.snarl_ps_lv_at_edge.not_in_vcf.high_lv_edge.txt.gz \\
        --kept_edges   chr1.independent_edge_list.txt \\
        --out          chr1.bp_allele_balance.txt
"""

import argparse
import gzip
import re
from itertools import chain

def open_maybe_gzip(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Find biallelic snarl edges and compute their allele balance against the GRCh38 edge."
)
parser.add_argument(
    "--edge_raw", required=True,
    help="Raw per-sample edge depth count file."
)
parser.add_argument(
    "--gfa", required=True,
    help="Chromosome GFA file; GRCh38 path nodes are read from the W line with column 2 == 'GRCh38'."
)
parser.add_argument(
    "--snarl_vcf", required=True,
    help="Step 3 output for snarls with a VCF-derived traversal."
)
parser.add_argument(
    "--snarl_notvcf", required=True,
    help="Step 3 output for snarls without a VCF-derived traversal."
)
parser.add_argument(
    "--kept_edges", required=True,
    help="Edge list kept after the bubble-puncture redundancy filter."
)
parser.add_argument(
    "--out", required=True,
    help="Output file: edge_id, allele_balance_per_sample."
)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Read raw edge depths
# ---------------------------------------------------------------------------
edge_dict = {}        # edge -> list of per-sample depth strings
node_edge_dict = {}    # (node1, node2) canonical -> original edge string

with open_maybe_gzip(args.edge_raw) as edge_raw_file:
    for line in edge_raw_file:
        parts = line.strip().split('\t')
        edge = parts[1]
        depth = parts[2:]
        edge_dict[edge] = depth
        nodes = re.split(r">|<", edge)[1:]
        if nodes[0] > nodes[1]:
            nodes[0], nodes[1] = nodes[1], nodes[0]
        node_edge_dict[(nodes[0], nodes[1])] = edge

n_samples = len(next(iter(edge_dict.values())))

edge_ave_depth = {
    edge: sum(int(d) for d in depths) / len(depths)
    for edge, depths in edge_dict.items()
} # '>69345828<69345829': 21.89069

# ---------------------------------------------------------------------------
# Read GRCh38 path nodes from the GFA
# ---------------------------------------------------------------------------
b38_node = set()
with open_maybe_gzip(args.gfa) as gfa_file:
    for line in gfa_file:
        if not line.startswith('W'):
            continue
        parts = line.rstrip('\n').split('\t')
        if parts[1] != "GRCh38":
            continue
        path = parts[6]
        b38_node.update(re.split(r">|<", path)[1:])


# ---------------------------------------------------------------------------
# Read kept edges (post bubble-puncture redundancy filter)
# ---------------------------------------------------------------------------
kept_edge_bp_set = set()
with open_maybe_gzip(args.kept_edges) as kept_edge_bp:
    for line in kept_edge_bp:
        kept_edge_bp_set.add(line.strip())


def edge_format_change(edge):
    """'30191785-30191786' -> '>30191785>30191786' (or whatever direction is on record)."""
    nodes = edge.split('-')
    if nodes[0] > nodes[1]:
        nodes[0], nodes[1] = nodes[1], nodes[0]
    return node_edge_dict[tuple(nodes)]


def is_desired_structure(snarl_line):
    """Return 4 for a diamond structure, 3 for a triangle, 0 otherwise."""
    snarl_info = snarl_line.strip().split('\t')
    # snarl ps lv at edge uniq_edge
    try:
        traversal = snarl_info[3].split(',')
        traversal_nodirection = [sorted(re.split(r">|<", i)[1:]) for i in traversal] # [[1,2,4],[1,2,3]]
        traversal_nodirection_uniq = list(map(list, set(tuple(s) for s in traversal_nodirection)))
        # if the snarl edges are in a diamond structure
        # >30192808>30192809>30192811,>30192808>30192810>30192811
        if (len(snarl_info[4].split(',')) == 4 and len(snarl_info[5].split(',')) == 4
                and len(traversal_nodirection_uniq) == 2
                and len(traversal_nodirection_uniq[0]) == 3 and len(traversal_nodirection_uniq[1]) == 3):
            return 4
        # if the snarl edges are in a triangle structure
        # >30192811>30192813,>30192811>30192812>30192813 
        elif (len(snarl_info[4].split(',')) == 3 and len(snarl_info[5].split(',')) == 3
                and len(traversal_nodirection_uniq) == 2):
            return 3
        else:
            return 0
    except Exception:
        return 0


def allele_balance(var_depth, b38_depth):
    return [
        0 if (var_depth[i] + b38_depth[i]) == 0 else var_depth[i] / (var_depth[i] + b38_depth[i])
        for i in range(n_samples)
    ]


# ---------------------------------------------------------------------------
# Find biallelic snarl edges and compute allele balance
# ---------------------------------------------------------------------------
with open_maybe_gzip(args.snarl_vcf) as snarl_invcf, \
     open_maybe_gzip(args.snarl_notvcf) as snarl_notinvcf, \
     open(args.out, 'w') as out_file:

    for line in chain(snarl_invcf, snarl_notinvcf):
        line_split = line.strip().split('\t')
        structure = is_desired_structure(line)

        if structure == 4:
            nodes = set(x for i in line_split[5].split(',') for x in i.split('-'))
            # not in GRCh38 path, skip
            if sum(n in b38_node for n in nodes) != 3:
                continue
            # if the snarl is at GRCh38 path
            mid_nodes = nodes - set(line_split[0].split('-'))
            mid_node1 = mid_nodes.pop()
            mid_node2 = mid_nodes.pop()
            start_node, end_node = line_split[0].split('-')
            # find the kept edge 
            if mid_node1 in b38_node:
                # edge_var1 = node_edge_dict[tuple(sorted((start_node, mid_node2)))] # >30192808>30192810
                edge_var1 = node_edge_dict.get(tuple(sorted((start_node, mid_node2))), None) # >30192808>30192810
                edge_var2 = node_edge_dict.get(tuple(sorted((mid_node2, end_node))), None)
                edge_b38_1 = node_edge_dict.get(tuple(sorted((start_node, mid_node1))), None)
                edge_b38_2 = node_edge_dict.get(tuple(sorted((mid_node1, end_node))), None)
            else:
                edge_var1 = node_edge_dict.get(tuple(sorted((start_node, mid_node1))), None)
                edge_var2 = node_edge_dict.get(tuple(sorted((mid_node1, end_node))), None)
                edge_b38_1 = node_edge_dict.get(tuple(sorted((start_node, mid_node2))), None)
                edge_b38_2 = node_edge_dict.get(tuple(sorted((mid_node2, end_node))), None)
            # if more than 1 edge been kept, cotinue 
            if sum([edge_var1 in kept_edge_bp_set, edge_var2 in kept_edge_bp_set]) != 1:
                continue
            
            if edge_var1 in kept_edge_bp_set:
                edge_var = edge_var1
            else:
                edge_var = edge_var2
            edge_var_depth = [int(i) for i in edge_dict[edge_var]]
            # choose the b38 edge with higher depth 
            edge_b38 = edge_b38_1 if edge_ave_depth[edge_b38_1] > edge_ave_depth[edge_b38_2] else edge_b38_2
            edge_b38_depth = [int(i) for i in edge_dict[edge_b38]]

            ab_line = '\t'.join(f"{x:.6f}" for x in allele_balance(edge_var_depth, edge_b38_depth))
            out_file.write(f"{edge_var}\t{ab_line}\n")

        elif structure == 3:
            nodes = set(x for i in line_split[5].split(',') for x in i.split('-'))
            # if the snarl is not at GRCh38 path
            if sum(n in b38_node for n in nodes) < 2:
                continue

            mid_node = (nodes - set(line_split[0].split('-'))).pop()
            start_node, end_node = line_split[0].split('-')
            # find the kept edge
            if mid_node in b38_node:
                edge_var = node_edge_dict.get(tuple(sorted((start_node, end_node))), None)
                edge_b38_1 = node_edge_dict.get(tuple(sorted((start_node, mid_node))), None)
                edge_b38_2 = node_edge_dict.get(tuple(sorted((mid_node, end_node))), None)

                if edge_var not in kept_edge_bp_set:
                    continue

                edge_var_depth = [int(i) for i in edge_dict[edge_var]]
                edge_b38 = edge_b38_1 if edge_ave_depth[edge_b38_1] > edge_ave_depth[edge_b38_2] else edge_b38_2
                edge_b38_depth = [int(i) for i in edge_dict[edge_b38]]
            else:
                edge_var_1 = node_edge_dict.get(tuple(sorted((start_node, mid_node))), None)
                edge_var_2 = node_edge_dict.get(tuple(sorted((mid_node, end_node))), None)
                edge_b38 = node_edge_dict.get(tuple(sorted((start_node, end_node))), None)

                if sum([edge_var_1 in kept_edge_bp_set, edge_var_2 in kept_edge_bp_set]) != 1:
                    continue

                edge_var = edge_var_1 if edge_var_1 in kept_edge_bp_set else edge_var_2
                edge_var_depth = [int(i) for i in edge_dict[edge_var]]
                edge_b38_depth = [int(i) for i in edge_dict[edge_b38]]

            ab_line = '\t'.join(f"{x:.6f}" for x in allele_balance(edge_var_depth, edge_b38_depth))
            out_file.write(f"{edge_var}\t{ab_line}\n")
                    



