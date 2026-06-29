"""
find allele traversal (AT) and edge for snarl not in vcf 

Input: 
    # path each chr 
    /vast/palmer/scratch/hall/wl474/projects/hprc_r2/edge_read_depth/chroms/chr18.gfa
    # all snarls after checking using vcf file 
    /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/chr1.snarl_ps_lv_at_edge.txt.gz

Output:
    /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/my_filter/chr1.snarl_ps_lv_at_edge.not_in_vcf.txt

Usage:
    python find.not_in_vcf.AT_edge.py chr1 
"""

"""
Replace biallelic edges by allele balance: step 2 - find allele traversal (AT)
and edges for snarls that had no traversal in the VCF (step 1 output), by
locating the snarl's start/end nodes along reference/haplotype paths (W lines)
in the GFA.

Input:
    --snarl_info : step 1 output (snarl_info_fromvcf.py): snarl, parent, level,
                   traversal, edges - one row per snarl ("None" if no VCF traversal) 

    --gfa        : chromosome GFA file; paths are read from "W" lines 

Output:
    --out : snarl, parent, level, traversal(s), edge(s) for snarls that had no
            VCF traversal (same columns as step 1; "None" if still not found) 
Usage:
    python snarl_info_fromgfa.py \\
        --snarl_info chr1.snarl_ps_lv_at_edge.txt.gz \\
        --gfa        chr1.gfa \\
        --out        chr1.snarl_ps_lv_at_edge.not_in_vcf.txt
"""
import argparse
import re
from collections import defaultdict

import pandas as pd

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Find allele traversals/edges from GFA paths for snarls absent from the VCF."
)
parser.add_argument(
    "--snarl_info", required=True,
    help="Step 1 output file (snarl_info_fromvcf.py): snarl, parent, level, traversal, edges."
)
parser.add_argument(
    "--gfa", required=True,
    help="Chromosome GFA file; paths are read from 'W' lines."
)
parser.add_argument(
    "--out", required=True,
    help="Output file: snarl, parent, level, traversal(s), edge(s)."
)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Read step 1 snarl info, keep only snarls with no VCF traversal (4th col = None)
# ---------------------------------------------------------------------------
snarl_df = pd.read_csv(args.snarl_info, sep="\t", header=None)
# snarl PS LV AT edges
# 45308515-45308517       45308511-45308541       10      >45308515>45308517,>45308515>45308516>45308517  45308515-45308517,45308515-45308516,45308516-45308517
snarl_novcf = snarl_df[snarl_df.iloc[:, 3].isna()].copy() # # select not in vcf snarls (no AT snarls), 4th column = None 
snarl_novcf.columns = ['snarl', 'ps', 'lv', 'at', 'edges']

# ---------------------------------------------------------------------------
# Read GFA paths (W lines), find traversal between each snarl's start/end node
# ---------------------------------------------------------------------------
snarl_novcf_traversal_dict = defaultdict(list)

with open(args.gfa, 'r') as gfa_file:
    # loop through each path to append snarl traversal 
    for line in gfa_file:
        # read in path each chr 
        if not line.startswith('W'):
            continue
        # record node, index and direction 
        path = line.strip().split('\t')[6] # W       GRCh38  0       chr9    0       138394717       >186851125>186851126>18685
        node_in_path = re.findall(r'[><]\d+', path)
        node_index_dir_dict = {}
        for n, node in enumerate(node_in_path):
            node_index_dir_dict[node[1:]] = [n, node[0]] # 186851125: 0, > 

        for snarl_id in snarl_novcf['snarl']: # 45308515-45308517
            start_node, end_node = snarl_id.split('-')
            # find path between start and end node of snarl >45308515<45308516>45308517
            if start_node not in node_index_dir_dict or end_node not in node_index_dir_dict:
                continue
            # find traversal between start and end node 
            start_index = node_index_dir_dict[start_node][0]
            end_index = node_index_dir_dict[end_node][0]
            if start_index > end_index:
                start_index, end_index = end_index, start_index

            snarl_traversal = "".join(node_in_path[start_index:end_index + 1])
            snarl_novcf_traversal_dict[snarl_id].append(snarl_traversal)

for snarl_id in snarl_novcf_traversal_dict:
    snarl_novcf_traversal_dict[snarl_id] = set(snarl_novcf_traversal_dict[snarl_id])
# add snarl traversal to snarl_novcf
snarl_novcf['at'] = snarl_novcf['snarl'].map(snarl_novcf_traversal_dict)

# ---------------------------------------------------------------------------
# Derive edges from each snarl's traversal(s) and write output
# ---------------------------------------------------------------------------
with open(args.out, "w") as out_f:
    for _, row in snarl_novcf.iterrows(): # sorted # snarl_id, parent_id  level traversal
        allele_traversal = row["at"]
        if pd.isna(allele_traversal): # no traversal 
            out_f.write(f"{row['snarl']}\t{row['ps']}\t{row['lv']}\tNone\tNone\n")
            continue
        # if traversal is not NaN
        edges = set()
        for traversal in allele_traversal:
            nodes = re.split(r'>|<', traversal)[1:] # split traversal into node by either > or <
            for node_i in range(len(nodes) - 1):
                n1, n2 = nodes[node_i], nodes[node_i + 1]
                edge = f"{n1}-{n2}" if n1 < n2 else f"{n2}-{n1}"
                edges.add(edge)
        # output snarl_id, parent_id, level, traversal, edge
        traversal_out = ",".join(allele_traversal)
        edge_out = ",".join(edges)
        out_f.write(f"{row['snarl']}\t{row['ps']}\t{row['lv']}\t{traversal_out}\t{edge_out}\n")
