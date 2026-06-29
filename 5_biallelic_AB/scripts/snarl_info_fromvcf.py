"""
Replace biallelic edges by allele balance: step 1 - snarl/traversal/edge mapping. (based on vcf file)

For one chromosome, builds the snarl hierarchy (parent/child, level from root),
maps each snarl to its allele traversals from the VCF (falling back to the
nearest ancestor snarl's traversal when a snarl has none directly), and derives
the edge set implied by each traversal.

Input:
    --vcf    : VCF recording allele traversals (AT INFO field), one chromosome 
    --snarls : snarl JSON (genome-wide; filtered here to nodes in --gfa) 
    --gfa    : chromosome GFA file; node IDs are read from "S" lines (column 2) 

Output:
    --out : snarl_id, parent_id, level, traversal(s), edge(s) - one row per snarl
            (traversal/edge are comma-separated; "None" if no traversal found) 

Usage:
    python snarl_info_fromvcf.py \\
        --vcf    hprc-v2.0-mc-grch38.raw.chr22.sorted.vcf.gz \\
        --snarls hprc-v2.0-mc-grch38.snarls.json \\
        --gfa    chr22.gfa \\
        --out    chr22.snarl_ps_lv_at_edge.txt
"""

import argparse
import json
import re
from collections import defaultdict, deque

import pandas as pd
from pysam import VariantFile


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Map snarls to allele traversals and edges for one chromosome."
)
parser.add_argument(
    "--vcf", required=True,
    help="VCF file recording allele traversals (AT INFO field) for this chromosome."
)
parser.add_argument(
    "--snarls", required=True,
    help="Snarl JSON file (genome-wide); filtered to snarls touching --gfa nodes."
)
parser.add_argument(
    "--gfa", required=True,
    help="Chromosome GFA file; node IDs are read from 'S' lines (column 2)."
)
parser.add_argument(
    "--out", required=True,
    help="Output file: snarl_id, parent_id, level, traversal(s), edge(s)."
)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Read VCF: snarl_id -> allele traversals
# ---------------------------------------------------------------------------
input_vcf = VariantFile(args.vcf)
# record snarl and their allele traversals from vcf 
snarl_traversal = {}
for snarl in input_vcf.fetch():
    snarl_traversal[snarl.id] = snarl.info["AT"]  # tuple of traversal strings
    # >178211212>178211214: (178211212>178211213>178211214,>178211212>178211214)

# ---------------------------------------------------------------------------
# Read node IDs for this chromosome from the GFA's "S" (segment) lines
# ---------------------------------------------------------------------------
chr_node = set()
with open(args.gfa, 'r') as gfa_file:
    for line in gfa_file:
        if line.startswith('S'):
            # S  <node_id>  <sequence>  ...
            chr_node.add(line.rstrip('\n').split('\t')[1])  # node id

# ---------------------------------------------------------------------------
# Read snarl JSON, filter to snarls touching this chromosome's nodes
# ---------------------------------------------------------------------------
rows = []
with open(args.snarls) as f:
    for line in f:
        # {'directed_acyclic_net_graph': True, 'end': {'node_id': '77308931'}, 'start': {'node_id': '77308930'}, 'start_end_reachable': True, 'type': 1}
        # {"directed_acyclic_net_graph": true, "end": {"node_id": "64560279"}, "end_self_reachable": true, "parent": {"end": {"node_id": "64580676"}, "start": {"node_id": "64556358"}}, "start": {"node_id": "64557435"}, "start_end_reachable": true, "start_self_reachable": true}
        snarl = json.loads(line)
        snarl_start = snarl["start"]["node_id"]
        snarl_end = snarl["end"]["node_id"]
        if snarl_start not in chr_node and snarl_end not in chr_node:
            continue
        try:
            parent_snarl_start = snarl["parent"]["start"]["node_id"]
            parent_snarl_end = snarl["parent"]["end"]["node_id"]
            # smaller node id first for patrent snarl 
            if parent_snarl_start > parent_snarl_end:
                parent_snarl_start, parent_snarl_end = parent_snarl_end, parent_snarl_start
            # smaller node id first for snarl
            if snarl_start > snarl_end:
                snarl_start, snarl_end = snarl_end, snarl_start
            row = {"snarl_id": f'{snarl_start}-{snarl_end}', "parent_id": f'{parent_snarl_start}-{parent_snarl_end}'}
        except KeyError:
            if snarl_start > snarl_end:
                snarl_start, snarl_end = snarl_end, snarl_start
            row = {"snarl_id": f'{snarl_start}-{snarl_end}', "parent_id": None}
        rows.append(row)

input_snarl_df = pd.DataFrame(rows)  # snarl_id, parent_id

# ---------------------------------------------------------------------------
# Compute snarl level by BFS from roots (parent-less snarls) downward
# ---------------------------------------------------------------------------
children_map = defaultdict(list)  # parent -> children
parents = {}  # snarl -> parent (or None)

for _, row in input_snarl_df.iterrows():
    snarl_id = row['snarl_id']
    parent_id = row['parent_id']
    parents[snarl_id] = parent_id
    if parent_id is not None:
        children_map[parent_id].append(snarl_id)

# Compute levels by BFS 
levels = {}
queue = deque()

# Start with roots (snarls with no parent)
for snarl in input_snarl_df['snarl_id']:
    if parents[snarl] is None:
        levels[snarl] = 0
        queue.append((snarl, 0)) # lvl = 0

while queue:
    current, lvl = queue.popleft()
    for child in children_map[current]:
        levels[child] = lvl + 1
        queue.append((child, lvl + 1))

# add level to input_snarl_df
input_snarl_df['level'] = input_snarl_df['snarl_id'].map(levels)

# ---------------------------------------------------------------------------
# Map each snarl to allele traversals: direct from VCF, or derived from the
# nearest ancestor snarl's traversal if this snarl has none directly
# ---------------------------------------------------------------------------
# find traversals for each snarl - based on vcf 
traversals = {}

for snarl in input_snarl_df['snarl_id']: # 111362662-111362666
    snarl_start, snarl_end = snarl.split("-")
    snarl_id_vcf = f'>{snarl_start}>{snarl_end}'
    snarl_id_vcf_rev = f'>{snarl_end}>{snarl_start}'

    if snarl_id_vcf in snarl_traversal or snarl_id_vcf_rev in snarl_traversal:
        try: 
            traversal = snarl_traversal[snarl_id_vcf]
        except KeyError:
            traversal = snarl_traversal[snarl_id_vcf_rev]
        traversals[snarl] = traversal
        continue

    # snarl not in VCF directly: walk up ancestors until one is found in the VCF
    ancestors = []
    current = snarl
    while current in parents and parents[current] is not None:
        ancestors.append(parents[current])
        current = parents[current]

    # loop through ancestors
    for ancestor in ancestors:
        anc_start, anc_end = ancestor.split("-")
        ancestor_id_vcf = f'>{anc_start}>{anc_end}'
        ancestor_id_vcf_rev = f'>{anc_end}>{anc_start}'
        
        if ancestor_id_vcf not in snarl_traversal and ancestor_id_vcf_rev not in snarl_traversal:
            continue

        # if ancestor is in vcf - find traversal from the ancestor traversal
        try: 
            traversal = snarl_traversal[ancestor_id_vcf]
        except KeyError:
            traversal = snarl_traversal[ancestor_id_vcf_rev]
        # extract this snarl's sub-traversal from the ancestor's traversal
        short_traversal = []
        for t in traversal:
            s_start, s_end = snarl_start, snarl_end
            # if snarl is in traversal, add it to short_traversal
            if s_start in t and s_end in t:
                start_index = t.find(s_start)
                end_index = t.find(s_end)
                if end_index <= start_index:
                    start_index, end_index = end_index, start_index
                    s_start, s_end = s_end, s_start
                # idividual snarl traversal
                short_traversal.append(t[start_index - 1:end_index + len(s_end)])
        # find unique traversals in short_traversal
        traversals[snarl] = set(short_traversal)
        break  # stop at the nearest ancestor found in the VCF
    # else: snarl is not in vcf and its recursive parent not in vcf, ignore in this script

# add traversal to input_snarl_df
input_snarl_df['traversal'] = input_snarl_df['snarl_id'].map(traversals)  # NaN if no traversal found


# ---------------------------------------------------------------------------
# Derive edges from each snarl's traversal(s); write output sorted by level
# (deepest first) so child edges are recorded before their ancestors
# ---------------------------------------------------------------------------
input_snarl_df_sort = input_snarl_df.sort_values(by=["level"], ascending=False)

edge_set = set()  # edges seen so far across all snarls

with open(args.out, "w") as out_f:
    for _, row in input_snarl_df_sort.iterrows():  # sorted # snarl_id, parent_id  level traversal
        allele_traversal = row["traversal"]
        if pd.isna(allele_traversal):
            out_f.write(f"{row['snarl_id']}\t{row['parent_id']}\t{row['level']}\tNone\tNone\n")
            continue

        edges = []
        for traversal in allele_traversal:
            # split traversal into node by either > or <
            nodes = re.split(r'>|<', traversal)[1:] # ['8177665', '8177666', '8177667']
            for node_i in range(len(nodes) - 1):
                n1, n2 = nodes[node_i], nodes[node_i + 1]
                edge = f"{n1}-{n2}" if n1 < n2 else f"{n2}-{n1}"
                if edge not in edges:
                    edges.append(edge)
                    edge_set.add(edge)

        traversal_out = ",".join(allele_traversal)
        edge_out = ",".join(edges)
        out_f.write(f"{row['snarl_id']}\t{row['parent_id']}\t{row['level']}\t{traversal_out}\t{edge_out}\n")


