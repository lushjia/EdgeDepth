"""
Replace biallelic edges by allele balance: step 3 - output the edges in each
snarl that are not already claimed by one of its (deeper-level) child snarls.

Input is processed from deepest level to shallowest (the level-descending sort
order produced by step 1/2). A running set of edges already seen at deeper
levels is used to determine, for each snarl, which of its edges are newly
introduced at this level rather than inherited from a child.

This is run twice, independently, on the two step 1/2 outputs (snarls with a
VCF-derived traversal, and snarls whose traversal was instead recovered from
GFA paths)

Input:
    --snarl_vcf     : step 1 output (snarl_info_fromvcf.py) # /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/$chr.snarl_ps_lv_at_edge.txt.gz
    --snarl_notvcf  : step 2 output (snarl_info_fromgfa.py) # /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/hprc_v2/edge/my_filter/$chr.snarl_ps_lv_at_edge.not_in_vcf.txt

Output:
    --out_vcf     : --snarl_vcf rows, each with an added column of edges first
                     introduced at this snarl's level
    --out_notvcf  : same, for --snarl_notvcf

Usage:
    python snarl_highlv_edge.py \\
        --snarl_vcf    chr1.snarl_ps_lv_at_edge.txt.gz \\
        --snarl_notvcf chr1.snarl_ps_lv_at_edge.not_in_vcf.txt \\
        --out_vcf      chr1.snarl_ps_lv_at_edge.high_lv_edge.txt \\
        --out_notvcf   chr1.snarl_ps_lv_at_edge.not_in_vcf.high_lv_edge.txt
"""

import argparse
import gzip

import pandas as pd

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Output edges in each snarl that are not already claimed by a child snarl."
)
parser.add_argument(
    "--snarl_vcf", required=True,
    help="Step 1 output file (snarl_info_fromvcf.py)."
)
parser.add_argument(
    "--snarl_notvcf", required=True,
    help="Step 2 output file (find.not_in_vcf.AT_edge.py)."
)
parser.add_argument(
    "--out_vcf", required=True,
    help="Output file for --snarl_vcf with the new high-level-edge column appended."
)
parser.add_argument(
    "--out_notvcf", required=True,
    help="Output file for --snarl_notvcf with the new high-level-edge column appended."
)
args = parser.parse_args()


def output_edge_high_lv(snarl_file, output_file):
    """Append, to each row, the subset of its edges not seen in any deeper-level row."""
    opener = gzip.open if snarl_file.endswith(".gz") else open
    edge_set = set()  # edges already seen at a deeper level

    with opener(snarl_file, 'rt') as snarl, open(output_file, "w") as out:
        for line in snarl:
            line_split = line.strip().split("\t")
            try:
                if pd.isna(line_split[3]):
                    continue
            except IndexError:
                continue

            snarl_edge = set()
            for edge in line_split[4].split(","):
                if edge not in edge_set:
                    snarl_edge.add(edge)
                    edge_set.add(edge)

            out.write(f"{line.strip()}\t{','.join(snarl_edge)}\n")

# snarl from vcf 
output_edge_high_lv(args.snarl_vcf, args.out_vcf)
# snarl from gfa 
output_edge_high_lv(args.snarl_notvcf, args.out_notvcf)



