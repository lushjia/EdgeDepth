"""
Filter edges by fitting a GenomeSTRiP Gaussian Mixture Model (GMM).

For each edge, the depth distribution across samples is fit to a GMM with
nclusters = ceil(max_depth / expected_depth) components.  An edge is filtered
(removed) if:
  - nclusters == 1  (max depth <= expected_depth, monomorphic signal), or
  - max(component weights) > weight_threshold  (one component dominates).

Edges passing the filter are written to --pass_out with columns:
    edge_id  nclusters  weights  means  stds
Edges failing the filter are written to --fail_out with columns:
    edge_id  reason  (reason: '1cluster' or '1weight')

Usage:
    python filter_edge_genomeSTRiP.split.py \\
        --input            chr1.chunk_gd \\
        --gmmscript        compute_gmm.R \\ # /gpfs/gibbs/pi/ycgh/lushjia/project/SV/AFGR/RNA/sample_430/edge_fix/monomorphic_check/genomestrip_gmm/compute_gmm.R edge1
        --pass_out         chr1.chunk_gd.pass.txt \\  # filter_genomeSTRiP.txt 
        --fail_out         chr1.chunk_gd.fail.txt \\ # filter_genomeSTRiP.filtered_id.txt
        --expected_depth   17 \\
        --min_outlier_samples 10 

"""

import argparse
import math
import os
import subprocess

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Filter edges using a GenomeSTRiP GMM on per-sample depth values."
)
parser.add_argument(
    "--input", required=True,
    help="Input file: tab-separated, columns chr, edge_id, sample_depth_1, ..., sample_depth_N."
)
parser.add_argument(
    "--gmmscript", required=True,
    help="Path to compute_gmm.R."
)
parser.add_argument(
    "--pass_out", required=True,
    help="Output file for edges passing the filter (edge_id, nclusters, weights, means, stds)."
)
parser.add_argument(
    "--fail_out", required=True,
    help="Output file for edges failing the filter (edge_id, reason)."
)
# parser.add_argument(
#     "--n_samples", type=int, default=430,
#     help="Number of samples (depth columns) to read per edge (default: 430)."
# )
parser.add_argument(
    "--expected_depth", type=float, required=True,
    help="Expected depth per haplotype copy, used to set nclusters and scale GMM "
         "output (default: 17)."
)
parser.add_argument(
    "--min_outlier_samples", type=int, default=10, 
    help="Filter edges where samples outside the dominant GMM component below this value "
         "(default: 10)."
)
args = parser.parse_args()

# Prefix for the per-edge temp file; derived from pass_out to be unique per job
prefix = os.path.basename(args.pass_out)


# ---------------------------------------------------------------------------
# Read input
# ---------------------------------------------------------------------------
vcf_file_df = pd.read_csv(args.input, sep="\t", comment='#', header=None)


# ---------------------------------------------------------------------------
# Fit GMM per edge and classify
# ---------------------------------------------------------------------------
filter_edge = []

with open(args.pass_out, 'w') as keep_f:
    # edge_id ncluster wieghts mean std 
    n_samples = vcf_file_df.shape[1] - 2 
    weight_threshold = 1 - args.min_outlier_samples / n_samples # reuire 10 samples with another genotype 
    for num in range(vcf_file_df.shape[0]):
        # edge depths in samples 
        edge_id = vcf_file_df.iloc[num].values[1]
        x = vcf_file_df.iloc[num].values[2:2 + n_samples]
        nclusters = math.ceil(np.max(x) / args.expected_depth)

        # max <= expected_depth, unable to fit — filter out
        if nclusters <= 1:
            filter_edge.append(f'{edge_id}\t1cluster')
            continue

        # if max > expected_depth, fit genomeSTRiP GMM
        # Write depths to temp file for R
        x_tmp = pd.DataFrame(x)
        x_tmp['expect'] = args.expected_depth
        x_tmp.to_csv(f"{prefix}.tmp", sep="\t", index=True, header=False)

        # Run GenomeSTRiP GMM
        result = subprocess.run(
            f'Rscript {args.gmmscript} edge1 {prefix}.tmp {nclusters}',
            shell=True, text=True, stdout=subprocess.PIPE
        )
        result_line = result.stdout.split("\n")
        weights  = [float(i) for i in result_line[1].split(" ")[:-1]]
        mean     = [float(i) * args.expected_depth for i in result_line[2].split(" ")[:-1]]
        variance = [float(i) * args.expected_depth for i in result_line[3].split(" ")[:-1]]
        std      = np.sqrt(variance)

        os.remove(f"{prefix}.tmp")

        # Dominant-component edges are effectively monomorphic — filter out
        if max(weights) > weight_threshold:
            filter_edge.append(f'{edge_id}\t1weight')
            continue

        join_weight = ";".join([str(i) for i in weights])
        join_mean   = ";".join([str(i) for i in mean])
        join_std    = ";".join([str(i) for i in std])
        keep_f.write(f'{edge_id}\t{nclusters}\t{join_weight}\t{join_mean}\t{join_std}\n')

# ---------------------------------------------------------------------------
# Write filtered edges
# ---------------------------------------------------------------------------
with open(args.fail_out, 'w') as filter_f:
    for line in filter_edge:
        filter_f.write(line + '\n')





