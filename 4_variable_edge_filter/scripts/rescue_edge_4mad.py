"""
Rescue edges that failed the GenomeSTRiP GMM filter based on a 4-MAD criterion.

Edges filtered by the GMM step (reason: '1cluster' or '1weight') are re-examined
using the normalized depth distribution.  An edge is rescued if:
  1. It has >= --min_outlier_samples samples outside the median ± 4*MAD range
     (indicating genuine depth variation across samples), AND
  2. median + 4*MAD > --mad_threshold
     (indicating the upper tail of the distribution has meaningful signal).

Outputs:
  --info_out    : all GMM-failed edges annotated with 4MAD columns
                  (cols: chr, edge_id, depths..., gmm_filter_reason, has_outliers, above_threshold)
  --rescued_out : edge IDs of edges that pass both rescue criteria (one per line)

Usage:
    python rescue_edge_4mad.py \\
        --depth      chr21.normalized_edge_depth.txt \\ 
        --fail       chr21.notpass_filter_genomeSTRiP.txt \\ 
        --info_out   chr21.notpass_filter_genomeSTRiP.4mad.info.txt \\
        --rescued_out chr21.rescued_edges.txt
"""

import argparse

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Rescue GMM-filtered edges that show genuine depth variation by 4-MAD criterion."
)
parser.add_argument(
    "--depth", required=True,
    help="Normalized edge depth file (tab-separated, columns: chr, edge_id, sample1, ...)."
)
parser.add_argument(
    "--fail", required=True,
    help="Edges that failed the GMM filter (tab-separated: edge_id, reason). "
         "Reason is '1cluster' or '1weight'."
)
parser.add_argument(
    "--info_out", required=True,
    help="Output annotation file: all GMM-failed edges with 4MAD summary columns appended."
)
parser.add_argument(
    "--rescued_out", required=True,
    help="Output file listing rescued edge IDs (one per line)."
)
parser.add_argument(
    "--mad_threshold", type=float, default=4.71,
    help="Rescue threshold: edge is rescued only if median + 4*MAD exceeds this value "
         "(default: 4.71)."
)
parser.add_argument(
    "--min_outlier_samples", type=int, default=10,
    help="Minimum number of samples outside median ± 4*MAD required to rescue an edge "
         "(default: 10)."
)
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Load normalized depth matrix
# ---------------------------------------------------------------------------
depth_df = pd.read_csv(args.depth, sep="\t", comment='#', header=None)
n_samples = depth_df.shape[1] - 2  # subtract chr and edge_id columns

# ---------------------------------------------------------------------------
# Load GMM-failed edges (edge_id -> reason)
# ---------------------------------------------------------------------------
fail_dict = {}
with open(args.fail, 'r') as f:
    for line in f:
        # >1>2 1cluster (edge_id, reason_been_filtered)
        parts = line.strip().split("\t")
        fail_dict[parts[0]] = parts[1]  # edge_id -> '1cluster' or '1weight'


# ---------------------------------------------------------------------------
# Extract failed edges from depth matrix
# ---------------------------------------------------------------------------
failed_df = depth_df[depth_df[1].isin(fail_dict)].copy()
# add a column to represent filter reason 
failed_df['gmm_filter_reason'] = failed_df[1].map(fail_dict)


# ---------------------------------------------------------------------------
# Apply 4-MAD rescue criteria row by row
# ---------------------------------------------------------------------------
has_outliers   = [] # boolean list to store if the node have >=10 samples in +-4MAD
above_threshold = [] # bollean list to store if the node have mean+4MAD above threshold 

for i in range(failed_df.shape[0]):
    x = failed_df.iloc[i, 2:2 + n_samples].values
    median = np.median(x)
    mad    = stats.median_abs_deviation(x)
    lower, upper = median - 4 * mad, median + 4 * mad

    n_outliers = int(np.sum((x < lower) | (x > upper)))
    has_outliers.append(n_outliers >= args.min_outlier_samples)
    above_threshold.append((median + 4 * mad) > args.mad_threshold)

failed_df['has_outliers']    = has_outliers
failed_df['above_threshold'] = above_threshold


# ---------------------------------------------------------------------------
# Write annotation file (all failed edges)
# ---------------------------------------------------------------------------
failed_df.to_csv(args.info_out, sep="\t", index=False, header=False)

# ---------------------------------------------------------------------------
# Write rescued edge IDs (pass both criteria)
# ---------------------------------------------------------------------------
rescued = failed_df[failed_df['has_outliers'] & failed_df['above_threshold']][1]
with open(args.rescued_out, 'w') as out:
    for edge_id in rescued:
        out.write(f"{edge_id}\n")



