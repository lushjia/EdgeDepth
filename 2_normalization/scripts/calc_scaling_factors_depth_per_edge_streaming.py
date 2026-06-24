#!/usr/bin/env python3
"""
Compute DESeq2-style size factors (median-of-ratios) for edge depth normalization.

Two-pass streaming algorithm — mathematically equivalent to the original script
but never holds the full N_edges x N_samples matrix in memory.

Pass 1 (rows): stream edges to compute per-edge log geometric mean.
               Memory: O(N_edges) float32 array for log means.

Pass 2 (columns, batched): for each batch of BATCH_SIZE samples, stream the file
               again and accumulate log ratios only for valid edges, then compute
               per-sample median.
               Memory: O(N_valid_edges x BATCH_SIZE) float32 per batch.

Total peak memory: O(N_edges + N_valid_edges x BATCH_SIZE) x 4 bytes.
Example: N=10M edges, BATCH_SIZE=20 → ~560 MB, vs ~16 GB for the full matrix.
"""
import sys
import math
import numpy as np

INPUT = sys.argv[1] # edge depth matrix: "all.hprc-v2.0-mc-grch38.depth_per_edge.chr1-22.txt" 
OUTPUT = sys.argv[2] # size factor per sample: "all.hprc-v2.0-mc-grch38.depth_per_edge.size_factors.txt"

# Tune this to control peak memory: larger = fewer file passes, more RAM.
# BATCH_SIZE=20 means ceil(N_samples/20) extra passes through the input file.
BATCH_SIZE = 40

# ---------------------------------------------------------------------------
# Pass 1 (row by row): compute per-edge log geometric mean
# ---------------------------------------------------------------------------
# Edges where any sample has count=0 are marked -inf and excluded, matching
# the original script's behavior (np.log(0) = -inf propagates through mean).
print("Pass 1: computing per-edge log geometric means...", file=sys.stderr)

samples = None
n_samples = None
log_means_list = []

with open(INPUT) as fh:
    for i, line in enumerate(fh):
        parts = line.rstrip("\n").split("\t")
        if i == 0:
            samples = parts[2:]  # first two columns are chr and edge_id
            n_samples = len(samples)
            continue

        counts = np.array(parts[2:], dtype=np.int32)
        if np.any(counts == 0):
            log_means_list.append(-math.inf)
        else:
            log_means_list.append(float(np.log(counts.astype(np.float32)).mean()))

log_means = np.array(log_means_list, dtype=np.float32)
del log_means_list

valid_mask = ~np.isinf(log_means)   # edges non-zero in every sample
n_total = len(log_means)
n_valid = int(valid_mask.sum())
valid_log_means = log_means[valid_mask]  # shape: (n_valid,)

print(f"  Total edges:  {n_total}", file=sys.stderr)
print(f"  Valid edges:  {n_valid} (non-zero in all samples)", file=sys.stderr)

if n_valid == 0:
    print("ERROR: no valid edges found.", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Pass 2 (row by row, column-batched): accumulate log ratios per sample batch,
# then compute per-sample median
# ---------------------------------------------------------------------------
# For each batch of BATCH_SIZE samples we do one additional pass through the
# file, accumulating only the n_valid rows that passed the filter in Pass 1.
size_factors = np.empty(n_samples, dtype=np.float64)

sample_batches = [
    list(range(i, min(i + BATCH_SIZE, n_samples)))
    for i in range(0, n_samples, BATCH_SIZE)
]
n_batches = len(sample_batches)

print(f"Pass 2: {n_batches} batch(es) of up to {BATCH_SIZE} samples each...", file=sys.stderr)

for batch_num, batch_cols in enumerate(sample_batches):
    b = len(batch_cols)
    print(f"  Batch {batch_num + 1}/{n_batches}: samples {batch_cols[0]}-{batch_cols[-1]}", file=sys.stderr)

    # log_ratios[valid_edge_idx, sample_in_batch] = log(count) - log_mean
    log_ratios = np.empty((n_valid, b), dtype=np.float32)

    with open(INPUT) as fh:
        valid_idx = 0
        for i, line in enumerate(fh):
            if i == 0:
                continue
            edge_idx = i - 1

            if not valid_mask[edge_idx]:
                continue

            parts = line.rstrip("\n").split("\t")
            # All counts for this edge are non-zero (guaranteed by valid_mask),
            # so np.log is safe.
            counts_batch = np.array([parts[2 + s] for s in batch_cols], dtype=np.float32)
            log_ratios[valid_idx] = np.log(counts_batch) - valid_log_means[valid_idx]
            valid_idx += 1

    # Median across valid edges for each sample in this batch (axis=0)
    medians = np.median(log_ratios, axis=0)   # shape: (b,)
    for j, s in enumerate(batch_cols):
        size_factors[s] = math.exp(medians[j])

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
with open(OUTPUT, "w") as out:
    out.write("Sample\tSize Factor\n")
    for s, sf in zip(samples, size_factors):
        out.write(f"{s}\t{sf:.6f}\n")

print(f"Done. Size factors written to {OUTPUT}", file=sys.stderr)
