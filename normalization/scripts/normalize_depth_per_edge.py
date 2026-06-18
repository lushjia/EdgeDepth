#!/usr/bin/env python3

# Usage: python normalize_depth_per_edge.py size_factors.txt depth_per_edge.txt depth_per_edge.normalzied.txt 

import sys 

size_factors = {}

size_factor = sys.argv[1] # "all.hprc-v2.0-mc-grch38.depth_per_edge.size_factors.txt"
input_raw_edge_depth = sys.argv[2] # "all.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.txt"
output_norm_edge_depth = sys.argv[3] # "all.hprc-v2.0-mc-grch38.depth_per_edge.added_chroms.normalized.txt"

with open(size_factor) as infile:
    for i, line in enumerate(infile):
        if i > 0:
            sample, size_factor = line.strip().split("\t")
            size_factor = float(size_factor)
            size_factors[sample] = size_factor

with open(input_raw_edge_depth) as infile:
    with open(output_norm_edge_depth, "w") as outfile:
        for i, line in enumerate(infile):
            if i == 0:
                samples = line.strip().split("\t")[2:]
                scaling_factors = [size_factors[sample] for sample in samples]
                outfile.write(line)
            else:
                cols = line.strip().split("\t")
                chrom = cols[0]
                edge_id = cols[1]
                outfile.write(f"{chrom}\t{edge_id}")
                for depth_per_edge, scaling_factor in zip(map(int, cols[2:]), scaling_factors):
                    normalized_depth_per_edge = depth_per_edge / scaling_factor
                    outfile.write(f"\t{normalized_depth_per_edge:.6f}")
                outfile.write("\n")
