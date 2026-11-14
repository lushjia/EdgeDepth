"""
Hard genotype: step 3 - merge one chromosome's per-sample genotype files into
an edge-by-sample genotype matrix.

All per-sample files from step 2 list the same edges in the same (sorted)
order, so they are read line by line in parallel; the edge id of every row is
checked across samples.

Input:
    --chrom   : chromosome name
    --samples : sample names, one per line (output column order)
    --indir   : directory with {chrom}.{sample}.edge_genotype.tsv files

Output:
    --out : edge_id, sample1, ..., sampleN (tab-separated; '.' = no call)

Usage:
    python merge_genotype_matrix.py \\
        --chrom   chr21 \\
        --samples samples.txt \\
        --indir   . \\
        --out     chr21.edge_genotype_matrix.tsv
"""

import argparse
import os
import sys
from contextlib import ExitStack


def main():
    parser = argparse.ArgumentParser(description="Merge per-sample edge genotypes into one matrix.")
    parser.add_argument("--chrom", required=True, help="Chromosome name.")
    parser.add_argument("--samples", required=True, help="Sample names, one per line (column order).")
    parser.add_argument("--indir", default=".", help="Directory with per-sample genotype files.")
    parser.add_argument("--out", required=True, help="Output matrix.")
    args = parser.parse_args()

    with open(args.samples) as fh:
        samples = [line.strip() for line in fh if line.strip()]
    paths = [os.path.join(args.indir, f"{args.chrom}.{s}.edge_genotype.tsv") for s in samples]
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit(f"ERROR: {len(missing)} per-sample file(s) missing, e.g. {missing[0]}")

    n_rows = 0
    with ExitStack() as stack, open(args.out, "w") as out:
        handles = [stack.enter_context(open(p)) for p in paths]
        for fh in handles:
            fh.readline()   # header
        out.write("edge_id\t" + "\t".join(samples) + "\n")
        for lines in zip(*handles):
            edge_id = None
            values = []
            for sample, line in zip(samples, lines):
                eid, value = line.rstrip("\n").split("\t")
                if edge_id is None:
                    edge_id = eid
                elif eid != edge_id:
                    sys.exit(f"ERROR: row {n_rows + 1}: {sample} has edge {eid}, expected {edge_id}")
                values.append(value)
            out.write(edge_id + "\t" + "\t".join(values) + "\n")
            n_rows += 1
        if any(fh.readline() for fh in handles):
            sys.exit("ERROR: per-sample files have different numbers of rows")

    print(f"{args.chrom}: wrote {n_rows} edge(s) x {len(samples)} sample(s) -> {args.out}")


if __name__ == "__main__":
    main()
