"""
Hard genotype: step 2 - genotype every edge of one chromosome for a batch of
samples, using the step 1 structure pickle.

For each sample:
  1. Local depth: depth at each GRCh38 articulation node, and per component
     the per-haplotype depth lambda (median over nearby articulation nodes / 2).
  2. Per component (see genotype_model.genotype_component):
       - allele (multinomial) model over the component DAG -> dosage per edge
       - copy-number (Poisson) model for cyclic edges and for allele-called
         edges whose depth deviates from lambda * dosage
  3. Recovered edges (step 1 recovery groups): copy-number model with lambda
     taken around the group's GRCh38 anchor node.
  4. Final genotype per edge:
       - allele dosage if the edge has an allele call (also when it was flagged
         as deviant and got a copy-number call; uncalled allele sites give 0),
       - otherwise the copy number,
       - '.' if no lambda could be estimated.

The depth file is read once per batch (all batch samples at once).

Input:
    --chrom       : chromosome name
    --structure   : step 1 output ({chrom}.structure.pkl)
    --depthfile   : per-chromosome raw edge-depth matrix (Chrom, Edge, sample...),
                    plain or .gz, header line optional
    --samples     : all sample names, one per line, in the depth-matrix column
                    order (column 3 onward)
    --batch       : samples to genotype in this run, one per line

Output (per sample, in --outdir):
    {chrom}.{sample}.edge_genotype.tsv      : edge_id, genotype (sorted by edge_id;
                                              identical rows for every sample)
    {chrom}.{sample}.edge_genotype_info.tsv : one row per call -
        edges      : comma-separated edge ids (allele site) or one edge id (CNV)
        CALL_TYPE  : Allele | CNV
        FORMAT     : GT:AD:DP:GL:GQ:PL
        SAMPLE     : GT is the allele genotype (e.g. 0/1, '.' uncalled) or the
                     copy number; AD the edge depth(s); GL log10 likelihoods

Usage:
    python genotype_edges.py \\
        --chrom     chr21 \\
        --structure chr21.structure.pkl \\
        --depthfile chr21.depth_per_edge.txt.gz \\
        --samples   samples.txt \\
        --batch     batch_000.txt \\
        --outdir    .
"""

import argparse
import os
import pickle
import sys
from array import array

from graph_utils import open_maybe_gzip, canonical_edge, parse_path
from genotype_model import (
    genotype_component,
    copy_number_info,
    compute_stem_total_depth,
    local_lambda,
)


def read_list(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def load_batch_depths(depthfile, column_indices):
    """
    One pass over the depth matrix. Returns (canonical edge ids, [array of
    depths per sample in column_indices order]).
    """
    edge_ids = []
    columns = [array('d') for _ in column_indices]
    n_fields = None
    with open_maybe_gzip(depthfile) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            edge = fields[1]
            if edge[:1] not in ('>', '<'):   # header line
                continue
            if n_fields is None:
                n_fields = len(fields)
                if max(column_indices) >= n_fields:
                    sys.exit(f"ERROR: {depthfile} has {n_fields} columns; "
                             f"sample list implies at least {max(column_indices) + 1}")
            edge_ids.append(canonical_edge(*parse_path(edge)))
            for col, idx in zip(columns, column_indices):
                col.append(float(fields[idx]))
    return edge_ids, columns


def format_cnv_row(edge_id, info):
    gl = ",".join(f"{ll:.4f}" for ll in info["log10_likelihoods"])
    pl = ",".join(str(p) for p in info["pl"])
    sample = f"{info['cn']}:{info['depth']}:{int(round(info['depth']))}:{gl}:{int(round(info['gq']))}:{pl}"
    return f"{edge_id}\tCNV\tGT:AD:DP:GL:GQ:PL\t{sample}\n"


def genotype_one_sample(chrom, sample, structure, edge_depths, outdir):
    stem_total_depth = compute_stem_total_depth(
        structure["articulation_nodes"], structure["outgoing_paths"],
        structure["incoming_paths"], structure["end_node"], edge_depths,
    )

    def lam(start, end):
        return local_lambda(structure["stem_positions"], structure["path_edges"],
                            structure["b38_node_order"], start, end, stem_total_depth)

    allele_genotype = {}   # canonical edge -> dosage
    allele_info = {}       # tuple(raw edge keys) -> GenotypeCall
    cnv_info = {}          # canonical edge -> copy-number call
    all_edges = set()      # every edge reported for this chromosome

    # components
    for entry in structure["components"]:
        stem_depth = lam(entry["start_component_node"], entry["end_component_node"])
        edge_genotype, edgegroup_info, component_cnv = genotype_component(entry, stem_depth, edge_depths)
        allele_genotype.update(edge_genotype)
        for edge_group, call in edgegroup_info.items():
            allele_info[tuple(e[2] for e in edge_group)] = call
        cnv_info.update(component_cnv)
        all_edges.update(entry["cyclic_edges"])

    # recovered edges
    for anchor, edges in structure["recovery_groups"]:
        stem_depth = lam(anchor, anchor)
        all_edges.update(edges)
        if stem_depth is None:
            continue
        for edge in edges:
            cnv_info[edge] = copy_number_info(edge_depths.get(edge, 0), stem_depth)

    all_edges.update(allele_genotype)

    # genotype: allele dosage > copy number > '.'
    prefix = os.path.join(outdir, f"{chrom}.{sample}")
    with open(f"{prefix}.edge_genotype.tsv", "w") as out:
        out.write("edge_id\tgenotype\n")
        for edge in sorted(all_edges):
            if edge in allele_genotype:
                value = allele_genotype[edge]
            elif edge in cnv_info:
                value = cnv_info[edge]["cn"]
            else:
                value = "."
            out.write(f"{edge}\t{value}\n")

    with open(f"{prefix}.edge_genotype_info.tsv", "w") as out:
        out.write("edges\tCALL_TYPE\tFORMAT\tSAMPLE\n")
        for edge_ids, call in allele_info.items():
            ad = ",".join(str(edge_depths.get(canonical_edge(*parse_path(e)), 0)) for e in edge_ids)
            gl = ",".join(f"{ll:.4f}" for ll in call.log10_likelihoods)
            pl = ",".join(str(p) for p in call.pl)
            sample_field = (f"{call.as_vcf_string()}:{ad}:{int(round(call.depth))}:"
                            f"{gl}:{int(round(call.gq))}:{pl}")
            out.write(f"{','.join(edge_ids)}\tAllele\tGT:AD:DP:GL:GQ:PL\t{sample_field}\n")
        for edge, info in cnv_info.items():
            out.write(format_cnv_row(edge, info))

    print(f"{chrom} / {sample}: {len(all_edges)} edge(s), {len(allele_genotype)} allele-called, "
          f"{len(cnv_info)} copy-number call(s)")


def main():
    parser = argparse.ArgumentParser(description="Genotype every edge of one chromosome for a batch of samples.")
    parser.add_argument("--chrom", required=True, help="Chromosome name.")
    parser.add_argument("--structure", required=True, help="Step 1 structure pickle.")
    parser.add_argument("--depthfile", required=True, help="Per-chromosome raw edge-depth matrix.")
    parser.add_argument("--samples", required=True,
                        help="All sample names in depth-matrix column order (one per line).")
    parser.add_argument("--batch", required=True, help="Samples to genotype in this run (one per line).")
    parser.add_argument("--outdir", default=".", help="Output directory.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    all_samples = read_list(args.samples)
    batch = read_list(args.batch)
    sample_column = {s: i + 2 for i, s in enumerate(all_samples)}   # Chrom, Edge, samples...
    unknown = [s for s in batch if s not in sample_column]
    if unknown:
        sys.exit(f"ERROR: sample(s) not in {args.samples}: {', '.join(unknown)}")

    with open(args.structure, "rb") as fh:
        structure = pickle.load(fh)

    edge_ids, columns = load_batch_depths(args.depthfile, [sample_column[s] for s in batch])
    print(f"{args.chrom}: {len(edge_ids)} edge(s) loaded for {len(batch)} sample(s)")

    for sample, column in zip(batch, columns):
        edge_depths = dict(zip(edge_ids, column))
        genotype_one_sample(args.chrom, sample, structure, edge_depths, args.outdir)


if __name__ == "__main__":
    main()
