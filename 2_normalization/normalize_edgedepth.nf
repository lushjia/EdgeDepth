#!/usr/bin/env nextflow

/*
 * Normalize pangenome edge depth across samples using median-of-ratios size factors.
 *
 * Steps:
 *   1. combine_depths    – paste per-sample depth columns into a single matrix
 *   2. calc_size_factors – compute DESeq2-style size factors (streaming, low memory)
 *   3. normalize_depth   – divide each sample's depths by its size factor
 *
 * Expected input: per-sample *.hprc-v2.0-mc-grch38.edge_depth.txt files produced
 * by the alignment_count/align_count_edgedepth.nf pipeline.
 */

nextflow.enable.dsl=2

params.edges         = "data/hprc-v2.0-mc-grch38.edges.txt"  // reference edge list (chr + edge_id columns)
params.depth_dir     = "results"                             // dir containing per-sample *.edge_depth.txt files
params.scripts_dir   = "scripts"                           // path to Python normalization scripts
params.outdir        = "results"


// ---------------------------------------------------------------------------
// Step 1: combine per-sample depth files into a single matrix
// ---------------------------------------------------------------------------
// Uses bash process substitution to prepend the sample name as a column header
// to each depth file before pasting alongside the edge reference columns.
// Sample names are derived from the depth filenames themselves (not a sample list),
// so the file set actually present in depth_dir is authoritative.
process combine_depths {
    publishDir params.outdir, mode: 'copy'

    input:
    path(depth_files)   // all per-sample depth files, collected into the work dir
    path(edges)

    output:
    path "all_sample.hprc-v2.0-mc-grch38.edge_depth.txt", emit: merged

    script:
    """
    cmd="paste ${edges}"
    for f in ${depth_files}; do
        sample=\$(basename "\${f}" .hprc-v2.0-mc-grch38.edge_depth.txt)
        cmd="\${cmd} <(cat \"\${f}\")"
    done
    eval "\${cmd}" > all_sample.hprc-v2.0-mc-grch38.edge_depth.txt
    """
}


// ---------------------------------------------------------------------------
// Step 2: compute per-sample size factors (DESeq2 median-of-ratios, streaming)
// ---------------------------------------------------------------------------
// The Python script expects input named all.hprc-v2.0-mc-grch38.depth_per_edge.chr1-22.txt;
// a symlink bridges the pipeline filename to the script's hardcoded constant.
process calc_size_factors {
    publishDir params.outdir, mode: 'copy'

    input:
    path(merged)
    path(script)

    output:
    path "all_sample.size_factors.txt", emit: size_factors

    script:
    """
    python3 ${script} ${merged} all_sample.size_factors.txt
    """
}


// ---------------------------------------------------------------------------
// Step 3: normalize edge depth per sample
// ---------------------------------------------------------------------------
// Symlinks bridge pipeline filenames to the names the normalize script expects.
// The normalized output is renamed to the final output filename.
process normalize_depth {
    publishDir params.outdir, mode: 'copy'

    input:
    path(merged)
    path(size_factors)   // staged with name matching what normalize script opens
    path(script)

    output:
    path "all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt"

    script:
    """
    python3 ${script} ${size_factors} ${merged} all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt
    """
}


workflow {
    // Collect all per-sample edge depth files from the depth directory
    ch_depth_files = Channel
        .fromPath("${params.depth_dir}/*.hprc-v2.0-mc-grch38.edge_depth.txt")
        .filter { it.name != "all_sample.hprc-v2.0-mc-grch38.edge_depth.txt" }
        .collect()

    combine_depths(
        ch_depth_files,
        file(params.edges)
    )

    calc_size_factors(
        combine_depths.out.merged,
        file("${params.scripts_dir}/calc_scaling_factors_depth_per_edge_streaming.py")
    )

    normalize_depth(
        combine_depths.out.merged,
        calc_size_factors.out.size_factors,
        file("${params.scripts_dir}/normalize_depth_per_edge.py")
    )
}
