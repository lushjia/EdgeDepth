#!/usr/bin/env nextflow

/*
 * Redundancy filtering of pangenome edges via the bubble-puncture method.
 *
 * Steps (run chromosome by chromosome):
 *   0. split_depth_by_chrom / split_norm_depth_by_chrom
 *                              – split the genome-wide raw and normalized depth
 *                              matrices into per-chromosome gzipped files, so every
 *                              downstream step reads only its own chromosome's data.
 *   1. step1_subgraph        – build graph, identify bridge edges, decompose into
 *                              biconnected-component GML files (one job per chr).
 *   2. step2_filter          – apply bubble-puncture filter to each subgraph GML
 *                              (one job per GML file, scattered from Step 1).
 *                              Large subgraphs (>15 MB, non-residual) use the faster
 *                              dp50 script; all others use the standard script.
 *   3. collect_edges         – cat all per-subgraph independent edge lists plus the
 *                              Step 1 bridge edges into one file per chromosome.
 *   4. filter_normalized_depth – subset this chromosome's normalized edge depth
 *                              matrix down to only the independent edges kept.
 */

nextflow.enable.dsl=2

params.chroms      = (1..22).collect { "chr${it}" }   // chromosomes to process
params.edges       = "data/hprc-v2.0-mc-grch38.edges.txt"   // edge reference (chr + edge_id)
params.depth       = "data/all_sample.hprc-v2.0-mc-grch38.edge_depth.txt"  // genome-wide raw depth matrix, split by chromosome below
params.norm_depth  = "data/all_sample.hprc-v2.0-mc-grch38.edge_depth_norm.txt"  // genome-wide normalized depth matrix, split by chromosome below
params.gfa_dir     = "data/gfa"                        // dir with one {chrom}.gfa per chromosome
params.scripts_dir = "scripts"                         // dir containing the Python filter scripts
params.outdir      = "results"
params.dp50_threshold = 15 * 1024 * 1024  // file size (bytes) above which dp50 script is used


// ---------------------------------------------------------------------------
// Step 0: split genome-wide depth matrices into per-chromosome gzipped files
// ---------------------------------------------------------------------------
process split_depth_by_chrom {
    tag "${chrom}"

    input:
    val(chrom)
    path(depth)

    output:
    tuple val(chrom), path("${chrom}.depth_per_edge.txt.gz"), emit: chr_depth

    script:
    """
    cat ${depth} | awk -v chr=${chrom} '\$1==chr' | gzip > ${chrom}.depth_per_edge.txt.gz
    """
}

process split_norm_depth_by_chrom {
    tag "${chrom}"

    input:
    val(chrom)
    path(norm_depth)

    output:
    tuple val(chrom), path("${chrom}.depth_per_edge.normalized.txt.gz"), emit: chr_norm_depth

    script:
    """
    cat ${norm_depth} | awk -v chr=${chrom} '\$1==chr' | gzip > ${chrom}.depth_per_edge.normalized.txt.gz
    """
}


// ---------------------------------------------------------------------------
// Step 1: build graph, find bridge edges, decompose into biconnected subgraphs
// ---------------------------------------------------------------------------
process step1_subgraph {
    tag "${chrom}"
    publishDir "${params.outdir}/${chrom}", mode: 'copy', pattern: "*.original_bridge_edge.txt"

    input:
    tuple val(chrom), path(gfa)
    path(edges)
    path(script)

    output:
    tuple val(chrom), path("${chrom}.original_bridge_edge.txt"), emit: bridge_edges
    tuple val(chrom), path("${chrom}_subgraphs/*.graph.gml"),    emit: subgraphs

    script:
    """
    python3 ${script} \\
        --chr          ${chrom} \\
        --edges        ${edges} \\
        --gfa          ${gfa} \\
        --bridge_out   ${chrom}.original_bridge_edge.txt \\
        --subgraph_dir ${chrom}_subgraphs
    """
}


// ---------------------------------------------------------------------------
// Step 2: bubble-puncture filter for one subgraph GML
// Large non-residual subgraphs (>15 MB) use the dp50 approximation.
// ---------------------------------------------------------------------------
process step2_filter {
    tag "${chrom} ${gml.name}"

    input:
    tuple val(chrom), path(gml), path(gfa), path(depth)
    path(script_regular)
    path(script_dp50)

    output:
    tuple val(chrom), path("*.independent_edge_list.txt"), emit: edge_lists

    script:
    def prefix     = gml.baseName.tokenize('.')[0]   // e.g. subgraph_5
    def use_dp50   = gml.size() > params.dp50_threshold && !gml.name.contains('subgraph_remaining')
    def run_script = use_dp50 ? script_dp50 : script_regular
    """
    python3 ${run_script} \\
        --gml   ${gml} \\
        --depth ${depth} \\
        --gfa   ${gfa} \\
        --out   ${chrom}.${prefix}.independent_edge_list.txt
    """
}


// ---------------------------------------------------------------------------
// Step 3: collect bridge edges (Step 1) + all subgraph results (Step 2)
//         into one file per chromosome
// ---------------------------------------------------------------------------
process collect_edges {
    tag "${chrom}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(chrom), path(bridge_edges), path(edge_lists)

    output:
    tuple val(chrom), path("${chrom}.independent_edge_list.txt"), emit: collected

    script:
    """
    cat ${bridge_edges} ${edge_lists} > ${chrom}.independent_edge_list.txt
    """
}


// ---------------------------------------------------------------------------
// Step 4: subset the genome-wide normalized depth matrix down to the
//         independent edges kept for this chromosome
// ---------------------------------------------------------------------------
process filter_normalized_depth {
    tag "${chrom}"
    publishDir params.outdir, mode: 'move'

    input:
    tuple val(chrom), path(independent_edges), path(norm_depth)

    output:
    path "${chrom}.independent_edge_list.normalized_depth.txt"

    script:
    """
    zcat ${norm_depth} | awk -F'\\t' '
        NR==FNR { keep[\$1]=1; next }
        \$2 in keep { print }
    ' ${independent_edges} - > ${chrom}.independent_edge_list.normalized_depth.txt
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    ch_chroms = Channel.fromList(params.chroms)

    // Step 0: split genome-wide depth matrices into per-chromosome gzipped files,
    // one split job per chromosome, so every downstream step reads only its own data
    split_depth_by_chrom(ch_chroms, file(params.depth))
    split_norm_depth_by_chrom(ch_chroms, file(params.norm_depth))

    // One channel item per chromosome, paired with its GFA file
    ch_step1_input = ch_chroms
        .map { chrom -> tuple(chrom, file("${params.gfa_dir}/${chrom}.gfa")) }

    step1_subgraph(
        ch_step1_input,
        file(params.edges),
        file("${params.scripts_dir}/filter_redundancy_bubble_puncture.step1_subgraph.py")
    )

    // Scatter: one channel item per GML file, re-attach GFA and this chromosome's
    // pre-split depth file
    ch_step2_input = step1_subgraph.out.subgraphs
        .transpose()   // [chrom, gml_list] -> one [chrom, gml] per file
        .map { chrom, gml ->
            tuple(chrom, gml, file("${params.gfa_dir}/${chrom}.gfa"))
        }
        .join(split_depth_by_chrom.out.chr_depth)   // [chrom, gml, gfa, chr_depth_gz]

    step2_filter(
        ch_step2_input,
        file("${params.scripts_dir}/filter_redundancy_bubble_puncture.step2_filter.py"),
        file("${params.scripts_dir}/filter_redundancy_bubble_puncture.step2_filter.dp50.py")
    )

    // Gather all per-subgraph edge lists back by chromosome, then join with bridge edges
    ch_edge_lists = step2_filter.out.edge_lists.groupTuple()   // [chrom, [file1, file2, ...]]
    ch_bridge     = step1_subgraph.out.bridge_edges            // [chrom, bridge_file]

    collect_edges(
        ch_bridge.join(ch_edge_lists)   // [chrom, bridge_file, [file1, file2, ...]]
    )

    // Join with this chromosome's pre-split normalized depth file
    ch_filter_norm_input = collect_edges.out.collected
        .join(split_norm_depth_by_chrom.out.chr_norm_depth)   // [chrom, independent_edge_list, chr_norm_depth_gz]

    filter_normalized_depth(ch_filter_norm_input)
}
