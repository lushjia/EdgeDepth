#!/usr/bin/env nextflow

/*
 * Redundancy filtering of pangenome edges via the bubble-puncture method.
 *
 * Steps (run chromosome by chromosome):
 *   1. step1_subgraph  – build graph, identify bridge edges, decompose into
 *                        biconnected-component GML files (one job per chr).
 *   2. step2_filter    – apply bubble-puncture filter to each subgraph GML
 *                        (one job per GML file, scattered from Step 1).
 *                        Large subgraphs (>15 MB, non-residual) use the faster
 *                        dp50 script; all others use the standard script.
 *   3. collect_edges   – cat all per-subgraph independent edge lists plus the
 *                        Step 1 bridge edges into one file per chromosome.
 */

nextflow.enable.dsl=2

params.chroms      = (1..22).collect { "chr${it}" }   // chromosomes to process
params.edges       = "data/hprc-v2.0-mc-grch38.edges.txt"   // edge reference (chr + edge_id)
params.depth       = "data/all_sample.hprc-v2.0-mc-grch38.edge_depth.txt"  // raw depth matrix
params.gfa_dir     = "data/gfa"                        // dir with one {chrom}.gfa per chromosome
params.scripts_dir = "scripts"                         // dir containing the Python filter scripts
params.outdir      = "results"

// file size threshold (bytes) above which dp50 script is used
def DP50_THRESHOLD = 15 * 1024 * 1024  // 15 MB


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
    tuple val(chrom), path(gml), path(gfa)
    path(depth)
    path(script_regular)
    path(script_dp50)

    output:
    tuple val(chrom), path("*.independent_edge_list.txt"), emit: edge_lists

    script:
    def prefix     = gml.baseName.tokenize('.')[0]   // e.g. subgraph_5
    def use_dp50   = gml.size() > DP50_THRESHOLD && !gml.name.contains('subgraph_remaining')
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
    path "${chrom}.independent_edge_list.txt"

    script:
    """
    cat ${bridge_edges} ${edge_lists} > ${chrom}.independent_edge_list.txt
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    // One channel item per chromosome, paired with its GFA file
    ch_step1_input = Channel
        .fromList(params.chroms)
        .map { chrom -> tuple(chrom, file("${params.gfa_dir}/${chrom}.gfa")) }

    step1_subgraph(
        ch_step1_input,
        file(params.edges),
        file("${params.scripts_dir}/filter_redundancy_bubble_puncture.step1_subgraph.py")
    )

    // Scatter: one channel item per GML file, re-attach GFA for that chromosome
    ch_step2_input = step1_subgraph.out.subgraphs
        .transpose()   // [chrom, gml_list] -> one [chrom, gml] per file
        .map { chrom, gml ->
            tuple(chrom, gml, file("${params.gfa_dir}/${chrom}.gfa"))
        }

    step2_filter(
        ch_step2_input,
        file(params.depth),
        file("${params.scripts_dir}/filter_redundancy_bubble_puncture.step2_filter.py"),
        file("${params.scripts_dir}/filter_redundancy_bubble_puncture.step2_filter.dp50.py")
    )

    // Gather all per-subgraph edge lists back by chromosome, then join with bridge edges
    ch_edge_lists = step2_filter.out.edge_lists.groupTuple()   // [chrom, [file1, file2, ...]]
    ch_bridge     = step1_subgraph.out.bridge_edges            // [chrom, bridge_file]

    collect_edges(
        ch_bridge.join(ch_edge_lists)   // [chrom, bridge_file, [file1, file2, ...]]
    )
}
