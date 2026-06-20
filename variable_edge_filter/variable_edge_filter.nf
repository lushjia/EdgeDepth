#!/usr/bin/env nextflow

/*
 * Variable edge filtering, run after redundancy_filter.nf.
 *
 * Input: per-chromosome normalized depth files for independent edges
 *        (output of redundancy_filter.nf's filter_normalized_depth step):
 *          {chrom}.independent_edge_list.normalized_depth.txt
 *
 * Steps (run chromosome by chromosome):
 *   1. filter0       – keep edges with depth >= min_depth in >= min_samples samples.
 *   2. filter1000     – drop edges whose mean depth across samples exceeds max_mean_depth
 *                       (likely mapping/collapse artifacts).
 *   3. split_chunks / gmm_filter / merge_gmm
 *                     – split each chromosome's edges into chunks (chunk_size lines),
 *                       fit a GenomeSTRiP GMM to each edge in parallel, then merge all
 *                       per-chunk pass/fail results back into one file per chromosome.
 *   4. rescue_4mad / collect_final
 *                     – re-examine GMM-failed edges with a 4-MAD criterion to rescue
 *                       genuinely variable edges, then combine GMM-pass + rescued edges
 *                       into the final kept edge list per chromosome.
 *   5. filter_final_depth
 *                     – subset this chromosome's normalized depth matrix down to
 *                       only the edges in the final kept edge list.
 */

nextflow.enable.dsl=2

params.norm_depth_dir = "data"     // dir with {chrom}.independent_edge_list.normalized_depth.txt
params.gmmscript      = "compute_gmm.R"   // path to GenomeSTRiP GMM-fitting R script
params.scripts_dir     = "scripts"  // dir containing the Python filter scripts
params.outdir          = "results"

// filter0: minimum depth and minimum sample count
params.min_depth   = 5
params.min_samples = 10

// filter1000: maximum allowed mean depth across samples
params.max_mean_depth = 1000

// GMM step
params.chunk_size              = 2000
params.expected_depth          = 17
params.gmm_min_outlier_samples = 10    // samples required outside the dominant GMM component

// 4-MAD rescue step
params.mad_threshold           = 4.71
params.mad_min_outlier_samples = 10    // samples required outside median +/- 4*MAD


// ---------------------------------------------------------------------------
// Step 1: keep edges with >= min_depth reads in >= min_samples samples
// ---------------------------------------------------------------------------
process filter0 {
    tag "${chrom}"

    input:
    tuple val(chrom), path(norm_depth)

    output:
    tuple val(chrom), path("${chrom}.independent.filter0.txt"), emit: filtered

    script:
    """
    awk '{s=0; for(i=3;i<=NF;i++) {if(\$i>=${params.min_depth}) {s+=1; if(s>=${params.min_samples}) {print \$0; next}}}}' \\
        ${norm_depth} > ${chrom}.independent.filter0.txt
    """
}


// ---------------------------------------------------------------------------
// Step 2: drop edges with abnormally high mean depth
// ---------------------------------------------------------------------------
process filter1000 {
    tag "${chrom}"

    input:
    tuple val(chrom), path(filter0_file)

    output:
    tuple val(chrom), path("${chrom}.independent.filter0.filter1000.txt"), emit: filtered

    script:
    """
    awk '{s=0; for(i=3;i<=NF;i++) s+=\$i; s /= (NF-2); if(s<=${params.max_mean_depth}){print \$0}}' \\
        ${filter0_file} > ${chrom}.independent.filter0.filter1000.txt
    """
}


// ---------------------------------------------------------------------------
// Step 3a: split each chromosome's edges into chunks for parallel GMM fitting
// ---------------------------------------------------------------------------
process split_chunks {
    tag "${chrom}"

    input:
    tuple val(chrom), path(filter1000_file)

    output:
    tuple val(chrom), path("${chrom}.chunk_*"), emit: chunks

    script:
    """
    split -l ${params.chunk_size} ${filter1000_file} ${chrom}.chunk_
    """
}


// ---------------------------------------------------------------------------
// Step 3b: fit GenomeSTRiP GMM to each edge in one chunk
// ---------------------------------------------------------------------------
process gmm_filter {
    tag "${chrom} ${chunk.name}"

    input:
    tuple val(chrom), path(chunk)
    path(gmmscript)
    path(script)

    output:
    tuple val(chrom), path("${chunk.name}.pass.txt"), emit: pass
    tuple val(chrom), path("${chunk.name}.fail.txt"), emit: fail

    script:
    """
    python3 ${script} \\
        --input              ${chunk} \\
        --gmmscript          ${gmmscript} \\
        --pass_out           ${chunk.name}.pass.txt \\
        --fail_out           ${chunk.name}.fail.txt \\
        --expected_depth     ${params.expected_depth} \\
        --min_outlier_samples ${params.gmm_min_outlier_samples}
    """
}


// ---------------------------------------------------------------------------
// Step 3c: merge all per-chunk pass/fail results into one file per chromosome
// ---------------------------------------------------------------------------
process merge_gmm {
    tag "${chrom}"

    input:
    tuple val(chrom), path(pass_files), path(fail_files)

    output:
    tuple val(chrom), path("${chrom}.filter_genomeSTRiP.pass.txt"), emit: pass
    tuple val(chrom), path("${chrom}.filter_genomeSTRiP.fail.txt"), emit: fail

    script:
    """
    cat ${pass_files} > ${chrom}.filter_genomeSTRiP.pass.txt
    cat ${fail_files} > ${chrom}.filter_genomeSTRiP.fail.txt
    """
}


// ---------------------------------------------------------------------------
// Step 4a: rescue GMM-failed edges that show genuine depth variation (4-MAD)
// ---------------------------------------------------------------------------
process rescue_4mad {
    tag "${chrom}"

    input:
    tuple val(chrom), path(filter1000_file), path(fail_file)
    path(script)

    output:
    tuple val(chrom), path("${chrom}.4mad.info.txt"),     emit: info
    tuple val(chrom), path("${chrom}.rescued_edges.txt"), emit: rescued

    script:
    """
    python3 ${script} \\
        --depth               ${filter1000_file} \\
        --fail                ${fail_file} \\
        --info_out            ${chrom}.4mad.info.txt \\
        --rescued_out         ${chrom}.rescued_edges.txt \\
        --mad_threshold       ${params.mad_threshold} \\
        --min_outlier_samples ${params.mad_min_outlier_samples}
    """
}


// ---------------------------------------------------------------------------
// Step 4b: combine GMM-pass edges + rescued edges into the final kept edge list
// ---------------------------------------------------------------------------
process collect_final {
    tag "${chrom}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(chrom), path(pass_file), path(rescued_file)

    output:
    tuple val(chrom), path("${chrom}.variable_edge_list.final.txt"), emit: final_edges

    script:
    """
    cut -f1 ${pass_file} > ${chrom}.variable_edge_list.final.txt
    cat ${rescued_file} >> ${chrom}.variable_edge_list.final.txt
    """
}


// ---------------------------------------------------------------------------
// Step 5: subset this chromosome's normalized depth matrix down to the
//         final kept edges
// ---------------------------------------------------------------------------
process filter_final_depth {
    tag "${chrom}"
    publishDir params.outdir, mode: 'move'

    input:
    tuple val(chrom), path(final_edges), path(norm_depth)

    output:
    path "${chrom}.variable_edge_list.normalized_depth.txt"

    script:
    """
    awk -F'\\t' '
        NR==FNR { keep[\$1]=1; next }
        \$2 in keep { print }
    ' ${final_edges} ${norm_depth} > ${chrom}.variable_edge_list.normalized_depth.txt
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    // One channel item per chromosome, derived from redundancy_filter.nf output
    ch_norm_depth = Channel
        .fromPath("${params.norm_depth_dir}/*.independent_edge_list.normalized_depth.txt")
        .map { f -> tuple(f.name.tokenize('.')[0], f) }

    filter0(ch_norm_depth)
    filter1000(filter0.out.filtered)
    split_chunks(filter1000.out.filtered)

    // Scatter: one channel item per chunk file
    ch_chunks = split_chunks.out.chunks.transpose()   // [chrom, chunk_file]

    gmm_filter(
        ch_chunks,
        file(params.gmmscript),
        file("${params.scripts_dir}/filter_edge_genomeSTRiP.split.py")
    )

    // Gather all chunk results back by chromosome
    ch_pass_grouped = gmm_filter.out.pass.groupTuple()   // [chrom, [pass_file, ...]]
    ch_fail_grouped = gmm_filter.out.fail.groupTuple()   // [chrom, [fail_file, ...]]

    merge_gmm(ch_pass_grouped.join(ch_fail_grouped))

    // Rescue GMM-failed edges using the pre-GMM depth matrix for this chromosome
    ch_rescue_input = filter1000.out.filtered.join(merge_gmm.out.fail)

    rescue_4mad(
        ch_rescue_input,
        file("${params.scripts_dir}/rescue_edge_4mad.py")
    )

    collect_final(
        merge_gmm.out.pass.join(rescue_4mad.out.rescued)
    )

    filter_final_depth(
        collect_final.out.final_edges.join(ch_norm_depth)
    )
}
