#!/usr/bin/env nextflow

/*
 * Optional: discrete (hard) genotypes for every edge, run after redundancy_filter.nf.
 *
 * Steps:
 *   1. prep_component_structure – per chromosome, sample-independent: split the
 *                                 graph into biconnected components anchored on
 *                                 GRCh38, find cyclic edges, precompute each
 *                                 component's DAG traversal, and group edges the
 *                                 components cannot reach (recovery groups).
 *   2. genotype_batch           – per chromosome x sample batch: allele
 *                                 (multinomial) model over each component DAG,
 *                                 copy-number (Poisson) model for cyclic, deviant
 *                                 and recovered edges. Allele dosage is preferred
 *                                 when an edge has both calls.
 *   3. merge_genotype_matrix    – per chromosome: merge per-sample genotypes into
 *                                 an edge-by-sample genotype matrix.
 */

nextflow.enable.dsl=2

params.chroms       = (1..22).collect { "chr${it}" }

params.gfa_dir      = "data/gfa"        // dir with one {chrom}.gfa per chromosome
params.edge_raw_dir = "data/edge_raw"   // dir with {chrom}.depth_per_edge.txt.gz (3_redundancy_filter output)
params.depth        = "data/all_sample.hprc-v2.0-mc-grch38.edge_depth.txt"  // genome-wide raw depth matrix; only its header (sample names) is read
params.batch_size   = 10                // samples per genotyping job

params.scripts_dir  = "${projectDir}/scripts"
params.outdir       = "results"


// Sample names = columns 3.. of the genome-wide depth matrix header
// (the per-chromosome depth files from 3_redundancy_filter have no header)
def read_sample_names(path) {
    def f = file(path)
    def stream = f.newInputStream()
    if (f.name.endsWith('.gz')) stream = new java.util.zip.GZIPInputStream(stream)
    def header = new BufferedReader(new InputStreamReader(stream)).readLine()
    stream.close()
    def fields = header.split('\t') as List
    return fields.size() > 2 ? fields[2..-1] : []
}


// ---------------------------------------------------------------------------
// Step 1: per-chromosome graph structure (topology only, no sample depths)
// ---------------------------------------------------------------------------
process prep_component_structure {
    tag "${chrom}"

    input:
    tuple val(chrom), path(gfa), path(edge_raw)
    path(script)
    path(graph_utils)

    output:
    tuple val(chrom), path("${chrom}.structure.pkl"), path(edge_raw), emit: structure

    script:
    """
    PYTHONPATH=\$PWD python3 ${script} \\
        --chrom     ${chrom} \\
        --gfa       ${gfa} \\
        --depthfile ${edge_raw} \\
        --out       ${chrom}.structure.pkl
    """
}


// ---------------------------------------------------------------------------
// Step 2: genotype every edge for one batch of samples
// ---------------------------------------------------------------------------
process genotype_batch {
    tag "${chrom}:batch${batch_id}"
    publishDir "${params.outdir}/genotype_info/${chrom}", mode: 'copy', pattern: '*.edge_genotype_info.tsv'

    input:
    tuple val(chrom), path(structure), path(edge_raw), val(batch_id), val(batch_samples)
    path(samples)
    path(script)
    path(libs)

    output:
    tuple val(chrom), path("${chrom}.*.edge_genotype.tsv"), emit: genotypes
    path("${chrom}.*.edge_genotype_info.tsv"), emit: info

    script:
    """
    printf '%s\\n' ${batch_samples.join(' ')} > batch.txt

    PYTHONPATH=\$PWD python3 ${script} \\
        --chrom     ${chrom} \\
        --structure ${structure} \\
        --depthfile ${edge_raw} \\
        --samples   ${samples} \\
        --batch     batch.txt \\
        --outdir    .
    """
}


// ---------------------------------------------------------------------------
// Step 3: merge per-sample genotypes into an edge-by-sample matrix
// ---------------------------------------------------------------------------
process merge_genotype_matrix {
    tag "${chrom}"
    publishDir params.outdir, mode: 'move'

    input:
    tuple val(chrom), path(genotype_files)
    path(samples)
    path(script)

    output:
    path("${chrom}.edge_genotype_matrix.tsv")

    script:
    """
    python3 ${script} \\
        --chrom   ${chrom} \\
        --samples ${samples} \\
        --indir   . \\
        --out     ${chrom}.edge_genotype_matrix.tsv
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    def chroms_list = (params.chroms instanceof List)
        ? params.chroms
        : params.chroms.toString().split(',').collect { it.trim() }

    def samples = read_sample_names(params.depth)
    if (!samples) error "No sample names found in the header of ${params.depth}"
    def batches = samples.collate(params.batch_size as int)
    log.info "${samples.size()} samples in ${batches.size()} batch(es) of up to ${params.batch_size}"

    ch_chroms = Channel.fromList(chroms_list)

    // sample names in depth-matrix column order, one per line
    ch_samples = Channel.fromList(samples)
        .collectFile(name: 'samples.txt', newLine: true, sort: false)
        .first()

    ch_batches = Channel.fromList(
        (0..<batches.size()).collect { i -> tuple(String.format('%03d', i), batches[i]) }
    )

    // Step 1
    ch_step1_input = ch_chroms.map { chrom ->
        tuple(chrom,
              file("${params.gfa_dir}/${chrom}.gfa"),
              file("${params.edge_raw_dir}/${chrom}.depth_per_edge.txt.gz"))
    }

    prep_component_structure(
        ch_step1_input,
        file("${params.scripts_dir}/prep_component_structure.py"),
        file("${params.scripts_dir}/graph_utils.py")
    )

    // Step 2: every chromosome x every sample batch
    genotype_batch(
        prep_component_structure.out.structure.combine(ch_batches),
        ch_samples,
        file("${params.scripts_dir}/genotype_edges.py"),
        [file("${params.scripts_dir}/graph_utils.py"), file("${params.scripts_dir}/genotype_model.py")]
    )

    // Step 3: wait for all batches of a chromosome, then merge
    ch_merge_input = genotype_batch.out.genotypes
        .groupTuple(size: batches.size())
        .map { chrom, file_lists -> tuple(chrom, file_lists.flatten()) }

    merge_genotype_matrix(
        ch_merge_input,
        ch_samples,
        file("${params.scripts_dir}/merge_genotype_matrix.py")
    )
}
