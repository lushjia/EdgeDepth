#!/usr/bin/env nextflow

/*
 * Replace biallelic edges by allele balance, run after redundancy_filter.nf and
 * variable_edge_filter.nf.
 *
 * Steps (run chromosome by chromosome):
 *   1. snarl_info_fromvcf    – map snarls to allele traversals/edges from the VCF.
 *   2. snarl_info_fromgfa    – for snarls with no VCF traversal, recover one from GFA paths.
 *   3. snarl_highlv_edge     – for both step 1 and step 2 outputs, keep only the
 *                              edges first introduced at each snarl's level (not
 *                              already claimed by a deeper child snarl).
 *   4. compute_edge_nonredundant_AB
 *                            – for triangle/diamond snarls on the GRCh38 path with
 *                              exactly one kept (post bubble-puncture) non-reference
 *                              edge, compute that edge's allele balance against the
 *                              GRCh38 edge. Produces an AB matrix (matrix 1):
 *                              edge_id, AB_sample1, ..., AB_sampleN.
 *   5. replace_edge_with_AB  – take the final variable-edge normalized depth matrix
 *                              from variable_edge_filter.nf (matrix 2: chr, edge_id,
 *                              sample1, ..., sampleN) and, for any edge also present
 *                              in matrix 1, replace its depth row with the AB row.
 */

nextflow.enable.dsl=2

params.chroms        = (1..22).collect { "chr${it}" }

params.vcf_dir        = "data/vcf"                  // dir with {prefix}{chrom}{suffix} VCFs
params.vcf_prefix      = "hprc-v2.0-mc-grch38.raw."
params.vcf_suffix      = ".sorted.vcf.gz"
params.gfa_dir         = "data/gfa"                  // dir with one {chrom}.gfa per chromosome
params.snarls          = "data/hprc-v2.0-mc-grch38.snarls.json"   // genome-wide snarl JSON

params.edge_raw_dir    = "data/edge_raw"             // dir with {chrom}.depth_per_edge.txt.gz (raw depth, split by chrom)
params.kept_edges_dir  = "data/kept_edges"           // dir with {chrom}.independent_edge_list.txt (redundancy_filter.nf output)
params.norm_depth_dir  = "data/variable_edge_filter_results"  // dir with {chrom}.variable_edge_list.normalized_depth.txt (variable_edge_filter.nf output)

params.scripts_dir     = "${projectDir}/scripts"
params.outdir          = "results"


// ---------------------------------------------------------------------------
// Step 1: snarl -> allele traversal/edge mapping from the VCF
// ---------------------------------------------------------------------------
process snarl_info_fromvcf {
    tag "${chrom}"

    input:
    tuple val(chrom), path(vcf), path(gfa)
    path(snarls)
    path(script)

    output:
    tuple val(chrom), path("${chrom}.snarl_ps_lv_at_edge.txt"), emit: snarl_vcf

    script:
    """
    python3 ${script} \\
        --vcf    ${vcf} \\
        --snarls ${snarls} \\
        --gfa    ${gfa} \\
        --out    ${chrom}.snarl_ps_lv_at_edge.txt
    """
}


// ---------------------------------------------------------------------------
// Step 2: recover allele traversal/edge from GFA paths for snarls absent from the VCF
// ---------------------------------------------------------------------------
process snarl_info_fromgfa {
    tag "${chrom}"

    input:
    tuple val(chrom), path(snarl_vcf), path(gfa)
    path(script)

    output:
    tuple val(chrom), path("${chrom}.snarl_ps_lv_at_edge.not_in_vcf.txt"), emit: snarl_notvcf

    script:
    """
    python3 ${script} \\
        --snarl_info ${snarl_vcf} \\
        --gfa        ${gfa} \\
        --out        ${chrom}.snarl_ps_lv_at_edge.not_in_vcf.txt
    """
}


// ---------------------------------------------------------------------------
// Step 3: keep only edges first introduced at each snarl's level
// ---------------------------------------------------------------------------
process snarl_highlv_edge {
    tag "${chrom}"

    input:
    tuple val(chrom), path(snarl_vcf), path(snarl_notvcf)
    path(script)

    output:
    tuple val(chrom),
          path("${chrom}.snarl_ps_lv_at_edge.high_lv_edge.txt"),
          path("${chrom}.snarl_ps_lv_at_edge.not_in_vcf.high_lv_edge.txt"), emit: high_lv_edges

    script:
    """
    python3 ${script} \\
        --snarl_vcf    ${snarl_vcf} \\
        --snarl_notvcf ${snarl_notvcf} \\
        --out_vcf      ${chrom}.snarl_ps_lv_at_edge.high_lv_edge.txt \\
        --out_notvcf   ${chrom}.snarl_ps_lv_at_edge.not_in_vcf.high_lv_edge.txt
    """
}


// ---------------------------------------------------------------------------
// Step 4: compute allele balance for qualifying biallelic edges (matrix 1)
// ---------------------------------------------------------------------------
process compute_edge_nonredundant_AB {
    tag "${chrom}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(chrom), path(edge_raw), path(gfa), path(snarl_vcf_hl), path(snarl_notvcf_hl), path(kept_edges)
    path(script)

    output:
    tuple val(chrom), path("${chrom}.bp_allele_balance.txt"), emit: ab_matrix

    script:
    """
    python3 ${script} \\
        --edge_raw     ${edge_raw} \\
        --gfa          ${gfa} \\
        --snarl_vcf    ${snarl_vcf_hl} \\
        --snarl_notvcf ${snarl_notvcf_hl} \\
        --kept_edges   ${kept_edges} \\
        --out          ${chrom}.bp_allele_balance.txt
    """
}


// ---------------------------------------------------------------------------
// Step 5: replace matrix 2 (variable-edge normalized depth) rows with the
// allele-balance row from matrix 1, for any edge present in both
// ---------------------------------------------------------------------------
process replace_edge_with_AB {
    tag "${chrom}"
    publishDir params.outdir, mode: 'move'

    input:
    tuple val(chrom), path(ab_matrix), path(norm_depth)

    output:
    path "${chrom}.variable_edge_list.normalized_depth.AB_replaced.txt"

    script:
    """
    awk -F'\\t' 'BEGIN{OFS="\\t"}
        NR==FNR {
            rest=\$0
            sub(/^[^\\t]*\\t/, "", rest)
            ab[\$1]=rest
            next
        }
        {
            if (\$2 in ab) print \$1, \$2, ab[\$2]
            else print \$0
        }
    ' ${ab_matrix} ${norm_depth} > ${chrom}.variable_edge_list.normalized_depth.AB_replaced.txt
    """
}


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    def chroms_list = (params.chroms instanceof List)
        ? params.chroms
        : params.chroms.toString().split(',').collect { it.trim() }

    ch_chroms = Channel.fromList(chroms_list)

    // Step 1
    ch_step1_input = ch_chroms.map { chrom ->
        tuple(chrom,
              file("${params.vcf_dir}/${params.vcf_prefix}${chrom}${params.vcf_suffix}"),
              file("${params.gfa_dir}/${chrom}.gfa"))
    }

    snarl_info_fromvcf(
        ch_step1_input,
        file(params.snarls),
        file("${params.scripts_dir}/snarl_info_fromvcf.py")
    )

    // Step 2
    ch_step2_input = snarl_info_fromvcf.out.snarl_vcf
        .map { chrom, f -> tuple(chrom, f, file("${params.gfa_dir}/${chrom}.gfa")) }

    snarl_info_fromgfa(
        ch_step2_input,
        file("${params.scripts_dir}/snarl_info_fromgfa.py")
    )

    // Step 3
    ch_step3_input = snarl_info_fromvcf.out.snarl_vcf
        .join(snarl_info_fromgfa.out.snarl_notvcf)

    snarl_highlv_edge(
        ch_step3_input,
        file("${params.scripts_dir}/snarl_highlv_edge.py")
    )

    // Step 4
    ch_step4_input = ch_chroms
        .map { chrom ->
            tuple(chrom,
                  file("${params.edge_raw_dir}/${chrom}.depth_per_edge.txt.gz"),
                  file("${params.gfa_dir}/${chrom}.gfa"),
                  file("${params.kept_edges_dir}/${chrom}.independent_edge_list.txt"))
        }
        .join(snarl_highlv_edge.out.high_lv_edges)
        .map { chrom, edge_raw, gfa, kept_edges, snarl_vcf_hl, snarl_notvcf_hl ->
            tuple(chrom, edge_raw, gfa, snarl_vcf_hl, snarl_notvcf_hl, kept_edges)
        }

    compute_edge_nonredundant_AB(
        ch_step4_input,
        file("${params.scripts_dir}/compute_edge_nonredundant_AB.py")
    )

    // Step 5
    ch_replace_input = compute_edge_nonredundant_AB.out.ab_matrix
        .map { chrom, ab ->
            tuple(chrom, ab, file("${params.norm_depth_dir}/${chrom}.variable_edge_list.normalized_depth.txt"))
        }

    replace_edge_with_AB(ch_replace_input)
}
