#!/usr/bin/env nextflow

/*
 * align WGS to pangenome reference and count edge depth per sample
 */

nextflow.enable.dsl=2


// 1. Define parameters
params.cram_list = "data/samples_list.txt"
params.b38_ref = "data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.gbz = "data/hprc-v2.0-mc-grch38.gbz"
params.hapl = "data/hprc-v2.0-mc-grch38.hapl"
params.edges = "data/hprc-v2.0-mc-grch38.edges.txt" // ensure output edges in same order for each sample
params.zjoin = "data/zjoin" // path to zjoin executable, used to join output edge depths with all edges

params.outdir = "results"

process cram2fastq {
    tag "sample ${cram.baseName}" 

    input:
    tuple file(cram), val(insert_size), val(std)
    path b38_ref 

    output:
    tuple val(cram.baseName), path("${cram.baseName}_1.fq"), path("${cram.baseName}_2.fq"), val(insert_size), val(std), emit: fastq_data

    script:
    """
    samtools sort --reference ${b38_ref} -n -@${task.cpus} ${cram} | samtools fastq --reference ${b38_ref} - -1 ${cram.baseName}_1.fq -2 ${cram.baseName}_2.fq -@${task.cpus} > /dev/null
    """
}

process fastq_align {
    tag "sample ${samp}"

    input:
    tuple val(samp), path(fq1), path(fq2), val(insert_size), val(std), path(gbz), path(hapl)

    output:
    tuple val(samp), path("${samp}.gbz"), path("${samp}.hprc-v2.0-mc-grch38.gam"),  emit: align_result

    script:
    """
    TMPDIR="tmp"
    mkdir -p \$TMPDIR
    # k-mer counting
    echo ${fq1} >> \$TMPDIR/${samp}.fastq.lst
    echo ${fq2} >> \$TMPDIR/${samp}.fastq.lst
    kmc -k29 -m128 -okff -t${task.cpus} -hp @\${TMPDIR}/${samp}.fastq.lst \${TMPDIR}/${samp} \${TMPDIR}
    # align to pangenome graph
    vg haplotypes -v 2 -t ${task.cpus} --include-reference --diploid-sampling -i ${hapl} -k \${TMPDIR}/${samp}.kff -g ${samp}.gbz ${gbz}
    vg giraffe -p -t ${task.cpus} -Z ${samp}.gbz -f ${fq1} -f ${fq2} --fragment-mean ${insert_size} --fragment-stdev ${std} -o gam > ${samp}.hprc-v2.0-mc-grch38.gam
    """
}


process count_edge_depth {
    tag "sample ${sample}"

    input:
    tuple val(sample), path(gbz), path(gam), path(edges), path(zjoin)

    output:
    path "${gam.baseName}.edge_depth.txt"

    script:
    """
    echo -e "Edge\\t${sample}" > ${gam.baseName}.depth_per_edge.txt
    vg pack --threads ${task.cpus} --xg ${gbz} --gam ${gam} --as-edge-table | awk -F'\\t' '
    BEGIN { OFS="\\t" }
    !/^from.id/ {
        node1 = \$1
        direct1 = (\$2 == 0) ? ">" : "<"
        node2 = \$3
        direct2 = (\$4 == 0) ? ">" : "<"
        depth = \$5

        if (node1 < node2) {
            strand1 = direct1 node1
            strand2 = direct2 node2
        } else {
            rev_direct2 = (direct2 == ">") ? "<" : ">"
            rev_direct1 = (direct1 == ">") ? "<" : ">"
            strand1 = rev_direct2 node2
            strand2 = rev_direct1 node1
        }

        print strand1 strand2, depth
    }
    ' >> ${gam.baseName}.depth_per_edge.txt

    # Join output edge depths with all edges, ensureing output edges in same order for each sample
    ./${zjoin} -a ${edges} -b ${gam.baseName}.depth_per_edge.txt -1 2 -2 1 -r -e 0 | cut -f4 > ${gam.baseName}.edge_depth.txt
    """
}

workflow {
    // Read the sample list file, one line per item
    Channel
        .fromPath(params.cram_list)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(file(row.cram), row.insert_size, row.std) }
        .set { ch_input_files }

    cram2fastq(ch_input_files, file(params.b38_ref)) 

    // Attach constant reference inputs to each fastq tuple
    ch_align_input = cram2fastq.out.fastq_data.map { sample, fq1, fq2, insert_size, std ->
        tuple(sample, fq1, fq2, insert_size, std, file(params.gbz), file(params.hapl))
    }

    fastq_align(ch_align_input)

    ch_count_input = fastq_align.out.align_result.map { sample, gbz, gam ->
        tuple(sample, gbz, gam, file(params.edges), file(params.zjoin))
    }

    count_edge_depth(ch_count_input)

}

