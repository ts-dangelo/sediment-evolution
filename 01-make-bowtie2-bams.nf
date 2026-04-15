#!/usr/bin/env nextflow
nextflow.enable.dsl=2

DIR_input  = "/home/tdangelo/epc/mmdb_analysis/qc_reads/" ## reads from 00
DIR_output = "/home/tdangelo/epc/mmdb_analysis/bams-bowtie2-to-AM-917-D13"
BT2_INDEX  = "/home/tdangelo/epc/mmdb_analysis/bt2_index" ## bowtie2 index of SAG AM-917-D13

params.dev        = false
params.num_inputs = 2

workflow {

    CH_samples = Channel
    .fromFilePairs(
        "${DIR_input}/*_QC_{1,2}.fastq.gz",
        size: 2,
        checkIfExists: true
    )
    .map { sample, reads ->
        tuple(sample, reads[0], reads[1])
    }

    CH_run = params.dev ? CH_samples.take(params.num_inputs) : CH_samples

    CH_sam = run_bowtie2(CH_run)

}

process run_bowtie2 {

    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    containerOptions "--bind ${BT2_INDEX}"

    publishDir "${DIR_output}", mode: 'copy'

    input:
    tuple val(base), path(read1), path(read2)

    output:
    path "${base}.bam"

    script:
    """
    bowtie2 \
        -x ${BT2_INDEX}/AM-917-D13 \
        -1 ${read1} \
        -2 ${read2} \
    | samtools view -b -o ${base}.bam

    """
}
