#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dev        = false
params.num_inputs = 2
params.min_length = 1000

// read produced from prior script - perform assemblies on raw mapped read and those after 95% ANI filtering
params.DIR_input_pysam  = "/home/tdangelo/epc/mmdb_analysis/bams-bowtie2-to-all_sags_concat/bams-pysam95-filtered/mapped_reads" 
params.DIR_input_raw    = "/home/tdangelo/epc/mmdb_analysis/bams-bowtie2-to-all_sags_concat/mapped_reads" 
params.DIR_output_pysam = "${params.DIR_input_pysam}/spades"
params.DIR_output_raw   = "${params.DIR_input_raw}/spades"

// Create output directories at startup
new File("${params.DIR_output_pysam}").mkdirs()
new File("${params.DIR_output_raw}").mkdirs()

workflow {

    // Channel from pysam95 filtered reads - tag with source
    CH_pysam = Channel
        .fromFilePairs(
            "${params.DIR_input_pysam}/*_mapped_R{1,2}.fastq.gz",
            size: 2,
            checkIfExists: true
        )
        .map { sample, reads ->
            tuple(sample, reads[0], reads[1], 'pysam', params.DIR_output_pysam)
        }

    // Channel from raw mapped reads - tag with source
    CH_raw = Channel
        .fromFilePairs(
            "${params.DIR_input_raw}/*_mapped_R{1,2}.fastq.gz",
            size: 2,
            checkIfExists: true
        )
        .map { sample, reads ->
            tuple(sample, reads[0], reads[1], 'raw', params.DIR_output_raw)
        }

    // Combine both channels
    CH_all = CH_pysam.mix(CH_raw)

    CH_run = params.dev ? CH_all.take(params.num_inputs) : CH_all

    SPADES(CH_run)
}

process SPADES {
    tag "${sample_id} (${source})"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'

    publishDir "${output_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2), val(source), val(output_dir)

    output:
    tuple val(sample_id), val(source), path("${sample_id}_${source}")

    script:
    """
    set -euo pipefail

    spades.py \
        -t ${task.cpus} \
        -1 ${reads_R1} \
        -2 ${reads_R2} \
        -k 21,33,55,77 \
        -o ${sample_id}_${source}

    mv ${sample_id}_${source}/contigs.fasta \
       ${sample_id}_${source}/${sample_id}_contigs.fasta

    seqkit seq -m ${params.min_length} \
        ${sample_id}_${source}/${sample_id}_contigs.fasta \
        > ${sample_id}_${source}/${sample_id}_contigs_filtered.fasta
    """
}
