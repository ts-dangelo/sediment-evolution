#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dev        = false
params.num_inputs = 2

params.DIR_input     = "/home/tdangelo/epc/mmdb_analysis/qc_reads/"
params.DIR_output    = "/home/tdangelo/epc/mmdb_analysis/bams-bowtie2-to-all_sags_concat"
params.GENOME_INDEX  = "/home/tdangelo/epc/mmdb_analysis/all-sags-input/all-misag-sags-concat"
params.sample_file   = "/home/tdangelo/epc/mmdb_analysis/bt2_95id_10x-AM-917-D13-mean_coverage_samples.txt"
params.filter_script = "/home/tdangelo/epc/code/filter_95_bam.py"

// Create output directories at startup
new File("${params.DIR_output}").mkdirs()
new File("${params.DIR_output}/mapped_reads").mkdirs()
new File("${params.DIR_output}/bams-pysam95-filtered").mkdirs()
new File("${params.DIR_output}/bams-pysam95-filtered/mapped_reads").mkdirs()

workflow {
    CH_sample_list = Channel
        .fromPath(params.sample_file)
        .splitText()
        .map { it.trim() }
        .filter { it }

    CH_all_samples = Channel
        .fromFilePairs(
            "${params.DIR_input}/*_QC_{1,2}.fastq.gz",
            size: 2,
            checkIfExists: true
        )
        .map { sample, reads ->
            tuple(sample, reads[0], reads[1])
        }

    CH_samples = CH_all_samples
        .join(CH_sample_list.map { s -> tuple(s, s) })
        .map { sample, r1, r2, _dummy -> tuple(sample, r1, r2) }

    CH_run = params.dev ? CH_samples.take(params.num_inputs) : CH_samples

    // Step 1: Bowtie2 alignment
    BOWTIE2_ALIGN(CH_run)

    // Step 2: Sort + index raw bams, extract proper pairs
    SAMTOOLS_SORT_INDEX(BOWTIE2_ALIGN.out)

    // Step 3: Convert mapped bams to paired fastq.gz
    BAM_TO_FASTQ(SAMTOOLS_SORT_INDEX.out.mapped_bam, "raw")

    // Step 4: Filter with pysam identity script
    PYSAM_FILTER(SAMTOOLS_SORT_INDEX.out.sorted_bam)

    // Step 5: Sort + index filtered bams, extract proper pairs
    SAMTOOLS_SORT_INDEX_FILTERED(PYSAM_FILTER.out)

    // Step 6: Convert filtered mapped bams to paired fastq.gz
    BAM_TO_FASTQ_FILTERED(SAMTOOLS_SORT_INDEX_FILTERED.out.mapped_bam, "filtered")
}

process BOWTIE2_ALIGN {
    tag "$sample_id"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    publishDir "${params.DIR_output}", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    set -euo pipefail
    bowtie2 \
        -x ${params.GENOME_INDEX} \
        -1 ${reads_R1} \
        -2 ${reads_R2} \
        --no-unal \
        --threads ${task.cpus} \
        | samtools view -bS - > ${sample_id}.bam
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    publishDir "${params.DIR_output}",              mode: 'copy', pattern: "*.sorted.bam*"
    publishDir "${params.DIR_output}/mapped_reads", mode: 'copy', pattern: "*_mapped.bam*"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"),     emit: sorted_bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: sorted_bai
    tuple val(sample_id), path("${sample_id}_mapped.bam"),     emit: mapped_bam
    tuple val(sample_id), path("${sample_id}_mapped.bam.bai"), emit: mapped_bai

    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    samtools view -@ ${task.cpus} -b -f 0x2 -F 0x4 -F 0x8 \
        ${sample_id}.sorted.bam > ${sample_id}_mapped.bam
    samtools index ${sample_id}_mapped.bam
    """
}

process BAM_TO_FASTQ {
    tag "$sample_id"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    publishDir "${params.DIR_output}/mapped_reads", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample_id), path(mapped_bam)
    val suffix

    output:
    tuple val(sample_id),
          path("${sample_id}_mapped_R1.fastq.gz"),
          path("${sample_id}_mapped_R2.fastq.gz")

    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -n ${mapped_bam} -o ${sample_id}_namesorted.bam
    samtools fastq -@ ${task.cpus} \
        -1 ${sample_id}_mapped_R1.fastq.gz \
        -2 ${sample_id}_mapped_R2.fastq.gz \
        -0 /dev/null \
        -s /dev/null \
        -n \
        ${sample_id}_namesorted.bam
    rm ${sample_id}_namesorted.bam
    """
}

process PYSAM_FILTER {
    tag "$sample_id"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    publishDir "${params.DIR_output}/bams-pysam95-filtered", mode: 'copy', pattern: "*_pysam95.bam"

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}_pysam95.bam")

    script:
    """
    set -euo pipefail
    python3 ${params.filter_script} \
        ${sorted_bam} \
        ${sample_id}_pysam95.bam \
        95
    """
}

process SAMTOOLS_SORT_INDEX_FILTERED {
    tag "$sample_id"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    publishDir "${params.DIR_output}/bams-pysam95-filtered",              mode: 'copy', pattern: "*.sorted.bam*"
    publishDir "${params.DIR_output}/bams-pysam95-filtered/mapped_reads", mode: 'copy', pattern: "*_mapped.bam*"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_pysam95.sorted.bam"),     emit: sorted_bam
    tuple val(sample_id), path("${sample_id}_pysam95.sorted.bam.bai"), emit: sorted_bai
    tuple val(sample_id), path("${sample_id}_pysam95_mapped.bam"),     emit: mapped_bam
    tuple val(sample_id), path("${sample_id}_pysam95_mapped.bam.bai"), emit: mapped_bai

    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -o ${sample_id}_pysam95.sorted.bam ${bam}
    samtools index ${sample_id}_pysam95.sorted.bam
    samtools view -@ ${task.cpus} -b -f 0x2 -F 0x4 -F 0x8 \
        ${sample_id}_pysam95.sorted.bam > ${sample_id}_pysam95_mapped.bam
    samtools index ${sample_id}_pysam95_mapped.bam
    """
}

process BAM_TO_FASTQ_FILTERED {
    tag "$sample_id"
    label 'charlie'
    container '/home/tdangelo/singularity_stuff/bowtie2.sif'
    publishDir "${params.DIR_output}/bams-pysam95-filtered/mapped_reads", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample_id), path(mapped_bam)
    val suffix

    output:
    tuple val(sample_id),
          path("${sample_id}_pysam95_mapped_R1.fastq.gz"),
          path("${sample_id}_pysam95_mapped_R2.fastq.gz")

    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -n ${mapped_bam} -o ${sample_id}_namesorted.bam
    samtools fastq -@ ${task.cpus} \
        -1 ${sample_id}_pysam95_mapped_R1.fastq.gz \
        -2 ${sample_id}_pysam95_mapped_R2.fastq.gz \
        -0 /dev/null \
        -s /dev/null \
        -n \
        ${sample_id}_namesorted.bam
    rm ${sample_id}_namesorted.bam
    """
}
