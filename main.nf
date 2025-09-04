#!/usr/bin/env nextflow

process FASTQC {
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "fastqc_${sample_id}_logs"
    
    script:
    """
    fastqc $reads -o fastqc_${sample_id}_logs
    """
}

process BWA_ALIGN {
    input:
    tuple val(sample_id), path(reads)
    path genome
    
    output:
    path "${sample_id}.bam"
    
    script:
    """
    bwa mem $genome $reads | samtools sort -o ${sample_id}.bam
    """
}

workflow {
    read_pairs = Channel.fromFilePairs(params.reads)
    FASTQC(read_pairs)
    BWA_ALIGN(read_pairs, params.genome)
}
