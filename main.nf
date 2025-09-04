#!/usr/bin/env nextflow

params.reads = "$baseDir/data/*_{R1,R2}.fastq.gz"
params.genome = "$baseDir/reference/genome.fa"
params.gtf = "$baseDir/reference/genes.gtf"
params.outdir = "$baseDir/results"

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.html"
    path "*.zip"
    
    script:
    """
    fastqc -t 2 -q ${reads}
    """
}

process TRIMMOMATIC {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*_trimmed_{R1,R2}.fastq.gz")
    
    script:
    """
    trimmomatic PE -threads 4 \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \
        ${sample_id}_trimmed_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process BWA_INDEX {
    input:
    path genome
    
    output:
    path "${genome}*"
    
    script:
    """
    bwa index ${genome}
    """
}

process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path genome
    path index
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    bwa mem -t 4 ${genome} ${reads} | \
    samtools sort -@ 4 -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}

process VARIANT_CALLING {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path genome
    
    output:
    path "${sample_id}.vcf.gz"
    
    script:
    """
    bcftools mpileup -Ou -f ${genome} ${bam} | \
    bcftools call -mv -Oz -o ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """
}

workflow {
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
    
    // Quality control
    FASTQC(read_pairs_ch)
    
    // Preprocessing
    trimmed_ch = TRIMMOMATIC(read_pairs_ch)
    
    // Alignment
    genome_ch = Channel.fromPath(params.genome)
    index_ch = BWA_INDEX(genome_ch)
    aligned_ch = BWA_MEM(trimmed_ch, genome_ch, index_ch)
    
    // Variant calling
    VARIANT_CALLING(aligned_ch, genome_ch)
}
