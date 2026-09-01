#!/usr/bin/env bash
# Read mapping to human reference genome (hs37d5) using BWA mem
# Paper: "we mapped both merged and unmerged reads to the human reference genome (hs37d5)
#         using bwa mem in single-end mode with default parameters"

REF=/path/to/hs37d5.fa
SAMPLE=JD206
THREADS=8

# Map merged (collapsed) reads — single-end mode
bwa mem -t $THREADS $REF ${SAMPLE}.collapsed.gz > ${SAMPLE}.collapsed.sam
bwa mem -t $THREADS $REF ${SAMPLE}.collapsed.truncated.gz > ${SAMPLE}.collapsed.truncated.sam

# Map unmerged reads — single-end mode
bwa mem -t $THREADS $REF ${SAMPLE}.singleton.truncated.gz > ${SAMPLE}.singleton.sam

# Convert, sort, and merge
samtools view -bS ${SAMPLE}.collapsed.sam         | samtools sort -o ${SAMPLE}.collapsed.bam
samtools view -bS ${SAMPLE}.collapsed.truncated.sam | samtools sort -o ${SAMPLE}.collapsed.truncated.bam
samtools view -bS ${SAMPLE}.singleton.sam          | samtools sort -o ${SAMPLE}.singleton.bam

samtools merge -f ${SAMPLE}_merged.bam \
    ${SAMPLE}.collapsed.bam \
    ${SAMPLE}.collapsed.truncated.bam \
    ${SAMPLE}.singleton.bam

# Mark duplicates
samtools sort ${SAMPLE}_merged.bam -o ${SAMPLE}_merged.sorted.bam
samtools rmdup -s ${SAMPLE}_merged.sorted.bam ${SAMPLE}_libmerged_rg_rmdup.bam
samtools index ${SAMPLE}_libmerged_rg_rmdup.bam

# Summary stats
samtools stats --threads $THREADS --reference $REF \
    ${SAMPLE}_libmerged_rg_rmdup.bam > ${SAMPLE}.stats.txt
