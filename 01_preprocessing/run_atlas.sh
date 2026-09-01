#!/usr/bin/env bash
# Soft-clipping of reads using ATLAS to minimize residual damage biases
# Paper: "we soft-clipped the reads to minimize residual damage biases using ATLAS"
# ATLAS version used: 2.0.0-rc.10

ATLAS=/path/to/atlas
BAM=JD206_libmerged_rg_rmdup.bam
REF=/path/to/hs37d5.fa
SAMPLE=JD206

# Step 1: Estimate errors / recalibration parameters
$ATLAS \
    bam=$BAM \
    fasta=$REF \
    out=${SAMPLE}_errors \
    poolRecal=all \
    task=estimateErrors

# Step 2: Soft-clip reads using estimated error model
$ATLAS \
    bam=$BAM \
    fasta=$REF \
    out=${SAMPLE}_filtered \
    task=softClip

# Step 3: BAM diagnostics (produces the mapped_reads_rg_* files)
$ATLAS \
    bam=${SAMPLE}_filtered_rg.bam \
    fasta=$REF \
    out=mapped_reads \
    task=BAMDiagnostics
