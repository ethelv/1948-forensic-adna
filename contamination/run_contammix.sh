#!/usr/bin/env bash
SAMPLE=JD206
BAM=${SAMPLE}_libmerged_rg_rmdup.bam
REF_MT=/path/to/mt_reference.fa      # mitochondrial reference (rCRS)
MALN=/path/to/alignment.maln         # multiple alignment of mt sequences

# Extract mitochondrial reads
samtools view -b $BAM MT > ${SAMPLE}_mt.bam
samtools index ${SAMPLE}_mt.bam

# Run contamMix
# Requires: R with contamMix package installed
Rscript $(which contamMix) \
    --samFn ${SAMPLE}_mt.bam \
    --malnFn $MALN \
    --figure ${SAMPLE}.contamMix_fig.pdf \
    --estimOut ${SAMPLE}.contamMix_estimates.txt
