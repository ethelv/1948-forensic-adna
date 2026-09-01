#!/usr/bin/env bash
# Contamination estimation using contamMix
# Paper: "no evidence of contamination (using contamMix), with a MAP contamination
#         rate estimate of 0.00048 (2.5%-97.5% quantiles: 0.000077 to 0.0112)"

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
