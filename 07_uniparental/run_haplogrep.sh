#!/usr/bin/env bash
# Mitochondrial DNA haplogroup inference using Haplogrep 3
# Paper: "The mitochondrial DNA haplogroup was inferred using Haplogrep 3"
# Result: haplogroup T1a2

SAMPLE=JD206
BAM=${SAMPLE}_libmerged_rg_rmdup.bam
REF_MT=/path/to/chrMT_rCRS.fa    # mitochondrial reference (rCRS, NC_012920)

# Step 1: Extract MT reads and call variants
samtools mpileup -uf $REF_MT $BAM | \
    bcftools call -mv -Ov -o ${SAMPLE}_mt.vcf

# Step 2: Classify haplogroup with Haplogrep 3
# Can also be run via the web interface: https://haplogrep.i-med.ac.at
haplogrep3 classify \
    --in ${SAMPLE}_mt.vcf \
    --format vcf \
    --out ${SAMPLE}_haplogrep.txt

# Result: T1a2
