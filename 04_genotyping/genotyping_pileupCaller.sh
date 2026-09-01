#!/usr/bin/env bash
# Pseudo-haploid genotype calling on the Human Origins SNP panel.
# One random allele is drawn per SNP (randomHaploid), mimicking a single ancient read.
# Output: EIGENSTRAT {geno,snp,ind} for the target genome, used for PCA and ADMIXTURE.
set -euo pipefail

# ---- edit for your setup ----
BAM=JD206_filtered_filtered.bam
REF=hs37d5.fa                        # GRCh37 + decoy
SNP=v54.1.p1_HO_public.snp           # AADR v54.1 Human Origins panel (597,573 SNPs)
SAMPLE=JD206
OUT=JD206_filtered_filtered

# 1. pileup at the panel sites (MAPQ >= 30, BAQ off, base quality >= 30)
samtools mpileup -R -B -q30 -Q30 -f "$REF" "$BAM" > pileup.txt

# 2. draw one random allele per covered SNP -> EIGENSTRAT
pileupCaller --randomHaploid --sampleNames "$SAMPLE" -f "$SNP" -e "$OUT" < pileup.txt

echo "wrote ${OUT}.geno / ${OUT}.snp / ${OUT}.ind"
