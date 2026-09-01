#!/usr/bin/env bash
set -euo pipefail

# ---- edit for your setup ----
BAM=JD206_libmerged_rg_rmdup.bam                 # full BAM (incl. unmapped)
DB=k2_pluspf_16gb                                # Kraken2 PlusPF DB (benlangmead.github.io/aws-indexes/k2)
N=1000000                                        # reads to subsample
THREADS=8
SEED=1948

# 1. extract unmapped reads and subsample to N (fixed seed for reproducibility)
samtools view -f 4 -b "$BAM" \
  | samtools fasta - \
  | seqtk sample -s"$SEED" - "$N" > unmapped_sample.fasta

# 2. classify with Kraken2
kraken2 --db "$DB" --threads "$THREADS" --use-names \
        --report kraken_report.txt --output kraken_out.txt \
        unmapped_sample.fasta

# 3. overall classified / unclassified fraction
awk '{c[$1]++} END{for (k in c) print k, c[k]}' <(cut -f1 kraken_out.txt)
