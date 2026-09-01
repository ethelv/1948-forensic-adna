#!/usr/bin/env bash
# Y chromosome haplogroup inference using YLeaf
# Paper: "The Y chromosome branch was inferred by YLeaf, accepting single-read calls
#         only when at least 90% of the bases agreed"
# Result: G-FT276712 (YFull) / G-FT18229 (FTDNA)

SAMPLE=JD206
BAM=${SAMPLE}_libmerged_rg_rmdup.bam
REF=/path/to/hs37d5.fa
YLEAF_POS=/path/to/yleaf/Position_files/hg19/   # position files for hg19

python yleaf.py \
    -bam $BAM \
    -ref $REF \
    -pos $YLEAF_POS \
    -r 1 \
    -q 30 \
    -base_majority 0.90 \
    -out ${SAMPLE}_yleaf

# -r 1            : minimum number of reads per position
# -q 30           : minimum mapping quality
# -base_majority  : minimum fraction of reads supporting a call (90% as per paper)

# Result: G-FT276712 (YFull nomenclature), present mostly in Arabian Peninsula and Levant
