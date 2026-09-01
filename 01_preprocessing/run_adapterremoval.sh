#!/usr/bin/env bash
# Read trimming and merging using AdapterRemoval v2
# Paper: "We trimmed and merged reads using AdapterRemoval v2"
# Input: paired-end FASTQ files (76bp reads)
# Output: merged (collapsed) reads + unmerged reads, used as input to BWA

AdapterRemoval \
    --file1 JD206_R1.fastq.gz \
    --file2 JD206_R2.fastq.gz \
    --basename JD206 \
    --collapse \
    --trimns \
    --trimqualities \
    --minquality 20 \
    --minlength 25 \
    --threads 8

# Output files used downstream:
#   JD206.collapsed.gz           → merged paired reads
#   JD206.collapsed.truncated.gz → merged + quality-trimmed
#   JD206.singleton.truncated.gz → unmerged reads
