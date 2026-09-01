#!/usr/bin/env bash
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
