#!/usr/bin/env bash
# Figure S1 – DNA damage patterns
# Generates Fragmisincorporation_plot.pdf and Length_plot.pdf
# Output is written to JD206_libmerged_rg_rmdup.mapDamage/

mapDamage \
    -i JD206_libmerged_rg_rmdup.bam \
    -r /path/to/hs37d5.fa \
    --single-stranded
