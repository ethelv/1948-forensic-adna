# Ancestry inference of a mid-20th-century femur

Analysis and plotting code for the paper. Input data (BAM/FASTQ, the AADR reference panel,
EIGENSTRAT intermediates, the Kraken2 database) are not distributed; each script has a short
editable path block.

## Contents

| folder | step | output |
|--------|------|--------|
| `preprocessing/` | AdapterRemoval v2, bwa mem (SE, hs37d5), ATLAS soft-clip | mapped, deduplicated BAM |
| `damage/`        | mapDamage2, single-stranded | **Figure S1** |
| `coverage_sex/`  | per-chromosome coverage & sex inference (`figure1_coverage.ipynb`) | **Figure 1** |
| `genotyping/`    | pileupCaller pseudo-haploid calls on the Human Origins panel | EIGENSTRAT genotypes |
| `pca/`           | SmartPCA projection onto the AADR panel + plot (`figure2_pca.ipynb`) | **Figure 2** |
| `admixture/`     | reference filtering + ADMIXTURE K=5-8 + plot (`plot_admixture.py`) | **Figure 3** (K6), **Figure S2** (K7,K8) |
| `contamination/` | contamMix on mitochondrial reads | MAP contamination 0.00048 |
| `ashkenazi_lowSNP_pca/` | project 544 Ashkenazi genomes onto the same PCA | supplementary figure |
| `unmapped_taxonomy/`    | Kraken2 taxonomy of the unmapped reads | supplementary result |

Read preprocessing and mapping were run via nf-core/eager 2.4.0 + Sentieon; the scripts in
`preprocessing/` document the corresponding commands.

## Reference data

- Reference genome: hs37d5 (GRCh37 + decoy).
- Reference panel: AADR v54.1 Human Origins (`v54.1.p1_HO_public`), 597,573 SNPs;
  446 West-Eurasian / North-African genomes used for projection.
- Kraken2 DB: PlusPF (benlangmead.github.io/aws-indexes/k2).

## References

Tool and data citations, as in the manuscript:

- AdapterRemoval v2 — Schubert M, et al. *BMC Research Notes* 9, 88 (2016).
- bwa — Li H, et al. *Bioinformatics* 25, 1754–1760 (2009).
- mapDamage2 — Jónsson H, et al. *Bioinformatics* 29, 1682–1684 (2013).
- ATLAS — Link V, et al. *bioRxiv* 105346 (2017). https://doi.org/10.1101/105346
- contamMix — Fu Q, et al. *Current Biology* 23, 553–559 (2013).
- pileupCaller (SequenceTools) — Schiffels S. *SequenceTools* (2024).
- SmartPCA (EIGENSOFT) — Patterson N, et al. *PLoS Genetics* 2, e190 (2006).
- ADMIXTURE — Alexander DH, et al. *Genome Research* 19, 1655–1664 (2009).
- AADR reference panel — Mallick S, et al. *Scientific Data* 11, 182 (2024).
- Ashkenazi reference genomes — Lencz T, et al. *Human Genetics* (2018). https://doi.org/10.1007/s00439-018-1886-z

Revision tools not cited in the manuscript: Kraken2, pong, seqtk.

## Software

samtools, bedtools, AdapterRemoval v2, bwa, ATLAS, mapDamage2, EIGENSOFT (mergeit, smartpca),
pileupCaller (SequenceTools), PLINK, ADMIXTURE, pong, contamMix, Kraken2, seqtk;
Python 3 (numpy, pandas, matplotlib, seaborn).
