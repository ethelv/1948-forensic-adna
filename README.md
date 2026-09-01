# Ancestry inference of a mid-20th-century femur

This repository holds the **custom analysis and plotting code** only. Steps performed
with standard tools at default settings are not included as scripts (see *Not in this
repo* below). Input data (BAM/FASTQ, the AADR reference panel, EIGENSTRAT intermediates,
the Kraken2 database) are not distributed; each script has a short editable path block.

## Contents

| folder | step | output |
|--------|------|--------|
| `02_coverage_sex/` | per-chromosome coverage & sex inference (`figure1_coverage.ipynb`) | **Figure 1** |
| `04_genotyping/`   | pileupCaller pseudo-haploid calls on the Human Origins panel | EIGENSTRAT genotypes |
| `05_pca/`          | SmartPCA projection onto the AADR panel + plot (`figure2_pca.ipynb`) | **Figure 2** |
| `06_admixture/`    | reference filtering + ADMIXTURE K=5-8 + plot (`plot_admixture.py`) | **Figure 3** (K6), **Figure S2** (K7,K8) |
| `08_ashkenazi_lowSNP_pca/` | project 544 Ashkenazi genomes onto the same PCA | supplementary figure |
| `09_unmapped_taxonomy/`    | Kraken2 taxonomy of the unmapped reads | supplementary result |

## Not in this repo (standard tools, default settings)

- **Read preprocessing & mapping:** AdapterRemoval v2 -> bwa mem (SE, hs37d5) -> dedup ->
  ATLAS soft-clip, run via nf-core/eager 2.4.0 + Sentieon.
- **DNA damage (Figure S1):** mapDamage2, single-stranded mode.
- **Uniparental markers & contamination:** Haplogrep3 (mtDNA -> T1a2), YLeaf
  (Y -> G-FT276712), contamMix (MAP contamination 0.00048).

## Reference data

- Reference genome: hs37d5 (GRCh37 + decoy).
- Reference panel: AADR **v54.1** Human Origins (`v54.1.p1_HO_public`), 597,573 SNPs;
  446 West-Eurasian / North-African genomes used for projection.
- Ashkenazi genomes: Lencz et al. 2018.
- Kraken2 DB: PlusPF (benlangmead.github.io/aws-indexes/k2).

## Software

samtools, bedtools, EIGENSOFT (mergeit, smartpca), pileupCaller (sequenceTools), PLINK,
ADMIXTURE, pong, Kraken2, seqtk; Python 3 (numpy, pandas, matplotlib, seaborn).
