# Ancestry inference of a mid-20th-century femur (sample 1948 / JD206)

Low-coverage whole-genome sequencing of a single skeletal sample was mapped to the
human reference, pseudo-haploid genotyped on the Human Origins SNP panel, and its
ancestry inferred by PCA and ADMIXTURE against the AADR reference panel, plus
uniparental markers (mtDNA, Y) and contamination estimation.

Each numbered folder is one analysis step. Scripts contain a short editable path block
at the top; **input data are not included** (BAM/FASTQ, the AADR v54.1 Human Origins
panel, the Kraken2 database, and intermediate EIGENSTRAT files).

## Pipeline

| folder | step | output |
|--------|------|--------|
| `01_preprocessing/` | AdapterRemoval -> bwa mem (SE, hs37d5) -> ATLAS soft-clip | mapped, deduplicated BAM |
| `02_coverage_sex/`  | per-chromosome coverage & sex inference | **Figure 1** |
| `03_damage/`        | mapDamage post-mortem damage profile | **Figure S1** |
| `04_genotyping/`    | pileupCaller pseudo-haploid calls on the HO panel | EIGENSTRAT genotypes |
| `05_pca/`           | SmartPCA projection onto the AADR panel + plot | **Figure 2** |
| `06_admixture/`     | ADMIXTURE K=5-8 + plot | **Figure 3** (K6), **Figure S2** (K7,K8) |
| `07_uniparental/`   | Haplogrep3 (mtDNA), YLeaf (Y), contamMix (contamination) | haplogroups, contamination |
| `08_ashkenazi_lowSNP_pca/` |project 544 Ashkenazi genomes onto the same PCA | supplementary figure |
| `09_unmapped_taxonomy/`    |Kraken2 taxonomy of the unmapped reads | supplementary result |

## Software

bwa, samtools, bedtools, AdapterRemoval v2, ATLAS, mapDamage2, EIGENSOFT (mergeit,
smartpca), pileupCaller (sequenceTools), PLINK, ADMIXTURE, pong, Haplogrep3, YLeaf,
contamMix, Kraken2, seqtk; Python 3 (numpy, pandas, matplotlib, seaborn).

## Reference data

- Reference genome: hs37d5 (GRCh37 + decoy).
- Reference panel: AADR **v54.1** Human Origins (`v54.1.p1_HO_public`), 597,573 SNPs;
  446 West-Eurasian / North-African genomes used for projection.
- Ashkenazi genomes: Lencz et al. 2018.
- Kraken2 DB: PlusPF (benlangmead.github.io/aws-indexes/k2).

