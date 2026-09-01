#!/usr/bin/env bash
# PCA ancestry inference using SmartPCA (EIGENSOFT)
# Paper: "We used principal components analysis (PCA) to infer the ancestry of the
#         target genome. We used SmartPCA with the target genome projected on the axes
#         formed by the AADR genomes. Least-squares projection was enabled and no
#         outlier filtering iterations were applied."
# Reference panel: 446 West Eurasian and North African genomes from AADR (v54.1)
# SNPs used: 597,573 (Human Origins panel); 36,098 covered by at least one read

# Step 1: Convert pseudo-haploid calls to EIGENSTRAT format
#         (see genotyping_pileupCaller.ipynb for the pileupCaller step)

# Step 2: Merge target genome with AADR reference panel
mergeit -p merge.par
# merge.par merges JD206 eigenstrat files with v54.1.p1_HO_public.{geno,snp,ind}

# Alternative merge order (target first):
# mergeit -p mergeit.par

# Step 3: Convert merged files to BED for filtering, then back to EIGENSTRAT
plink --file merged --make-bed --out merged

# Step 4: Run SmartPCA
#   - lsqproject: YES  → project target onto axes without influencing them
#   - numoutlieriter: 0 → no outlier removal
smartpca -p smartpca.par > merged_pca.log

# Step 5: Plot results
# See figure2_pca.ipynb
