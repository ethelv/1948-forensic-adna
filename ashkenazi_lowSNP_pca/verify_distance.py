#!/usr/bin/env python
"""
Verify that sample 1948 (JD206) is farther, in Euclidean distance over (PC1, PC2),
from the Ashkenazi cluster centroid than EVERY Ashkenazi genome.

Distance in (PC1,PC2) is invariant to smartpca's arbitrary per-axis sign, so no
axis re-orientation is needed. Uses the canonical 544-sample projection.
"""
import pandas as pd, numpy as np

EVEC = "./merged_aj_pca.evec"

cols = ["Sample"] + [f"PC{i}" for i in range(1, 11)] + ["Group"]
d = pd.read_csv(EVEC, delim_whitespace=True, comment="#", header=None, names=cols)

# target = sample 1948 (group "Unknown" in the EIGENSTRAT ind)
t = d.loc[d.Group == "Unknown", ["PC1", "PC2"]].values[0]


def report(name, sub):
    C = sub[["PC1", "PC2"]].mean().values                      # cluster centroid
    di = np.sqrt(((sub[["PC1", "PC2"]].values - C) ** 2).sum(axis=1))   # each AJ -> centroid
    d1948 = np.sqrt(((t - C) ** 2).sum())                      # sample 1948 -> centroid
    n_ge = int((di >= d1948).sum())
    print(f"--- Ashkenazi cluster = {name} (n={len(sub)}) ---")
    print(f"  centroid (PC1,PC2)                    = ({C[0]:+.4f}, {C[1]:+.4f})")
    print(f"  max Ashkenazi->centroid distance      = {di.max():.4f}  ({sub.iloc[di.argmax()].Sample})")
    print(f"  mean Ashkenazi->centroid distance     = {di.mean():.4f}")
    print(f"  sample 1948 -> centroid distance      = {d1948:.4f}")
    print(f"  Ashkenazim as far or farther than 1948= {n_ge} / {len(sub)}")
    print(f"  1948 farther than ALL Ashkenazim      = {d1948 > di.max()}  "
          f"(x{d1948/di.max():.2f} the farthest; +{d1948-di.max():.4f})\n")


aj = d[d.Group == "Lencz_AJ"]              # 544 Lencz AJ, 35k SNPs, pseudo-haploid
ho = d[d.Group == "Jew_Ashkenazi.HO"]      # 7 HO Ashkenazi, full SNPs

report("Lencz AJ 35k-SNP (n=544)", aj)
report("HO Ashkenazi full-SNP (n=7)", ho)
report("all Ashkenazi combined", pd.concat([aj, ho]))
