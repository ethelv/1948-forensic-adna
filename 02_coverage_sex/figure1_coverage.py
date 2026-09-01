#!/usr/bin/env python
"""
Figure 1 - per-chromosome coverage and sex inference.

Coverage from the mapped, deduplicated BAM, keeping reads with MAPQ >= 10.
Both panels use mappable sequence only: reference assembly gaps / heterochromatin
(centromeres, 1q/9q/16q blocks, acrocentric p-arms, Yq12) and the pseudoautosomal
regions (PAR1/PAR2 on X & Y) are excluded; chrY is restricted to its assembled
male-specific euchromatin (0-28.82 Mb).

Panel A - mean depth +/- SD (1 Mb windows) per chromosome. Autosomes uniform;
          chrX ~ 1/2 the autosomal mean; chrY low -> male.
Panel B - coverage profile, x-axis = % of each chromosome's mappable length.

Inputs (not shipped): coverage_mapq{MAPQ}.bed and coverage.bed (all-reads), produced by
    samtools view -b -q {MAPQ} in.bam | bedtools genomecov -ibam - -bga
Reference build: hs37d5 (GRCh37 + decoy).
"""
import os, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# ---- paths (edit for your setup) ----
DATA = "."                                          # dir holding the BAM / BED files
BAM  = f"{DATA}/JD206_libmerged_rg_rmdup.mapped.bam"
MAPQ = 10
BED     = f"{DATA}/coverage_mapq{MAPQ}.bed"          # MAPQ-filtered per-interval depth (bedgraph)
BED_ALL = f"{DATA}/coverage.bed"                     # all-reads depth -> detect assembly gaps
CHROMS  = [str(i) for i in range(1, 23)] + ["X", "Y"]
WIN = 1_000_000

PAR = {"X": [(60001, 2699520), (154931044, 155260560)],
       "Y": [(10001, 2649520), (59034050, 59373566)]}
Y_EUCHR_END = 28_819_151          # GRCh37 Yq12 N-gap start (Skaletsky 2003, Nature 423:825)
EFF = {"Y": Y_EUCHR_END}
GAP_MIN, MIN_MAPPABLE, MIN_NONPAR = 500_000, 200_000, WIN // 10

# ---- one-time MAPQ-filtered coverage ----
if not os.path.exists(BED):
    subprocess.run(f"samtools view -b -q {MAPQ} {BAM} | bedtools genomecov -ibam - -bga > {BED}",
                   shell=True, check=True)

cov = pd.read_csv(BED, sep="\t", header=None, names=["chrom", "start", "end", "depth"],
                  dtype={"chrom": str, "start": np.int64, "end": np.int64, "depth": np.float64})
cov = cov[cov.chrom.isin(CHROMS)]

# assembly gaps / heterochromatin = long depth-0 runs in the all-reads coverage
_a = pd.read_csv(BED_ALL, sep="\t", header=None, names=["chrom", "start", "end", "depth"],
                 dtype={"chrom": str, "start": np.int64, "end": np.int64, "depth": np.float64})
_a = _a[_a.chrom.isin(CHROMS)]
GAPS = {c: [(int(r.start), int(r.end)) for r in
            _a[(_a.chrom == c) & (_a.depth == 0) & (_a.end - _a.start >= GAP_MIN)].itertuples()]
        for c in CHROMS}


def _cum(chrom):
    """C(x) = integral of depth over [0, x], and chromosome length L."""
    sub = cov[cov.chrom == chrom]
    s, e, d = sub.start.to_numpy(), sub.end.to_numpy(), sub.depth.to_numpy()
    ce = np.cumsum(d * (e - s)); cs = ce - d * (e - s)
    def C(x):
        x = np.asarray(x, float); i = np.clip(np.searchsorted(e, x, side="right"), 0, len(e) - 1)
        return cs[i] + d[i] * (x - s[i])
    return C, e[-1]


def _excl(C, intervals, a, b):
    """Depth-integral and length of `intervals` overlapping each window [a, b)."""
    ii = np.zeros_like(a, float); ll = np.zeros_like(a, float)
    for p0, p1 in intervals:
        o0, o1 = np.clip(a, p0, p1), np.clip(b, p0, p1)
        ii += C(o1) - C(o0); ll += (o1 - o0)
    return ii, ll


# ---- Panel A: mean depth + SD over mappable sequence (gaps + PAR excluded, Y over MSY euchromatin) ----
rows = []
for chrom in CHROMS:
    C, L = _cum(chrom)
    eff = EFF.get(chrom, L)
    exclude = GAPS[chrom] + PAR.get(chrom, [])
    ex_i, ex_l = _excl(C, exclude, np.array([0.0]), np.array([float(eff)]))
    mean_depth = (C(np.array([float(eff)]))[0] - ex_i[0]) / (eff - ex_l[0])
    edges = np.append(np.arange(0, eff, WIN, dtype=np.int64), eff)
    a, b = edges[:-1].astype(float), edges[1:].astype(float)
    ei, el = _excl(C, exclude, a, b)
    with np.errstate(invalid="ignore", divide="ignore"):
        wm = np.where((b - a) - el >= MIN_MAPPABLE, ((C(b) - C(a)) - ei) / ((b - a) - el), np.nan)
    rows.append((chrom, mean_depth, float(np.nanstd(wm))))
stats = pd.DataFrame(rows, columns=["chrom", "mean_depth", "sd_windows"])
auto_mean = stats.loc[~stats.chrom.isin(["X", "Y"]), "mean_depth"].mean()


# ---- Panel B: PAR-masked profile, x-axis = % of mappable length (Y -> MSY euchromatin) ----
def profile_pct(chrom):
    C, L = _cum(chrom)
    eff = EFF.get(chrom, L)
    edges = np.append(np.arange(0, eff, WIN, dtype=np.int64), eff)
    a, b = edges[:-1].astype(float), edges[1:].astype(float)
    ei, el = _excl(C, PAR.get(chrom, []), a, b)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where((b - a) - el >= MIN_NONPAR, ((C(b) - C(a)) - ei) / ((b - a) - el), np.nan)
    return (a + b) / 2 / eff * 100.0, m
pct = {c: profile_pct(c) for c in ["1", "X", "Y"]}

si = stats.set_index("chrom")
print(f"autosomal {auto_mean:.4f}; X {si.loc['X','mean_depth']:.4f} "
      f"(X/auto={si.loc['X','mean_depth']/auto_mean:.2f}); Y {si.loc['Y','mean_depth']:.4f} "
      f"(Y/auto={si.loc['Y','mean_depth']/auto_mean:.2f})")

# ---- plot ----
order = [str(i) for i in range(1, 23)] + ["X", "Y"]
s = stats.set_index("chrom").loc[order].reset_index()
colors = np.where(s["chrom"].isin(["X", "Y"]), "lightsalmon", "#4C72B0")
fig = plt.figure(figsize=(11, 7.5))
gs = gridspec.GridSpec(2, 1, hspace=0.42)

axA = fig.add_subplot(gs[0]); x = np.arange(len(s))
axA.bar(x, s.mean_depth, yerr=[np.minimum(s.sd_windows, s.mean_depth), s.sd_windows],
        color=colors, edgecolor="black", linewidth=0.4, capsize=2.5, error_kw=dict(lw=0.8))
axA.set_xticks(x); axA.set_xticklabels(s["chrom"], fontsize=9)
axA.set_xlabel("Chromosome", fontsize=11); axA.set_ylabel("Mean depth", fontsize=11)
axA.text(-0.055, 1.02, "A", transform=axA.transAxes, fontsize=13, fontweight="bold")
axA.legend(handles=[Patch(color="#4C72B0", label="Autosomes"),
                    Patch(color="lightsalmon", label="Sex chromosomes")], fontsize=9, frameon=False)

axB = fig.add_subplot(gs[1])
for chrom, col in [("1", "#4C72B0"), ("X", "seagreen"), ("Y", "crimson")]:
    p, m = pct[chrom]; axB.plot(p, m, color=col, lw=1.4, marker="o", ms=3, label="chr" + chrom)
axB.set_xlim(0, 100); axB.set_ylim(bottom=0)
axB.set_xlabel("Position along chromosome (% of mappable length)", fontsize=11)
axB.set_ylabel(f"Mean depth in {WIN//10**6} Mb window", fontsize=11)
axB.text(-0.055, 1.02, "B", transform=axB.transAxes, fontsize=13, fontweight="bold")
axB.legend(fontsize=9, frameon=False)

for ext in ("png", "pdf"):
    fig.savefig(f"Figure1_coverage.{ext}", dpi=300, bbox_inches="tight")
