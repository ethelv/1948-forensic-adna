#!/usr/bin/env python
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

AJDIR = "."

# ---- load projection + restrict to the manuscript's population set ----
cols = ["Sample"] + [f"PC{i}" for i in range(1, 11)] + ["Group"]
df = pd.read_csv(f"{AJDIR}/merged_aj_pca.evec", delim_whitespace=True, comment="#",
                 header=None, names=cols)
pops = [l.strip() for l in open(f"{AJDIR}/../pops_list") if l.strip()]
keep = set(pops) | {"Unknown", "Lencz_AJ"}
df = df[df.Group.isin(keep)].copy()

# axis orientation already matches the paper (JD206 ~ +0.078, -0.051); flip only if needed
def mpc(pop, pc):
    s = df.loc[df.Group == pop, pc]; return s.mean() if len(s) else np.nan
if mpc("BedouinB.HO", "PC1") < mpc("French.HO", "PC1"): df["PC1"] = -df["PC1"]
if mpc("Moroccan.HO", "PC2") > mpc("French.HO", "PC2"): df["PC2"] = -df["PC2"]

# ---- group scheme identical to Figure 2 ----
european = ['French.HO','Basque.HO','Bulgarian.HO','Hungarian.HO','Lithuanian.HO','Mordovian.HO',
 'Russian.HO','Ukrainian.HO','Belarusian.HO','Estonian.HO','Czech.HO','Icelandic.HO','Scottish.HO',
 'English.HO','Spanish.HO','Spanish_North.HO','Finnish.HO','Maltese.HO','Croatian.HO','Italian_North.HO',
 'Italian_South.HO','Sardinian.HO','Norwegian.HO','Albanian.HO','Romanian.HO','Sicilian.HO','Orcadian.HO','Cypriot.HO']
caucasus = ['Adygei.HO','Abkhasian.HO','Armenian.HO','Balkar.HO','Chechen.HO','Lezgin.HO','Kumyk.HO',
 'Georgian.HO','Ossetian.HO','Iran_Zoroastrian.HO','Turkish.HO','Iranian.HO','Jew_Iranian.HO','Assyrian.HO',
 'Jew_Turkish.HO','Jew_Georgian.HO','Syrian.HO','Lebanese_Christian.HO','Jew_Iraqi.HO']
grey_groups = set(european + caucasus)

# Ashkenazi HO is drawn separately (the reference we compare the Lencz cloud to)
non_grey = [g for g in df.Group.unique()
            if g not in grey_groups and g not in ("Unknown", "Lencz_AJ", "Jew_Ashkenazi.HO")]

palette = sns.color_palette("tab20", n_colors=len(non_grey))
markers = ["o","s","D","^","v","P","X","<",">","*","h","H","d","p","8"]
style = {g: {"color": palette[i % len(palette)], "marker": markers[i % len(markers)]}
         for i, g in enumerate(non_grey)}

# ---- figure (Figure-2 layout) ----
plt.figure(figsize=(11, 8.5))

# grey Europe/Caucasus backdrop
g = df[df.Group.isin(grey_groups)]
sns.scatterplot(x="PC1", y="PC2", data=g, color="lightgrey", s=20, alpha=0.7, label="Europe/Caucasus")

# coloured reference groups
for grp in non_grey:
    d = df[df.Group == grp]; st = style[grp]
    plt.scatter(d.PC1, d.PC2, color=st["color"], marker=st["marker"], s=40, alpha=0.6,
                label=grp.replace(".HO", "").replace("_", " "))

# Lencz AJ cloud (the focal distribution) — behind the anchors
aj = df[df.Group == "Lencz_AJ"]
plt.scatter(aj.PC1, aj.PC2, color="#ff7f0e", marker="o", s=22, alpha=0.5, zorder=4,
            edgecolor="none", label=f"Lencz AJ (35k SNPs, n={len(aj)})")

# full-SNP HO Ashkenazi reference
ho = df[df.Group == "Jew_Ashkenazi.HO"]
plt.scatter(ho.PC1, ho.PC2, color="#08306b", marker="*", s=200, zorder=6,
            edgecolor="white", linewidth=0.8, label=f"Ashkenazi (HO, full SNPs, n={len(ho)})")

# sample 1948 — true point + arrow
s48 = df[df.Group == "Unknown"]
if not s48.empty:
    px, py = s48.PC1.iloc[0], s48.PC2.iloc[0]
    plt.scatter(px, py, color="red", marker="o", s=45, edgecolor="black", linewidth=0.8,
                zorder=7, label="sample 1948")
    plt.annotate("sample 1948", xy=(px, py), xycoords="data",
                 xytext=(0.093, -0.075), textcoords="data",
                 fontsize=13, fontweight="bold", color="red", ha="center", va="top",
                 arrowprops=dict(arrowstyle="-|>", color="red", lw=2.2, shrinkB=4), zorder=8)

plt.xlabel("Principal Component 1", fontsize=17)
plt.ylabel("Principal Component 2", fontsize=17)
plt.ylim(-0.12, 0.055); plt.xlim(-0.05, 0.105); plt.grid(False)
plt.xticks(fontsize=14); plt.yticks(fontsize=14)

# legend outside right, Lencz AJ / Ashkenazi / sample 1948 last
ax = plt.gca(); h, l = ax.get_legend_handles_labels()
tail = [x for x in l if x.startswith("Lencz AJ") or x.startswith("Ashkenazi") or x == "sample 1948"]
order = [x for x in l if x not in tail] + tail
hd = dict(zip(l, h))
ax.legend([hd[x] for x in order], order, title="Group", bbox_to_anchor=(1.02, 1),
          loc="upper left", markerscale=1.3, fontsize=9, title_fontsize=11)

plt.tight_layout()
for ext in ("png", "pdf"):
    plt.savefig(f"{AJDIR}/FigureSX_AJ_lowSNP_PCA.{ext}", dpi=300, bbox_inches="tight")
print(f"saved FigureSX_AJ_lowSNP_PCA.png/.pdf  (Lencz AJ n={len(aj)}, Ashkenazi HO n={len(ho)})")
