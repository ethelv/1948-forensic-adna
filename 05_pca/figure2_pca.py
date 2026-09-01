#!/usr/bin/env python
"""
Figure 2 - PCA ancestry inference.

Plots the SmartPCA projection (merged_pca.evec) of the target genome onto axes
formed by the AADR West-Eurasian / North-African reference panel. Europe/Caucasus
populations are greyed; Middle-Eastern / North-African / Jewish groups are coloured.
The target (sample 1948) is drawn at its true position and indicated with an arrow.

Input (not shipped): merged_pca.evec (SmartPCA output; see run_pca.sh) and pops_list.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

EVEC = "merged_pca.evec"
POPS = "pops_list"

# ---- load projection, restrict to the reference population set ----
cols = ["Sample"] + [f"PC{i}" for i in range(1, 11)] + ["Group"]
df = pd.read_csv(EVEC, delim_whitespace=True, comment="#", header=None).reset_index()
df.columns = cols
pops = list(pd.read_csv(POPS, header=None)[0]) + ["Unknown"]
df = df[df.Group.isin(pops)].copy()
df["Group"] = df["Group"].replace("Unknown", "sample_1948")

# ---- greyed background populations (Europe + Caucasus) ----
european = ['French.HO','Basque.HO','Bulgarian.HO','Hungarian.HO','Lithuanian.HO','Mordovian.HO',
 'Russian.HO','Ukrainian.HO','Belarusian.HO','Estonian.HO','Czech.HO','Icelandic.HO','Scottish.HO',
 'English.HO','Spanish.HO','Spanish_North.HO','Finnish.HO','Maltese.HO','Croatian.HO','Italian_North.HO',
 'Italian_South.HO','Sardinian.HO','Norwegian.HO','Albanian.HO','Romanian.HO','Sicilian.HO','Orcadian.HO','Cypriot.HO']
caucasus = ['Adygei.HO','Abkhasian.HO','Armenian.HO','Balkar.HO','Chechen.HO','Lezgin.HO','Kumyk.HO',
 'Georgian.HO','Ossetian.HO','Iran_Zoroastrian.HO','Turkish.HO','Iranian.HO','Jew_Iranian.HO','Assyrian.HO',
 'Jew_Turkish.HO','Jew_Georgian.HO','Syrian.HO','Lebanese_Christian.HO','Jew_Iraqi.HO']
grey_groups = set(european + caucasus)

special = "sample_1948"
non_grey = [g for g in df.Group.unique() if g not in grey_groups and g not in ("Unknown", special)]

palette = sns.color_palette("tab20", n_colors=len(non_grey))
markers = ["o","s","D","^","v","P","X","<",">","*","h","H","d","p","8"]
style = {g: {"color": palette[i % len(palette)], "marker": markers[i % len(markers)]}
         for i, g in enumerate(non_grey)}

# ---- plot ----
plt.figure(figsize=(10, 8))
sns.scatterplot(x="PC1", y="PC2", data=df[df.Group.isin(grey_groups)],
                color="lightgrey", s=20, alpha=0.7, label="Europe/Caucasus")
for grp in non_grey:
    d = df[df.Group == grp]; st = style[grp]
    plt.scatter(d.PC1, d.PC2, color=st["color"], marker=st["marker"], s=40, alpha=0.6,
                label=grp.replace(".HO", "").replace("_", " "))

# target: true data point + arrow (do NOT paste an oversized marker over it)
h = df[df.Group == special]
if not h.empty:
    px, py = h.PC1.iloc[0], h.PC2.iloc[0]
    plt.scatter(px, py, color="red", marker="o", s=45, edgecolor="black", linewidth=0.8,
                zorder=6, label="sample 1948")
    plt.annotate("sample 1948", xy=(px, py), xytext=(0.093, -0.075), textcoords="data",
                 fontsize=13, fontweight="bold", color="red", ha="center", va="top",
                 arrowprops=dict(arrowstyle="-|>", color="red", lw=2.2, shrinkB=4), zorder=7)

plt.xlabel("Principal Component 1", fontsize=17)
plt.ylabel("Principal Component 2", fontsize=17)
plt.ylim(-0.12, 0.055); plt.xlim(-0.05, 0.105); plt.grid(False)
plt.xticks(fontsize=14); plt.yticks(fontsize=14)

ax = plt.gca(); handles, labels = ax.get_legend_handles_labels()
seen, hh, ll = set(), [], []
for hnd, lab in zip(handles, labels):
    if lab not in seen:
        seen.add(lab); hh.append(hnd); ll.append(lab)
ax.legend(hh, ll, title="Group", bbox_to_anchor=(1.05, 1), loc="upper left",
          markerscale=1.5, fontsize=10, title_fontsize=12)

plt.tight_layout()
plt.savefig("pca_plot_high_res.png", dpi=300, bbox_inches="tight")
