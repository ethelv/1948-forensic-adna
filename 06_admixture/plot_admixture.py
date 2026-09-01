"""
Recreate the paper's pong-style ADMIXTURE figures (Figure 3 = K6, Figure S2 = K7,K8)
in matplotlib, but drawing the target genome (JD206 / "1948 sample") as a WIDE,
clearly-visible band instead of a 1-px bar (reviewer #2).

Faithful to the pong run that made the paper figures:
  pong -m filemap -i ind2pop_cleaned.txt -n pop_order.txt -s 1
  - default 9-colour palette (K_max = 8 <= 9)
  - population order from pop_order.txt
  - individuals sorted within each population (membership gradient)
  - major-mode representative run per K: K6r1 (Fig 3); K7r1, K8r2 (Fig S2)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import os

ADMIX_DIR = './admixture_runs'
OUT_DIR   = '.'
os.chdir(ADMIX_DIR)

# pong default palette (static/pong.js:  var colors = [...])
PALETTE = ['#E04B4B', '#6094C3', '#63BC6A', '#A76BB2', '#F0934E',
           '#FEFB54', '#B37855', '#EF91CA', '#A4A4A4']

TARGET_POP = 'Unknown'

# Named palette entries
RED, BLUE, GREEN, PURPLE, ORANGE, YELLOW, BROWN, PINK, GREY = PALETTE

# Map each ADMIXTURE cluster to a fixed colour by the reference population that
# defines it, so the colours match the published pong figure (Iran_Zoroastrian =
# blue, Druze = green, BedouinB = orange, North Africans + target = red, etc.).
SIG_COLOR = {
    'Iran_Zoroastrian.HO': BLUE,
    'Druze.HO': GREEN,
    'BedouinB.HO': ORANGE,
    'Algerian.HO': RED, 'Moroccan.HO': RED, 'Tunisian.HO': RED, 'Libyan.HO': RED,
    'Palestinian.HO': BROWN,
    'Jew_Yemenite.HO': YELLOW, 'BedouinA.HO': YELLOW, 'Saudi.HO': YELLOW,
    'Jew_Ashkenazi.HO': PURPLE, 'Italian_South.HO': PURPLE,
    'Jew_Georgian.HO': PURPLE, 'Turkish.HO': PURPLE,
    'Jew_Libyan.HO': PINK, 'Jew_Tunisian.HO': PINK, 'Jew_Moroccan.HO': PINK,
}

# ── Inputs shared by all panels ──────────────────────────────────────────────
ind_pop = pd.read_csv('ind2pop_cleaned.txt', header=None)[0].values   # pop per Q row
pop_order = pd.read_csv('pop_order.txt', header=None)[0].tolist()
ref_pops = [p for p in pop_order if p != TARGET_POP]
target_row = int(np.where(ind_pop == TARGET_POP)[0][0])

# Layout constants for the wide target band
GAP          = 6.0
TARGET_WIDTH = 16.0


def cluster_color_perm(Q):
    """Assign each cluster the palette colour of its defining (signature)
    population, so the panel matches the published pong colours. Any cluster
    whose signature population is not in SIG_COLOR falls back to a spare colour."""
    K = Q.shape[1]
    pop_mean = {p: Q[np.where(ind_pop == p)[0]].mean(axis=0) for p in ref_pops}
    M = np.array([pop_mean[p] for p in ref_pops])          # pops x K
    ccolor, used = {}, set()
    for k in range(K):
        sig_pop = ref_pops[int(np.argmax(M[:, k]))]
        c = SIG_COLOR.get(sig_pop)
        if c is None or c in used:
            c = next(col for col in PALETTE if col not in used)
        ccolor[k] = c
        used.add(c)
    return ccolor


def sort_within_pop(rows, Q):
    """Order individuals within a population by their membership gradient
    (descending in the population's dominant cluster), like pong/distruct."""
    if len(rows) <= 1:
        return list(rows)
    dom = np.argmax(Q[rows].mean(axis=0))
    return list(rows[np.argsort(-Q[rows, dom])])


def draw_panel(ax, Q, ccolor, top_panel=False, show_xlabels=False):
    K = Q.shape[1]

    # Build the ordered list of reference individuals + per-pop label centres
    ordered_rows, boundaries, label_pos, label_txt = [], [], [], []
    x = 0
    for p in ref_pops:
        rows = np.where(ind_pop == p)[0]
        rows = sort_within_pop(rows, Q)
        if x > 0:
            boundaries.append(x - 0.5)
        label_pos.append(x + (len(rows) - 1) / 2.0)
        label_txt.append(p.replace('.HO', '').replace('_', ' '))
        ordered_rows.extend(rows)
        x += len(rows)
    n_ref = len(ordered_rows)
    ordered_rows = np.array(ordered_rows)

    target_center = (n_ref - 1) + GAP + TARGET_WIDTH / 2.0
    x_right = target_center + TARGET_WIDTH / 2.0 + 2.0

    # Reference individuals (unit width)
    ref_mat = Q[ordered_rows]
    bottom = np.zeros(n_ref)
    for k in range(K):
        ax.bar(np.arange(n_ref), ref_mat[:, k], bottom=bottom,
               color=ccolor[k], width=1.0, linewidth=0)
        bottom += ref_mat[:, k]

    # Target: single wide band after the gap
    tgt = Q[target_row]
    b = 0.0
    for k in range(K):
        ax.bar(target_center, tgt[k], bottom=b, color=ccolor[k],
               width=TARGET_WIDTH, linewidth=0, zorder=5)
        b += tgt[k]
    ax.add_patch(Rectangle((target_center - TARGET_WIDTH / 2.0, 0),
                           TARGET_WIDTH, 1.0, fill=False,
                           edgecolor='black', linewidth=1.6, zorder=6))

    for vl in boundaries:
        ax.axvline(x=vl, color='black', linewidth=0.8)
    # frame the reference block
    ax.add_patch(Rectangle((-0.5, 0), n_ref, 1.0, fill=False,
                           edgecolor='black', linewidth=1.2, zorder=6))

    ax.set_xlim(-0.5, x_right)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.tick_params(axis='y', labelsize=9)
    ax.set_xticks([])

    if top_panel:
        ax.annotate('1948 sample\n(target)',
                    xy=(target_center, 1.01), xytext=(target_center, 1.30),
                    xycoords='data', textcoords='data',
                    ha='center', va='bottom', fontsize=12, fontweight='bold',
                    color='red',
                    arrowprops=dict(arrowstyle='-|>', color='red', lw=2.5),
                    annotation_clip=False)

    if show_xlabels:
        ax.set_xticks(list(label_pos) + [target_center])
        ax.set_xticklabels(label_txt + ['1948 sample'],
                           rotation=45, ha='right', fontsize=7.5)
        for lab in ax.get_xticklabels():
            if lab.get_text() == '1948 sample':
                lab.set_color('red')
                lab.set_fontweight('bold')


def load_Q(fname):
    return pd.read_csv(fname, sep=r'\s+', header=None).values


def make_figure(runs, output_file):
    n = len(runs)
    fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(16, 2.8 * n + 2.5))
    if n == 1:
        axes = [axes]
    for i, (ax, fname) in enumerate(zip(axes, runs)):
        Q = load_Q(fname)
        ccolor = cluster_color_perm(Q)
        draw_panel(ax, Q, ccolor,
                   top_panel=(i == 0), show_xlabels=(i == n - 1))
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print('Saved', output_file)


# Figure 3: K=6 major mode (K6r1)
make_figure(['merged_pruned_subset_keep.K6.r1.Q'],
            os.path.join(OUT_DIR, 'Figure_3_admixture_K6.png'))

# Figure S2: K=7 major mode (K7r1) then K=8 major mode (K8r2)
make_figure(['merged_pruned_subset_keep.K7.r1.Q',
             'merged_pruned_subset_keep.K8.r2.Q'],
            os.path.join(OUT_DIR, 'Figure_S2_admixture_K7_K8.png'))
