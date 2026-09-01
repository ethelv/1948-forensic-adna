#!/usr/bin/env python
"""
Build a pseudo-haploid EIGENSTRAT .geno for the 544 Lencz AJ genomes over the
35,968 HO SNPs that were covered in JD206.

For each SNP x sample: draw ONE random allele (mimicking a single ancient read).
 - Sites present in the AJ cohort VCF: draw from the called diploid GT.
 - Sites absent from the VCF: the 30x genome is hom-ref -> drawn allele = b37 ref.
Encoding (EIGENSTRAT, count of allele a1 = JD206.snp col5):
   drawn == a1 -> 2 ; drawn == a2 -> 0 ; missing / third allele -> 9
All 35,968 sites are on the b37 forward strand (verified: 0 strand flips),
so VCF alleles match HO alleles directly.
"""
import random, sys

random.seed(1948)

AJDIR = "."

# --- site table (order = AJ.snp order) ---
snp_order = []            # list of (chr,pos)
a1a2 = {}                 # (chr,pos) -> (a1,a2)
for line in open(f"{AJDIR}/AJ.snp"):
    rs, c, gpos, pos, a1, a2 = line.split()
    snp_order.append((c, pos))
    a1a2[(c, pos)] = (a1, a2)

b37ref = {}
with open(f"{AJDIR}/sites_full.tsv") as fh:
    next(fh)
    for line in fh:
        c, pos, rs, a1, a2, rb, flip = line.rstrip("\n").split("\t")
        b37ref[(c, pos)] = rb

# --- sample order actually in the extracted subset ---
samples = [x.strip() for x in open(f"{AJDIR}/aj_subset_samples.txt") if x.strip()]
nsam = len(samples)
sys.stderr.write(f"samples: {nsam}\n")

# --- parse extracted genotypes for present sites ---
# aj_gts.tsv: CHROM  POS  REF  ALT  GT1 GT2 ... GTn   (n = nsam, same order as samples)
present = {}   # (chr,pos) -> (REF, [ALT...], [gt strings])
with open(f"{AJDIR}/aj_gts.tsv") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        c, pos, ref, alt = f[0], f[1], f[2], f[3]
        gts = f[4:]
        present[(c, pos)] = (ref, alt.split(","), gts)
sys.stderr.write(f"present sites in VCF: {len(present)} / {len(snp_order)}\n")

def draw_allele(gt, ref, alts):
    """Return one randomly chosen allele base, or None if missing/symbolic."""
    sep = "/" if "/" in gt else "|"
    parts = gt.replace("|", "/").split("/")
    if len(parts) < 2:
        parts = parts * 2
    a, b = parts[0], parts[1]
    if a == "." or b == ".":
        return None
    idx = random.choice((a, b))
    ai = int(idx)
    base = ref if ai == 0 else alts[ai - 1]
    if base in ("<NON_REF>", "*") or base.startswith("<"):
        return None
    return base

geno_rows = []
counts = {"present": 0, "absent": 0}
miss_total = 0
for key in snp_order:
    a1, a2 = a1a2[key]
    row = []
    if key in present:
        counts["present"] += 1
        ref, alts, gts = present[key]
        for gt in gts:
            base = draw_allele(gt, ref, alts)
            if base == a1:
                row.append("2")
            elif base == a2:
                row.append("0")
            else:
                row.append("9")
    else:
        counts["absent"] += 1
        rb = b37ref[key]                 # hom-ref for every sample
        val = "2" if rb == a1 else ("0" if rb == a2 else "9")
        row = [val] * nsam
    miss_total += row.count("9")
    geno_rows.append("".join(row))

with open(f"{AJDIR}/AJ.geno", "w") as out:
    out.write("\n".join(geno_rows) + "\n")

# rewrite AJ.ind to match extracted sample order
with open(f"{AJDIR}/AJ.ind", "w") as out:
    for s in samples:
        out.write(f"{s}\tU\tLencz_AJ\n")

tot = len(snp_order) * nsam
sys.stderr.write(f"SNPs present={counts['present']} absent={counts['absent']}\n")
sys.stderr.write(f"missing genotypes: {miss_total}/{tot} = {miss_total/tot:.4f}\n")
sys.stderr.write("wrote AJ.geno, AJ.ind\n")
