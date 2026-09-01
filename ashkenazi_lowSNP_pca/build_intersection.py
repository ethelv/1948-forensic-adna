#!/usr/bin/env python
import sys
BASE = ".."
AJ = "."

aj_snp = [l.split()[0] for l in open(f"{AJ}/AJ.snp")]
target = set(aj_snp)
aj_rows = {}
with open(f"{AJ}/AJ.geno") as g:
    for rs, row in zip(aj_snp, g):
        aj_rows[rs] = row.strip()
w = len(next(iter(aj_rows.values())))
sys.stderr.write(f"target SNPs: {len(target)}  AJ width: {w}\n")

n = 0
with open(f"{BASE}/merged.snp") as snp, open(f"{BASE}/merged.geno") as geno, \
     open(f"{AJ}/merged_aj.snp", "w") as osnp, open(f"{AJ}/merged_aj.geno", "w") as ogeno:
    for sline, gline in zip(snp, geno):
        rs = sline.split()[0]
        if rs in target:
            ogeno.write(gline.rstrip("\n") + aj_rows[rs] + "\n")
            osnp.write(sline)
            n += 1
sys.stderr.write(f"wrote {n} intersection rows to merged_aj.geno / merged_aj.snp\n")
