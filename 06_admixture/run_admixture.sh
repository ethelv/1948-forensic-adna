#!/usr/bin/env bash
set -euo pipefail              # safer shell defaults

# --------- edit only these lines ---------
PREFIX=merged_pruned_subset_keep    # base name of your BED/FAM/BIM trio
CORES=8                        # threads for ADMIXTURE  (-j flag)
CONV=1e-3                      # convergence criterion   (-C flag)
Ks="5 6 7 8"                   # K range
REPS=3                         # how many replicates per K
SEEDS=(123 456 789)            # supply exactly $REPS unique ints
# -----------------------------------------

mkdir -p admixture_runs/logs

for K in $Ks; do
  for r in $(seq 1 $REPS); do
    SEED=${SEEDS[$((r-1))]}
    echo "▶  K=$K  replicate=$r  seed=$SEED"

    admixture -j"$CORES" -C "$CONV" --seed="$SEED" \
              "${PREFIX}.bed" "$K" | \
        tee "admixture_runs/logs/${PREFIX}.K${K}.r${r}.log"

    mv "${PREFIX}.${K}.Q" "admixture_runs/${PREFIX}.K${K}.r${r}.Q"
    mv "${PREFIX}.${K}.P" "admixture_runs/${PREFIX}.K${K}.r${r}.P"
  done
done
