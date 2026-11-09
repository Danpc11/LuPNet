#!/usr/bin/env bash
# SLURM array: the three cohorts of the manuscript re-simulated at the barometric pressure of Mexico City
# (PB = 585 mmHg, 2,240 m), same seeds, so the only difference from big2 / adv2 / altF is the altitude.
#
#   mkdir -p logs && sbatch scripts/altitude_cdmx.slurm
#   # when the array finishes:
#   for t in big2_cdmx adv2_cdmx altF_cdmx; do python3 -c "import sys; sys.path.insert(0,'src'); from cohort_io import pack; pack('results/ph_cohorts','$t')"; done
#   python3 src/ph_cohort.py analyze --tag big2_cdmx
#
# 20 + 6 + 10 = 36 shards of 10 cases, about 5 min each on one core.
#SBATCH --job-name=lupnet-cdmx
#SBATCH --array=0-35%36
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:40:00
#SBATCH --output=logs/cdmx_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# source .venv/bin/activate

PB=585
t=$SLURM_ARRAY_TASK_ID
if   [ "$t" -lt 20 ]; then TAG=big2_cdmx; ARGS="--n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85"; s=$t;           n=200
elif [ "$t" -lt 26 ]; then TAG=adv2_cdmx; ARGS="--n 60  --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45"; s=$(( t - 20 )); n=60
else                       TAG=altF_cdmx; ARGS="--n 100 --seed 7000 --pseed 40000 --vmode fibrosis";    s=$(( t - 26 )); n=100
fi
i0=$(( s * 10 )); i1=$(( i0 + 10 )); [ "$i1" -gt "$n" ] && i1=$n
echo "$TAG cases $i0-$((i1 - 1)) PB=$PB on $(hostname)"
# shellcheck disable=SC2086
python3 src/ph_cohort.py run --tag "$TAG" $ARGS --PB $PB --i0 "$i0" --i1 "$i1"
