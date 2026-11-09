#!/usr/bin/env bash
# SLURM array: the calibrated cohort re-simulated with the UIP and the fibrotic NSIP pattern, same seeds and
# lesions as big2, so the only difference between the three cohorts is the histological pattern.
#
#   mkdir -p logs && sbatch scripts/patterns.slurm
#   # when the array finishes:
#   for t in big2_uip big2_nsip; do python3 -c "import sys; sys.path.insert(0,'src'); from cohort_io import pack; pack('results/ph_cohorts','$t')"; done
#   python3 src/ph_cohort.py analyze --tag big2_uip; python3 src/ph_cohort.py analyze --tag big2_nsip
#
# 2 cohorts x 20 shards of 10 cases, about 5-8 min each on one core. Add --PB 585 to ARGS for Mexico City.
#SBATCH --job-name=lupnet-pattern
#SBATCH --array=0-39%40
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/pattern_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# source .venv/bin/activate

PATS=(uip nsip)
c=$(( SLURM_ARRAY_TASK_ID / 20 )); s=$(( SLURM_ARRAY_TASK_ID % 20 ))
i0=$(( s * 10 )); i1=$(( i0 + 10 ))
ARGS="--n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85"
echo "big2_${PATS[$c]} cases $i0-$((i1 - 1)) on $(hostname)"
# shellcheck disable=SC2086
python3 src/ph_cohort.py run --tag "big2_${PATS[$c]}" $ARGS --pattern "${PATS[$c]}" --i0 "$i0" --i1 "$i1"
