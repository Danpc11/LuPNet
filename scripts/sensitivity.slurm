#!/usr/bin/env bash
# SLURM array that generates the four sensitivity cohorts, one task per shard.
#
#   sbatch scripts/sensitivity.slurm
#   # when the array finishes:
#   python3 src/sensitivity.py --ref big2 --tags tau1.0 tau3.0 pv0.50 pv0.85 \
#       --labels 'T_les 1 y' 'T_les 3 y' '70% to 50%' '70% to 85%'
#
# 4 cohorts x 20 shards = 80 tasks of 10 cases, about 5 min each. Adjust the partition,
# the account and the module or environment lines to the cluster before submitting.
#SBATCH --job-name=lupnet-sens
#SBATCH --array=0-79%40
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:40:00
#SBATCH --output=logs/sens_%A_%a.out
set -euo pipefail

N=200
SHARDS=20
COMMON="--n $N --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85"
TAGS=(tau1.0 tau3.0 pv0.50 pv0.85)
EXTRA=("--tau 1.0" "--tau 3.0" "--p_vasc 0.50" "--p_vasc 0.85")

# One thread per task: the array already uses the cores.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# source ~/venvs/lupnet/bin/activate      # or: module load python/3.12

c=$(( SLURM_ARRAY_TASK_ID / SHARDS ))
s=$(( SLURM_ARRAY_TASK_ID % SHARDS ))
i0=$(( s * N / SHARDS )); i1=$(( (s + 1) * N / SHARDS ))

echo "cohort ${TAGS[$c]} cases $i0-$((i1 - 1)) on $(hostname)"
# shellcheck disable=SC2086
python3 src/ph_cohort.py run --tag "${TAGS[$c]}" $COMMON ${EXTRA[$c]} --i0 "$i0" --i1 "$i1"
