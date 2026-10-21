#!/usr/bin/env bash
# Generate the four sensitivity cohorts across several cores, then build Table S16 and
# Figure S2. Run from the root of the repository.
#
#   bash scripts/sensitivity_parallel.sh                  # uses every core it can see
#   JOBS=16 SHARDS=25 bash scripts/sensitivity_parallel.sh
#   N=60 bash scripts/sensitivity_parallel.sh             # a smaller, faster pilot
#
# One simulated case takes about 25 s of one core, so the four cohorts of 200 cases are
# roughly 5.5 core-hours: about 20 min on 16 cores. Every task is independent and every
# case is written atomically, so the run is resumable: re-run the script and it skips
# the cases that are already on disk.
set -euo pipefail

N=${N:-200}                                   # cases per cohort
SHARDS=${SHARDS:-20}                          # tasks per cohort
JOBS=${JOBS:-$( (command -v nproc >/dev/null && nproc) || echo 4 )}
REF=${REF:-big2}                              # the calibrated cohort, already in the repo

# One thread per worker: numpy's own threads would fight the workers for the same cores.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

COMMON="--n $N --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85"
# tag<TAB>extra arguments. Same seeds everywhere, so a difference between cohorts comes
# from the parameter and not from a different draw of cases.
SPECS=$'tau1.0\t--tau 1.0\ntau3.0\t--tau 3.0\npv0.50\t--p_vasc 0.50\npv0.85\t--p_vasc 0.85'

TASKS=$(mktemp)
trap 'rm -f "$TASKS"' EXIT
while IFS=$'\t' read -r tag extra; do
    for s in $(seq 0 $((SHARDS - 1))); do
        i0=$(( s * N / SHARDS )); i1=$(( (s + 1) * N / SHARDS ))
        [ "$i0" -lt "$i1" ] && printf 'python3 src/ph_cohort.py run --tag %s %s %s --i0 %d --i1 %d\n' \
            "$tag" "$COMMON" "$extra" "$i0" "$i1" >> "$TASKS"
    done
done <<< "$SPECS"

echo "$(wc -l < "$TASKS") tasks on $JOBS workers, $N cases per cohort"
xargs -P "$JOBS" -I CMD bash -c 'CMD' < "$TASKS"

for tag in tau1.0 tau3.0 pv0.50 pv0.85; do
    have=$(ls results/ph_cohorts/${tag}_*.tsv 2>/dev/null | wc -l)
    [ "$have" -eq "$N" ] || { echo "cohort $tag incomplete: $have of $N cases; re-run to resume" >&2; exit 1; }
done

python3 src/sensitivity.py --ref "$REF" --tags tau1.0 tau3.0 pv0.50 pv0.85 \
    --labels 'T_les 1 y' 'T_les 3 y' '70% to 50%' '70% to 85%'
