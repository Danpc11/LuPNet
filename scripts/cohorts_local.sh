#!/usr/bin/env bash
# Generate the pattern cohorts (UIP, fibrotic NSIP) and the Mexico City altitude cohorts on the cores of
# one machine, without a scheduler. Run from the root of the repository.
#
#   bash scripts/cohorts_local.sh                       # everything, on every core it can see
#   JOBS=96 bash scripts/cohorts_local.sh               # leave some cores free
#   WHAT=patterns bash scripts/cohorts_local.sh         # only big2_uip and big2_nsip
#   WHAT=altitude bash scripts/cohorts_local.sh         # only big2_cdmx, adv2_cdmx and altF_cdmx
#   WHAT=advanced bash scripts/cohorts_local.sh         # only adv2_uip and adv2_nsip (advanced disease by pattern)
#   PB=585 WHAT=patterns bash scripts/cohorts_local.sh  # pattern cohorts at Mexico City altitude
#   LIMIT=4 bash scripts/cohorts_local.sh               # pilot: the first 4 cases of each cohort
#   nohup bash scripts/cohorts_local.sh > logs/cohorts.log 2>&1 &     # detached; follow with tail -f
#
# Same seeds as the manuscript cohorts, so a difference between cohorts comes from the pattern or the
# altitude and not from a different draw of cases. One case takes 25 to 60 s of one core: the 760 cases
# are about 8 core-hours, some 6 min on 96 cores. Every case is written atomically and finished cases
# are skipped, so the run is resumable: run the script again.
set -euo pipefail

JOBS=${JOBS:-$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4 )}
WHAT=${WHAT:-all}
PB=${PB:-760}                                 # barometric pressure of the pattern cohorts
PY=${PY:-python3}
LIMIT=${LIMIT:-0}                             # >0: only the first LIMIT cases of each cohort (pilot)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p logs results/ph_cohorts

BIG="--n 200 --seed 5000 --pseed 20000 --lo 0.25 --hi 0.85"
ADV="--n 60 --seed 6000 --pseed 30000 --lo 0.12 --hi 0.45"
ALT="--n 100 --seed 7000 --pseed 40000 --vmode fibrosis"
SUF=""; [ "$PB" != "760" ] && SUF="_pb$PB"
# tag<TAB>cases<TAB>arguments
SPECS=""
if [ "$WHAT" = all ] || [ "$WHAT" = patterns ]; then
    SPECS+="big2_uip$SUF"$'\t'"200"$'\t'"$BIG --pattern uip --PB $PB"$'\n'
    SPECS+="big2_nsip$SUF"$'\t'"200"$'\t'"$BIG --pattern nsip --PB $PB"$'\n'
fi
if [ "$WHAT" = all ] || [ "$WHAT" = altitude ]; then
    SPECS+="big2_cdmx"$'\t'"200"$'\t'"$BIG --PB 585"$'\n'
    SPECS+="adv2_cdmx"$'\t'"60"$'\t'"$ADV --PB 585"$'\n'
    SPECS+="altF_cdmx"$'\t'"100"$'\t'"$ALT --PB 585"$'\n'
fi
if [ "$WHAT" = all ] || [ "$WHAT" = advanced ]; then
    SPECS+="adv2_uip$SUF"$'\t'"60"$'\t'"$ADV --pattern uip --PB $PB"$'\n'
    SPECS+="adv2_nsip$SUF"$'\t'"60"$'\t'"$ADV --pattern nsip --PB $PB"$'\n'
fi
[ -n "$SPECS" ] || { echo "WHAT must be all, patterns, altitude or advanced" >&2; exit 1; }

TASKS=$(mktemp); trap 'rm -f "$TASKS"' EXIT
while IFS=$'\t' read -r tag n args; do
    [ -n "$tag" ] || continue
    [ "$LIMIT" -gt 0 ] && [ "$LIMIT" -lt "$n" ] && n=$LIMIT
    for i in $(seq 0 $((n - 1))); do          # one case per task: the workers stay balanced
        printf '%s src/ph_cohort.py run --tag %s %s --i0 %d --i1 %d\n' "$PY" "$tag" "$args" "$i" "$((i + 1))" >> "$TASKS"
    done
done <<< "$SPECS"

echo "$(wc -l < "$TASKS") cases on $JOBS workers"
xargs -P "$JOBS" -I CMD bash -c 'CMD' < "$TASKS" > logs/cohorts_local.out 2>&1 || true

status=0
while IFS=$'\t' read -r tag n args; do
    [ -n "$tag" ] || continue
    [ "$LIMIT" -gt 0 ] && [ "$LIMIT" -lt "$n" ] && n=$LIMIT
    have=$(ls results/ph_cohorts/${tag}_[0-9]*.tsv 2>/dev/null | wc -l | tr -d ' ')
    if [ "$have" -eq "$n" ]; then
        "$PY" -c "import sys; sys.path.insert(0,'src'); from cohort_io import pack; pack('results/ph_cohorts','$tag')"
        echo "== $tag: $have cases, packed to results/ph_cohorts/$tag.tsv.gz"
        [ "$LIMIT" -gt 0 ] || "$PY" src/ph_cohort.py analyze --tag "$tag" || true
    else
        echo "== $tag incomplete: $have of $n cases; run the script again to resume (see logs/cohorts_local.out)" >&2; status=1
    fi
done <<< "$SPECS"
exit $status
