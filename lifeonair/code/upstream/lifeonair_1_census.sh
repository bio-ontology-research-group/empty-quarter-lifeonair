#!/bin/bash
# lifeonair step 1 -- CRUDE KO CENSUS over the current dRep representatives.
# Counts representatives carrying the trace-gas marker KOs in the per-rep eggNOG KO table:
#   hydB  = K06281  (H2 uptake hydrogenase, large subunit)
#   coxL  = K03520  (form-I CO dehydrogenase, large subunit)
# This is a FIRST PASS only: KO presence over-counts, because it does not distinguish
# high-affinity (atmospheric) enzyme forms. Step 2 types them phylogenetically and refines
# the guild down to bona-fide high-affinity carriers. Markers are taken from
# scripts/82_metabolism_marker_census.sh (not guessed).
#
#   bash scripts/lifeonair_1_census.sh
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
A=$M/06_atmospheric/lifeonair_clean
KO=$M/04_functional/repko_clean/results/rep_ko_v2.tsv          # rep_id <TAB> comma-separated KOs
R=$A/results; mkdir -p "$R"
[ -s "$KO" ] || { echo "ERROR: missing per-rep KO table $KO"; exit 2; }

awk -F'\t' -v out="$R" '
  { h = ($2 ~ /(^|,)K06281(,|$)/); c = ($2 ~ /(^|,)K03520(,|$)/); tot++
    if (h) { nH++; print $1 > (out "/hydB_KO_reps.txt") }
    if (c) { nC++; print $1 > (out "/coxL_formI_KO_reps.txt") }
    if (h && c) nB++
    if (h || c) nA++ }
  END {
    printf "== crude-KO census (pre-typing) ==\n"
    printf " total reps: %d\n", tot
    printf " hydB K06281 (H2): %d (%.1f%%)\n", nH, 100*nH/tot
    printf " coxL form-I K03520 (CO): %d (%.1f%%)\n", nC, 100*nC/tot
    printf " both: %d (%.1f%%)\n", nB, 100*nB/tot
    printf " any (H2 or CO): %d (%.1f%%)\n", nA, 100*nA/tot
  }' "$KO"

echo "naming sanity (expect __single_/__self_ = current merged catalogue, NOT plain __SemiBin):"
head -1 "$R/hydB_KO_reps.txt"; head -1 "$R/coxL_formI_KO_reps.txt"
echo "wrote: $R/hydB_KO_reps.txt  $R/coxL_formI_KO_reps.txt"
