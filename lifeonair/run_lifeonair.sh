#!/bin/bash
# Life on air (P3): regenerate every cited table, figure and the manuscript from the frozen inputs.
# Run from the repository root after `bash reproducibility/fetch_large_files.sh` (per-sample read counts).
# Cluster stages (typing, read screening, gRodon, benchmarks) are in lifeonair/code/upstream/README.md.
#   bash lifeonair/run_lifeonair.sh
set -euo pipefail
if [[ "${1:-}" == "--verify" ]]; then
  python3 lifeonair/reproducibility/verify_p3.py
  exit
fi
# Statistical bootstrap stages should run in a CPU allocation (IBEX or unimatrix01).
R="micromamba run -n eq-lifeonair"
L=lifeonair/code
if [[ ! -d lifeonair/results/reads/per_sample ]]; then
  tar -xzf lifeonair/results/reads/per_sample_counts.tar.gz -C lifeonair/results/reads
fi
python3 lifeonair/reproducibility/verify_p3.py
python3 $L/r3_aggregate.py lifeonair/results/reads lifeonair/results/reads/orphan_sub_reps.txt guild_genome_equivalents_v2.tsv db/marker_labels_v2.tsv
python3 $L/r3_aggregate.py lifeonair/results/reads lifeonair/results/reads/excl_v2_noflagged.txt guild_genome_equivalents_v2_noflagged.tsv db/marker_labels_v2.tsv
python3 $L/r3_aggregate.py lifeonair/results/reads lifeonair/results/review_r3/mag_exclusion_union.txt guild_genome_equivalents_v2_r3noflagged.tsv db/marker_labels_v2.tsv
python3 $L/r3_aggregate.py lifeonair/results/reads lifeonair/results/reads/excl_v2_strict1l.txt guild_genome_equivalents_v2_strict1l.tsv db/marker_labels_v2.tsv
python3 $L/r3_aggregate.py lifeonair/results/reads lifeonair/results/reads/orphan_sub_reps.txt guild_genome_equivalents_v2_motifonly.tsv db/marker_labels_v2_motifonly.tsv
python3 $L/coxl_deconv.py lifeonair/results/benchmark/classifier_benchmark_confusion.tsv clu50
python3 $L/apply_sub_fix.py
python3 $L/operon_summary.py lifeonair/results/operon lifeonair/results/inputs/coxL_class_v2.tsv
# env_site/env_site_trip are frozen inputs. Their upstream derivation additionally
# requires the checksum-pinned amplicon cache; see reproducibility/README.md.
if [[ "${P3_REBUILD_ENV:-0}" == "1" ]]; then
  $R Rscript $L/prep_environment.R
fi
$R Rscript $L/h123_v2.R
$R Rscript $L/h1_mechanism_v2.R
$R Rscript $L/a5_genome_traits.R
$R Rscript $L/r3_mag_flag_sensitivity.R
python3 $L/r3_evidence_checks.py
$R Rscript $L/lineages_v2.R
$R Rscript $L/fig_v2.R
python3 $L/make_tables.py
python3 lifeonair/reproducibility/verify_p3.py
(cd lifeonair/manuscript && latexmk -pdf -interaction=nonstopmode main.tex >/dev/null)
echo "life on air: done"
