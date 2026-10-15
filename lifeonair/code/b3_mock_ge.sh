#!/bin/bash
# Life on air (P3): mock-community calibration of the genome-equivalent (GE) estimator (review round 2, R2-05).
#
# Communities are built from the representative genomes with the CoverM abundance profiles of real bulk
# (Deep) and root-adjacent (PlantRoot) samples from the same site visits (b3_mock_sim.py select). Each array
# task simulates 4 M reads for one community (b3_mock_sim.py simulate; truth = cell-weighted marker copies
# per genome) and screens them exactly as r2_reads_array.sh screens real samples (100-nt reads, DIAMOND blastx
# --id 50 --evalue 1e-5 --max-target-seqs 1, best hit with alignment >= 25 aa, counts at 50/60/70%), against
# (a) the production database and (b) the database without any EQ marker or SCG protein of the community's
# own genomes. The aggregation step runs r3_aggregate.py (v2 labels) on both and joins the truth.
#
#   sbatch b3_mock_ge.sh select                          # choose visits (no array)
#   sbatch --array=1-20 b3_mock_ge.sh                    # one task per community
#   sbatch --dependency=afterok:<array> b3_mock_ge.sh aggregate
#
#SBATCH --job-name=loa_mock
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r2/logs/mock_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r2/logs/mock_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
W=$R/06_atmospheric/lifeonair_review_r2
M=$W/mock
DB=$R/06_atmospheric/lifeonair_reads/db
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
DIAMOND=$CONDA/eggnog/bin/diamond
T=${SLURM_CPUS_PER_TASK:-8}
NREADS=4000000
source "$CONDA/../bin/activate" bakta   # any env with python3
mkdir -p "$M"
PY="python3 $W/scripts/b3_mock_sim.py"

screen() {   # reads.fa, diamond db, out dir, sample name
    local TMP=${SLURM_TMPDIR:-/tmp}/mock_${SLURM_JOB_ID}_$4; mkdir -p "$TMP" "$3"
    "$DIAMOND" blastx -q "$1" -d "$2" -p "$T" --id 50 --evalue 1e-5 --max-target-seqs 1 --block-size 4 \
        --index-chunks 1 --tmpdir "$TMP" --outfmt 6 qseqid sseqid pident length evalue bitscore --quiet -o "$TMP/hits.tsv"
    awk -F'\t' '$4 >= 25 && !seen[$1]++ {
            c50[$2]++; if ($3 >= 60) c60[$2]++; if ($3 >= 70) c70[$2]++ }
         END { for (s in c50) printf "%s\t%d\t%d\t%d\n", s, c50[s], c60[s] + 0, c70[s] + 0 }' \
        "$TMP/hits.tsv" | sort > "$3/subject_counts.tsv"
    printf "sample\treads_screened\tR1\tseed\ttruncate_nt\tmin_aln_aa\n%s\t%s\tsimulated\t%s\t100\t25\n" \
        "$4" "$(grep -c '^>' "$1")" "$SLURM_ARRAY_TASK_ID" > "$3/meta.tsv"
    rm -rf "$TMP"
}

if [ "${1:-}" = select ]; then
    $PY select "$M" 10 11
    cat "$M/communities.tsv"
    exit 0
fi

if [ "${1:-}" = aggregate ]; then
    for V in full novel; do
        A=$M/agg_$V; mkdir -p "$A/db"
        ln -sf "$W/inputs/marker_labels_v2.tsv" "$A/db/marker_labels_v2.tsv"
        ln -sf "$DB/subject_lengths.tsv" "$A/db/subject_lengths.tsv"
        python3 "$W/scripts/r3_aggregate.py" "$A" - mock_ge.tsv db/marker_labels_v2.tsv
    done
    $PY compare "$M"
    cp "$M/mock_ge_vs_truth.tsv" "$M/communities.tsv" "$W/results/"
    sha256sum "$W/results/mock_ge_vs_truth.tsv" "$W/results/communities.tsv" > "$W/results/mock_ge.sha256"
    echo "[$(date -Iseconds)] done"
    exit 0
fi

C=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR == i + 1 {print $1}' "$M/communities.tsv")
[ -n "$C" ] || { echo "no community for task $SLURM_ARRAY_TASK_ID"; exit 0; }
echo "[$(date -Iseconds)] $C simulate"
$PY simulate "$M" "$C" "$NREADS" "$SLURM_ARRAY_TASK_ID"
echo "[$(date -Iseconds)] $C screen (production DB)"
screen "$M/comm/$C.reads.fa" "$DB/loa_reads" "$M/agg_full/per_sample/$C" "$C"
echo "[$(date -Iseconds)] $C screen (own genomes withheld)"
$PY novel_db "$M" "$C"
"$DIAMOND" makedb --in "$M/comm/$C.novel.faa" -d "$M/comm/$C.novel" -p "$T" --quiet
screen "$M/comm/$C.reads.fa" "$M/comm/$C.novel" "$M/agg_novel/per_sample/$C" "$C"
rm -f "$M/comm/$C.reads.fa" "$M/comm/$C.novel.faa" "$M/comm/$C.novel.dmnd"
echo "[$(date -Iseconds)] $C done"
