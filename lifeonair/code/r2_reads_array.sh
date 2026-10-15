#!/bin/bash
# Life on air (P3), read-based guild abundance -- stage R2: one array task per paired-read sample.
# Subsamples 10 M R1 reads (seqtk, fixed seed), truncates every read to 100 nt so that the
# 101-nt (trip 4) and 150-nt libraries are screened identically, DIAMOND blastx against the
# stage-R1 competitive database, keeps each read's single best hit, and writes per-subject read
# counts at three identity floors (50/60/70%) with alignment length >= 25 aa.
#
#   sbatch --array=1-573%60 r2_reads_array.sh
#
#SBATCH --job-name=loa_reads
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/reads_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/reads_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
A=$R/06_atmospheric/lifeonair_reads
DB=$A/db
MAN=$R/rebuild_20260729/picrust2_current_20260803/coverm_current/prepared/sample_input_manifest.tsv
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
DIAMOND=$CONDA/eggnog/bin/diamond
SEQTK=$CONDA/eco_probes/bin/seqtk
T=${SLURM_CPUS_PER_TASK:-16}
NSUB=10000000
SEED=11

[ -s "$DB/DB_DONE" ] || { echo "ERROR: database not built"; exit 1; }
LINE=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR == i + 1' "$MAN")
[ -n "$LINE" ] || { echo "no manifest row for task $SLURM_ARRAY_TASK_ID"; exit 0; }
S=$(cut -f1 <<< "$LINE"); R1=$(cut -f6 <<< "$LINE")
OUT=$A/per_sample_L100/$S
mkdir -p "$OUT"
[ -s "$OUT/DONE" ] && { echo "$S already done"; exit 0; }

TMP=${SLURM_TMPDIR:-/tmp}/loa_${SLURM_JOB_ID}_${S}
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

echo "[$(date -Iseconds)] $S subsample $NSUB reads (seed $SEED) from $R1"
"$SEQTK" sample -s "$SEED" "$R1" "$NSUB" | "$SEQTK" trimfq -L 100 - | "$SEQTK" seq -A -L 75 - > "$TMP/r1.fa"
NREADS=$(grep -c '^>' "$TMP/r1.fa")

echo "[$(date -Iseconds)] $S blastx ($NREADS reads)"
"$DIAMOND" blastx -q "$TMP/r1.fa" -d "$DB/loa_reads" -p "$T" --id 50 --evalue 1e-5 \
    --max-target-seqs 1 --block-size 4 --index-chunks 1 --tmpdir "$TMP" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore --quiet -o "$TMP/hits.tsv"

# one best hit per read (DIAMOND may report several HSPs for the top subject)
awk -F'\t' '$4 >= 25 && !seen[$1]++ {
        c50[$2]++; if ($3 >= 60) c60[$2]++; if ($3 >= 70) c70[$2]++ }
     END { for (s in c50) printf "%s\t%d\t%d\t%d\n", s, c50[s], c60[s] + 0, c70[s] + 0 }' \
    "$TMP/hits.tsv" | sort > "$OUT/subject_counts.tsv"
printf "sample\treads_screened\tR1\tseed\ttruncate_nt\tmin_aln_aa\n%s\t%s\t%s\t%s\t100\t25\n" "$S" "$NREADS" "$R1" "$SEED" > "$OUT/meta.tsv"
date -Iseconds > "$OUT/DONE"
echo "[$(date -Iseconds)] $S done: $(wc -l < "$OUT/subject_counts.tsv") subjects hit"
