#!/bin/bash
#SBATCH --job-name=gapseq_pilot
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --array=0-14%8
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/gapseq_traits_pilot_clean/logs/sim_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/gapseq_traits_pilot_clean/logs/sim_%A_%a.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
R=/ibex/scratch/projects/c2014/rob/metagenomics
W=$R/06_atmospheric/gapseq_traits_pilot_clean
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" "$W/scripts/metabolic_selftest.py"
for offset in 0 1 2 3; do
 ix=$((SLURM_ARRAY_TASK_ID * 4 + offset))
 "$PY" "$W/scripts/gapseq_pilot_simulate.py" "$R" "$ix"
done
