#!/bin/bash
#SBATCH --job-name=loa_model_audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/metabolic_pilot_clean/logs/inventory_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/metabolic_pilot_clean/logs/inventory_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
W=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/metabolic_pilot_clean
M=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/gapsmith_models
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" "$W/scripts/metabolic_inventory.py" --inputs "$W/inputs" --models "$M" --out "$W/results"
"$PY" "$W/scripts/metabolic_inspect.py" "$W" "$M"
