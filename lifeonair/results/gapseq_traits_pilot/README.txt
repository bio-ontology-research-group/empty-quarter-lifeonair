Gapseq-first small predictive pilot

Read REPORT.txt for findings and the decision. analysis_plan.txt specifies the
panel, two metabolic predictors, fixed baseline/ridge penalty, blocked validation
and scale-up criteria before phenotype/association inspection. Gapsmith outputs
are not used. The previous inventory is reused only for current-ID metadata.

Full workbench:
/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/gapseq_traits_pilot_clean
R = /ibex/scratch/projects/c2014/rob/metagenomics
W = R/06_atmospheric/gapseq_traits_pilot_clean
Python = /ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
(COBRApy 0.31.1, GLPK; numpy for predictions). CPU analyses run through Slurm.

Input files in W/inputs: model_inventory.tsv, samples.tsv, rep_growth_v2.tsv,
coverm_rep_abundance_long.tsv.gz, guild_genome_equivalents_v2_r3noflagged.tsv.
The first contains current taxonomy, MAG quality, contamination flags and marker
groups; legacy gapsmith availability fields are ignored. Copy code/gapseq_pilot*
and code/metabolic_simulate.py plus metabolic_selftest.py into W/scripts.

Execution order, using absolute paths and the environment Python:
  python W/scripts/gapseq_pilot_inventory.py R
  sbatch W/scripts/gapseq_pilot.sbatch
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python W/scripts/gapseq_pilot_field.py R
Submit the inventory and field commands as CPU Slurm jobs, and wait for each
previous stage to finish successfully. The supplied simulation array has capacity
for 60 MAGs; tasks beyond the selected panel exit without model simulation.
No source model is modified. job_history.tsv records original run identifiers.
Input/model/script hashes and per-model runtime records accompany the outputs.

Retrieve compact results to this directory. Local checks and figure regeneration:
  python3 lifeonair/code/gapseq_pilot_verify.py
  python3 lifeonair/code/gapseq_pilot_plot.py
These commands run from the repository root. Plotting needs numpy/matplotlib.
Verification recomputes metrics from predictions, not metabolic fluxes or fits.
Run sha256sum -c checksums.sha256 from this directory to check frozen artifacts.
Raw per-MAG shards remain in the IBEX workbench and are ignored here.

Genus validation tests held-out genera; order can be shared across its train/test
sets. Site validation tests geographic repeatability for represented genomes,
not prediction for unseen genomes. Site and genus holdouts are separate tests,
not a joint double holdout. All imputation/scaling/encoding is trained inside
each fold. Site-fold errors give every MAG equal total weight. The within-order
shuffle is a screening reference and is not reported as a formal p-value.

ATP production is a feasible upper bound at specified resource caps, not growth,
measured substrate use, regulation or in situ flux. Root-adjacent and bulk
abundances are restricted to study sites 1-60 and are relative CoverM percentages; no panel renormalisation is used.
Both-zero visits are excluded from a MAG's response. Missing CoverM samples are
excluded rather than interpreted as nondetection. This pilot uses an existing
exploratory field dataset and is not independent external biological validation.
