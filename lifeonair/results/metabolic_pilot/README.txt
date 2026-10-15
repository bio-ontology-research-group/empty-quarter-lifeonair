Metabolic feasibility pilot for lifeonair, 12 October 2026

Start with REPORT.txt and pilot_summary.pdf. These are exploratory QC outputs,
not new manuscript claims. No field association model was fitted: the planned
prerequisite of defensible metabolic phenotypes was not met.

Reproduce compact summaries from this repository (Python standard library):
  python3 lifeonair/code/metabolic_summarise.py --from-frozen
  python3 lifeonair/code/metabolic_compare.py
  python3 lifeonair/code/metabolic_plot.py
The plotting command also requires numpy and matplotlib (used 3.10.1).
From this directory, sha256sum -c checksums.sha256 verifies frozen inputs/results
and source scripts. Figure files are frozen artifacts; PDF metadata can change
when regenerated. These checks recompute summaries, not metabolic fluxes.

Full simulation workspace on IBEX:
  /ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/metabolic_pilot_clean
Its scripts/, inputs/, results/ and logs/ contain the run. Define R as
/ibex/scratch/projects/c2014/rob/metagenomics and W as the workspace above.
Use the metabolism-env Python under
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
(COBRApy 0.31.1, Python 3.10.20, GLPK; per-model runtime_manifest.json).

Copy metabolic_*.py and metabolic_*.sbatch from ../../code/ to W/scripts/.
Preserve and verify W/inputs/ against upstream_inputs.sha256 before rerunning.
The inputs are the current representative quality table, GTDB R232 taxonomy,
curated hydrogenase classification, CoxL guild list and R3 exclusion-ID union.
Source reconstructions/proteomes remain under R/06_darkmatter/; do not alter them.
The SHA256 ledger input_hashes.tsv pins both gapsmith XMLs, added-reaction ledger
and Prodigal proteome for every pilot MAG. gapseq_input_hashes.tsv pins secondary
XMLs/genomes, the source-genome manifest and curated marker FASTA. Its hashes
were collected after the read-only run; filled-model hashes can additionally be
cross-checked against those captured during execution in gapseq_identity.tsv.

Submit CPU jobs (account pi-hohndor, using a valid partition/time allocation):
  sbatch W/scripts/metabolic_inventory.sbatch
  sbatch W/scripts/metabolic_simulate.sbatch
Replace W with the absolute workspace path. The inventory must finish first.
The array runs 18 tasks with four MAGs per task, 1 CPU/task, 4 GB, and also runs
metabolic_selftest.py. The final array job was 53631046; the earlier partial
runs are superseded. job_history.tsv and slurm_accounting.tsv record completion.

The secondary CPU job runs:
  python W/scripts/metabolic_gapseq_check.py R
followed after completion by:
  python W/scripts/metabolic_gapseq_map.py R
and the read-only provenance helper:
  python W/scripts/metabolic_provenance.py R
Substitute the environment Python and absolute R/W paths. Original Slurm jobs
53631047 and 53631420 completed successfully. The first requires COBRApy; the
mapper imports the shared simulation helpers. Models are not reconstructed.

Retrieve W/results/simulations/ to this directory, then run
  python3 lifeonair/code/metabolic_summarise.py
from the repository root. It checks 72 model runtimes against the simulation
script SHA256, closes 24 three-group blocks, aggregates 216 audits and 10,584
scenario rows, and checks 3,240 H2 monotonicity comparisons. Models with no
internal H2 reaction legitimately have no per-model hydrogen_reactions shard.
The compact aggregate tables and runtime manifest are tracked; duplicate shards
remain on IBEX and are ignored locally. The secondary comparison checks 25 exact
sequence pairs, 750 unique scenarios and closed/inorganic ATP controls.

Limitations: exact marker matching can miss alternative predicted starts; formula
coverage is insufficient for a complete mass/charge audit; ATP feasibility is
not growth, flux in situ, enzyme affinity or measured substrate use. The pilot
is order-matched but not a random catalogue sample. gapseq coverage is restricted
to the available single-sample subset. Model media and biomass formulations
were inherited; no conclusion relies on a negative growth prediction alone.

Reviews and validation logs are under ../../reviews/metabolic_pilot/. The
364-check lifeonair verification passes. The full repository harness retains
30 pre-existing P4 failures; see its log. Full external reproduction also needs
release of the upstream model/genome/proteome bundle; this directory alone is
not a substitute for that deposition.
