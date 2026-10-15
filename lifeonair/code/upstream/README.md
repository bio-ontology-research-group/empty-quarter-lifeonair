# Upstream IBEX scripts (copies)

Copies of the cluster scripts that produced the frozen inputs in `lifeonair/results/inputs/` and
`lifeonair/results/reads/`. They run on IBEX against the workbench
(`/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/`) and are kept here so the full derivation is
readable from a clean clone. Order:

1. `lifeonair_1_census.sh`, `lifeonair_2_type_guild.sh` (workbench `lifeonair_clean/`): extract K06281/K03520
   large subunits from the 5,438 representatives, type them against HydDB and the CoxL references.
   The CoxL motif rule in `lifeonair_2_type_guild.sh` is inverted; `../fix_coxL_v2.py` corrects it (R1-01).
2. `../r1_build_db.sh`, `../r2_reads_array.sh` (workbench `lifeonair_reads/`): build the competitive DIAMOND
   database and screen 10 M reads per sample (Slurm array 53524441); `../r3_aggregate.py` aggregates locally.
3. `../g1_growth_rate.sh` with `run_grodon.R`: gRodon2 doubling times; `../s1_sub_reconcile.sh` re-estimates
   the `_sub` representatives (53614016).
4. `genome_size.sh` (workbench `lifeonair_review_r1/`, job 53611804): completeness-corrected genome size.
5. `../b1_classifier_benchmark.sh` (53613713), `../b2_disjoint_benchmark.sh`, `../b3_mock_ge.sh`: benchmarks.
`lifeonair_3_abundance.py` and `lifeonair_4_energetics.py` belong to the superseded catalogue census.
