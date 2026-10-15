# Consultation: adjacent and alternative hypotheses for the "Life on air" (P3) paper after round-1 corrections

You reviewed this manuscript in round 1 (your review: `REPO/lifeonair/reviews/round1/review.md` and
`findings.tsv`; REPO = /home/leechuck/Documents/papers/empty-quarter-mags). Work read-only as before (you may
run read-only R/Python one-liners; R env `micromamba run -n eq-lifeonair`). IBEX is mounted read-only at
`/home/leechuck/mnt/ibex-metagenomics-ro/`.

## What changed since your review (verified by the authors)

1. CoxL motif bug confirmed and fixed (`lifeonair/code/fix_coxL_v2.py`): form I = AYxCSFR (all three Cordero
   references carry it), AYRGAGR = form II; no-motif proteins form I only if >55% id / >=90% cov to Cordero refs.
   285 form-I proteins in 246 reps (was 1,901 / 1,170). Group 2a now its own column; two orphan `_sub` reps'
   proteins excluded. New table: `lifeonair/results/reads/guild_genome_equivalents_v2.tsv` (id_floor 60 is primary).
2. Corrected headline numbers (sites 1-60 only, 546 samples; visit medians, log2, pseudocount 0.01):
   - bulk medians GE: 1l 1.08, 1h 0.14, 2a 0.03, form-I CoxL 0.105, other NiFe (no 2a) 0.08, other CoxL 2.29
   - root vs bulk (106 visits / 47 sites): 1l -2.27; 1h -0.33; 2a -0.34; form I -0.36; controls -0.29;
     excess beyond controls: 1l -1.83 (p_site<1e-4), 1h -0.01, form I -0.05
   - surface vs bulk (68 visits / 36 sites): 1l -1.05, 1h -0.63, 2a -0.37, form I -0.51, controls -0.16;
     excess 1l -0.83, 1h -0.42, form I -0.24 (p 0.008)
3. Biomass/DNA yield (`lifeonair/results/inputs/sample_ledger_dna.tsv`, extract concentration; available for
   267 usable samples, trips 1-3 and few of trip 5, none of trip 4): root-adjacent extracts ~9x bulk (log2 +3.2);
   prokaryotic SCG reads per million are similar across positions (plant DNA does not explain it). DNA-scaled
   per-gram proxy, root vs bulk: 1l +0.4 (n.s.), 1h +1.5, 2a +2.3, form I +2.5.
   Exploratory mixed models (log2 GE ~ log2 DNA + position + trip + (1|site)): slope for 1l -0.32 (SE 0.03),
   within bulk only -0.23 (SE 0.05); all other markers -0.02 to -0.05. Across bulk sites, log DNA rises with
   MAP (rho 0.41, p 0.01) and site-median 1l falls with log DNA (rho -0.45, p 0.005).
4. The PI decided: H1 is reframed as relative exclusion (1l carriers do not take part in the larger root-zone
   community) rather than loss of capacity; H2 is narrowed (the root contrast is 1l-specific; CO capacity is no
   longer claimed to be position-independent). The PI asked for ADJACENT OR ALTERNATIVE HYPOTHESES that these
   data can test and that would be of ISME Journal quality.

## Candidate hypothesis from the authors (evaluate critically)

"Biomass-decoupled scavengers": group 1l carriers' per-gram abundance scales sub-proportionally with total
microbial biomass (1l GE ~ biomass^-0.3, other markers ~ biomass^0), so their share falls wherever productivity
is higher: beside roots, at wetter sites, at greener sites. If true, a single mechanism (an atmospheric-energy-
limited population that does not track organic-carbon-driven growth) would unify the root contrast, the aridity
gradient and the NDVI association, and the aridity association would be partly a biomass-dilution effect rather
than selection by aridity. Weigh measurement error in DNA yield (regression dilution, kit/trip effects,
relic DNA), circularity (GE denominators), and what would falsify it.

## Data available beyond the P3 folder (read-only)

- Amplicon companion repo `/home/leechuck/Documents/papers/empty-quarter-amplicon/repository/cache/`:
  `metadata.tsv`, `xrf_site_compartment_panel.tsv`, `xrf_per_compartment.tsv` (XRF elements, per site x
  position), `smap_by_sample.parquet`, `ndvi_by_sample.parquet`, `topography.tsv`, PMA/relic analyses
  (`REPO/results/relic_climate/`, incl. `priorfree/arm_b_measured_pma.tsv`), 16S feature tables.
- MAG catalogue: GTDB R232 taxonomy, CoverM abundances, Bakta/eggNOG annotation on IBEX (`04_*`, `03_annotation`),
  gRodon growth (`lifeonair/results/inputs/rep_growth.tsv`), operon context, hydB/coxL trees.
- IBEX has other analyses that may be relevant (`08_dormancy`, `09_kraken2`, `10_thermal`, `05_metabolism`,
  `06_atmospheric/carbon_fixation`, `06_atmospheric/autotrophy_fba`); read their README.md files.
- Literature: `REPO/lifeonair/literature/` (MANIFEST.tsv; digest `search/round1/DIGEST.md`).

## Task

Propose 6-10 adjacent or alternative hypotheses (including refinements of H1-H3 and the biomass hypothesis
above) that (a) are testable with the existing data or with modest new computation on existing data, (b) would
be of ISME Journal quality and novel relative to the literature (cite the specific papers they extend or
contradict), and (c) do not overclaim capacity as activity. Examples of directions to consider, not to limit you:
genomic traits distinguishing 1l from 1h carriers (genome size, RuBisCO/carbon fixation, CO co-occurrence,
growth potential, dormancy/sporulation genes); the biomass-dilution hypothesis; why the surface depletes 1l,
1h and form I together (UV, temperature, desiccation); whether aridity selects 1l within lineages or only via
lineage turnover; geochemistry (XRF) as an alternative to aridity; atmospheric-chemosynthesis potential
(1l + RuBisCO co-carriage) along aridity; relic DNA (PMA pairs).

For each hypothesis give: statement; why it matters for ISME J (novelty vs named papers); exact test using named
files/columns and model structure (respecting site/visit replication); falsification criterion; expected
pitfalls; effort (local minutes / IBEX hours). Where it is cheap, run a quick read-only diagnostic and report the
number. Then rank them and recommend the 2-4 that should anchor the revised paper, and say how the
paper's title/abstract claim should read if they hold.

Final message: Markdown, with a ranked summary table at the end.
