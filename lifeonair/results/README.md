# P3 "Life on air": analyses for the ISME Journal paper

## Current state (12 October 2026, after review round 3)

**Read this section first; everything below "Historical record" is pre-review and superseded where it
conflicts.** Corrections are registered in `notes/CORRECTIONS.md` (entry 2026-10-12); the review and the
point-by-point response are in `lifeonair/reviews/round1/`.

Framing (decided with the PI): the former H1–H3 are now exploratory questions. H1 = *relative* depletion
of group 1l beside roots (copies per genome, not per gram); H2 narrowed to "the root effect is specific to
group 1l, the surface is a broad filter"; the bulk-soil aridity association is not robust to spatial blocking.

Current folders and scripts (run from the repository root):
- `review_r3/` — seven recovered expedition-4 candidates and sequences, 387-ID exclusion union, focused genus/trait/CBB refits, paired GE contrasts and residual-homology audit. See `reviews/round3/response_round3.md`; the all-reference component holdout is not strictly identity-disjoint.
- `h123_v2/` — `code/h123_v2.R`: position contrasts, joint models, sensitivity, spatial (block bootstrap,
  Moran's I, partials), temporal within/between with ±0.5 equivalence, PERMANOVA, A1 DNA scaling,
  A2 turnover decomposition, A3 surface exposures, A4 CBB co-carriage, A6 F/M/C decomposition.
- `h1_mechanism_v2/` — `code/h1_mechanism_v2.R`: NDVI, gRodon growth, excess vs growth shift, genus models.
- `a5_genome_traits/` — `code/a5_genome_traits.R`: genome size, GC, doubling time by carrier status.
- `benchmark/` — `code/b1_classifier_benchmark.sh` (IBEX 53613713): held-out classifier benchmark.
- `decontam/` — `code/decontam_batch.R`: expedition×kit-stratified frequency test.
- `operon/*_v2.tsv` — `code/operon_summary.py` v2 mode.
- `inputs/` — corrected inputs: `atmo_classification_v2.tsv`, `rep_growth_v2.tsv`, `guild_coxL_reps_v2.txt`
  (`code/apply_sub_fix.py` after `code/s1_sub_reconcile.sh`, IBEX 53614016), `rep_genome_size.tsv`
  (IBEX 53611804), `era5_daily_1995_2025.tsv` (`code/fetch_era5_2025.py`).
- Figures: `code/fig_v2.R`; tables: `code/make_tables.py`. Workbench record:
  `/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1/README.md`.

Key current numbers: root vs bulk group 1l −2.27 log2 (106 visits, 47 sites), excess beyond controls −1.83
(CI −2.62 to −1.00); 1h/2a/form I ≈0. Group 1l GE −0.32 per doubling of DNA; root effect −2.12 → −1.36
after DNA adjustment. Genus turnover carries 97% of the carrier-share change. Excess at zero growth shift
−1.57 (−2.22 to −0.93). Surface: all four groups depleted (1l −0.83, 1h −0.42, 2a −0.22, form I −0.24).
Bulk 1l vs ERA5 MAP rho −0.37, block CIs include 0; surface rho −0.48 to −0.54, robust.

## Historical record (9–10 October 2026)

Working record for the split-off atmospheric trace-gas paper (`plans/manuscript-split-2026-10.md`).
Nothing here is yet a registered manuscript claim in `reproducibility/claims.tsv`. Every number below is a
genomic **capacity** measure, not activity or flux. Genome equivalents (GE) are copies per
single-copy-gene genome equivalent. They are relative capacity, not cell fractions: group 1l genomes
carry about 1.9 copies on average and group 1h genomes about 1.1.

## Hypotheses (decided 10 October 2026)

The original P3 hypothesis was a stable, compartment-independent baseline despite strong temporal
turnover. It is contradicted by the data:
- Trip explains about 2% of community composition (`stability/`).
- Compartment has a strong effect on hydB.

The replacement hypotheses were set after an independent second assessment of the same data and literature:

- **H1.** Root-adjacent soil selectively loses group 1l H2 scavengers. Root-adjacent and bulk "Deep" samples
  are both taken at 5-10 cm (data-paper methods), so this is a same-depth contrast. The 1l loss must exceed
  the loss of control markers (other [NiFe] large subunits, non-form-I CoxL).
- **H2.** H2 and CO scavenging capacity are decoupled across microhabitats: high-affinity hydB varies,
  form-I CoxL does not.
- **H3.** Long-term aridity sets the spatial 1l baseline, but between-trip change at a site does not track
  antecedent precipitation (space-for-time substitution fails).

Dropped:
- Functional redundancy with compensatory dynamics. Carrier genera fluctuate in step (phi 0.115 vs null 0.048).
- A persistence advantage of carriers. It is not present within genus or family.
- Rain-pulse responses.
- The RQ/D3 16S capacity index as primary evidence.

## Results (all 573 paired-read samples screened; 550 with metadata and >= 10 SCG genes)

**H1** (`h123/h1_*`; 107 site x trip pairs at 48 sites):
- Group 1l GE is a median log2 −2.23 (about 4.7-fold) lower at roots than in bulk sand at the same depth
  (site bootstrap 95% CI −2.82 to −1.57; lower in 89% of pairs; site-level Wilcoxon p 8e-9).
- Excess over the control markers: −1.97 (CI −2.53 to −1.05, p 3e-8).
- Group 1h shows no excess (−0.06, p 0.90). Form-I CoxL is slightly less depleted than the controls (+0.17).
- Marker x compartment interaction: LRT chi2 308, df 4, p 2e-65. Only the 1l term departs from the reference (−1.74, SE 0.12).
- Holds in all 5 trips (−1.6 to −3.4).
- All ten major 1l-carrier orders are depleted (Acidimicrobiales −2.8, Solirubrobacterales −1.7, Rubrobacterales −2.3, Gaiellales −2.4).

**H2** (`h123/h2_*`):
- Raw-read ratio of high-affinity hydB to form-I coxL reads (no shared normaliser):
  - roots vs bulk: log2 −1.60 (CI −1.94 to −1.15)
  - surface vs bulk: log2 −0.95 (CI −1.45 to −0.48)
- Surface vs bulk: 1l −0.97, form-I CoxL −0.02 (p 0.18).
- Compartment variance share: 1l 20.5%, 1h 5.5%, form-I CoxL 3.1%, controls 6–7%.
- CO capacity is spread over about twice as many effective genera as H2 capacity (about 20 vs 10–11 per sample).

**H3** (`h123/h3_*`; site medians, bulk 5-10 cm, 51 sites):

Spatial:
- 1l vs ERA5 1995–2024 mean annual precipitation: Spearman rho −0.51 (p 1e-4).
- 1l vs NASA POWER 2001–2020: rho −0.43 (p 0.0015).
- Spatial-block bootstrap CI −0.78 to −0.18. Leave-one-block-out range −0.61 to −0.37.
- Partial on latitude/longitude/elevation: −0.48 (p 4e-4).
- Mixed model across compartments: 1l −0.49 log2 per SD log MAP (p 4e-5); 1h −0.14 (p 0.09); form-I CoxL −0.04 (p 0.41).
- **Caveat:** after partialling out the first community-composition axis, rho falls to −0.20 (p 0.16).
  Aridity and the main composition axis covary, and 1l carriers are a large part of that axis.

Temporal:
- Within-site, between-trip change does not track antecedent precipitation in bulk soil (P365 rho 0.02).
- The 95% CI of the temporal slope (−1.57 to 1.70 per log unit) excludes the spatial slope (−4.15).
- Surface 1l rises after wet years (rho 0.38), but this is aliased with trip: trip-adjusted p 0.31.
  Form-I CoxL keeps a trip-adjusted surface response (p 0.010).

**Sensitivity** (`h123/sensitivity.tsv`): the H1 and H3 effects are stable at identity 50/60/70%, with
contamination-flagged reps removed, without trip 5 (excess −1.61), and without trip 4 (101-nt libraries;
excess −1.74).

**Weak 1l calls** (Alejandra's open item B9: weakest 1l call 37.8% identity / 20.5% coverage):
- Restricting the 1l label to the 1,002 of 1,124 EQ 1l proteins typed at >= 50% identity and >= 80% coverage
  gives an H1 log2 fold change of −2.19 (all calls: −2.23). Weak calls contribute 0.07 GE in Deep.
- Registered as sensitivity variant `strict_1l_calls` (reads to the 122 weak proteins in `reads/weak_1l_proteins.txt`
  discarded, `reads/guild_genome_equivalents_strict1l.tsv`): ratio −2.19, excess −1.97.

**Operon context** (`operon/`; loci with >= 5 genes on both sides):
- 1h: small subunit adjacent 84%, maturation genes 96%.
- 1l: 73% and 88%. Recurrent DUF948 and NifU neighbours.
- Form-I CoxL: 73% in complete coxMSL operons.

**Contamination** (`decontam/`):
- 376 reps flagged by blank genus. No rep shows the decontam frequency signature (3,453 tested; minimum score 0.19).
- Guild CoverM abundance: 18.2% → 15.9% of the whole community after removing flagged reps.

## Layout and provenance

| Path | Content | Produced by |
|---|---|---|
| `inputs/` | Compact frozen inputs: R232 taxonomy, guild typing (IBEX `06_atmospheric/lifeonair_clean`, Aug 2026), CoverM long table, samples, environment (ERA5, NASA POWER, NDVI/SMAP), DNA concentrations | `lifeonair/code/prep_environment.R`, `fetch_nasa_power_precip.py`; sources in `inputs/env_sources.sha256` |
| `reads/` | Read-based GE (`guild_genome_equivalents*.tsv`), DB labels/lengths; per-sample counts on IBEX (`reproducibility/large_files.sha256`) | `lifeonair/code/r1_build_db.sh` (job 53522301), `r2_reads_array.sh` (array 53524441; 53522544 superseded), `r3_aggregate.py` |
| `decontam/` | MAG contamination screen | `lifeonair/code/decontam_screen.R` |
| `stability/` | Trip/compartment mixed models, PERMANOVA, synchrony | `lifeonair/code/stability_turnover.R` |
| `h123/` | Registered H1-H3 analyses and sensitivity | `lifeonair/code/h123_analysis.R` |
| `h1_mechanism/` | H1 mechanism proxies (NDVI, gRodon growth rate) | `lifeonair/code/h1_mechanism.R`; growth rates from `lifeonair/code/g1_growth_rate.sh` (job 53541094) into `inputs/rep_growth.tsv` |
| `operon/` | Operon context | `lifeonair/code/c1_operon_context.sh` (job 53524544), `operon_summary.py` |
| `trees/` | hydB and coxL trees | `lifeonair/code/t1_trees.sh` (job 53535825) |
| `../figures/` | Figure 1 (tree rooted on group 2, H1, H2, H3), Figure 2 (mechanism proxies), supplementary CoxL tree | `lifeonair/code/fig_h123.R`, `fig_mechanism.R` |

R environment: `lifeonair/code/environment.yml` (micromamba env `eq-lifeonair`). The IBEX working
directory is `06_atmospheric/lifeonair_reads`, linked as
`EmptyQuarter_Data/processed/sequencing-results/catalog/analyses/lifeonair_reads`.

## H1 mechanism proxies (`h1_mechanism/`, 10 October 2026)

Soil organic carbon was never measured: P1 has XRF elements, pH and weather only, and the RQ26 confounder
audit lists organic carbon as unmeasured. Li et al. 2023's organic-carbon mechanism is therefore tested
through two proxies.

**M1, plant input (NDVI, 60-day window, site x trip).** Mixed models with log MAP and site/trip random effects:
- 1l falls with NDVI in bulk 5-10 cm (−0.64 log2 per SD, SE 0.16) and at the surface (−1.23, SE 0.15).
  This is independent of MAP; site NDVI vs MAP rho is 0.57.
- Root-adjacent 1l is flat (−0.02).
- The root excess loss is smaller at greener sites (site-median rho +0.36, p 0.014).
- Reading: vegetation depresses 1l throughout the soil, and the root contrast shrinks where bulk soil is
  already plant-influenced.

**M2, copiotrophy (gRodon2 2.7.3 predicted minimum doubling time).** Job 53541094; 4,817 of 5,438 reps
estimated. The first run, 53541016, matched only old ribosomal-protein names and is superseded.
- **a. Carriers grow slowly.** Within order: 1l +0.12 log doubling time (SE 0.03; medians 4.2 vs 3.3 h);
  1h +0.11; form-I CoxL +0.06. Oligotrophy is shared by 1h carriers, so it cannot by itself explain why only
  1l is lost.
- **b. Root-adjacent communities are more copiotrophic.** The abundance-weighted log doubling time is
  −0.25 lower than in bulk soil (faster in 79% of 106 pairs; site Wilcoxon p 1.5e-5). The surface shift is
  −0.29.
- **c. The 1l loss scales with that shift, but only partly.**
  - Excess loss vs shift: Spearman rho 0.22 (p 0.023); lmer slope 1.55 log2 per log unit (SE 0.62).
  - The predicted excess loss at zero shift is still −1.59. Copiotrophy accounts for a minor part of
    the −1.97 excess.
- **d. Carrying 1l predicts root depletion beyond growth rate and phylum.** 288 genera with ≥0.05% mean
  bulk abundance; genus root log2 fold change regressed on median log doubling time and the carrier
  fractions:
  - growth slope −0.37 (p 6e-4): slow growers are depleted;
  - fraction carrying 1l −0.75 (p 3e-5);
  - 1h (+0.07) and form-I (−0.02) carriage are not associated with depletion;
  - R² 0.30.

**Interpretation for H1.** Two findings are consistent with Li et al.'s plant-carbon/copiotroph mechanism:
root-adjacent soil shifts to faster growers, and greener sites have less 1l in bulk soil. But a large
1l-specific loss remains after growth rate, phylum and the community shift are accounted for. 1h carriers
are equally slow-growing and are not lost. H1 therefore stands as a group-1l-specific effect. Organic carbon
availability is a likely contributor, not the whole explanation. A remaining candidate is the higher
H2/O2 or moisture regime beside roots, which would favour 1h over 1l; this is untested.

## Literature read on 10 October 2026 (`../literature/`)

**Li et al. 2023, ISME J (fertile islands)** is the closest precedent for H1.
- Mu Us Desert, semiarid (180 mm). 60 metagenomes under and outside Hedysarum shrubs, 0-60 cm.
- H2- and CO-oxidation genes are lower under shrubs, most strongly at 0-10 cm. H2 genes track soil and
  dissolved organic carbon (rho −0.48, −0.43); CO genes track dissolved organic carbon only (−0.31).
- Added soluble carbon suppresses measured H2 and CO uptake. The trace-gas oxidisers are oligotrophs.
- Consequence: "plants deplete trace-gas oxidisers" is not new and H1 must cite it. What H1 adds:
  - specific to group 1l, with 1h unchanged (Li did not separate 1h from 1l; their typing names group 1d/2a);
  - a loss beyond control genes, with copies per genome rather than TPM;
  - root-adjacent soil against bulk soil at the same depth;
  - hyperarid sand, 48 sites and 5 trips;
  - form-I CO capacity not depleted. Li found a weaker CO decline, so this extends H2 rather than contradicting it.
- Li's organic-carbon mechanism is the obvious explanation for H1. Soil organic carbon was never measured;
  see the proxy tests under "H1 mechanism proxies" above.

**Tribbia et al. 2026, Appl Environ Microbiol (Namib fog/aridity transect)** supports the direction of H3.
- 10 sites, 1-5 cm, one time point, activity measurements and amplicons, no gene typing.
- Biomass, ATP and respiration fall with aridity, but atmospheric H2 oxidation per 16S copy is highest in
  the most arid fog zone. Rubrobacteriaceae are enriched there.
- It neither tests nor contradicts the H3 temporal result.
- Practical point for a gas-uptake add-on: their assays used 2 g of soil stored at 4 °C, not frozen.

Neither paper bears on P4 (the dark proteome) or P2 (ecology).

## Open items

- Manuscript draft exists (`../manuscript/main.tex`); register its headline numbers in `reproducibility/verify.py` and `claims.tsv`.
- H3 composition-axis confound: decide whether a causal-mediation framing (aridity → lineages → 1l) is defensible.
- Relic DNA: only 9 PMA pairs exist; state as a limitation.
- Optional: gas-chromatography uptake on archived soil (2–3 g at −20 °C is workable for potential rates).
