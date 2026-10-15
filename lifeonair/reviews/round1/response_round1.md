> Round-3 update: all seven expedition-4 blank candidates and sequences are now retained and the expanded exclusion union has been tested. See ../round3/response_round3.md.

# Response to review round 1 (Astra, gpt-6-astra, 12 October 2026)

Decisions taken with the PI at the agreed stop point: H1 is reframed as **relative** exclusion of group 1l
(not loss per gram), and H2 is **narrowed** to "the root effect is specific to group 1l; the surface is a broad
filter". The former three confirmatory hypotheses are replaced by exploratory questions, with that history
stated in the Introduction. New code is listed per finding; all corrected numbers are registered in
`notes/CORRECTIONS.md` (entry 2026-10-12). Status: **A** accepted and fixed, **P** partially addressed,
**R** rebutted / kept with justification.

| ID | Status | Response |
|---|---|---|
| R1-01 | A | Motif rule inverted. `fix_coxL_v2.py`: AYXCSFR = form I, AYRGAGR = form II, motif-less proteins form I only at >55% id / ≥90% cov to Cordero refs. 1,901→284 form-I proteins (245 reps after R1-19). Form-I loci now show 96% coxMSL co-occurrence (was 73%), which supports the correction. Reads relabelled (`marker_labels_v2.tsv`), GE recomputed (`guild_genome_equivalents_v2.tsv`). Bulk form-I median 0.616→0.105. All CO claims rewritten; "CO capacity independent of position" withdrawn. |
| R1-02 | A | Group 2a separated into its own column (`r3_aggregate.py` with v2 labels), removed from the other-[NiFe] control, cited Islam 2020 and Grinter 2023 (Huc), and analysed as a third high-affinity group throughout. 2a shows no root excess (−0.02) but surface depletion (−0.22). |
| R1-03 | A | Pure ERA5 daily precipitation fetched to 2025-10-31 (`fetch_era5_2025.py`); the original Open-Meteo default ("best_match", ERA5+IFS) was not reproducible on re-query and is now a sensitivity product. Temporal analysis covers all 5 expeditions, 31 months (103 bulk visits at 35 sites). |
| R1-04 | A | Space-for-time and "independent of rain" claims withdrawn. Within/between decomposition with an explicit ±0.5 log2 equivalence bound: bulk 1l is neither associated nor equivalent; surface and root-adjacent 1l are lower after wetter 90-day periods (q 0.013, 0.020). Text now says the data neither show nor exclude a bulk-soil response. |
| R1-05 | A | Introduction states that initial expectations were not supported and that the questions were formulated after inspecting the data; "as H1 predicted" language removed; BH families declared. |
| R1-06 | A | Sites 62/64 and metadata-less samples excluded: 546 samples, sites 1–60 (59 with data), 106 root/bulk visits at 47 sites. |
| R1-07 | A (reframed) | Title, abstract and results now say "relative depletion / fewer copies per genome". New analysis A1: 1l GE scales with log2 DNA concentration (−0.32/doubling; within-visit −0.30, between −0.35, single-kit −0.36, bulk-only −0.23; 1h/2a/controls ≈0); root DNA 10× bulk; root effect −2.12 → −1.36 after DNA adjustment. Stated explicitly that DNA is not calibrated biomass; Azarbad 2022 cited; qPCR/spike-in named as the test. |
| R1-08 | A | Causal wording removed ("relatively depleted beside roots"); Limitations list moisture, oxygen, plant identity and microsite as unseparated by depth matching. |
| R1-09 | A | Mediation language removed; the composition-PC1 partial (→0.01) is reported as removal of the association, without causal interpretation. |
| R1-10 | A | "minor part" removed. Reported as associations: growth shift slope, excess at zero shift with site-bootstrap CI, and NDVI coefficients; López 2023 cited. |
| R1-11 | A | Replaced by within/between mixed model with expedition fixed, site random, site bootstrap (400), equivalence bound. |
| R1-12 | A | Site-visit random effects added to joint marker×position, NDVI (site:trip), DNA-scaling, H3 mixed models. |
| R1-13 | A | Block bootstrap for k = 4/6/8/10, Moran's I of residuals, partials for geography, XRF PCs, Ca/Cl/S and composition PC1. Bulk 1l association is **not** robust to blocking and the manuscript now says so; surface associations are robust (1l retains residual autocorrelation, reported). |
| R1-14 | P | Three precipitation products; Wang 2025 SaRa cited for ERA5 uncertainty in Saudi Arabia; MAP treated as a site ranking. Other water-balance drivers (PET, dew, groundwater) not available; stated. Air temperature and SMAP soil moisture added for the surface exposure models. |
| R1-15 | A | PERMANOVA with expedition permuted within site (`how(blocks=site)`), dispersion test reported; persistence language removed. |
| R1-16 | A | Synchrony restricted to bulk soil (8 sites ≥4 trips); weak and no longer cited as evidence of persistence (dropped from the manuscript). |
| R1-17 | A | Retention/invariance language removed; only "changed no more than controls" with CIs and q-values. |
| R1-18 | A | BH within declared families (position contrasts, spatial correlations, temporal slopes, surface exposures); q-values reported. |
| R1-19 | A | `s1_sub_reconcile.sh` (IBEX 53614016) mapped all CDS of the 1,010 `_sub` reps to retained contigs by sequence MD5; 220 database subjects lie on removed contigs. Reads to these remain correctly labelled community genes, so GE are unchanged; carrier lists (`atmo_classification_v2.tsv`, form-I reps), growth (`rep_growth_v2.tsv`, 891 re-estimated) and trait analyses use retained contigs only. 1l carriers 599→597, 2a 136→134, form-I 246→245. |
| R1-20 | A | Held-out benchmark (`b1_classifier_benchmark.sh`, IBEX 53613713): leave-genome-out and leave-50%-cluster-out, 5 folds. At 60% floor, 50%-cluster holdout: 1l recall 0.89 / precision 0.999, 1h 0.99/0.98, 2a 0.94/0.999, form I 0.82/0.77; 2 of 1.5 M negative fragments called 1l. Table S1. |
| R1-21 | A | Exact-cancellation claim corrected ("up to the pseudocount and separate replicate medians"); pseudocount 0/0.001/0.01/0.1 reported (excess −1.82/−1.82/−1.83/−1.32). |
| R1-22 | P | Genus model rerun with 2a, nonzero-ratio genera and ≥2-rep genera, carriers from retained contigs; site bootstrap of the visit-level intercept. Finer phylogenetic control (order random effect at genus level) not done; acknowledged that carriage may mark lineage ecology rather than the enzyme. |
| R1-23 | P | Batch-stratified decontam (`decontam_batch.R`, expedition×kit): 5 of 3,449 testable reps flagged, 1 a carrier. Correction in R2: expedition 4 also has a shotgun no-template control and six amplicon blanks; the latter candidate list remains unavailable. See response_round2.md R2-08. DNA concentrations are missing for trip 4. |
| R1-24 | P | The PMA pairs are 16S amplicons from two sites with no bulk soil; mapping to MAGs at 97% identity cannot resolve hydrogenase carriage. Stated as a limitation rather than analysed. |
| R1-25 | A | Novelty narrowed to the group-resolved, replicated contrast; Li 2023 subgroups (1d, 2a) stated; López 2023 cited. |
| R1-26 | A | Moisture paragraph rewritten citing both Jordaan 2020 (relative decline, uptake stimulated) and Imminger 2024 (1l transcription up); direction stated as unresolved. |
| R1-27 | P | All new code is in `lifeonair/code/`; IBEX workbench scripts copied; public data archive and claims.tsv registration remain TODO (claims.tsv is coordination-only per CLAUDE.md). |
| R1-28 | A | Sensitivity table now covers root and surface excess, 1h/form-I excess, and both MAP products; temporal analyses are reported separately with all products. |
| R1-29 | A | "Operon" replaced by "neighbourhood co-occurrence" with an explicit statement that order and strand are not tested. |
| R1-30 | A | Figures rebuilt (`fig_v2.R`): larger text, CIs (Fig 1d, 2c, 2d), temporal panel (Fig S2), 2a coloured. |
| R1-31 | A | Replaced with gain in marginal R² from adding position to an expedition-only model. |
| R1-32 | A | README/SUMMARY regenerated from v2 outputs; estimands labelled (visit-median ratio vs mixed-model coefficient). |
| R1-33 | A | "predicted growth potential" used throughout; the hydration study is described as a desert soil without misattributed location. |

Additional analyses suggested in the consultation (`round1_consult/`): A1 biomass scaling, A2 turnover vs
within-genus (root and aridity), A3 surface exposures, A4 CBB co-carriage (not significant within order; reported
as descriptive), A5 genome traits (`a5_genome_traits.R`), A6 copy number (F/M/C decomposition).
