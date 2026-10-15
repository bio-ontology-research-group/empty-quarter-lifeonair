> Round-3 update: the expedition-4 candidates have been recovered and MAG sensitivities completed; Wang full text was checked. The R2-03 explanation based only on half-length coverage was incorrect: residual matches also meet that coverage threshold. See ../round3/response_round3.md for corrections and current dispositions. The text below records the round-2 response.

# Response to round 2 (Astra)

Revision completed by Codex after takeover from Claude, 12 October 2026. The central relative-exclusion framing remains unchanged. Items below distinguish corrections from remaining evidence/release gaps; partial items are not claimed resolved. Opus will independently review the revised version in round 3.

## R2-01: Partially addressed

The two topography rows for site 19 are collapsed to site means before joining, with unique-key assertions. Spatial and mixed models were regenerated (51 bulk, 45 surface, 59 root sites). Both source records remain frozen; their individual provenance has not been adjudicated. The averaged fields are slope and TWI, not the elevation used in the reported geographic partial correlations.

## R2-02: Addressed

Abstract, aridity Results, Discussion and Figure 2 caption distinguish ERA5 support from POWER/blended sensitivity. Bulk ERA5 rho is -0.401; surface POWER/blended intervals span zero.

## R2-03: Addressed with bounded scope

The original benchmark is labelled catalogue-cluster holdout with external references retained. A new all-reference component holdout is in benchmark_r2/disjoint_benchmark_summary.tsv and b2_disjoint_benchmark.sh (IBEX 53616234). At 70% linkage, maximum held-out identity reaches 83.3% because linkage also requires half-length coverage. At 50%, single-linkage connects almost whole guilds and recall collapses. No strict all-pairs identity ceiling is claimed.

## R2-04: Addressed as sensitivity

coxl_deconv.py inverts the two-class assignment matrix from the clu50 benchmark. h123_v2/sensitivity.tsv reports corrected form-I contrasts and the motif-only variant; neither changes the group 1l conclusion. This is a bounded error sensitivity, not validated functional activity.

## R2-05: Addressed with bounded scope

Twenty simulated communities from ten paired visits test the complete marker/SCG estimator (jobs 53616264, 53616320, 53616321). Table S2 reports calibration bias, correlations and ratio error with full and community-genome-removed references. The other-NiFe control is dominated by off-target hits in these mocks; the other-CoxL-only sensitivity is now reported. The simulations cover catalogue genomes with uniform coverage.

## R2-06: Addressed

apply_sub_fix.py propagates corrected CoxL labels and affinity flags to atmo_classification_v2.tsv. operon_summary.py removes orphan/removed-contig loci: complete-context n=445/224/53/125. Representative carriage is distinguished from community read references in the script documentation.

## R2-07: Addressed

coxL_motif_region.py inspected aligned active-site regions: the 59 motif-less calls have missing aligned regions. Motif variants from Imaura are recognised. The motif-only GE sensitivity treats these calls as unresolved. The aggregator v2-detection bug for the motif-only filename was also fixed, preserving the separate 2a column.

## R2-08: Partially addressed

The expedition-4 shotgun control SNTC is inspected and described: no focal-guild reads, with other-NiFe hits from Escherichia. Six expedition-4 amplicon extraction blanks are acknowledged; their candidate list remains unavailable. This corrects the R1 response statement that blanks existed only for expedition 5. Complete assay/batch contamination coverage remains unavailable.

## R2-09: Partially addressed

The exclusion list now includes all five batch-frequency flags, with 383 entries in its union. Read-based sensitivities were regenerated and remain consistent (group 1l excess -1.76). Genus-level and other MAG-based models have not all been repeated with this union; their contamination dependence remains a limitation.

## R2-10: Addressed

Genus models now include site-resampling through ratio construction, prevalence thresholds, high-completeness carriage and order fixed effects. Results are in m2d_genus_models.tsv. The group 1l estimate remains negative; order adjustment weakens the growth association and also yields a group 2a coefficient.

## R2-11: Addressed with bounded scope

Coverage by sample and position and 80%/90% restrictions plus fastest/slowest missing-trait assignments are in m2_growth_coverage_by_* and m2c_coverage_sensitivity.tsv. The intercept persists; the slope is sensitive to the fastest assignment. Statistical uncertainty of individual gRodon predictions is not estimated.

## R2-12: Addressed

Absent genera are zero-filled before site averaging; closure is asserted. Observed and turnover-only MAP correlations are -0.430 and -0.433.

## R2-13: Addressed

DNA adjustment is described as conditional attenuation (36%), with microbial, plant and extraction contributions made explicit. The compositional-expansion exclusion was removed.

## R2-14: Addressed

Single-kit sensitivity retains kit adjustment; bulk and DNA/aridity models include site and visit intercepts. Model tables were regenerated.

## R2-15: Addressed with bounded scope

NASA POWER daily data now supplies matched 90/365-day exposures for all five expeditions. Within-site product agreement is weak for 90 days (rho 0.193), stronger for 365 days (0.612). ERA5 surface/root associations are qualified by the POWER results on identical cohorts. Local gauge validation remains unavailable.

## R2-16: Partially addressed

ERA5-series-cell bootstrap and leave-one-expedition-out estimates are in H3_temporal_robustness.tsv. Bulk annual estimates change sign on omission of an expedition. The manuscript reports this and sparse sampling. The bootstrap currently retains finite estimates without a full per-resample convergence/singularity audit; it is conditional uncertainty, not uncertainty in the rainfall product itself.

## R2-17: Addressed by narrowing interpretation

The estimand is explicitly log2(GE+0.01), scaled by the IQR of within-site log(1+precipitation). The +/-0.5 bound is exploratory (sqrt(2) on the shifted scale), and selected marker/exposure models have 0.001/0.1 pseudocount sensitivity. Broad biological equivalence claims were removed.

## R2-18: Addressed with design qualification

PERMANOVA and dispersion now use site/expedition/position mean profiles and permutations within site x position. Position-specific analyses are also reported. Bulk expedition R2=3.7%, restricted p=0.036, dispersion p=0.26. The combined analysis permutes within position strata; its cross-position dependence needs to be read with the separate-position analyses.

## R2-19: Addressed

The within-genus term is described as a weighted community-share contribution, with the missing-side convention stated. Genus-turnover share is 97.4% (site-bootstrap 96.7-98.3%), conditional on the recovered catalogue.

## R2-20: Addressed

Methods now name the weighted Gaussian mixed model of clipped logits, its sqrt(abundance) weights and clipping thresholds; the result is a logit slope. Sensitivity slopes are -0.383 and -0.210.

## R2-21: Addressed by narrowing interpretation

The CBB result is a candidate three-KO co-occurrence screen. The >=90% complete subset is reported (28% vs 9%); pathway/sequence-form validation remains a further requirement for a functional interpretation.

## R2-22: Addressed

Methods specify the four BH families and their sizes. DNA/genus/NDVI/trait/co-occurrence models are labelled exploratory with unadjusted p values.

## R2-23: Addressed

Abstract, Results and Discussion use no detectable additional depletion rather than equality. Raw contrasts and excess intervals remain distinct.

## R2-24: Partially addressed

Baril and Kropp full-text evidence is incorporated with group/substrate specificity; Peng and Browne qualify lineage/coverage interpretations. Wang comparison is limited to its available publisher abstract; the full text remains inaccessible. The round-2 Flash digest and downloaded sources are retained.

## R2-25: Partially addressed; open before submission

Added selected TABLE claim checks, independent cohort/contrast recomputation, frozen input hashes, explicit conda lock, runtime versions, deterministic mechanism seed, upstream scripts and an entry point. The frozen-evidence verifier passes. Public deposition, unauthenticated fetching, portable upstream sequence derivation and complete numeric-prose coverage remain open. See reproducibility/README.md. No public reproducibility claim is made.

## R2-26: Addressed pending independent review

Figures were regenerated at their manuscript width, titles shortened, Figure 2c cohorts/intervals labelled, tree clades labelled with SH-like support, and the sensitivity table split into readable panels. Rendered PDF inspected; the third review should check final-size readability.

## R2-27: Addressed

DNA control slopes, copies-per-carrier contrast/cohort, surface CoxL Moran p and predicted-growth wording corrected. Aridity site counts and temporal PERMANOVA values now match regenerated outputs.

## R2-28: Addressed

Temperature is identified as POWER surface skin temperature; expedition explains 96% of its variation. The positive coefficient means weaker surface depletion. Residual variation and leave-expedition-out range are reported, with observational wording.

