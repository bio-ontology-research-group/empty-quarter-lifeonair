# Response to Gemini Pro, round 3 — 12 October 2026

Reviewer: Gemini 3.1 Pro High through agy, replacing the failed Opus attempt at Robert's request. The review recommends Minor Revision but labels three individual items major; we retain both statements without reclassifying the findings. Public release is still a submission gate. The reviewer did not refit the numerical analyses, and its initial response omitted explicit dispositions for R2-01–03; the focused addendum supplies these dispositions.

## R3-01 — Addressed

Recovered the authoritative September-9 control replay in the companion data-paper repository. Seven expedition-4 candidates from six blanks and 95 linked biological profiles are retained in `results/review_r3/trip4_candidate_calls.tsv`, with all seven sequences, a source checksum and genus-to-MAG mapping. The three named genera map to 14 MAGs, four beyond the previous exclusion union. Four unnamed genera cannot be mapped by genus and have not been assigned to an entire family. This is a conservative cross-assay flag, not proof of MAG contamination. The manuscript now states the assay, thresholds, scope and limitations. `code/import_trip4_controls.py` records the import; normal reanalysis uses the frozen exports without requiring the companion repository.

## R3-02 — Addressed for the requested models

The union contains 387 IDs including the two already excluded catalogue orphans. `code/r3_mag_flag_sensitivity.R` removes the union before constructing predictors and refits both primary genus models, all three within-order trait models and all three CBB models. Results are in `results/review_r3/mag_sensitivity_models.tsv`; cohort counts and runtime versions accompany them. Baseline genus, trait and CBB coefficients reproduce the existing outputs to numerical precision. These are focused model refits without fresh bootstrap intervals; model SEs are clearly labelled.

The retained cohort contains 5,051 genomes, 4,430 with predicted growth, 281 genera, 594 group-1l carriers and 395 group-1h carriers. Genus group-1l coefficients are -0.770 (SE 0.179) with phylum and -0.964 (0.215) with order. Group-1l trait coefficients are -0.0161 on log genome size, +0.661 GC percentage points and +0.117 on log predicted doubling time. CBB log-odds are +0.786 (0.401), +1.036 (0.425) after completeness adjustment and +1.495 (0.574) among >=90%-complete genomes. Thus the interpretation is retained; we do not claim unchanged significance for every exploratory coefficient.

Read GE were also reaggregated using the expanded union. In 106 paired visits at 47 sites the root/bulk group-1l median changes from -2.275 to -2.285 log2 units, and excess depletion from -1.827 to -1.760. `code/r3_evidence_checks.py` independently reproduces the baseline contrasts and saves all exclusion contrasts. This extension updates paired medians; it does not claim a new full spatial/temporal bootstrap analysis with the additional four exclusions.

## R3-03 — Addressed by qualification

The existing round-2 analysis already resamples ERA5 grid cells. Methods now states that uncertainty is conditional on those grids, does not propagate local rainfall measurement error, and retains finite bootstrap estimates without a full convergence audit. Product dependence and residual spatial structure remain explicit. No unconditional rainfall effect or demonstrated equivalence is asserted.

## R3-04 — Addressed by primary full-text verification

The author-uploaded Wang et al. full text became accessible via ResearchGate. Sections 2.1–2.3 and 4.3 support an across-lineage association of trace-gas oxidation with smaller genomes, but the study covers northern China and uses a different growth-prediction method. The manuscript's erroneous “global soil survey” wording is corrected, and the comparison is bounded accordingly. Source extracts and provenance are in `literature/wang_2025_*extract.txt` and `literature/round3_source_check.md`. The PDF download remained blocked; we do not label the extracted text a downloaded PDF.

## R3-05 — Open before submission

The local reconstruction and frozen-evidence route is documented and executable. A public archive of the required sequences/counts/inputs, a tested unauthenticated fetch, portable upstream derivation, ENA release confirmation and final author/funding metadata remain required. The named repository is private. No public accession or public reproducibility claim has been fabricated. This author-facing release work remains outside the completed scientific revision loop.

## Additional benchmark correction

Our independent final audit found that 34 of 6,478 comp70 held-out proteins still have >70% identity matches at >=half-length coverage (maximum 83.3%). The round-2 explanation attributing this solely to the coverage threshold was incorrect. Graph discovery used E<=1e-5; verification used E<=10 and found edges absent from that graph. Methods and the table caption now describe a graph-component holdout with residual homology, not a guaranteed identity-disjoint benchmark. The comp50 audit likewise finds 27 proteins above 50% (maximum 66.7%). Counts are independently regenerated in `results/review_r3/benchmark_residual_homology.json`. This limits benchmark interpretation; it does not alter the primary GE or paired ecological estimates.

Final implementation check: the CBB sensitivity now uses the primary model's named-order convention and excludes unclassified orders; explicit baseline guards pass for all three CBB fits. This corrects a small intermediate discrepancy before final reporting.

## Review addendum and correction

Gemini accepts R2-01–03 and R3-01–04, with benchmark scope and public release explicitly limited. It initially alleged a baseline group-1l genome-size coefficient of +0.08. This was not present in the manuscript or data: -0.0189549 is the coefficient and 0.0847315 is the p-value. After receiving the exact rows, Gemini explicitly withdrew the allegation and reported no remaining scientific blocker from its addendum (`clarification_output.txt`). The initial addendum remains intact for auditability. Its acceptance of R2-01 does not establish provenance of the duplicate topography records: that source-adjudication limitation is retained.
