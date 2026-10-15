# P3 "Life on air" (ISME Journal, manuscript draft)

Atmospheric H2 and CO scavenging capacity in Empty Quarter soils. Split from the MAG/dark-proteome
paper on 10 October 2026 (source provenance in `../MIGRATION.md`). The paper concerns encoded capacity; heavy computation runs on IBEX.

Start with `results/README.md` (top section: current state after review round 3), then `reviews/` for the
critical-review loop (Astra rounds 1–2, Gemini Pro round 3, and point-by-point responses).

| Folder | Content |
|---|---|
| `code/` | Analysis scripts (IBEX Slurm jobs, R and Python); R environment `code/environment.yml`. `make_lifeonair_figure.py` is Alejandra's earlier catalogue-census figure. |
| `results/` | Compact inputs and outputs; per-sample read counts are fetched by `../reproducibility/fetch_large_files.sh` |
| `figures/` | Figure 1 (`fig1_guilds`), Figure 2 (`fig2_scaling`), Figs S1–S2, all from `code/fig_v2.R`; older `fig_h123`/`fig_mechanism` are superseded |
| `literature/` | Novelty-check literature (further PDFs in `../resources/literature/trace-gas-deserts/`) |
| `manuscript/` | The ISME J manuscript `main.tex` with `references.bib` and generated supplementary table rows (`make_tables.py`); text removed from P4 on 10 Oct 2026 |
| `reviews/` | Review loop: `round*/` independent review, findings table, response letter |
| `plans/` | Earlier atmospheric-metabolism plan |

Scripts are run from the repository root, for example
`micromamba run -n eq-lifeonair Rscript lifeonair/code/h123_v2.R`.

The proteome paper and its historical atmospheric summaries are maintained in the separate `empty-quarter-proteome` repository.

## Verification and review status

Round 2 revisions and the response are in `reviews/round2/`. Run
`bash lifeonair/run_lifeonair.sh --verify` for frozen-input and selected claim checks;
see `reproducibility/README.md` for model reanalysis, the explicit environment lock,
and the remaining public-archive requirements. The scientific conclusions concern
relative group 1l representation; climate associations depend on rainfall product
and spatial uncertainty. Full public reproduction remains a submission requirement.

Round 3 was reviewed by Gemini 3.1 Pro High through agy at Robert's request, after
Opus failed before producing findings. The review recommends Minor Revision;
`reviews/round3/response_round3.md` records the final revisions and remaining public
release requirement. New control candidates, sensitivity fits and benchmark audit
are in `results/review_r3/`. The three scientific review/edit rounds are complete;
this does not mean the draft is ready for submission without the release tasks.

A metabolic-model feasibility pilot was completed on 12 October 2026:
[`results/metabolic_pilot/REPORT.txt`](results/metabolic_pilot/REPORT.txt).
It covers 72 matched gapsmith models and 25 exact-genome gapseq comparisons.
Hydrogenase representation and reconstruction sensitivity prevent mechanistic
interpretation at present; the conditional field regression was not pursued.
The manuscript's biological conclusions do not incorporate these pilot outputs.

The subsequent [gapseq-only predictive pilot](results/gapseq_traits_pilot/REPORT.txt)
used 36 exact-genome reconstructions and the paper's 106 paired root/bulk visits.
Energy checks passed, but two carbon-energy traits worsened held-out-genus
prediction while providing only a small held-out-site gain. The prespecified
scale-up criteria were not met; the larger eligible cohort was not analysed.
