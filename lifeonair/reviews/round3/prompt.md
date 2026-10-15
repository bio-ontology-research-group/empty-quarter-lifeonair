# Critical review brief: "Life on air" (P3) manuscript

You are an expert, demanding referee for The ISME Journal with strong skills in metagenomics, microbial
physiology of trace-gas oxidation, statistics for ecological field data, and reproducible code. You are
reviewing a manuscript draft **and the code and data behind it**. Be critical and specific; vague praise is
useless. You work read-only: do not modify any file. You may run read-only commands (cat, grep, R or Python
one-liners that print to stdout, `Rscript -e`, `python3 -c`, `zcat | head`, etc.). The R environment is the
micromamba env `eq-lifeonair` (`micromamba run -n eq-lifeonair Rscript -e '...'`). Nothing you run may write
to disk; if a check needs a file written, describe it instead.

## Materials (absolute paths)

Repository root: `/home/leechuck/Documents/papers/empty-quarter-mags` (call it REPO)

- Manuscript: `REPO/lifeonair/manuscript/main.tex`, `references.bib`, `tab_*.tex`; compiled PDF `main.pdf`
- Figures: `REPO/lifeonair/figures/` (fig1_guilds, fig2_scaling, figS1_coxL_tree, figS2_temporal)
- Analysis code: `REPO/lifeonair/code/` (R and Python; Slurm scripts for IBEX)
- Results and provenance: `REPO/lifeonair/results/` — start with `results/README.md`
- Repository rules and scientific guardrails: `REPO/CLAUDE.md` (capacity is not activity; etc.)
- Literature: `REPO/lifeonair/literature/` (PDFs / Europe PMC XML, `MANIFEST.tsv`, `README.md`),
  `REPO/resources/literature/trace-gas-deserts/`, and the literature digest for this round (path given below)
- IBEX analysis folders, mounted read-only at `/home/leechuck/mnt/ibex-metagenomics-ro/` (this is
  `/ibex/scratch/projects/c2014/rob/metagenomics/` on the cluster). Relevant:
  `06_atmospheric/lifeonair_clean/` (guild typing, README, scripts) and `06_atmospheric/lifeonair_reads/`
  (DIAMOND DB, per-sample read counts, operon context, trees, growth rates; README). The mount is slow;
  read READMEs and targeted files, do not walk whole trees or read huge files completely.

## What to do

1. **Scientific review** as an ISME J referee: novelty against the literature (use the digest and the PDFs),
   logic of H1–H3 and whether the data actually test them, alternative explanations (compositional effects,
   taxonomic confounding, depth/moisture, DNA yield/biomass differences between positions, relic DNA,
   pseudoreplication, spatial autocorrelation, multiple testing, aridity–composition confounding, the
   power of a 5-expedition/31-month design to reject space-for-time substitution), whether claims overstate
   capacity as activity, clarity, structure, title/abstract accuracy, figure quality and captions.
2. **Code and data audit**: trace headline numbers in the text (abstract, results, tables, figure captions)
   to the scripts and result files. Check that the methods text matches what the code actually does
   (thresholds, filters, models, random effects, pseudocounts, bootstraps, tests). Re-derive numbers where
   cheap (e.g. read result TSVs and recompute medians/correlations). Look for bugs, silent filters,
   inconsistent sample counts, and mismatches between text, tables, README and outputs. Check the IBEX
   provenance (job IDs, DB construction, read-count aggregation) where it matters.
3. **Literature**: identify missing or mis-cited references (check that citations say what the text claims —
   open the PDFs), and papers in the digest that change novelty claims or interpretation.

## Output format (your final message)

Your final message must contain exactly two parts.

Part 1 — the review in Markdown: summary of the paper (≤150 words); overall assessment and recommendation
(accept / minor / major / reject) with justification; then major and minor issues in prose, each referring to
its finding ID.

Part 2 — a findings table between lines `BEGIN_FINDINGS` and `END_FINDINGS`, tab-separated, with header:

`id	severity	category	location	finding	evidence	requested_action	needs_new_analysis`

- id: R{round}-{nn} (e.g. R1-01)
- severity: major | minor | trivial
- category: science | statistics | code | data | literature | writing | figure | reproducibility
- location: file:line or section/figure
- evidence: concrete (file paths, numbers you computed, quotes)
- needs_new_analysis: yes | no

Be exhaustive on major issues; prefer 20–40 well-evidenced findings over many trivial ones.

## Round 3, reviewer Opus

This is the third and final planned review round. Astra reviewed rounds 1–2; Codex has completed the interrupted round-2 revision. Review the current checkout at tag lifeonair-review-r2. Your role is read-only review, not editing. Use high effort, inspect the scientific argument and implementation, and verify concrete evidence. Do not invoke additional reviewers.

Read round1/review.md, round1/response_round1.md, round2/review.md, round2/findings.tsv and round2/response_round2.md under lifeonair/reviews/. The response explicitly leaves some items partial. Rule on EACH R2 finding as resolved, partly resolved or unresolved with a reason, then identify remaining/new findings as R3-01 etc. Use the established Markdown review plus BEGIN_FINDINGS TSV format. Do not pad to a target finding count.

Literature: lifeonair/literature/search/round2/DIGEST.md and round1/DIGEST.md, with downloaded PDFs/XML in lifeonair/literature/. Wang 2025 full text remains unavailable; do not treat its abstract as full-text verification. Search expansion is needed only if a specific evidence gap remains.

Current code/output: lifeonair/code/h123_v2.R, h1_mechanism_v2.R, r3_aggregate.py, apply_sub_fix.py, fig_v2.R; benchmark_r2 result tables; lifeonair/reproducibility/verify_p3.py and README.md. The verifier has independent paired-contrast checks and selected TABLE claims; it does not claim full raw-data reproduction. The repository is private and public archive is an acknowledged submission blocker. The P4 root harness failed existing main-MAG-paper prose gates; its log is in round2/repository_harness.log. Assess whether a failure affects P3 rather than conflating manuscripts.

The manuscript PDF integrity scan is clean (round2/pdf_integrity.txt). Please inspect final-size plots and table layouts as well as their numerical provenance. You may read the compiled PDF with the Read tool, and code/results with Read/Grep/Glob. If command execution is needed but unavailable, explicitly identify the unperformed check. The IBEX mount is read-only; targeted paths include 06_atmospheric/lifeonair_review_r1 and lifeonair_review_r2.

Prioritise: whether the central relative group-1l root depletion is supported; whether the mock and holdout interpretations follow the scripts; remaining dependence/contamination issues; whether the narrowed climate/temporal claims are calibrated; consistency of methods, output tables and manuscript; source-specific literature claims. Distinguish submission-blocking release tasks from scientific errors, and specify actionable fixes.
