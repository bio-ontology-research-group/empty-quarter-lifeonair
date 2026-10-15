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
- Figures: `REPO/lifeonair/figures/` (fig_h123, fig_mechanism, figS_coxL_tree)
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
