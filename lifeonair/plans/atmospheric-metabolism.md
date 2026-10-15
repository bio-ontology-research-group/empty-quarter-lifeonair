> **⚠ SUPERSEDED (2026-08-03).** Pre-pivot / pre-merge plan, kept for history. Current plan → `../notes/PLAN.md`; current narrative → `../notes/PAPER_STRUCTURE.md`; current numbers → `../notes/CATALOGUE_NUMBERS.md`.

# Atmospheric trace-gas metabolism — phylogenetic confirmation & integration

> **Historical planning record.** Status and headline language below may be superseded. Use
> `manuscript/main.tex`, `handover.md`, and `reproducibility/claims.tsv` for current evidence.

**Goal.** Turn the first-pass "life on air" marker census (coxL/hydB present by KO) into a
*defensible, citable* claim that the dominant EQ MAGs encode **high-affinity, atmospheric**
trace-gas-oxidising enzymes — the energy basis for life in the most carbon-starved biome on Earth.
KO presence cannot distinguish atmospheric (high-affinity) from non-atmospheric enzymes; clade/form
assignment can. This is the load-bearing analysis for the MAGs-paper headline.

## Why KO is not enough
- `K03520` (coxL) labels *all* Mo-CODH large subunits — but only **form I** CoxL is a bona fide
  CO oxidiser; **form II** is divergent and not confirmed to oxidise CO.
- `K06281` (hydB) labels generic [NiFe]-hydrogenase large subunits across ≥29 classes. Only
  **groups 1h and 1l** are the high-affinity, O₂-tolerant enzymes that scavenge **atmospheric H₂**.
  (KO `K18016/K18017/K18023` returned 0 hits, confirming KO granularity is insufficient.)

## Inputs already on disk
| Resource | Path (IBEX `R=/ibex/scratch/projects/c2014/rob/metagenomics`) |
|---|---|
| Marker hits (mag, protein, KO) | `R/06_metabolism/marker_hits.tsv` — coxL 13,675 / 5,801 MAGs; hydB 6,602 / 3,665 MAGs |
| Per-MAG proteins / genes / gff | `R/06_mobilome/mag_index.tsv` (mag_id, faa, fna, gff3) |
| GTDB-Tk r226 taxonomy | joined by `mag_id` (= `sample__bin`) |
| gapsmith FBA models (draft + gap-filled SBML) | `R/06_darkmatter/gapsmith_models/<mag>/<mag>-{draft,filled}.xml` |
| gapseq drafts | `R/05_metabolism/gapseq_out` |
| Dark-matter / structurally-novel protein catalog + dossiers | `R/06_darkmatter/`, workbench `results/darkmatter_function_assignment/` |
| 16S amplicon community (1,237 samples) | companion repo `~/Documents/papers/empty-quarter-amplicon/` |

## Reference methods (field-standard, citable)
- **CoxL form I** — Cordero et al. 2019, *ISME J* (PMC6794299). Operational definition: **>55% identity
  & 90% query coverage** to three RefSeq form-I references — *Oligotropha carboxidovorans*
  `WP_013913730.1`, *Mycobacterium smegmatis* `WP_003892166.1`, *Natronorubrum bangense*
  `WP_006067999.1`. Active-site motif **AYRGAGR** distinguishes form I; we derive the motif column
  *empirically* from the references (avoids literature direction ambiguity) and report each query's state.
- **[NiFe]-hydrogenase class** — Greening lab **HydDB** (Søndergaard et al. 2016, *Sci Rep* srep34212;
  repo `github.com/GreeningLab/HydDB`). Run `hydclassify.sh` (DIAMOND best-hit vs `HydDB_all_hydrogenases.faa`
  + HMM `HydDB_all_MM2022.hmm`) → per-protein class. **Group 1h + 1l = high-affinity atmospheric H₂ uptake.**

## Pipeline (scripts in workbench `scripts/ibex/`, all via Slurm)
- **82a_atmo_refs.sh** — clone HydDB, build DIAMOND db; fetch CoxL form-I refs (NCBI efetch), build db. `REFS_DONE`.
- **82b_atmo_extract.sh** — extract coxL/hydB(/coxM/coxS) query proteins from per-MAG faa → `queries/*.faa` (`mag_id|protein_id`), lengths.
- **82c_atmo_classify.sh** — hydB→HydDB class (flag 1h/1l); coxL→DIAMOND vs form-I refs + Cordero threshold + empirical motif; mafft+FastTree(refs+dereplicated queries) for the figure. → `atmo_classification.tsv`.
- **82d_atmo_integrate.(py|sh)** — Stage D below.

## Stage D — integration (the "mechanism, not catalog" layer)
1. **Per-MAG guild table** joined to GTDB tax + compartment + trip: confirmed atmospheric H₂ oxidiser
   (1h/1l), confirmed CO oxidiser (form I), both, or neither. Length-filter fragments.
2. **gapsmith/gapseq cross-check** — does the MAG's FBA model contain the CO/H₂-oxidation reactions, and
   are they feasible (flux-carrying) in the gap-filled model? Confirms metabolic plausibility independent
   of sequence. Reconcile id token (`<sample>__{single,self}_SemiBin_N` ↔ `sample__bin`).
3. **Genomic context (money shot)** — from gff, extract the operon neighbourhood (±N genes) of each
   confirmed coxL/hydB. Check for (a) maturation/accessory genes (cox operon `coxMSG`/`coxDEFG`, hyp genes),
   (b) **co-located dark-matter / structurally-novel proteins** → candidate novel-fold members of the gas
   guild. A structurally-novel family recurrently embedded in confirmed gas-uptake operons is the figure
   that makes this *Nature*, not a catalog.
4. **16S abundance weighting** — map confirmed-guild taxonomy to amplicon relative abundance across the
   1,237 samples → ecosystem-level prevalence ("X% of the community can live on air").

## Deliverables → `empty-quarter-mags/results/` + `figures/`
- `atmo_classification.tsv`, `atmo_guild_per_mag.tsv`, `atmo_guild_by_phylum.tsv`
- Figs: CoxL form-I tree, HydDB class distribution (1h/1l highlighted), guild-by-phylum, guild abundance vs 16S, novel-fold-in-operon panel.

## Caveats to state in the paper
- Identity/HMM clade assignment is sequence-based; flag fragments (<X aa) as low-confidence.
- 1h/1l assignment confirms *capacity*, not in-situ activity — gapseq feasibility + abundance strengthen but do not prove it.
