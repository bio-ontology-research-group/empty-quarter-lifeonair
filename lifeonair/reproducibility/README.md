# P3 evidence verification and reproduction

From the repository root, `bash lifeonair/run_lifeonair.sh --verify` performs offline checks of frozen evidence. It independently rebuilds the 546-sample/59-site cohort and all paired-position median contrasts from sample metadata and frozen GE. It checks selected reported model estimates against their source tables, annotation/carrier closure, unique environment keys, input hashes, table columns and citation keys. `claims_p3.tsv` labels model estimates TABLE: these checks do not refit models or validate every numerical sentence. The proteome registry is in its separate empty-quarter-proteome repository.

`inputs.sha256` freezes local inputs. `environment-linux-64.lock` is the explicit Linux conda package lock from the environment that ran the round-2 analyses; create it with `micromamba create -n eq-lifeonair -f lifeonair/reproducibility/environment-linux-64.lock`. `R-session.txt` records runtime versions. Python scripts use Python 3 standard libraries unless their headers state otherwise. LaTeX build requires latexmk and a TeX distribution.

For statistical reanalysis, run `bash lifeonair/run_lifeonair.sh` from the repository root in a CPU allocation: the site bootstraps take hours. The script regenerates GE, annotation summaries, models, tables, figures and PDF. It requires the per-sample counts retrieved by `bash reproducibility/fetch_large_files.sh`. The default uses frozen environmental tables. Rebuilding those tables with `P3_REBUILD_ENV=1` additionally requires the companion amplicon cache at the path in `code/prep_environment.R`; its source hashes are in `results/inputs/env_sources.sha256`.

Cluster input derivations and job IDs are described in `code/upstream/README.md`, the benchmark scripts and `reviews/round2/response_round2.md`. Some upstream scripts retain historical workstation/IBEX paths and require their documented reference databases and source sequences. This is an institutional reconstruction route, not a public raw-read-to-paper workflow.

## Open release requirements

The manuscript is a draft. A DOI archive, reference-sequence deposition, ENA release/accession confirmation and final publication metadata remain required before submission. Derived counts are available through the checksum-verified public release route. The repository and derived-count download are public. Full raw-read-to-paper reproducibility remains distinct from a passing frozen-evidence check. R2-25 remains partially addressed.

## Round-3 sensitivity

The runner also regenerates GE with the 387-ID union, fits the targeted genus/trait/CBB models and recomputes paired contrasts and residual benchmark homology. These targeted models report model SEs, without fresh bootstrap intervals. `code/import_trip4_controls.py` is an optional provenance import requiring the companion data-paper repository; ordinary reanalysis uses the committed candidate sequences, calls and union. Their source checksum is in `results/review_r3/trip4_provenance.json`. The verifier independently recomputes the added paired medians, benchmark audit and exclusion-set closure. It checks selected R3 model estimates as TABLE claims.
