# Life on Air: reproducibility package

The paper tests how relative atmospheric hydrogen and carbon-monoxide scavenging capacity varies across Rub’ al-Khali soil positions. Its conclusions concern DNA-based gene representation and encoded capacity.

## Verify the evidence and build the paper

Requirements: Python 3 standard library, curl, latexmk and a TeX distribution including natbib and authblk.

```bash
git clone https://github.com/bio-ontology-research-group/empty-quarter-lifeonair.git
cd empty-quarter-lifeonair
bash reproducibility/fetch_large_files.sh
bash reproducibility/run_all.sh
```

The download is a public, SHA-256-verified release asset containing per-sample derived read counts. The offline verifier independently recomputes cohorts and paired contrasts, checks frozen inputs and selected model-output claims, and validates tables/citations. It does not refit the statistical models. The manuscript is built at `lifeonair/manuscript/main.pdf`. No other repository checkout or institutional credentials are needed for the verification commands above.

## Statistical reanalysis

```bash
micromamba create -n eq-lifeonair -f lifeonair/reproducibility/environment-linux-64.lock
bash lifeonair/run_lifeonair.sh
```

Run the full statistical reanalysis in a CPU allocation: the site bootstraps take hours. It regenerates models, tables and figures from the frozen inputs. The default uses retained environmental tables; rebuilding environmental inputs and raw-read screening requires the separately documented source data and databases. See `lifeonair/reproducibility/README.md` and `lifeonair/code/upstream/README.md`.

## Layout and availability

`lifeonair/` retains the historical relative paths for its code, results, figures, manuscript and reviews. This repository is self-contained; the directory name is not a dependency on the old combined repository. `reproducibility/` provides the common public-download and verification entry points. The source and checksums are documented in `MIGRATION.md`.

ENA release confirmation, genome/reference-sequence deposition and an archival DOI remain separate submission requirements. Publishing derived counts and frozen analysis inputs does not constitute full raw-read-to-paper reproduction.

Manuscript editing: https://www.overleaf.com/project/6ad0d1b4338ca7cd7d291188, local `~/Documents/papers/empty-quarter-lifeonair/`. This package is checked out at `~/Public/software/empty-quarter-lifeonair/`.

Proteome companion package: https://github.com/bio-ontology-research-group/empty-quarter-proteome.
