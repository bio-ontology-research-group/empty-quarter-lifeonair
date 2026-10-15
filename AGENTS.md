# Paper-specific package guidance

This is the public `empty-quarter-lifeonair` reproducibility package in the `bio-ontology-research-group` organisation. Its dedicated manuscript editing checkout is `~/Documents/papers/empty-quarter-lifeonair/`; the package checkout is `~/Public/software/empty-quarter-lifeonair/`. Keep the manuscript snapshot synchronised with that dedicated checkout. The original combined empty-quarter-mags repository is retained as historical provenance.

Robert explicitly requires bulk data to remain on IBEX. Do not upload large datasets to GitHub, including releases or LFS. Commit code, compact evidence, checksum manifests and institutional fetch instructions. Preserve original source data. CPU analysis belongs in IBEX Slurm allocations, never on login nodes; use ORIX for GPU work.

Distinguish frozen-output verification from upstream raw-data reanalysis. Keep unresolved access, provenance and deposition requirements explicit.
