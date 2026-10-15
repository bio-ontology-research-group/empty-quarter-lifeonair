# Literature behind the P3 (life on air) novelty check, 9 October 2026

PDFs were fetched with `fetch_papers.py` (Unpaywall, then EuropePMC); `MANIFEST.tsv` records the
outcome per paper. Where the repository already held a copy in `../../resources/literature/trace-gas-deserts/`, that copy is
used instead of a duplicate (column "Local copy"). Paywalled or blocked items without a copy are kept as
citation + DOI only. On 10 Oct 2026 two blocked open-access papers were added as Europe PMC full-text
XML (PDF download returned HTTP 403). Later the same day Bay 2021 (Nat Microbiol) was downloaded from the
KAUST network and the southern-Africa paper added as Europe PMC XML. The MMBR review, first listed as Greening 2023, is Ray et al. 2023 and was
already in trace-gas-deserts. Li et al. 2023 and Tribbia et al. 2026 were supplied by Robert as PDFs (10 Oct 2026). Every listed paper now has a local copy.

| Key | Citation | DOI / PMCID | Relevance to P3 | Local copy |
|---|---|---|---|---|
| bay_2021_trace-gas-oxidizers-widespread | Bay et al., Nat Microbiol 2021 | 10.1038/s41564-020-00811-w | Trace-gas oxidisers are widespread across soils; most soil taxa encode H2/CO uptake. Sets the "already known" baseline. | this folder (downloaded from the KAUST network, 10 Oct 2026) |
| ji_2017_atmospheric-chemosynthesis-antarctica | Ji et al., Nature 2017 | 10.1038/nature25014 | Atmospheric chemosynthesis in Antarctic desert soil; measured H2/CO uptake rates. | this folder |
| ortiz_2021_antarctic-multiple-energy-group1l | Ortiz et al., PNAS 2021 | 10.1073/pnas.2025322118 | Group 1l [NiFe]-hydrogenase drives H2 oxidation in Antarctic desert soils; encoded by nine phyla. Our 1l dominance is precedented here (cold desert). | ../../resources/literature/trace-gas-deserts/Ortiz_2021_antarctic_energy.pdf |
| jordaan_2020_hot-desert-H2-hydration | Jordaan et al., mSystems 2020 | 10.1128/mSystems.01131-20 / PMC7677003 | Hot deserts (Namib, Gobi, Mojave): H2 oxidisers abundant, uptake stimulated by hydration. Closest hot-desert precedent; single time point. | ../../resources/literature/trace-gas-deserts/Jordaan_2020_hydrogen_desert.pdf |
| leung_2020_desert-energetics | Leung et al., mSystems 2020 | 10.1128/mSystems.00495-19 | Australian desert genome-resolved metagenomics; 1h and 1l dominant; hydration changes expression (1l down). | ../../resources/literature/trace-gas-deserts/Leung_2020_energetic_basis.pdf |
| bay_2021_aridity-gradient-chemosynthesis | Bay et al., ISME J 2021 | 10.1038/s41396-021-01001-0 | Aridity gradient; cell-specific H2 uptake rises with aridity alongside high-affinity hydrogenase abundance. | this folder |
| cordero_2019_atmospheric-CO-widespread | Cordero et al., ISME J 2019 | 10.1038/s41396-019-0479-8 | Form-I CoxL definition used in our typing; CO oxidation widespread. | ../../resources/literature/trace-gas-deserts/Cordero_2019_atmospheric_CO_widespread_ISME.pdf |
| ray_2022_chemosynthesis-cold-deserts | Ray et al., ISME J 2022 | 10.1038/s41396-022-01298-5 / PMC9561532 | Cold-desert MAG survey; most coxL MAGs also carry 1h/1l. Comparable guild-core statistic. | ../../resources/literature/trace-gas-deserts/Ray_2022_atmospheric_chemosynthesis_colddeserts.pdf |
| polloquere_2024_andean-hot-spring-trace-gas | ISME J 2024 | 10.1093/ismejo/wrae062 / PMC11094475 | Spatial (not temporal) shift toward trace-gas oxidisers in an Andean cold desert. | this folder, full text XML (Europe PMC): polloquere_2024_andean-hot-spring-trace-gas.fulltext.xml |
| fertile-islands_2023_trace-gas-oligotrophs | Li, Yang, Wei et al., ISME J 2023;17:1257-1266 | 10.1038/s41396-023-01437-6 / PMC10356767 | Mu Us Desert (semiarid, 180 mm), 60 metagenomes under vs outside Hedysarum shrubs at 0-60 cm: H2- and CO-oxidation genes (TPM) lower under shrubs, strongest at 0-10 cm; H2 genes track SOC and DOC (rho -0.48/-0.43), CO genes only DOC (-0.31); soluble C addition suppresses measured H2/CO uptake; trace-gas oxidisers are oligotrophs. No 1h/1l split. Closest precedent for H1. | this folder |
| reji_2026_moisture-threshold-H2-uptake | Reji, Bertagni, Paulot et al., Nat Commun 2026 | 10.1038/s41467-025-67208-3 / PMC12804893 | Activity down to -70 to -100 MPa; abundance of high-affinity H2 oxidisers predicts uptake; arid regions matter more for the global sink. Supports using genome equivalents as a capacity proxy and raises the stakes of a hot-desert dataset. | ../../resources/literature/trace-gas-deserts/2025_soil_moisture_threshold_H2_uptake_NatComm.pdf |
| ray_2023_clearing-the-air-review | Ray, Tribbia, Cowan & Ferrari, MMBR 2023 | 10.1128/mmbr.00048-23 | Review of atmospheric chemosynthesis; lists open questions (temporal dynamics, hot deserts). | ../../resources/literature/trace-gas-deserts/Ray_2023_atmospheric_chemosynthesis_review.pdf |
| warming_2025_trace-gas-oxidizers-preprint | bioRxiv 2025 | 10.1101/2025.11.18.689010 | Decadal warming: H2/CO uptake rises with cell numbers, per-cell activity falls. Closest "long-term" study; not repeated sampling. | this folder |
| southern-africa_2024_dryland-functional-redundancy | FEMS Microbiol Ecol 2024 | 10.1093/femsec/fiae157 | Dryland coxL spread over many genera, read as functional redundancy; single time point. | this folder, full text XML (Europe PMC): southern-africa_2024_dryland-functional-redundancy.fulltext.xml |
| antarctic_2026_trace-gas-temperature | 2026 | PMC12998438 | Antarctic trace-gas oxidisers active over a wide temperature range. | this folder, full text XML (Europe PMC): antarctic_2026_trace-gas-temperature.fulltext.xml |
| greening_2016_hydrogenase-distribution | Greening et al., ISME J 2016 | 10.1038/ismej.2015.153 | Hydrogenase classification framework and survey. | ../../resources/literature/trace-gas-deserts/Greening_2016_hydrogenase_survey.pdf |
| sondergaard_2016_hyddb | Sondergaard et al., Sci Rep 2016 | 10.1038/srep34212 | HydDB, used for our typing. | ../../resources/literature/trace-gas-deserts/Sondergaard_2016_hyddb.pdf |
| tribbia_2026_namib-fog-aridity-trace-gas | Tribbia, Lebre, ... Cowan, Ferrari, Appl Environ Microbiol 2026;92(8) | 10.1128/aem.00265-26 | Namib 200 km aridity/fog transect, 10 sites, 1-5 cm, single time point: 16S abundance, ATP and respiration fall with aridity, but per-16S atmospheric H2 oxidation is highest in the most arid fog zone (dry soils oxidise below 530 ppbv at 2 of 10 sites); Rubrobacteriaceae enriched there. Activity, amplicon only, no gene typing. Supports the H3 direction; assay on 2 g soil stored at 4 C. | this folder |

## Gap found
No study found pairs repeated multi-year sampling of the same sites with metagenome-based tracking of
trace-gas oxidiser capacity and composition. No trace-gas (or MAG) study from the Arabian Peninsula or the
Rub' al Khali was found; the only Empty Quarter microbiome record is our own amplicon knowledge base
(Zenodo 10.5281/zenodo.23134168).

## Expanded search, review loop round 1 (12 October 2026)

Eight search directions (hot-desert chemosynthesis, enzyme biochemistry, rhizosphere, aridity and
space-for-time, capacity vs activity, methods, environmental drivers, skeptical referee) were run with
Gemini 3.8 Flash through `agy` (`search/run_flash.sh`, prompts in `search/round1/prompts/`). All 122
proposed DOIs were checked against Crossref (`search/process_hits.py`); five mismatched DOIs were
corrected by hand (`search/round1_fix/`). Open-access copies were fetched via Unpaywall, bioRxiv and
Europe PMC; failures were retried from the KAUST network on ws (`search/fetch_ws.py`). Elsevier, ASM,
PNAS, Wiley and Annual Reviews return HTTP 403 to scripted downloads, so those items (about 40) are kept
as citation + DOI only in `MANIFEST.tsv`. The digest for the reviewer is `search/round1/DIGEST.md`;
its one-line relevance notes are Flash summaries and are unverified until the paper is read.
