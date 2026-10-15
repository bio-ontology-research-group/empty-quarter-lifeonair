**I would build the revised paper around differential biomass scaling, lineage turnover, and the distinct root and surface contrasts. Validated carbon-fixation co-carriage is the strongest candidate for a fourth anchor.** The data support a more interesting question than simple depletion: why does group 1l occupy a smaller share of expanding communities, while other atmospheric-gas markers behave differently?

The biomass hypothesis is promising, but “biomass-decoupled” and “do not take part in the larger root-zone community” are too strong. If carrier abundance \(A\) scales with microbial biomass \(B\) as \(A\propto B^\alpha\), then, assuming stable marker copy number,

\[
GE\propto B^{\alpha-1}.
\]

A GE slope of −0.3 therefore implies \(\alpha\approx0.7\): carriers increase sub-proportionally with biomass. An approximately constant carrier population would predict a GE slope near −1. The proposed atmospheric-energy limitation is one explanation for sub-proportional scaling; the existing data can discriminate several alternatives, but establish associations and genomic potential.

I inspected the corrected files and ran read-only Python diagnostics. The cohort reproduced as **546 samples at 59 sites**, including **267 samples with unambiguous positive DNA concentrations**. The R launcher required prohibited temporary-file writes, so the new regression diagnostics below use ordinary least squares with fixed effects; their coefficients are exploratory, without new inferential \(p\)-values.

For paths below, `P` means [`lifeonair/results/`](/home/leechuck/Documents/papers/empty-quarter-mags/lifeonair/results), `C` means the [amplicon cache](/home/leechuck/Documents/papers/empty-quarter-amplicon/repository/cache), and `I` means the [IBEX mount](/home/leechuck/mnt/ibex-metagenomics-ro).

Across the proposed analyses, retain sites 1–60, `id_floor=60`, and `scg_genes_detected>=10`. Use site–trip–compartment medians for paired contrasts; bootstrap sites while retaining their visits. Sample-level models need site and site–visit effects, with trip treated as fixed. Environmental comparisons need spatial sensitivity analyses, and genomic comparisons need phylogenetic control. These remain exploratory hypotheses; cross-validation assesses prediction, rather than retrospectively making them confirmatory.

1. **Group-1l representation scales sub-proportionally with community biomass, with an additional root-associated component.**

   This is the strongest refinement of H1 and the authors’ candidate hypothesis. Its distinctive contribution would be estimating a subgroup-specific scaling relationship across microsites and the landscape, then testing whether it predicts both. Li et al. already connected desert fertile islands, organic carbon and trace-gas oxidisers; Azarbad et al. showed that quantitative and relative rhizosphere profiles can imply different abundance patterns; Tribbia et al. measured low microbial abundance alongside high biomass-normalised hydrogen oxidation in hyperarid Namib soils. The novelty lies in the quantitative, group-resolved connection among these observations. [Li et al., 2023](https://doi.org/10.1038/s41396-023-01437-6), [Azarbad et al., 2022](https://doi.org/10.3389/fmicb.2021.798023), [Tribbia et al., 2026](https://journals.asm.org/doi/10.1128/aem.00265-26)

   Join `P/reads/guild_genome_equivalents_v2.tsv` to `P/inputs/sample_ledger_dna.tsv` through `sample`/`sample_id`, retaining `dna_concentration`, `dna_kit`, `site`, `trip` and `compartment`. Fit each marker using:

   `log2(GE + 0.01) ~ logDNA_within_visit + mean_logDNA_visit + compartment + trip + identifiable_kit_effects + (1|site) + (1|site:trip)`

   Then test `logDNA_within_visit × compartment`. For bulk samples, separately decompose DNA into the site mean and visit deviation. Add `map_era5_mm` and `ndvi_W60d` in small, prespecified models on the **same complete-case cohort**. Compare predictions and coefficients before and after DNA adjustment; attenuation is descriptive evidence about shared variation.

   My sample-level diagnostic with **site–visit fixed effects and position** gave DNA slopes of **−0.281 for 1l**, **−0.0055 for 1h**, **−0.0096 for 2a**, and **−0.0084 for form I**. Restricting to visits with one recorded extraction kit gave **−0.330 for 1l**, across 212 samples. The 1l root coefficient changed from **−2.472 to −1.683** after DNA adjustment. This supports partial explanatory value, with a substantial residual root association.

   Three safeguards are essential:

   - DNA concentration needs extraction mass, elution volume and recovery information before becoming a per-gram biomass estimate. Similar SCG reads per million strengthen the comparison, but also depend on genome size, ploidy and mapping.
   - Independent concentration measurement means a negative GE–DNA slope is not mechanically inevitable. Nevertheless, GE measures a community share, and multiplying GE by DNA then regressing that product on DNA builds the predictor into the response. The resulting slope of \(1+\beta\) supplies no independent confirmation.
   - Classical independent DNA measurement error attenuates slopes towards zero. Kit effects, extraction inhibition and relic DNA can bias them in either direction. Report an errors-in-variables sensitivity over plausible reliability values; biological replicates alone do not identify extraction-error variance.

   A particularly useful check already passed provisionally: slopes for **un-pseudocounted 1l/1h, 1l/2a and 1l/form-I ratios** were **−0.320, −0.319 and −0.312**, respectively, under visit and position fixed effects. These ratios cancel the common GE denominator at sample level, subject to rounding.

   Reject a *single biomass-only explanation* if matched-cohort predictions retain a substantial root effect or require different slopes across positions. Reject subgroup-specific scaling if calibrated or same-kit analyses give equivalent slopes across markers. Broad intervals remain inconclusive. Effort: **60–120 local minutes** for models and sensitivity analyses; no IBEX work initially.

2. **The aridity association arises mainly through turnover of carrier lineages, with a separately testable within-lineage component.**

   This refines H3 into an ecological assembly question. Bay et al. already reported increasing group-1l representation along a desert aridity gradient. Peng et al. subsequently found extensive hydrogenase transfer and selection in Antarctic communities, making the distinction between lineage sorting and trait-associated ecological sorting especially relevant. The advance would be quantifying where the association resides in the taxonomic hierarchy. [Bay et al., 2021](https://doi.org/10.1038/s41396-021-01001-0), [Peng et al., 2026](https://www.nature.com/articles/s41467-026-76323-8)

   Use `P/inputs/atmo_classification.tsv` (`rep_id`, `assignment`, `protein_id`), the R232 summaries (`user_genome`, `classification`), and `coverm_rep_abundance_long.tsv.gz`, whose three headerless columns are sample, representative and abundance. Restrict taxonomy to the canonical representative set: the R232 input files contain additional genomes.

   For each genus, calculate its community abundance \(w_g\) and the abundance fraction belonging to carriers \(q_g\). Overall carrier representation is \(\sum_g w_gq_g\). Decompose paired root/bulk differences into changes in genus abundance and changes in within-genus carrier representation, explicitly handling genera absent from one position.

   For sufficiently observed mixed genera, model within-genus carrier fractions using a fractional-response mixed model:

   `carrier_fraction ~ logMAP + compartment + trip + completeness_summary + (1|site) + (1|site:trip) + (1 + logMAP|genus)`

   CoverM percentages are continuous abundance estimates, so do not treat them as binomial read counts. Add matched comparisons among close relatives and leave-family-out analyses. Within-genus sorting is still sorting among genomes; the term “within-lineage” must specify that resolution.

   The feasibility check is encouraging: after excluding the two flagged orphan representatives, **82 R232 genera contain carriers and non-carriers**, comprising **547 representatives, including 180 carriers**. Both types were observed in **71 genera across bulk visits**. However, catalogue-mixed genera account for a median **23.6% of bulk mapped 1l-carrier abundance**, limiting how broadly the conditional result can generalise.

   A substantial, reproducible conditional aridity coefficient across independent clades would contradict turnover-only explanation. Conversely, precise near-zero within-genus effects plus successful prediction from genus weights would favour turnover. Pitfalls include missing genes in incomplete MAGs, catalogue recovery bias, ambiguous read assignment, and overadjustment using community axes containing the focal organisms. Effort: **60–180 local minutes**; **2–8 IBEX hours** if a suitable species tree or close-relative matching needs construction.

3. **Root proximity and surface exposure impose different ecological filters on atmospheric-gas guilds.**

   The revised H2 should predict a strong group-1l-specific root contrast and a broader surface contrast involving 1l, 1h and form I. This creates a mechanistic comparison within the same transect: organic-resource-associated community expansion near roots versus exposure-associated sorting at the surface.

   Li et al. provide the root-resource precedent, while Imminger et al. demonstrate rapid desert-community resuscitation and hydrogenase responses during hydration. Reji et al. show that hydrogen uptake depends on moisture and soil physical conditions. A useful extension would connect the distinct **DNA-level habitat distributions** to independently annotated stress traits. [Imminger et al., 2024](https://doi.org/10.1038/s41467-024-46920-6), [Reji et al., 2026](https://www.nature.com/articles/s41467-025-67208-3)

   First fit a joint marker-by-position model using the four focal markers and the two reference marker groups, allowing correlated responses within sample. Estimate both root/bulk and surface/bulk contrasts, including their differences across markers. Treat “other CoxL” and “other NiFe” as reference groups with their own ecology.

   Next use `P/inputs/env_site_trip.tsv` (`smap_sm_W30d`, `ndvi_W60d`), `C/per_trip_site_temperature.tsv` (`T_d30`, `T_d90`, `T_d365`), and `I/10_thermal/results/mechanism_copy_number_by_family.tsv`. Test paired surface/bulk contrasts against one exposure predictor at a time:

   `surface_bulk_log2FC ~ temperature_or_moisture + trip + (1|site)`

   Independently ask whether surface-enriched genomes possess validated DNA-repair, heat-protection or compatible-solute traits, accounting for phylogeny and completeness. Compare those relationships with root enrichment.

   The authors’ corrected contrasts already support the descriptive premise: surface excesses of **−0.83, −0.42 and −0.24 log₂** for 1l, 1h and form I, versus root excesses of **−1.83, −0.01 and −0.05**.

   The shared surface-filter hypothesis weakens if the broad contrast disappears after lineage, copy-number or biomass checks. The specific stress interpretation weakens if neither exposure interactions nor validated trait associations occur. Satellite moisture and air temperature have coarse spatial support; existing observations offer limited separation of UV, temperature, desiccation and phototrophic competition. Effort: **60–180 local minutes**, plus **2–8 IBEX hours** for focused trait validation.

4. **Aridity preferentially enriches the subset of group-1l carriers with carbon-fixation potential.**

   This could separate hydrogen-supported maintenance from a candidate hydrogen-supported carbon-acquisition strategy. Atmospheric chemosynthesis itself is established: Ji, Bay and Ray and colleagues already combined genomic and experimental evidence. The new question is whether carbon-fixation co-carriage explains the distinctive biomass and habitat associations of group 1l within this hot desert. [Ji et al., 2017](https://doi.org/10.1038/nature25014), [Ray et al., 2022](https://doi.org/10.1038/s41396-022-01298-5)

   A quick screen of the **current, retained-contig-aware** `I/04_functional/repko_clean/results/rep_ko_v2.tsv` found both `K01601` and `K00855` in:

   - **150/597 group-1l carriers: 25.1%**
   - **43/410 group-1h carriers: 10.5%**
   - **35/136 group-2a carriers: 25.7%**

   These are candidate RuBisCO–phosphoribulokinase co-occurrences. They justify follow-up. KEGG identifies K00855 as phosphoribulokinase, with additional orthology represented by K28462 that should enter the annotation audit. [KEGG phosphoribulokinase](https://www.kegg.jp/entry/2.7.1.19)

   Validate RuBisCO form, catalytic residues, phosphoribulokinase identity and the remaining CBB pathway, allowing recognised pathway variants. Check contig membership, contamination and phototrophic machinery. A two-KO rule alone is insufficient for declaring a complete pathway.

   Partition CoverM abundance into validated `1l+CBB`, `1l without detected CBB`, and corresponding 1h/2a categories. Test both total community representation and the CBB-positive fraction **within each hydrogenase guild**:

   `guild_subtype_abundance ~ subtype × (logMAP + compartment + logDNA) + trip + (1|site) + (1|site:trip)`

   Use the full cohort for environmental tests and the DNA subset separately. At genome level, test CBB co-carriage against hydrogenase group with phylogeny and completeness.

   Reject preferential chemosynthetic-potential enrichment if validation removes the co-carriage difference, or if CBB-positive and CBB-negative 1l carriers show equivalent habitat and aridity associations. Missing pathways in incomplete MAGs require an “undetermined” category.

   The old `carbon_fixation/` table contains **6,010 older MAGs with zero direct identifier matches** to current representatives. Its accompanying FBA analysis reports no biomass-objective autotrophic growth and uses gap-filled models. Neither is ready to supply the revised headline. Effort: **30–60 local minutes** for screening and joins; **4–12 IBEX hours** for sequence validation and pathway reconstruction.

5. **A combination of genomic traits, beyond predicted growth potential, distinguishes 1l from 1h carriers and predicts root enrichment.**

   The defensible hypothesis is that carbon-acquisition flexibility, maintenance/resuscitation machinery, or stress protection distinguishes the groups despite similar predicted maximal growth. A generic “scavengers are oligotrophs” finding would add little: López et al. established growth potential as a strong rhizosphere predictor, and **Wang et al. (2025), absent from the supplied manifest, already linked aridity, smaller genomes, slower predicted growth and trace-gas oxidation**. [López et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10432406/), [Wang et al., 2025](https://onlinelibrary.wiley.com/doi/10.1111/gcb.70498)

   Join `rep_growth.tsv` (`rep`, `d_hours`, `n_ribo`, `lower`, `upper`), corrected carrier assignments, current KO content, and `I/10_thermal/index/genome_index.tsv` (`genome_id`, `nucleotide`, `completeness`, `contamination`). Obtain genome length from retained nucleotide sequences. Start with a small trait set: genome size, growth potential, validated CBB, corrected form-I co-carriage, resuscitation, and coherent sporulation/compatible-solute modules.

   Compare 1l-only, 1h-only, dual carriers and non-carriers with phylogenetic models. Then predict per-genome paired root/bulk enrichment:

   `root_log2FC ~ carrier_group + log(d_hours) + genome_size + selected_traits + completeness + trip + (1|site) + phylogenetic_genome_effect`

   Use abundance/detection-aware uncertainty and leave-clade-out prediction. Avoid the previous genus-summary regression with many zero changes.

   Existing diagnostics already narrow the possibilities: carrier median predicted doubling times are **4.20 h for 1l and 4.15 h for 1h**. Corrected form-I co-carriage is also similar: **64/597 versus 46/410**, or **10.7% versus 11.2%**. Growth rate or CO carriage alone therefore looks unpromising as the distinction.

   Reject a proposed trait explanation if differences disappear within clades or fail to improve held-out-clade prediction. The thermal panel needs semantic curation: its aggregate “ectoine” column includes broad aminotransferase and acetyltransferase families. Likewise, SpoVG presence alone is insufficient to classify sporulation. `08_dormancy` is empty, but thermal results contain useful starting annotations. Effort: **90–240 local minutes**, plus **4–12 IBEX hours** for focused annotation repair.

6. **Changes in hydrogenase copies per carrier amplify the habitat contrast.**

   This alternative matters because GE measures marker copies per genome equivalent. It offers an ecological extension of subgroup distribution studies: communities may differ in both carrier prevalence and genomic investment in hydrogen metabolism. [Greening et al., 2016](https://pubmed.ncbi.nlm.nih.gov/26405831/)

   In the corrected carrier set, **1,121 annotated group-1l proteins occur in 597 representatives**; **300 representatives have multiple proteins**, with a maximum of seven. For 1h, **445 proteins occur in 410 representatives**, with 32 multiple-protein representatives.

   For sample \(s\), use CoverM abundance \(a_{is}\) and validated copy count \(c_i\) to calculate carrier representation \(F_s\), mean copies per carrier \(C_s\), and projected marker abundance \(M_s=F_sC_s\). Fit the same habitat and environmental models to all three. Validate copy counts using `atmo_classification.tsv`, `P/operon/marker_context.tsv` (`contig`, `position`, `genes_to_contig_end`), and sequence lengths; adjacent fragments or duplicated assembly segments need reconciliation.

   The quick check already argues against a dosage-only explanation. Across **105 root/bulk visits at 47 sites**, median changes were **−2.17 log₂ for mapped-community carrier representation**, **−2.30 for copy-weighted representation**, and **−0.20 for mean copies per carrier**. The separately computed medians need not add exactly. Surface values were **−0.71, −0.89 and −0.18**, respectively.

   Thus, dosage may amplify the contrast, while carrier representation supplies most of the observed MAG-based difference. Reject a substantial dosage explanation if validated copy number stays nearly constant and copy-normalised associations persist. Catalogue incompleteness and shared mapping remain important. Effort: **30–90 local minutes**, plus **2–6 IBEX hours** for validating suspect loci. This is an early validation requirement even if it remains supplementary.

7. **Mineralogy and salt-associated geochemistry explain part of the apparent aridity association.**

   This is a serious competing explanation for H3 because precipitation, substrate and salt accumulation covary geographically. Monteiro et al. showed the value of geochemically defined desert gradients; the southern African dryland study found substantial salinity-associated community structure. The P3 extension would identify whether atmospheric-gas traits follow climate or substrate after accounting for their covariance. [Monteiro et al., 2022](https://doi.org/10.3389/fmicb.2021.783767), [Functional redundancy in southern African drylands, 2024](https://doi.org/10.1093/femsec/fiae157)

   Use **`C/xrf_site_compartment_panel.tsv`**, which contains measurements by site and compartment. `xrf_per_compartment.tsv` contains previously computed correlations and is unsuitable as the measurement table.

   Start with bulk-site medians and a small predictor set: `logMAP`, an independently constructed lithology axis, and a salt-associated axis. Use compositional transformations or justified element ratios, explicitly handling below-detection zeros. Include `C/topography.tsv` (`elev_m`, `slope_deg`, `twi`) only in targeted alternatives. Compare models with spatially blocked prediction. For root/bulk contrasts, use paired compartment differences in XRF rather than treating site means as visit measurements.

   Across **51 bulk sites**, I reproduced **ρ(1l, MAP)=−0.512**. Associations with individual XRF elements were weaker: **Fe +0.018, Al +0.091, Ca −0.110, Cl −0.202, S −0.199**. MAP correlated with Ca, Cl and S at **+0.558, +0.542 and +0.566**. Partial rank correlation between 1l and MAP controlling those three elements remained **−0.528**.

   These checks currently favour retaining the climate association. Reject geochemistry as its principal explanation if prespecified geochemical adjustment leaves the association and spatial prediction intact. XRF measures total elements, with limited inference about bioavailable nickel, soluble salts or osmotic stress. Effort: **45–120 local minutes**; no IBEX computation required.

8. **Preservation of extracellular or damaged-cell DNA strengthens the apparent dry-site enrichment of 1l-associated lineages.**

   This directly challenges the biological interpretation of H3 and the DNA-yield scaling. Carini et al. demonstrated that relic-DNA removal can affect detection of temporal community dynamics. The useful extension would establish whether atmospheric-gas carrier associations are especially sensitive to preservation. [Carini et al., 2020](https://pubmed.ncbi.nlm.nih.gov/31964728/)

   The available PMA evidence has narrow support: **nine pairs from two sites during trip 5, comprising six rhizosphere and three surface pairs, with no bulk pairs**. Existing QC reports a median pair–pair correlation of **0.397**. The table `results/relic_climate/priorfree/arm_b_measured_pma.tsv` projects PMA-associated ASV categories across the transect; its “living” and “relic” labels require more cautious terminology.

   Begin with `C/test6_disconfirmation/qc_per_asv_pair_ratios.parquet` (`asv_id`, `pair_id`, `site`, `comp_code`, `T`, `UT`, `log_ratio`) and `per_eq_asv_viability.tsv` (`eq_asv`, `weighted_median_ratio`, `max_pid`). Revalidate the chain from PMA sequences to EQ ASVs to MAGs using `C/mag_16s/asv_to_mag_match.tsv`. The current PMA-to-EQ mapping accepts **97% identity**, which is too broad to assume identical hydrogenase carriage.

   Among unambiguously linked carrier/non-carrier groups, compare treatment-associated changes within actual pairs, treating pair and ASV as repeated factors and the two sites explicitly. Then perform sensitivity bounds on transect carrier abundance under different preservation weights. Report how much observed carrier abundance has a defensible PMA link.

   The existing projected-ASV result is a reason to investigate: precipitation-associated diversity has **ρ=−0.315** for the PMA-depleted group and **−0.012** for the PMA-retained group. These are ASV-diversity results, not carrier-abundance results.

   Consistent retention of well-linked carriers and persistence of associations across plausible preservation bounds would weaken the relic explanation. Poor linkage or wide uncertainty leaves it unresolved. PMA retention reflects membrane exclusion and treatment accessibility; relative T/UT ratios also carry compositional effects. Effort: **60–180 local minutes**, plus **1–4 IBEX hours** if representative 16S recovery is needed. This is supporting evidence with limited scope for becoming an anchor.

9. **Spatial habitat sorting coexists with weak short-term change in genomic guild composition.**

   This is the defensible temporal refinement of H3. Jordaan et al. and Imminger et al. observed strong physiological responses to hydration; persistent genomic representation can coexist with such responses. The contribution would be a bounded estimate of compositional change over the sampled interval. [Jordaan et al., 2020](https://doi.org/10.1128/mSystems.01131-20), [Imminger et al., 2024](https://doi.org/10.1038/s41467-024-46920-6)

   Use bulk visit medians and `env_site_trip.tsv` (`P90_era5`, `P365_era5`, `smap_sm_W30d`, `trip_date`). For each exposure \(x\), fit:

   `log2(GE + 0.01) ~ site_mean_x + within_site_x + trip + (1|site)`

   Add site-specific temporal slopes only if supported by the number of visits. Repeat for carrier representation and validated guild subtypes. Choose one primary exposure window before examining the new models; assess the others as sensitivity analyses.

   I confirmed **71 bulk visits at 27 repeatedly sampled sites**, spanning **17 March 2023–20 August 2024**, for the existing P365 analysis. Trip 5 remains missing. Four shared expeditions provide limited temporal replication.

   Define “weak” through an ecological equivalence bound, for example a change smaller than **0.5 log₂ units across the observed within-site exposure IQR**, if the authors consider that meaningful. Evidence for weak response requires the interval to lie inside the bound. A reproducible larger within-site effect falsifies it; an interval spanning both weak and substantial changes is inconclusive.

   Compare spatial and temporal slopes only for commensurate exposures and scaling. MAP and antecedent P365 answer different questions. Effort: **45–120 local minutes**; existing data suffice for the restricted interval.

I would execute the dosage check and annotation audit first, then develop **1, 2 and 3 as the core argument**, with **4 as the conditional fourth anchor**. The current CBB screen gives that fourth direction a concrete starting point. Trait analysis should explain these contrasts, while geochemistry and relic DNA provide substantive alternative-explanation tests. The temporal result belongs in a secondary position unless its confidence bounds become informative.

If those anchors hold, a suitable title is:

**“Distinct biomass scaling and habitat associations of atmospheric hydrogenase guilds in hyperarid desert soils”**

That title requires a defensible biomass calibration. With extract concentration alone, use:

**“Distinct habitat associations of atmospheric hydrogenase guilds across a hyperarid desert”**

An abstract claim could then read:

> Across 546 metagenomes from 59 Rub’ al-Khali sites, group-1l hydrogenase genes showed a pronounced root-associated reduction in community representation, while surface soils showed broader reductions across atmospheric hydrogen and carbon monoxide oxidation markers. Group-1l representation declined with DNA yield, consistent with sub-proportional scaling during community expansion. Genome-resolved analyses distinguished contributions from lineage turnover, marker copy number and carbon-fixation potential. Together, these findings identify contrasting ecological distributions among genomic strategies for atmospheric energy acquisition.

After completing the tests, replace the generic genome-resolved sentence with the actual supported direction. Claims about atmospheric-energy limitation or carbon flux would require additional physiological evidence.

The ranking below concerns contribution to the paper’s argument; validation tasks can deserve earlier execution.

| Rank | Hypothesis | Evidence now | Recommended role | Estimated additional effort |
|---:|---|---|---|---|
| 1 | Sub-proportional biomass scaling with residual root association | 1l slope −0.281 under visit fixed effects; other markers near zero | Core anchor, with proxy qualification | 60–120 local min |
| 2 | Lineage turnover versus within-lineage sorting | 82 mixed genera; conditional analysis covers a limited abundance fraction | Core anchor | 60–180 local min; 2–8 IBEX h if tree needed |
| 3 | Distinct root and surface filters | Corrected marker contrasts strongly motivate it | Core anchor; stress mechanism conditional | 60–180 local min; 2–8 IBEX h |
| 4 | Preferential enrichment of 1l+CBB potential | Candidate co-carriage: 25.1% of 1l versus 10.5% of 1h carriers | Conditional fourth anchor | 30–60 local min; 4–12 IBEX h |
| 5 | Genomic traits beyond growth potential | Similar growth predictions and CO co-carriage narrow the possibilities | Mechanistic support | 90–240 local min; 4–12 IBEX h |
| 6 | Copy-number amplification of contrasts | Root carrier change −2.17; copies/carrier change −0.20 log₂ | Early validation; likely secondary result | 30–90 local min; 2–6 IBEX h |
| 7 | Geochemistry as an alternative to aridity | Initial XRF adjustment preserves MAP association | Competing-explanation test | 45–120 local min |
| 8 | Relic-DNA preservation bias | Relevant companion pattern, but only two PMA sites and no bulk pairs | Bounded sensitivity analysis | 60–180 local min; 1–4 IBEX h if needed |
| 9 | Weak short-term genomic change despite spatial sorting | 71 visits, 27 sites, four expeditions | Secondary H3 refinement | 45–120 local min |