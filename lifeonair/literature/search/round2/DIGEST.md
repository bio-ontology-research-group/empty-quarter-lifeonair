# Round-2 literature digest (P3 life on air), 12 October 2026

Scope: Astra finding R2-24 and related round-2 gaps. Four Flash directions (`gemini-3.8-flash-high` via `agy
--sandbox`, prompts in `prompts/`, raw outputs in `out/`): d1 genome traits along aridity, d2 carbon regulation
of high-affinity hydrogenases/CODH, d3 validation of read-based gene quantification, d4 CoxL form I/II. d2 failed
on the first run (headless permission denial) and was rerun with "web search only" appended.
`process_hits.py round2`: 51 unique DOIs, 50 Crossref-verified (Dobbek 2002 failed the title check and is
dropped), 30 new, 20 downloaded; `fetch_ws.py` on ws-ts recovered 1 more (Nishimura 2010). All rows are
in `candidates.tsv` and `../../MANIFEST.tsv`; PDFs/XML are in `lifeonair/literature/`. Citation-only, with no
local copy (publisher 403 from KAUST, no PMC PDF): Berney 2014, Cunliffe 2011, King 2003, Lalonde 2016,
Liu 2023 iMeta, Lorite 2000, Myers 2020, Pelzmann 2009, Weber 2010.

Verification standard: "verified" findings below were read from the local full text (PDF via pdftotext, or
Europe PMC XML) unless marked *abstract only*. Flash relevance notes were **not** relied on; two of them
are wrong (see "Flash errors").

Proposed BibTeX (Crossref) for the papers recommended for citation: `new_refs.bib` (14 entries, keys given
below). `references.bib` was not edited.

---

## (a) Wang et al. 2025, Global Change Biology — focused comparison

**Citation.** Wang X, Wang W, Deng L, Li T, Lei S, Zhang L, Liao L, Song Z, Liu G, Zhang C (2025) Shifts in
genome size and energy utilization strategies sustain microbial functions along an aridity gradient. *Glob
Change Biol* 31(9). DOI 10.1111/gcb.70498. **Already in references.bib as `wang2025genomesize`** (same DOI)
but not cited in main.tex.

**Local copy.** Full text NOT obtained. Not open access (Unpaywall `closed`, no PMC). From ws-ts both Wiley
PDF routes (`pdfdirect`, `pdf`) return a Cloudflare "Just a moment..." challenge page (5.6 kB HTML), and
the browser extension was not connected. Saved instead: the Europe PMC record with the full abstract,
`lifeonair/literature/wang_2025_shifts-genome-size-energy.epmc_abstract.json` (PMID 40958540). MANIFEST row
updated to `ABSTRACT ONLY`. **Everything below is from the abstract**: the statistical design (phylogenetic
control, completeness correction, how "trace-gas oxidizer" was defined, which hydrogenase groups) cannot be
checked until someone with Wiley access downloads the PDF.

**What the abstract states.** 374 MAGs from 13 phyla along a semi-humid forest → arid desert gradient (China).
(i) Communities in more arid soils have smaller genomes, lower GC and fewer 16S rRNA copies. (ii) Reliance on
trace-gas oxidation (H2, CO, CH4, H2S) increases, and organic-energy use decreases, with aridity; ex situ
H2/CO/CH4 consumption was higher in desert soils. (iii) The shift is correlated with declining soil organic C;
trace-gas oxidizers (specialised and multi-gas) increase as organic C falls, non-oxidizers decline.
(iv) **Trace-gas oxidizers have smaller genomes, fewer rRNA operons and slower predicted growth**
(oligotrophs); copiotrophic non-oxidizers have larger genomes and faster growth.

**Comparison with our genome-trait results** (`results/a5_genome_traits/`; within-order models over 5,438
representatives):

| Trait | Wang 2025 (abstract) | Ours | Agreement |
|---|---|---|---|
| Predicted growth | oxidizers slower | 1l +0.12, 1h +0.12, 2a +0.20 log doubling time vs non-carriers within order; form I +0.02 (ns) | agrees for hydrogenase carriers; not for form-I CoxL |
| Genome size | oxidizers smaller | 1l −1.9% (q 0.11, ns); 1h +7%, 2a +11%, form I +17% (all q<1e-5) | **disagrees**: no oxidizer group is smaller within order, and three are larger |
| GC | arid communities lower GC | 1l +0.67 pp (q 0.001), form I +0.85 pp; 1h, 2a ns | different level (community vs carrier); ours is the opposite sign for 1l |
| Aridity | oxidizers increase with aridity | bulk 1l vs MAP rho −0.37, not robust to spatial blocks; surface all groups robust (ERA5) | consistent in direction only |

Likely reasons for the genome-size difference, all of which need the full text to confirm: (1) our contrast is
**within order**, whereas the abstract's comparison is across MAGs, where small-genome Actinomycetota
dominate the oxidizers (in our catalogue the 1l-only carriers have a median of 3.78 Mb vs 3.97 Mb for
non-carriers *without* phylogenetic control, so the raw medians agree with Wang); (2) we separate hydrogenase
groups and CoxL form, whereas "trace-gas oxidizer" in Wang pools them and includes CH4/H2S; (3) our genomes
are completeness-corrected, theirs unknown. The GC difference is at different levels (community mean vs
carrier vs non-carrier) and is not a contradiction. Chuckran 2023 (below) links high GC with small genomes
under low C:N and high pH, which fits our 1l result better than Wang's community-level low GC.

**Manuscript placement.** Results §"Genera carrying group 1l ..." para. 3 (main.tex ~351–356, "Group 1l
carriers also differed in genome traits ...") and the Discussion sentence at ~422 ("Group 1l carriers are
mostly Actinomycetota (59%) and Chloroflexota (15%) with high GC content ..."). Suggested Discussion addition:
> Along a Chinese aridity gradient, trace-gas oxidizers had smaller genomes and slower predicted growth than
> non-oxidizers \citep{wang2025genomesize}. Within orders, our carriers had slower predicted growth but not
> smaller genomes (group 1l), or larger ones (groups 1h and 2a and form-I CoxL), so the small genomes of
> oxidizer communities may reflect which lineages carry the enzymes rather than a genome-size cost of
> carrying them.

Hedge it as an abstract-level comparison until the full text has been checked.

---

## (b) Baril & Constant 2023 and Kropp et al. 2026 — what they actually show

### Baril X, Constant P (2023) FEMS Microbiol Ecol 99(12): fiad159

**Title.** Carbon amendments in soil microcosms induce uneven response on H2 oxidation activity and microbial
community composition. DOI 10.1093/femsec/fiad159. Local copy: `lifeonair/literature/baril_2023_carbon-amendments-microcosms-induce.fulltext.xml`
(Europe PMC PMC10716739). Key `baril2023carbon`.

Verified from the full text:
- **Design** (Methods, "Soil microcosms"): a single (temperate) soil amended with sucrose, starch or cellulose
  at 0, 0.1, 1, 3 and 5% C-eq per g dry soil. H2 uptake was measured, plus ddPCR and amplicon sequencing of
  16S rRNA, ITS and **hhyL**.
- **The marker is group 1h only** (Methods, "Abundance of taxonomic and functional genes": "hhyL gene encoding
  the large subunit of group 5/1h NiFe-hydrogenase"). Groups 1l and 2a were not measured.
- **Results** (Abstract; Results "Carbon preference traits"; Conclusion):
  - Only 0.1% cellulose stimulated activity.
  - Sucrose abated activity by 50% at 1% C-eq. hhyL copy number and richness fell, and at 3–5% sucrose hhyL
    could not be amplified for sequencing.
  - The response of individual hhyL ASVs was uneven: only 6 hhyL ASVs (1%) responded to amendments, versus
    12% of bacterial and 8% of fungal ASVs. The authors conclude that carbon preference is "a difficult
    trait to predict".
- **Mechanism is inferred, not shown** (Discussion, "Direct and indirect effects..."). Catabolite repression
  is "presumed" and argued from isolates (Ralstonia, M. smegmatis, Streptomyces). The authors offer three
  alternatives: fermentation, catabolic repression and microbial interactions.

What it supports: labile carbon at high doses depletes a group 1h scavenger guild and its activity in one
temperate soil. What it does **not** support: a group-1l-specific mechanism. If anything, it shows group 1h
can be carbon-sensitive, while in our root contrast group 1h did not change beyond the controls. It also
shows that carbon responses differ between genotypes within a group.

### Kropp A, Archer JD, Jespersen M, Watts TD, Solari J, Huang C, Schittenhelm RB, Grinter R, Greening C (2026) mSystems 11(4)

**Title.** Atmospheric hydrogen consumption is regulated by glycerol-mediated catabolite repression in
mycobacteria. DOI 10.1128/msystems.01678-25. Local copy: `lifeonair/literature/kropp_2026_atmospheric-hydrogen-consumption-regulated.fulltext.xml`
(PMC13098282). Key `kropp2026glycerol`.

Verified from the full text:
- **Enzymes** (Introduction): M. smegmatis has "Huc (group 2a) and Hhy (group 1h)". The regulated enzyme in
  this paper is **Huc = group 2a**. Hhy is "active during long-term persistence", while Huc "enables M.
  smegmatis to grow mixotrophically and is most active during the transition from growth to dormancy".
- **Glycerol specificity** (Results, "GylR is not required for growth and Huc activity with alternative
  organic substrates"). WT and ΔgylR were grown on glycerol, glucose, acetate or succinate. "During growth
  with glucose, acetate, or succinate as the sole carbon source, no difference in Huc activity was observed
  between the WT and gylR mutant strains, in contrast to the substantially increased Huc activity observed by
  the gylR mutant during growth with glycerol" (Fig. 2B).
- **Indirect repression** (Abstract; Conclusion): GylR, Crp1 and Crp2 do not bind the huc operon. Repression is
  mediated by unidentified transcription factors downstream of GylR.
- **Explicit limits** (Conclusion): the authors ask whether other catabolite sensors regulate Huc on other
  carbon sources, and "how are Hhy and CO dehydrogenase induced during starvation, given they are
  unresponsive to GylR?". Hhy (1h) and CODH are therefore **not** GylR-regulated.

What it supports: in one actinobacterium, a group 2a hydrogenase is repressed specifically by glycerol
availability. It says nothing about group 1l, and it shows that the 1h and CODH responses run through
different (unknown) regulation.

### Implication for the manuscript (R2-24)

main.tex currently cites neither paper. The round-1 digest (`round1/DIGEST.md` lines 29, 43) summarised them
as a "direct mechanism for root-adjacent group 1l loss". That is wrong: Baril measured group 1h and Kropp
studied group 2a. If carbon regulation is discussed (Discussion ~414–420, after the Li 2023 sentence
"Li et al. measured organic carbon, carbon additions and uptake and attributed the depletion to oligotrophic
oxidisers being outcompeted where carbon is available"), the evidence should be presented as substrate- and
group-specific. Suggested text:
> Carbon responses of high-affinity hydrogenases are substrate- and group-specific: sucrose suppressed
> hydrogen uptake and group 1h gene copies in soil microcosms, whereas cellulose at low dose stimulated uptake
> \citep{baril2023carbon}, and in Mycobacterium smegmatis glycerol, but not glucose, acetate or succinate,
> repressed the group 2a hydrogenase Huc, while the group 1h enzyme was not under this control
> \citep{kropp2026glycerol}. No comparable experiment exists for group 1l, so these studies do not explain why
> only group 1l was depleted beside roots.

This also suits the existing sentence "Something that distinguishes group 1l carriers, or the habitats they
occupy, is therefore needed to explain the pattern."

---

## Other papers by direction

### d1 — genome traits and aridity

| Key (new_refs) | Paper | Local copy | Verified finding | Manuscript sentence |
|---|---|---|---|---|
| `leung2025aridity` | Leung PM, Bay SK, ... Greening C (2025) Aridity drives global convergence of desert microbiomes and biogeochemical activities. bioRxiv 10.1101/2025.11.29.691269 (preprint) | `leung_2025_aridity-drives-global-convergence.pdf` | Completeness-normalised genome size of MAGs from nonpolar deserts averages 3.7 Mb, smaller than other biomes (Fig. 2A). RuBisCO abundance correlates most strongly with high-affinity [NiFe]-hydrogenases (R² 0.80, Fig. 3G); two-thirds of RuBisCO-encoding MAGs co-encode uptake hydrogenases. Uses the Greening marker DB with 50% identity (60% for CoxL) thresholds. | Context for our 1l-only median of 3.78 Mb (trait section ~351) and **directly relevant to the CBB co-carriage sentence** ("Calvin–Benson–Bassham genes were encoded by 18% of group 1l-only and 8% of group 1h-only carriers ...", ~354). Also supports the identity thresholds in Methods (~165). Cite as a preprint. |
| `liu2023genomesize` | Liu S, et al. (2023) Warmer and drier ecosystems select for smaller bacterial genomes in global soils. iMeta 2: e70. 10.1002/imt2.70 | none (Wiley 403) | *Abstract only*: across 237 global sites (16S profiling), aridity is a major regulator of community genome size; warmer and drier ecosystems select for smaller genomes. | Same as Wang 2025: a community-level, cross-lineage pattern. Optional citation next to Wang. |
| `chuckran2023edaphic` | Chuckran PF, et al. (2023) Edaphic controls on genome size and GC content of bacteria in soil microbial communities. Soil Biol Biochem 178: 108935. 10.1016/j.soilbio.2022.108935 | `chuckran_2023_edaphic-controls-genome-size.pdf` (preprint version) | Smaller genomes with high GC are associated with high pH and low soil C:N, interpreted as selection against costly AT-rich (more C/N per base) nucleotides. | Possible support for interpreting the high GC of 1l carriers (~352, ~422); optional. |
| `dragone2022antarctic` | Dragone NB, et al. (2022) Elevational constraints on the composition and genomic attributes of microbial communities in Antarctic soils. mSystems 7: e01330-21. 10.1128/msystems.01330-21 | `dragone_2022_..._genomic.fulltext.xml` | At higher, harsher elevations, gRodon doubling times are longer (P = 0.042) and trace-gas genes are over-represented. | Context for the predicted-growth sentence (~337, "carriers ... had longer predicted doubling times"); optional. |
| — | Simonsen AK (2022) ISME J 16: 423 (10.1038/s41396-021-01082-x) | `simonsen_2021_...pdf` | Genome streamlining along aridity, heat, salinity and acidity, **within one species** (Bradyrhizobium diazoefficiens). | Not recommended; different scale. |
| — | Chuckran 2025 PNAS (10.1073/pnas.2413032122) | PDF | 18O-SIP shows ribosomal-protein codon bias is the strongest predictor of in situ growth after rewetting. | Optional support that gRodon predictions track realised growth after wetting. The manuscript correctly says "predicted growth potential". |

### d2 — carbon regulation (besides Baril and Kropp)

| Key | Paper | Local copy | Verified finding | Manuscript sentence |
|---|---|---|---|---|
| `berney2014three` | Berney M, Greening C, Hards K, Collins D, Cook GM (2014) Three different [NiFe] hydrogenases confer metabolic flexibility in the obligate aerobe Mycobacterium smegmatis. Environ Microbiol 16: 318. 10.1111/1462-2920.12320 | none (403) | *Abstract only*: Hyd1 (= Huc, group 2a) and Hyd2 (= Hhy, group 1h) **both respond to carbon and oxygen limitation**. | With Kropp: carbon limitation induces both 1h and 2a in one organism, so carbon availability alone does not predict a group-1l-only response. |
| `cordero2019huchhy` | Cordero PRF, et al. (2019) Two uptake hydrogenases differentially interact with the aerobic respiratory chain during mycobacterial growth and persistence. J Biol Chem 294: 18980. 10.1074/jbc.RA119.011076 | `cordero_2019_uptake-hydrogenases-differentially-interact.pdf` | Huc (2a) is active in late exponential/early stationary phase and supports mixotrophic growth; Hhy (1h) is most active during long-term persistence when carbon is depleted. | Physiology sentence in the Introduction (~55–60) and Discussion: groups differ in growth stage. Optional. |
| — | Greening 2015 PNAS (already cited) | — | — | — |
| `peng2026hgt` | Peng (2026) Widespread horizontal transfer and strong selection enhance microbial adaptation in Antarctic soils. Nat Commun. 10.1038/s41467-026-76323-8 | `peng_2026_widespread-horizontal-transfer-strong.pdf` | Group 1l and 1h uptake hydrogenases were the most horizontally transferred genes in these polar desert communities (24 HGT events, 42 genes, 8 phyla); about a quarter of CODH and [NiFe]-hydrogenases are predicted to be horizontally acquired. The authors call group 1l "uncharacterized" and "hypothesized" to be O2-tolerant and mixotrophic. | **Bears on** "the depletion was carried by whole genera, so the enzyme may mark a lineage ecology rather than be the target of selection" (~423–424) and on the genus-carriage model: frequent HGT of 1l means carriage is not fixed by lineage, which weakens a pure lineage-marker reading. Also supports stating that group 1l physiology is uncharacterised. |
| — | Tveit 2021 (XML), Meier 2021 (XML), Firrincieli 2026 (XML) | local | Tveit: Methylocapsa gorgona co-oxidises CH4, CO and H2 to grow on air. Meier: biocrust MAGs show widespread mixotrophy and desiccation survival. | Not needed. |

### d3 — validation of read-based quantification (R2-03/04/05)

| Key | Paper | Local copy | Verified finding | Manuscript sentence |
|---|---|---|---|---|
| `nayfach2015microbecensus` | Nayfach S, Pollard KS (2015) Average genome size estimation improves comparative metagenomics ... Genome Biol 16: 51. 10.1186/s13059-015-0611-7 | `nayfach_2015_average-genome-size-estimation.pdf` | Single-copy-gene genome equivalents (MicrobeCensus), validated on simulated communities. | Methods "The result, genome equivalents (GE), is the number of marker copies per genome" (~170). Cite as the precedent for SCG-based GE and its validation; this supports R2-05's request to qualify GE as an estimated marker-to-SCG ratio. |
| `orellana2017rocker` | Orellana LH, Rodriguez-R LM, Konstantinidis KT (2017) ROCker ... Nucleic Acids Res 45: e14. 10.1093/nar/gkw900 | `orellana_2016_rocker-accurate-detection-quantification.pdf` | Fixed e-value/bitscore cutoffs produce false positives among homologous enzymes; position-specific bitscore models improve precision. | Methods (DIAMOND thresholds, ~163–168) and the CoxL precision issue (R2-04: form-I precision 0.77). Supports a stricter, position-aware CoxL read classifier. |
| `browne2020gcbias` | Browne PD, et al. (2020) GC bias affects genomic and metagenomic reconstructions, underrepresenting GC-poor organisms. GigaScience 9: giaa008. 10.1093/gigascience/giaa008 | `browne_2020_bias-affects-genomic-metagenomic.pdf` | MiSeq and NextSeq workflows show strong GC bias outside 45–65% GC, giving falsely low coverage of GC-rich and especially GC-poor sequences; HiSeq and PacBio have a different profile. | New Limitation: group 1l carriers are GC-rich (median 71%), so platform- or library-dependent GC bias could shift 1l GE relative to the community denominator. It would cancel only if libraries were processed identically across positions. Worth one sentence, or a check that platform and library are balanced across positions. |
| — | Tamames J, et al. (2019) BMC Genomics 20: 960 (10.1186/s12864-019-6289-6) | PDF | Read-based functional annotation recovers more functions but "makes a substantial number of over-predictions"; database completeness drives read-assignment performance. | Supports R2-04; optional. |
| — | Carr & Borenstein 2014 PLoS One (10.1371/journal.pone.0105776) | PDF | Many short-read translated matches are tied or mis-annotated across homologous families. | Optional. |
| — | Manor & Borenstein 2015 MUSiCC (10.1186/s13059-015-0610-8) | PDF | Normalising to total reads biases gene abundance; universal single-copy genes with correction give accurate gene abundances. | Optional, alongside Nayfach. |
| — | Seppey 2020 LEMMI, Meyer 2022 CAMI II, Beghini 2021 bioBakery 3, Zhou 2024 GCompip | local | Classifier and abundance benchmarks on mock and simulated communities. | Reference designs for the mock-community GE benchmark requested in R2-05; not needed in the text. |

### d4 — CoxL form I/II (R2-07, R1-01)

| Key | Paper | Local copy | Verified finding | Manuscript sentence |
|---|---|---|---|---|
| `imaura2026coxl` | Imaura (2026) Diversity and geographical distribution of potential carbon monoxide oxidizers using molybdenum-containing enzymes in the ocean. mSphere. 10.1128/msphere.00062-26 | `imaura_2026_..._potential.fulltext.xml` | Calls form I by **phylogeny + active-site motif + coxMSL gene order**. Notes that earlier surveys "may have overestimated" CO oxidizers because KEGG/COG homology also catches form II. Six form-I motif variants among 707 references: AYACSFR, AYRCSFR, SYRCSFR, AYRCSLR, AYSCSFR, SYACSFR. | Methods form-I rule (~136–139, "form I when it carried the form-I active-site motif AYXCSFR"). Two of these variants (SYxCSFR, AYRCSLR) do **not** match AYXCSFR, so a strict AYXCSFR regex would misclassify them; confirmed: `fix_coxL_v2.py:27` uses `AY.CSFR`, so SYxCSFR and AYRCSLR proteins would be called form II or unresolved. Their use of coxMSL order is the direct answer to R1-29/R2-07 (motif-less calls): order and strand plus phylogeny can resolve them. |
| `cunliffe2011roseobacter` | Cunliffe M (2011) Correlating carbon monoxide oxidation with cox genes in the abundant Marine Roseobacter Clade. ISME J 5: 685. 10.1038/ismej.2010.170 | none (403) | *Abstract only*: 18/29 Roseobacter genomes have only form II, 10 have both, none form I only; form-I maturation genes are absent from form-II-only genomes; only strains with both forms oxidised CO. | Support for treating form II as non-CO-oxidising in the Methods typing paragraph (~136). Cite with King 2007. |
| — | Hogendoorn C, et al. (2020) Front Microbiol 11: 951 (10.3389/fmicb.2020.00951) | PDF | Kyrpidia spormannii strains encode only a form-II CODH (AYRGAGR) and are described as H2- and CO-utilising. The paper itself states "no firm evidence is reported for CO oxidation by Form II enzymes" and cites Cunliffe 2011. | An optional caveat that form-II function is unresolved. **Flash's claim that the paper "demonstrates" CO growth via form II is overstated.** |
| — | Nishimura 2010 Fish Sci (10.1007/s12562-010-0277-8) | `nishimura_2010_...pdf` (via ws) | Purified an active Mo-CODH from Aeropyrum pernix, 2.1 µmol CO/min/mg at 95 °C. **Flash's claim that it "lacks the canonical form I cysteine motif" is not stated in the text read.** | Not recommended. |
| — | Lorite 2000 AEM (10.1128/aem.66.5.1871-1876.2000) | none | *Abstract only*: Bradyrhizobium japonicum grows chemolithoautotrophically on CO with a purified Mo-CODH. The abstract says nothing about form II. | Not recommended without checking the coxL sequence. |
| — | King 2003, Weber 2010, Pelzmann 2009, Lalonde 2016, Myers 2020 | none | Not verified. | — |

---

## Flash errors found (do not propagate)

1. d1, Kropp 2026: "via group 1l [NiFe]-hydrogenase (Huc)". **Wrong**: Huc is group 2a (Kropp Introduction).
2. d1/d2 synthesis: catabolite repression "explains group 1l depletion". Not shown for group 1l by any paper.
3. d4: Hogendoorn 2020 and Lorite 2000 are presented as proof that form II oxidises CO. Hogendoorn
   explicitly says there is no firm evidence; the Lorite abstract does not mention form II.
4. d4: Nishimura 2010 "lacks canonical motif": not found in the text.
5. d1 synthesis: "group 1l scavengers exhibit streamlined genome sizes". Our within-order 1l effect is
   −1.9% and not significant; do not describe 1l as streamlined.
6. Round-1 digest entries for Baril 2023 and Kropp 2026 (`round1/DIGEST.md` lines 29, 43) misstate them as
   group-1l mechanisms; corrected above.

---

## Recommended manuscript changes (for the parent; main.tex not edited)

1. **Trait section and Discussion**: add the Wang 2025 comparison (abstract-level, hedged). Our within-order
   result (slower growth, not smaller genomes) differs from the cross-MAG pattern (`wang2025genomesize`,
   optionally `liu2023genomesize`, `leung2025aridity`). Do not call group 1l carriers streamlined.
2. **Carbon paragraph** (after the Li 2023 sentence, ~415): add the substrate- and group-specific evidence from
   `baril2023carbon` (1h, sucrose) and `kropp2026glycerol` (2a Huc, glycerol only; Hhy and CODH not under
   GylR), optionally `berney2014three` (1h and 2a both induced by carbon limitation). State explicitly that
   no experiment addresses group 1l.
3. **Lineage-marker sentence** (~423): qualify it with `peng2026hgt` (1l among the most horizontally
   transferred genes in polar desert soils), so carriage is not purely lineage-fixed.
4. **CBB sentence** (~354): cite `leung2025aridity`; RuBisCO co-occurs with high-affinity hydrogenases in
   desert MAGs. Our 18% vs 8% is consistent in direction, but it is a three-KO candidate screen (R2-21).
5. **Methods, GE** (~170): cite `nayfach2015microbecensus` for SCG genome equivalents. Call GE an estimated
   marker-to-SCG ratio (R2-05).
6. **Methods/Limitations, classifier** (~163–168; R2-04): cite `orellana2017rocker` for homologue
   false-positive risk under fixed thresholds.
7. **Limitations**: add one sentence on GC-dependent coverage bias (`browne2020gcbias`), because 1l carriers are
   GC-rich (median 71%).
8. **CoxL typing** (~136–139): cite `cunliffe2011roseobacter` (only form-I-containing strains oxidised CO) and
   `imaura2026coxl` (phylogeny + motif + coxMSL order). Check that the AYXCSFR rule admits the SYxCSFR and
   AYRCSLR variants reported by Imaura; this matters for the motif-less and unresolved calls in R2-07.

BibTeX to add from `new_refs.bib`: baril2023carbon, kropp2026glycerol, berney2014three, cordero2019huchhy,
leung2025aridity, liu2023genomesize, chuckran2023edaphic, dragone2022antarctic, peng2026hgt, imaura2026coxl,
cunliffe2011roseobacter, browne2020gcbias, nayfach2015microbecensus, orellana2017rocker (wang2025genomesize is
already in references.bib).
