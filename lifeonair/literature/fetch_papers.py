#!/usr/bin/env python3
"""Download OA PDFs for the literature cited in the 2026-10-09 P3 life-on-air novelty
check. Resolves each DOI via Unpaywall; falls back to
EuropePMC render-PDF for PMC-only items. Writes MANIFEST.tsv (outcome per paper).
Paywalled items with no OA copy are recorded as such (citation kept regardless)."""
import json, os, sys, time, urllib.request, urllib.parse

EMAIL = "robert.hoehndorf@kaust.edu.sa"
HERE = os.path.dirname(os.path.abspath(__file__))

# key, doi, pmcid, citation
PAPERS = [
    ("bay_2021_trace-gas-oxidizers-widespread", "10.1038/s41564-020-00811-w", None,
     "Bay et al. Trace gas oxidizers are widespread and active members of soil microbial communities. Nat Microbiol 2021."),
    ("ji_2017_atmospheric-chemosynthesis-antarctica", "10.1038/nature25014", None,
     "Ji et al. Atmospheric trace gases support primary production in Antarctic desert surface soil. Nature 2017."),
    ("ortiz_2021_antarctic-multiple-energy-group1l", "10.1073/pnas.2025322118", None,
     "Ortiz et al. Multiple energy sources and metabolic strategies sustain microbial diversity in Antarctic desert soils. PNAS 2021."),
    ("jordaan_2020_hot-desert-H2-hydration", "10.1128/mSystems.01131-20", "PMC7677003",
     "Jordaan et al. Hydrogen-oxidizing bacteria are abundant in desert soils and strongly stimulated by hydration. mSystems 2020."),
    ("leung_2020_desert-energetics", "10.1128/mSystems.00495-19", None,
     "Leung et al. Energetic basis of microbial growth and persistence in desert ecosystems. mSystems 2020."),
    ("bay_2021_aridity-gradient-chemosynthesis", "10.1038/s41396-021-01001-0", None,
     "Bay et al. Chemosynthetic and photosynthetic bacteria contribute differentially to primary production across a steep desert aridity gradient. ISME J 2021."),
    ("cordero_2019_atmospheric-CO-widespread", "10.1038/s41396-019-0479-8", None,
     "Cordero et al. Atmospheric carbon monoxide oxidation is a widespread mechanism supporting microbial survival. ISME J 2019."),
    ("ray_2022_chemosynthesis-cold-deserts", "10.1038/s41396-022-01298-5", "PMC9561532",
     "Ray et al. Atmospheric chemosynthesis is phylogenetically and geographically widespread and contributes significantly to carbon fixation throughout cold deserts. ISME J 2022."),
    ("polloquere_2024_andean-hot-spring-trace-gas", "10.1093/ismejo/wrae062", "PMC11094475",
     "Prevalence of trace gas-oxidizing soil bacteria increases with radial distance from Polloquere hot spring within a high-elevation Andean cold desert. ISME J 2024."),
    ("fertile-islands_2023_trace-gas-oligotrophs", None, "PMC10356767",
     "Reduced trace gas oxidizers as a response to organic carbon availability linked to oligotrophs in desert fertile islands. ISME J 2023;17:1257."),
    ("reji_2026_moisture-threshold-H2-uptake", "10.1038/s41467-025-67208-3", "PMC12804893",
     "Reji, Bertagni, Paulot et al. Global implications of a low soil moisture threshold for microbial hydrogen uptake. Nat Commun 2026."),
    ("greening_2023_clearing-the-air-review", "10.1128/mmbr.00048-23", None,
     "Clearing the air: unraveling past and guiding future research in atmospheric chemosynthesis. Microbiol Mol Biol Rev 2023."),
    ("warming_2025_trace-gas-oxidizers-preprint", "10.1101/2025.11.18.689010", None,
     "Soil trace gas oxidizers divergently respond to short- and long-term warming. bioRxiv 2025 (preprint)."),
    ("southern-africa_2024_dryland-functional-redundancy", "10.1093/femsec/fiae157", None,
     "Functional redundancy buffers the effect of poly-extreme environmental conditions on southern African dryland soil microbial communities. FEMS Microbiol Ecol 2024."),
    ("antarctic_2026_trace-gas-temperature", None, "PMC12998438",
     "Resilient Antarctic soil bacteria consume trace gases across wide temperature ranges. 2026."),
    ("greening_2016_hydrogenase-distribution", "10.1038/ismej.2015.153", None,
     "Greening et al. Genomic and metagenomic surveys of hydrogenase distribution indicate H2 is a widely utilised energy source for microbial growth and survival. ISME J 2016."),
    ("sondergaard_2016_hyddb", "10.1038/srep34212", None,
     "Sondergaard, Pedersen & Greening. HydDB: a web tool for hydrogenase classification and analysis. Sci Rep 2016."),
]


def get(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (paper-archive)"})
    return urllib.request.urlopen(req, timeout=timeout)


def unpaywall_pdf(doi):
    url = f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email={EMAIL}"
    try:
        data = json.load(get(url))
    except Exception as e:
        return None, f"unpaywall-error:{e}"
    loc = data.get("best_oa_location") or {}
    pdf = loc.get("url_for_pdf")
    return (pdf, "oa") if pdf else (None, "closed")


def download(pdf_url, dest):
    try:
        with get(pdf_url) as r:
            blob = r.read()
        if blob[:4] != b"%PDF":
            return False, f"not-pdf({len(blob)}b)"
        open(dest, "wb").write(blob)
        return True, f"{len(blob)//1024}KB"
    except Exception as e:
        return False, f"dl-error:{e}"


def main():
    rows = []
    for key, doi, pmcid, cite in PAPERS:
        dest = os.path.join(HERE, key + ".pdf")
        if os.path.exists(dest) and os.path.getsize(dest) > 10000:
            rows.append((key, "cached", doi or pmcid, cite)); print("cached", key); continue
        pdf_url, status = (None, "no-id")
        if doi:
            pdf_url, status = unpaywall_pdf(doi)
        if not pdf_url and pmcid:
            pdf_url, status = (f"https://europepmc.org/articles/{pmcid}?pdf=render", "europepmc")
        outcome = status
        if pdf_url:
            ok, info = download(pdf_url, dest)
            outcome = f"OK {info}" if ok else f"FAIL({status}):{info}"
        print(f"{outcome:24s} {key}")
        rows.append((key, outcome, doi or pmcid, cite))
        time.sleep(0.4)
    with open(os.path.join(HERE, "MANIFEST.tsv"), "w") as o:
        o.write("key\toutcome\tdoi_or_pmcid\tcitation\n")
        for r in rows:
            o.write("\t".join(str(x) for x in r) + "\n")
    okn = sum(1 for r in rows if str(r[1]).startswith(("OK", "cached")))
    print(f"\n{okn}/{len(rows)} PDFs archived -> {HERE}")


if __name__ == "__main__":
    main()
