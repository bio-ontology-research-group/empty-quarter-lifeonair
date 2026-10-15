#!/usr/bin/env python3
"""Merge Flash search-direction outputs, verify every DOI against Crossref, deduplicate against the
existing P3 literature, and download open-access copies.

Usage: python3 process_hits.py <round_dir>
Writes <round_dir>/candidates.tsv (all parsed rows with verification status) and appends downloaded
papers to ../MANIFEST.tsv (key, outcome, doi, citation). PDFs go to ../ (lifeonair/literature/).
Failed downloads are kept as citation + DOI and retried later from the KAUST network (fetch_ws.py).
"""
import csv, glob, json, os, re, sys, time, unicodedata, urllib.parse, urllib.request
from difflib import SequenceMatcher

EMAIL = "robert.hoehndorf@kaust.edu.sa"
HERE = os.path.dirname(os.path.abspath(__file__))
LIT = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(LIT, "..", ".."))
COLS = ["doi", "first_author", "year", "title", "journal", "manuscript_link", "stance", "relevance", "oa_url"]


def get(url, timeout=40, accept=None):
    h = {"User-Agent": f"Mozilla/5.0 (literature-archive; mailto:{EMAIL})"}
    if accept:
        h["Accept"] = accept
    return urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=timeout)


def norm_doi(d):
    d = d.strip().strip(".").lower()
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:\s*)", "", d)
    return d


def norm_title(t):
    t = unicodedata.normalize("NFKD", t).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9 ]+", " ", t).split()


def known_dois():
    s = set()
    for f in [os.path.join(LIT, "MANIFEST.tsv"), os.path.join(LIT, "..", "manuscript", "references.bib"),
              os.path.join(ROOT, "resources", "literature", "trace-gas-deserts", "downloads_manifest.tsv")]:
        if os.path.exists(f):
            s |= {norm_doi(m) for m in re.findall(r"10\.\d{4,9}/[^\s\t\"}<>,]+", open(f).read())}
    return s


def parse_outputs(rdir):
    rows = []
    for f in sorted(glob.glob(os.path.join(rdir, "out", "*.txt"))):
        txt = open(f, errors="ignore").read()
        m = re.search(r"BEGIN_TSV\s*\n(.*?)\n\s*END_TSV", txt, re.S)
        if not m:
            print(f"no TSV block in {f}", file=sys.stderr)
            continue
        for line in m.group(1).splitlines():
            parts = [p.strip() for p in line.strip().strip("|").split("\t")]
            if len(parts) < 4 or parts[0].lower() == "doi" or not parts[0].startswith(("10.", "http", "doi")):
                continue
            r = dict(zip(COLS, parts + ["NA"] * (len(COLS) - len(parts))))
            r["doi"] = norm_doi(r["doi"])
            r["direction"] = os.path.basename(f)[:-4]
            rows.append(r)
    return rows


def crossref(doi):
    try:
        msg = json.load(get(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"))["message"]
    except Exception as e:
        return None, f"crossref:{getattr(e, 'code', e)}"
    title = (msg.get("title") or [""])[0]
    year = (msg.get("issued", {}).get("date-parts") or [[None]])[0][0]
    au = msg.get("author") or [{}]
    return {"title": title, "year": year, "journal": (msg.get("container-title") or [""])[0],
            "first_author": au[0].get("family", ""), "authors": [a.get("family", "") for a in au]}, "ok"


def slug(cr, doi):
    a = re.sub(r"[^a-z]", "", unicodedata.normalize("NFKD", cr["first_author"]).encode("ascii", "ignore").decode().lower()) or "anon"
    words = [w for w in norm_title(cr["title"]) if len(w) > 3 and w not in {"with", "from", "that", "their", "into", "soil", "soils"}][:4]
    return f"{a}_{cr['year']}_{'-'.join(words)}"


def is_pdf_to(url, dest):
    try:
        with get(url, timeout=60) as r:
            blob = r.read()
    except Exception as e:
        return False, f"dl:{getattr(e, 'code', e)}"
    if blob[:4] != b"%PDF":
        return False, f"not-pdf({len(blob)}b)"
    open(dest, "wb").write(blob)
    return True, f"OK {len(blob)//1024}KB"


def fetch(doi, oa_hint, dest_base):
    urls = []
    pmcid = None
    try:
        up = json.load(get(f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email={EMAIL}"))
        for loc in [up.get("best_oa_location")] + (up.get("oa_locations") or []):
            if loc and loc.get("url_for_pdf"):
                urls.append(loc["url_for_pdf"])
    except Exception:
        pass
    if oa_hint and oa_hint.startswith("http") and oa_hint.lower().endswith(".pdf"):
        urls.append(oa_hint)
    if doi.startswith("10.1101/"):
        urls.append(f"https://www.biorxiv.org/content/{doi}.full.pdf")
    try:
        q = json.load(get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22{urllib.parse.quote(doi)}%22&format=json&resultType=lite"))
        res = q.get("resultList", {}).get("result", [])
        if res and res[0].get("pmcid"):
            pmcid = res[0]["pmcid"]
            urls.append(f"https://europepmc.org/articles/{pmcid}?pdf=render")
    except Exception:
        pass
    last = "no-oa-url"
    for u in dict.fromkeys(urls):
        ok, msg = is_pdf_to(u, dest_base + ".pdf")
        if ok:
            return msg, "pdf"
        last = msg
    if pmcid:
        try:
            xml = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML").read()
            if b"<article" in xml[:5000]:
                open(dest_base + ".fulltext.xml", "wb").write(xml)
                return f"OK fulltext XML (Europe PMC {pmcid})", "xml"
        except Exception as e:
            last = f"xml:{getattr(e, 'code', e)}"
    return f"FAIL({last})", None


def main():
    rdir = os.path.abspath(sys.argv[1])
    rows = parse_outputs(rdir)
    known = known_dois()
    manifest = os.path.join(LIT, "MANIFEST.tsv")
    have = {norm_doi(l.split("\t")[2]) for l in open(manifest).read().splitlines()[1:] if l.count("\t") >= 2}
    seen, out = {}, []
    for r in rows:
        d = r["doi"]
        if d in seen:
            seen[d]["direction"] += "," + r["direction"]
            continue
        seen[d] = r
        cr, st = crossref(d)
        time.sleep(0.2)
        if cr is None:
            r["verify"] = st
        else:
            sim = SequenceMatcher(None, " ".join(norm_title(cr["title"])), " ".join(norm_title(r["title"]))).ratio()
            r["verify"] = "ok" if sim >= 0.6 else f"title-mismatch({sim:.2f})"
            r["cr_title"], r["cr_year"], r["cr_journal"], r["cr_first_author"] = cr["title"], cr["year"], cr["journal"], cr["first_author"]
            r["key"] = slug(cr, d)
        r["dup"] = "known" if (d in known or d in have) else "new"
        out.append(r)
    # download verified new papers
    with open(manifest, "a") as mf:
        for r in out:
            r["download"] = ""
            if r["verify"] != "ok" or r["dup"] != "new":
                continue
            msg, kind = fetch(r["doi"], r.get("oa_url", ""), os.path.join(LIT, r["key"]))
            r["download"] = msg
            cite = " ".join(f"{r['cr_first_author']} et al. {r['cr_title']}. {r['cr_journal']} {r['cr_year']}.".split())
            mf.write(f"{r['key']}\t{msg}\t{r['doi']}\t{cite}\n")
            print(r["key"], msg, flush=True)
            time.sleep(0.5)
    cols = ["direction", "doi", "verify", "dup", "download", "key", "cr_first_author", "cr_year", "cr_title",
            "cr_journal", "manuscript_link", "stance", "relevance", "title", "first_author", "year"]
    with open(os.path.join(rdir, "candidates.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in out:
            w.writerow({c: " ".join(str(r.get(c, "")).split()) for c in cols})
    n = len(out)
    print(f"{n} unique DOIs; verified {sum(r['verify']=='ok' for r in out)}; new {sum(r['verify']=='ok' and r['dup']=='new' for r in out)}; "
          f"downloaded {sum(r['download'].startswith('OK') for r in out)}", file=sys.stderr)


if __name__ == "__main__":
    main()
