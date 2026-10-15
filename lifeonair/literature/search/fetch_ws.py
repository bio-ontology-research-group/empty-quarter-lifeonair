#!/usr/bin/env python3
"""Retry paywalled downloads from the KAUST network (run on ws: `ssh ws-ts`).
Input: TSV key<TAB>doi. Tries publisher PDF URL patterns; writes <outdir>/<key>.pdf and prints key<TAB>outcome."""
import os, sys, urllib.request, http.cookiejar

UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36"
opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))


def candidates(doi):
    d = doi.lower()
    suf = doi.split("/", 1)[1]
    u = []
    if d.startswith("10.1038/"):
        u.append(f"https://www.nature.com/articles/{suf}.pdf")
    if d.startswith(("10.1007/", "10.1186/")):
        u.append(f"https://link.springer.com/content/pdf/{doi}.pdf")
    if d.startswith(("10.1002/", "10.1111/", "10.1046/", "10.1029/")):
        u += [f"https://onlinelibrary.wiley.com/doi/pdfdirect/{doi}", f"https://onlinelibrary.wiley.com/doi/pdf/{doi}"]
    if d.startswith("10.1073/"):
        u.append(f"https://www.pnas.org/doi/pdf/{doi}")
    if d.startswith("10.1128/"):
        u.append(f"https://journals.asm.org/doi/pdf/{doi}")
    if d.startswith("10.1126/"):
        u.append(f"https://www.science.org/doi/pdf/{doi}")
    if d.startswith("10.1146/"):
        u.append(f"https://www.annualreviews.org/doi/pdf/{doi}")
    if d.startswith("10.1021/"):
        u.append(f"https://pubs.acs.org/doi/pdf/{doi}")
    if d.startswith("10.1016/"):
        pii = None
        try:
            r = opener.open(urllib.request.Request(f"https://doi.org/{doi}", headers={"User-Agent": UA}), timeout=40)
            import re
            m = re.search(r"pii/(S?[0-9X]{16,17})", r.geturl() + r.read(200000).decode("latin1"))
            pii = m.group(1) if m else None
        except Exception:
            pass
        if pii:
            u.append(f"https://www.sciencedirect.com/science/article/pii/{pii}/pdfft?isDTMRedir=true&download=true")
    u.append(f"https://doi.org/{doi}")
    return u


def main():
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for line in open(sys.argv[1]):
        key, doi = line.rstrip("\n").split("\t")[:2]
        if not doi.startswith("10."):
            print(f"{key}\tskip-non-doi", flush=True)
            continue
        res = "FAIL(ws)"
        for url in candidates(doi):
            try:
                r = opener.open(urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/pdf,*/*"}), timeout=60)
                blob = r.read()
                if blob[:4] == b"%PDF":
                    open(os.path.join(out, key + ".pdf"), "wb").write(blob)
                    res = f"OK {len(blob)//1024}KB (KAUST network via ws)"
                    break
                res = f"FAIL(ws:not-pdf)"
            except Exception as e:
                res = f"FAIL(ws:{getattr(e, 'code', type(e).__name__)})"
        print(f"{key}\t{res}", flush=True)


if __name__ == "__main__":
    main()
