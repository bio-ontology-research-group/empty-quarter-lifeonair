#!/usr/bin/env python3
"""Life on air (P3): daily ERA5 precipitation 1995-01-01 to 2025-10-31 (review round 1, R1-03).

The frozen climate record (amplicon cache climate_historical_1995_2024.parquet, Open-Meteo archive, ERA5
0.25 deg) ends on 2024-12-31, so antecedent precipitation was missing for the October 2025 expedition. It
was also fetched with Open-Meteo's default "best_match" model, which blends ERA5 with ECMWF IFS from 2017, and
re-querying that default on 12 Oct 2026 returns different 2024 totals (site 1: 22.9 mm vs 102.7 mm frozen).
This script fetches pure ERA5 (models=era5) for the whole period at the site coordinates, so long-term MAP
and antecedent sums come from one consistent reanalysis. Only site coordinates are sent.

    python3 lifeonair/code/fetch_era5_2025.py
Writes lifeonair/results/inputs/era5_daily_1995_2025.tsv (site, date, P) and era5_daily_1995_2025.urls.tsv
(request URL per site). Resumable: sites already in the output are skipped.
"""
import csv
import json
import time
import urllib.request

IN = "lifeonair/results/inputs"
URL = ("https://archive-api.open-meteo.com/v1/archive?latitude={lat:.5f}&longitude={lon:.5f}"
       "&start_date=1995-01-01&end_date=2025-10-31&daily=precipitation_sum&models=era5&timezone=UTC")

import os
sites = list(csv.DictReader(open(f"{IN}/env_site.tsv"), delimiter="\t"))
OUT = f"{IN}/era5_daily_1995_2025.tsv"
done = {l.split("\t")[0] for l in open(OUT)} if os.path.exists(OUT) else set()
new = not done
with open(OUT, "a") as out, open(f"{IN}/era5_daily_1995_2025.urls.tsv", "a") as uo:
    if new:
        out.write("site\tdate\tP\n")
        uo.write("site\turl\n")
    for s in sites:
        if s["site"] in done:
            continue
        url = URL.format(lat=float(s["lat"]), lon=float(s["lon"]))
        for attempt in range(5):
            try:
                d = json.load(urllib.request.urlopen(url, timeout=60))
                break
            except Exception:
                time.sleep(30 * (attempt + 1))
        else:
            raise SystemExit(f"failed: site {s['site']}")
        for day, p in zip(d["daily"]["time"], d["daily"]["precipitation_sum"]):
            out.write(f"{s['site']}\t{day}\t{'' if p is None else p}\n")
        uo.write(f"{s['site']}\t{url}\n")
        out.flush(); uo.flush()
        time.sleep(5.0)
print("wrote", f"{IN}/era5_daily_1995_2025.tsv")
