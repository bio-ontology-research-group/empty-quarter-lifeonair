#!/usr/bin/env python3
"""Life on air (P3): NASA POWER daily precipitation 2021-01-01 to 2025-10-31 (review round 2, R2-15).

An exposure product independent of ERA5 for the within-site temporal models: NASA POWER PRECTOTCORR
(MERRA-2 reanalysis precipitation bias-corrected with GPCP gauge-satellite data, 0.5 x 0.625 deg), fetched at
the site coordinates. Antecedent 90- and 365-day sums are then computed for the same visits and windows as
ERA5 (prep_environment.R). Only site coordinates are sent.

    python3 lifeonair/code/fetch_power_daily.py
Writes lifeonair/results/inputs/power_daily_2021_2025.tsv (site, date, P). Resumable.
"""
import csv
import json
import os
import time
import urllib.request

IN = "lifeonair/results/inputs"
URL = ("https://power.larc.nasa.gov/api/temporal/daily/point?parameters=PRECTOTCORR&community=AG"
       "&longitude={lon:.5f}&latitude={lat:.5f}&start=20210101&end=20251031&format=JSON")
OUT = f"{IN}/power_daily_2021_2025.tsv"
sites = list(csv.DictReader(open(f"{IN}/env_site.tsv"), delimiter="\t"))
done = {l.split("\t")[0] for l in open(OUT)} if os.path.exists(OUT) else set()
with open(OUT, "a") as out:
    if not done:
        out.write("site\tdate\tP\n")
    for s in sites:
        if s["site"] in done:
            continue
        url = URL.format(lat=float(s["lat"]), lon=float(s["lon"]))
        for attempt in range(5):
            try:
                j = json.load(urllib.request.urlopen(url, timeout=120))
                break
            except Exception as e:  # noqa: BLE001 - retry transient API errors
                print("retry", s["site"], e); time.sleep(10 * (attempt + 1))
        p = j["properties"]["parameter"]["PRECTOTCORR"]
        for d, v in sorted(p.items()):
            out.write(f"{s['site']}\t{d[:4]}-{d[4:6]}-{d[6:]}\t{'' if v < -900 else v}\n")
        out.flush(); print("site", s["site"], len(p)); time.sleep(1)
