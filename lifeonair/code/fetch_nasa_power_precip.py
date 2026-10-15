#!/usr/bin/env python3
"""Life on air (P3), H3: second, independent precipitation product for the 60 sites.

NASA POWER climatology (MERRA-2 based, corrected precipitation PRECTOTCORR, community AG), mean
annual precipitation = sum over months of mean mm/day x days (values are reported to 0.01 mm/day). Only site coordinates are sent. Output is frozen as
lifeonair/results/inputs/nasa_power_precip_climatology.tsv with the request URL per site.

    python3 lifeonair/code/fetch_nasa_power_precip.py
"""
import csv
import json
import time
import urllib.request

IN = "lifeonair/results/inputs"
URL = ("https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=PRECTOTCORR"
       "&community=AG&longitude={lon:.5f}&latitude={lat:.5f}&format=JSON")

DAYS = {"JAN": 31, "FEB": 28.25, "MAR": 31, "APR": 30, "MAY": 31, "JUN": 30, "JUL": 31,
        "AUG": 31, "SEP": 30, "OCT": 31, "NOV": 30, "DEC": 31}
sites = list(csv.DictReader(open(f"{IN}/env_site.tsv"), delimiter="\t"))
with open(f"{IN}/nasa_power_precip_climatology.tsv", "w") as out:
    out.write("site\tmap_power_mm\tclimatology_period\turl\n")
    for s in sites:
        url = URL.format(lon=float(s["lon"]), lat=float(s["lat"]))
        for attempt in range(4):
            try:
                d = json.load(urllib.request.urlopen(url, timeout=60))
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
        else:
            raise SystemExit(f"failed: site {s['site']}")
        p = d["properties"]["parameter"]["PRECTOTCORR"]
        ann = sum(p[m] * n for m, n in DAYS.items())   # monthly mm/day x days; ANN is rounded harder
        period = d["header"]["range"]
        out.write(f"{s['site']}\t{ann:.2f}\t{period}\t{url}\n")
        time.sleep(0.5)
print("wrote", f"{IN}/nasa_power_precip_climatology.tsv")
