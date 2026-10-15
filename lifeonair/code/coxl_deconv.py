#!/usr/bin/env python3
"""Life on air (P3): form-I CoxL genome equivalents corrected for cross-assignment (review round 2, R2-04).

The held-out benchmark gives, per identity floor, the fraction of fragments from true form-I and true other
CoxL that the classifier assigns to each class (a = FI->FI, b = other->FI, c = FI->other, d = other->other).
Observed values then satisfy O_FI = a T_FI + b T_oth and O_oth = c T_FI + d T_oth; solving the 2 x 2 system
per sample gives error-corrected T_FI and T_oth (negative T_FI set to 0). Unassigned fragments rescale both
classes, so the corrected values estimate what a complete classifier would report.

    python3 lifeonair/code/coxl_deconv.py <confusion.tsv> <scheme> [OUT_NAME]
e.g. lifeonair/results/benchmark/classifier_benchmark_confusion.tsv clu50
Writes lifeonair/results/reads/<OUT_NAME> (default guild_genome_equivalents_v2_coxLcorr.tsv).
"""
import csv
import sys

conf, scheme = sys.argv[1], sys.argv[2]
out_name = sys.argv[3] if len(sys.argv) > 3 else "guild_genome_equivalents_v2_coxLcorr.tsv"
R = "lifeonair/results/reads"
rate = {}
for r in csv.DictReader(open(conf), delimiter="\t"):
    if r["scheme"] == scheme and r["true"] in ("coxL_formI", "coxL_other"):
        rate[(r["id_floor"], r["true"], r["assigned"])] = float(r["fraction"])
rows = list(csv.DictReader(open(f"{R}/guild_genome_equivalents_v2.tsv"), delimiter="\t"))
with open(f"{R}/{out_name}", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader()
    for r in rows:
        fl = r["id_floor"]
        a = rate.get((fl, "coxL_formI", "coxL_formI"), 0); c = rate.get((fl, "coxL_formI", "coxL_other"), 0)
        b = rate.get((fl, "coxL_other", "coxL_formI"), 0); d = rate.get((fl, "coxL_other", "coxL_other"), 0)
        of, oo = float(r["coxL_formI"]), float(r["coxL_other"])
        det = a * d - b * c
        tf = (d * of - b * oo) / det; to = (a * oo - c * of) / det
        r["coxL_formI"] = f"{max(tf, 0):.6f}"; r["coxL_other"] = f"{to:.6f}"
        w.writerow(r)
for fl in ("50", "60", "70"):
    print(fl, {k[1][5:] + ">" + k[2][5:]: v for k, v in rate.items() if k[0] == fl})
