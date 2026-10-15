#!/usr/bin/env python3
"""Life on air (P3): summarise the mock-community GE calibration (review round 2, R2-05).

    python3 lifeonair/code/b3_mock_summary.py   (from the repository root)

Input lifeonair/results/benchmark_r2/mock_ge_vs_truth.tsv (b3_mock_ge.sh). Outputs in the same folder:
mock_ge_accuracy.tsv   per database variant, identity floor and group: median and range of log2(estimated/true)
                       GE over communities, Spearman rho and Pearson r of log2 estimated vs true
mock_ge_root_bulk.tsv  per visit, variant, floor and group: true and estimated log2(root/bulk) GE, plus the
                       group 1l excess over the mean of the two controls (other [NiFe], other CoxL), as in the
                       manuscript's excess-change definition
"""
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

D = Path("lifeonair/results/benchmark_r2")
rows = list(csv.DictReader(open(D / "mock_ge_vs_truth.tsv"), delimiter="\t"))
l2 = lambda x: math.log2(x) if x > 0 else float("nan")


def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2
        i = j + 1
    return r


def pearson(a, b):
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    sa = math.sqrt(sum((x - ma) ** 2 for x in a)); sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb) if sa and sb else float("nan")


acc = defaultdict(list)
for r in rows:
    t, e = float(r["true_ge"]), float(r["est_ge"])
    acc[(r["db"], r["id_floor"], r["group"])].append((t, e, r["compartment"]))
with open(D / "mock_ge_accuracy.tsv", "w") as o:
    o.write("db\tid_floor\tgroup\tn\tmedian_true\tmedian_est\tmedian_log2_bias\tmin_log2_bias\tmax_log2_bias\t"
            "median_bias_bulk\tmedian_bias_root\tspearman_log\tpearson_log\n")
    for k, v in sorted(acc.items()):
        ok = [(t, e, c) for t, e, c in v if t > 0 and e > 0]
        b = [l2(e / t) for t, e, _ in ok]
        bb = [l2(e / t) for t, e, c in ok if c == "Deep"]; br = [l2(e / t) for t, e, c in ok if c == "PlantRoot"]
        lt = [l2(t) for t, _, _ in ok]; le = [l2(e) for _, e, _ in ok]
        f = lambda x: f"{x:.4f}"
        o.write("\t".join([*k, str(len(ok)), f(statistics.median(t for t, _, _ in v)), f(statistics.median(e for _, e, _ in v)),
                           f(statistics.median(b)) if b else "NA", f(min(b)) if b else "NA", f(max(b)) if b else "NA",
                           f(statistics.median(bb)) if bb else "NA", f(statistics.median(br)) if br else "NA",
                           f(pearson(rank(lt), rank(le))) if len(ok) > 2 else "NA",
                           f(pearson(lt, le)) if len(ok) > 2 else "NA"]) + "\n")

vis = defaultdict(dict)
for r in rows:
    vis[(r["db"], r["id_floor"], r["site"], r["trip"])][(r["compartment"], r["group"])] = (float(r["true_ge"]), float(r["est_ge"]))
with open(D / "mock_ge_root_bulk.tsv", "w") as o:
    o.write("db\tid_floor\tsite\ttrip\tgroup\ttrue_log2_root_bulk\test_log2_root_bulk\terror\n")
    for k, v in sorted(vis.items()):
        lr = {}
        for g in ("hydB_1l", "hydB_1h", "hydB_2a", "hydB_other_nife", "coxL_formI", "coxL_other"):
            if ("Deep", g) in v and ("PlantRoot", g) in v:
                (tb, eb), (tr, er) = v[("Deep", g)], v[("PlantRoot", g)]
                lr[g] = (l2(tr / tb) if tb > 0 and tr > 0 else float("nan"), l2(er / eb) if eb > 0 and er > 0 else float("nan"))
        if all(g in lr for g in ("hydB_1l", "hydB_other_nife", "coxL_other")):
            lr["excess_hydB_1l"] = tuple(lr["hydB_1l"][i] - (lr["hydB_other_nife"][i] + lr["coxL_other"][i]) / 2 for i in (0, 1))
        for g, (t, e) in lr.items():
            o.write("\t".join([*k, g, f"{t:.4f}", f"{e:.4f}", f"{e - t:.4f}"]) + "\n")
for p in ("mock_ge_accuracy.tsv", "mock_ge_root_bulk.tsv"):
    print(f"== {p}"); print(open(D / p).read()[:6000])
