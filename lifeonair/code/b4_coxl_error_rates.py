#!/usr/bin/env python3
"""Life on air (P3): CoxL form-I / other-CoxL misassignment rates for per-sample deconvolution (review round 2, R2-04).

    python3 lifeonair/code/b4_coxl_error_rates.py   (from the repository root)

Reads the held-out benchmark confusion tables (round 1: genome and clu50 schemes; round 2: comp50 and comp70,
sequence-disjoint over EQ and reference proteins) and writes, per scheme and identity floor, the 2x2 matrix
P(assigned a | true t) for t, a in {form I, other CoxL}, the rates of loss (unassigned or another label) and
the rate at which non-marker fragments are assigned to each CoxL class. With observed GE vector o and true
vector x, o = A x with A[a][t] = P(assigned a | true t); the corrected estimate is x = A^-1 o (per-length
normalisation is common to both classes, so the fragment-level rates carry over to RPK-based GE).
"""
import csv
from pathlib import Path

RES = Path("lifeonair/results")
SRC = [(RES / "benchmark/classifier_benchmark_confusion.tsv", ("genome", "clu50")),
       (RES / "benchmark_r2/disjoint_benchmark_confusion.tsv", ("comp50", "comp70"))]
OUT = RES / "benchmark_r2/coxl_error_rates.tsv"
C = ("coxL_formI", "coxL_other")

rows = []
for path, schemes in SRC:
    if not path.exists():
        print("missing", path); continue
    conf = list(csv.DictReader(open(path), delimiter="\t"))
    for s in schemes:
        for floor in ("50", "60", "70"):
            sub = [r for r in conf if r["scheme"] == s and r["id_floor"] == floor]
            if not sub:
                continue
            n = {t: next((int(r["reads_true"]) for r in sub if r["true"] == t), 0) for t in C + ("none",)}
            k = lambda t, a: sum(int(r["reads"]) for r in sub if r["true"] == t and r["assigned"] == a)
            A = {(a, t): k(t, a) / n[t] if n[t] else float("nan") for a in C for t in C}
            det = A[(C[0], C[0])] * A[(C[1], C[1])] - A[(C[0], C[1])] * A[(C[1], C[0])]
            rows.append([s, floor, n["coxL_formI"], n["coxL_other"], n["none"],
                         A[("coxL_formI", "coxL_formI")], A[("coxL_other", "coxL_formI")],
                         A[("coxL_formI", "coxL_other")], A[("coxL_other", "coxL_other")],
                         k("none", "coxL_formI") / n["none"] if n["none"] else float("nan"),
                         k("none", "coxL_other") / n["none"] if n["none"] else float("nan"), det])
with open(OUT, "w") as o:
    o.write("scheme\tid_floor\tfragments_formI\tfragments_other\tfragments_negative\tp_formI_given_formI\t"
            "p_other_given_formI\tp_formI_given_other\tp_other_given_other\tp_formI_given_negative\t"
            "p_other_given_negative\tdeterminant\n")
    for r in rows:
        o.write("\t".join(f"{x:.6f}" if isinstance(x, float) else str(x) for x in r) + "\n")
print(open(OUT).read())
