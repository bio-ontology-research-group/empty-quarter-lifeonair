#!/usr/bin/env python3
"""Life on air (P3): supplementary LaTeX table bodies from the result tables (revised after review rounds 1-2).

    python3 lifeonair/code/make_tables.py   (from the repository root)

Writes lifeonair/manuscript/tab_benchmark.tex, tab_mock.tex, tab_sensitivity.tex and tab_lineages.tex (rows
only; the table environments are in main.tex). Inputs: results/benchmark/classifier_benchmark_summary.tsv,
results/benchmark_r2/disjoint_benchmark_summary.tsv (70%-identity components, round 2),
results/benchmark_r2/mock_ge_accuracy.tsv and mock_ge_root_bulk.tsv, results/h123_v2/sensitivity.tsv and
results/h123_v2/h1_lineages_1l_carriers.tsv.
"""
import csv
from pathlib import Path

RES = Path("lifeonair/results")
OUT = Path("lifeonair/manuscript")
LABEL = {"primary": "Primary (60\\% identity)", "id50": "50\\% identity", "id70": "70\\% identity",
         "noflagged": "Contamination-flagged removed", "strict1l": "Weak group 1l calls removed",
         "without_trip5": "Without expedition 5", "without_trip4": "Without expedition 4",
         "motifonly": "Motif-less form-I CoxL unresolved", "coxLcorr": "CoxL misassignment-corrected"}
GROUP = {"hydB_1l": "{[NiFe]} group 1l", "hydB_1h": "{[NiFe]} group 1h", "hydB_2a": "{[NiFe]} group 2a",
         "other_nife": "Other {[NiFe]}", "coxL_formI": "CoxL form I", "coxL_other": "Other CoxL"}


def sci(p: float) -> str:
    m, e = f"{p:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))

b = {(r["scheme"], r["group"]): r for r in rd(RES / "benchmark/classifier_benchmark_summary.tsv") if r["id_floor"] == "60"}
b.update({(r["scheme"], r["group"]): r for r in rd(RES / "benchmark_r2/disjoint_benchmark_summary.tsv")
          if r["id_floor"] == "60" and r["scheme"] == "comp70"})
with open(OUT / "tab_benchmark.tex", "w") as o:
    for g, lab in GROUP.items():
        o.write(" & ".join([lab] + [f"{float(b[(s, g)][k]):.2f}" for s in ("genome", "clu50", "comp70")
                                    for k in ("recall", "precision")]) + "\\\\\n")

m = {(r["db"], r["group"]): r for r in rd(RES / "benchmark_r2/mock_ge_accuracy.tsv") if r["id_floor"] == "60"}
err = {}
for r in rd(RES / "benchmark_r2/mock_ge_root_bulk.tsv"):
    if r["id_floor"] == "60" and r["error"] != "nan":
        err.setdefault((r["db"], r["group"]), []).append(float(r["error"]))
with open(OUT / "tab_mock.tex", "w") as o:
    for g0, lab in GROUP.items():
        g = "hydB_other_nife" if g0 == "other_nife" else g0
        row = [lab, f"{float(m[('full', g)]['median_true']):.3f}"]
        for db in ("full", "novel"):
            e = err.get((db, g), [])
            row += [f"{float(m[(db, g)]['median_log2_bias']):.2f}", f"{float(m[(db, g)]['pearson_log']):.2f}",
                    f"{sum(e) / len(e):.2f} ({len(e)})" if e else "--"]
        o.write(" & ".join(row) + "\\\\\n")

for secondary in (False, True):
    name = "tab_sensitivity_secondary.tex" if secondary else "tab_sensitivity.tex"
    with open(OUT / name, "w") as o:
        for r in rd(RES / "h123_v2/sensitivity.tsv"):
            row = [LABEL.get(r["variant"], r["variant"])]
            if secondary:
                row += [f"{float(r[k]):.2f}" for k in ("root_excess_1h", "root_excess_formI", "surface_excess_1l", "rho_bulk_1l_ERA5", "rho_bulk_1l_POWER")]
            else:
                row += [r["n_samples"], f"{float(r['root_1l']):.2f}", f"{float(r['root_excess_1l']):.2f}", sci(float(r["root_excess_1l_p"])), f"{float(r['root_excess_1l_vs_coxLother']):.2f}"]
            o.write(" & ".join(row) + "\\\\\n")

with open(OUT / "tab_lineages.tex", "w") as o:
    for r in rd(RES / "h123_v2/h1_lineages_1l_carriers.tsv"):
        phylum, order = r["order"].split(";")
        o.write(" & ".join([phylum, order, f"{float(r['deep_mean_pct_mapped']):.2f}",
                            f"{float(r['median_log2fc_plantroot_vs_deep']):.2f}",
                            f"{100 * float(r['frac_cells_below0']):.0f}\\%"]) + "\\\\\n")
print("tables written to", OUT)
