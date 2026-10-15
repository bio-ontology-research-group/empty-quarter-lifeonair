#!/usr/bin/env python3
"""Life on air (P3): carrier lists and trait inputs restricted to retained contigs (review round 1, R1-19).

s1_sub_reconcile.sh (IBEX job 53614016) mapped every CDS of the 1,010 "_sub" representatives to whether its
contig is retained in the representative. Proteins on removed contigs are genuine community genes, so reads
assigned to them keep their label and the read-based genome equivalents are unchanged; but they must not make
the representative a carrier, nor enter its gRodon estimate. This script writes:
  inputs/atmo_classification_v2.tsv  typed hydB/coxL rows minus proteins on removed contigs and orphan reps;
                                     CoxL rows relabelled with the corrected form calls (R2-06)
  inputs/guild_coxL_reps_v2.txt      reps with >= 1 retained form-I CoxL (corrected calls, fix_coxL_v2.py)
  inputs/rep_growth_v2.tsv           rep_growth.tsv with _sub rows replaced by the retained-CDS estimates

    python3 lifeonair/code/apply_sub_fix.py
"""
import csv
import gzip

IN = "lifeonair/results/inputs"
ret = {}
for r in csv.DictReader(gzip.open(f"{IN}/sub_protein_retention.tsv.gz", "rt"), delimiter="\t"):
    ret[(r["rep"], r["protein_id"])] = r["retained"] == "1"
orph = set(open("lifeonair/results/reads/orphan_sub_reps.txt").read().split())
keep = lambda rep, pid: rep not in orph and ret.get((rep, pid), True)

rows = list(csv.DictReader(open(f"{IN}/atmo_classification.tsv"), delimiter="\t"))
kept = [r for r in rows if keep(r["rep_id"], r["protein_id"])]
# CoxL rows carry the corrected form calls of fix_coxL_v2.py (review round 2, R2-06): "<call>|<motif>", with
# call form_I / form_II / unresolved; is_high_affinity = 1 only for form I
cxc = {tuple(r["query"].split("|")): r for r in csv.DictReader(open(f"{IN}/coxL_class_v2.tsv"), delimiter="\t")}
for r in kept:
    if r["marker"] == "coxL":
        c = cxc[(r["rep_id"], r["protein_id"])]
        r["assignment"] = f"{c['form_call']}|{c['motif']}"; r["is_high_affinity"] = str(int(c["form_call"] == "form_I"))
with open(f"{IN}/atmo_classification_v2.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(kept)
print(f"atmo_classification: {len(rows)} -> {len(kept)} rows")
for g in ("Group_1l", "Group_1h", "Group_2a"):
    b = {r["rep_id"] for r in rows if r["marker"] == "hydB" and r["assignment"].endswith(g)}
    a = {r["rep_id"] for r in kept if r["marker"] == "hydB" and r["assignment"].endswith(g)}
    print(f"  {g} carriers {len(b)} -> {len(a)}")

cx = [r for r in csv.DictReader(open(f"{IN}/coxL_class_v2.tsv"), delimiter="\t") if r["form_call"] == "form_I"]
reps = sorted({r["query"].split("|")[0] for r in cx if keep(*r["query"].split("|"))})
old = open(f"{IN}/guild_coxL_reps_v2.txt").read().split()
open(f"{IN}/guild_coxL_reps_v2.txt", "w").write("\n".join(reps) + "\n")
print(f"form-I carriers {len(old)} -> {len(reps)} (proteins {sum(keep(*r['query'].split('|')) for r in cx)})")

fix = {r["rep"]: r for r in csv.DictReader(open(f"{IN}/rep_growth_sub_fixed.tsv"), delimiter="\t")}
g = list(csv.DictReader(open(f"{IN}/rep_growth.tsv"), delimiter="\t"))
with open(f"{IN}/rep_growth_v2.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(g[0]), delimiter="\t", lineterminator="\n"); w.writeheader()
    for r in g:
        w.writerow(fix.get(r["rep"], r))
print(f"growth: {len(fix)} _sub reps re-estimated")
