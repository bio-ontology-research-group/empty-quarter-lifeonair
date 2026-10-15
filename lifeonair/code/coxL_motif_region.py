#!/usr/bin/env python3
"""Life on air (P3): active-site region of form-I CoxL calls without a motif (review round 2, R2-07).

fix_coxL_v2.py calls a CoxL protein form I without the AY.CSFR motif when it shares > 55% identity over >= 90%
query coverage with the Cordero et al. 2019 form-I references. Absence of the motif string can mean that the
region is missing (fragment) or that it is present but diverged. Each such protein is aligned locally
(BLOSUM62, gap open -10, extend -0.5) to the three references; the residues aligned to the reference motif are
reported and classified:
  region_missing  the best alignment does not span the reference motif
  divergent       the region is aligned but does not read AY.CSFR (or AYRGAGR)
Also writes a label table in which these proteins are "coxL_unresolved", for a GE sensitivity
(r3_aggregate.py ... marker_labels_v2_motifonly.tsv).

    python3 lifeonair/code/coxL_motif_region.py <IBEX lifeonair_clean/queries/coxL.faa> <lifeonair_clean/classify/coxL_reps_plus_refs.faa>
Writes lifeonair/results/inputs/coxL_motifless_formI_region.tsv, lifeonair/results/reads/db/marker_labels_v2_motifonly.tsv
"""
import csv
import re
import sys

from Bio import Align
from Bio.Align import substitution_matrices

IN, DB = "lifeonair/results/inputs", "lifeonair/results/reads/db"
REFS = ("WP_013913730.1", "WP_003892166.1", "WP_006067999.1")


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


seqs = {k: v for f in sys.argv[1:] for k, v in read_fasta(f)}
al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-10, extend_gap_score=-0.5)
targets = [r for r in csv.DictReader(open(f"{IN}/coxL_class_v2.tsv"), delimiter="\t")
           if r["form_call"] == "form_I" and r["motif"] == "none"]
out = []
for r in targets:
    q = seqs[r["query"]]
    best = None
    for ref in REFS:
        s = seqs[ref]; m = re.search(r"AY.CSFR", s)
        a = al.align(s, q)[0]
        # map reference motif positions to query residues through the aligned blocks
        res = []
        for i in range(m.start(), m.end()):
            hit = "-"
            for (rs, re_), (qs, qe) in zip(*a.aligned):
                if rs <= i < re_:
                    hit = q[qs + i - rs]; break
            res.append(hit)
        cov = sum(c != "-" for c in res)
        if best is None or a.score > best[0]:
            best = (a.score, ref, "".join(res), cov, a.aligned[1][0][0], a.aligned[1][-1][1])
    score, ref, region, cov, qs, qe = best
    cls = "region_missing" if cov < 5 else ("divergent" if not re.fullmatch(r"AY.CSFR", region) else "motif_aligned")
    out.append({"query": r["query"], "length": len(q), "pident_formI": r["pident_formI"], "qcov_formI": r["qcov_formI"],
                "best_ref": ref, "aligned_motif_region": region, "region_residues_aligned": cov,
                "query_aln_start": qs + 1, "query_aln_end": qe, "class": cls})
with open(f"{IN}/coxL_motifless_formI_region.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
drop = {"EQ|coxL|" + o["query"] for o in out}
rows = list(csv.DictReader(open(f"{DB}/marker_labels_v2.tsv"), delimiter="\t"))
n = 0
with open(f"{DB}/marker_labels_v2_motifonly.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader()
    for r in rows:
        if r["subject"] in drop:
            r["label"] = "coxL_unresolved"; n += 1
        w.writerow(r)
from collections import Counter
print(Counter(o["class"] for o in out), "relabelled", n)
for o in out:
    if o["length"] >= 600:
        print(o["query"], o["length"], o["aligned_motif_region"], o["class"])
