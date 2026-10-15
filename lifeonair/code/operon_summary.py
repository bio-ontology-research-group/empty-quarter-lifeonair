#!/usr/bin/env python3
"""Life on air (P3), analysis 5 summary: operon features of typed hydB/coxL by class.

Reports all loci and, separately, loci with >= 5 genes on both sides on the contig (complete
context), because absent partners near a contig end mostly reflect assembly fragmentation.
Also lists the most frequent neighbours within +/-5 genes of high-affinity hydB.

    python3 lifeonair/code/operon_summary.py lifeonair/results/operon [lifeonair/results/inputs/coxL_class_v2.tsv]

With the corrected CoxL calls (review round 1, R1-01) as second argument, CoxL loci are relabelled form I /
form II / unresolved, group 2a is reported separately (R1-02), and the outputs get a _v2 suffix. In that mode
loci on contigs removed from "_sub" representatives and loci of the two orphan representatives are dropped
(review round 2, R2-06), so denominators equal the retained typed proteins.
"""
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path


V2 = {}
DROP = set()


def klass(r):
    c = r["class"]
    if r["marker"] == "hydB":
        if V2:
            return ("hydB_1h" if c.endswith("Group_1h") else "hydB_1l" if c.endswith("Group_1l")
                    else "hydB_2a" if c.endswith("Group_2a") else "hydB_other")
        return "hydB_1h" if c.endswith("Group_1h") else "hydB_1l" if c.endswith("Group_1l") else "hydB_other"
    if V2:
        return {"form_I": "coxL_formI", "form_II": "coxL_formII"}.get(V2[f"{r['rep']}|{r['protein_id']}"], "coxL_unresolved")
    return "coxL_formI" if c == "form_I" else "coxL_other"


def pct(v, f):
    return f"{100 * sum(map(f, v)) / len(v):.1f}" if v else "nan"


def main(d: Path, sfx: str = "") -> None:
    rows = [r for r in csv.DictReader(open(d / "marker_context.tsv"), delimiter="\t")
            if (r["rep"], r["protein_id"]) not in DROP and r["rep"] not in DROP]
    groups = defaultdict(list)
    for r in rows:
        k = klass(r)
        groups[(k, "all")].append(r)
        if int(r["genes_to_contig_end"]) >= 5:
            groups[(k, "complete_context")].append(r)
    with open(d / f"operon_summary{sfx}.tsv", "w") as o:
        o.write("class\tsubset\tn\tsmall_subunit_adjacent_pct\tany_maturation_pct\tcoxM_pct\tcoxS_pct\t"
                "coxMSL_pct\tcox_accessory_pct\tany_hypothetical_pct\n")
        for (k, s), v in sorted(groups.items()):
            o.write("\t".join([k, s, str(len(v)),
                               pct(v, lambda r: r["small_subunit_adjacent"] == "1"),
                               pct(v, lambda r: int(r["n_maturation"]) > 0),
                               pct(v, lambda r: r["coxM"] == "1"), pct(v, lambda r: r["coxS"] == "1"),
                               pct(v, lambda r: r["coxM"] == "1" and r["coxS"] == "1"),
                               pct(v, lambda r: int(r["n_cox_accessory"]) > 0),
                               pct(v, lambda r: int(r["n_hypothetical"]) > 0)]) + "\n")
    ctx = {(r["rep"], r["protein_id"]): klass(r) for r in rows}
    nb = Counter(); nloc = Counter()
    for r in csv.DictReader(gzip.open(d / "neighbours.tsv.gz", "rt"), delimiter="\t"):
        k = ctx.get((r["rep"], r["protein_id"]))
        if k is None:
            continue
        off = int(r["offset"])
        if off == 0:
            nloc[k] += 1
            continue
        if abs(off) > 5 or k not in ("hydB_1h", "hydB_1l"):
            continue
        name = r["eggnog_name"] if r["eggnog_name"] not in ("-", "") else (r["bakta_gene"] or r["bakta_product"])
        nb[(k, name)] += 1
    with open(d / f"hydB_ha_top_neighbours{sfx}.tsv", "w") as o:
        o.write("class\tneighbour\tcount\tper_locus\n")
        for (k, name), c in nb.most_common():
            if c >= 20:
                o.write(f"{k}\t{name}\t{c}\t{c / nloc[k]:.3f}\n")
    print(open(d / f"operon_summary{sfx}.tsv").read())


if __name__ == "__main__":
    if len(sys.argv) > 2:
        V2.update({r["query"]: r["form_call"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")})
        inp = Path(sys.argv[2]).parent
        DROP.update((r["rep"], r["protein_id"]) for r in csv.DictReader(
            gzip.open(inp / "sub_protein_retention.tsv.gz", "rt"), delimiter="\t") if r["retained"] == "0")
        DROP.update(open(inp.parent / "reads/orphan_sub_reps.txt").read().split())
    main(Path(sys.argv[1]), "_v2" if len(sys.argv) > 2 else "")
