#!/usr/bin/env python3
"""Life on air (P3), read-based guild abundance -- stage R3: per-sample genome-equivalent table.

For each sample and identity floor (50/60/70%), reads per kilobase (RPK, subject length in nt)
are summed per label; each marker is divided by the median RPK of the 14 single-copy ribosomal
protein genes, giving copies per genome equivalent. For single-copy markers this approximates the
fraction of genomes carrying the capacity.

    python3 lifeonair/code/r3_aggregate.py lifeonair/results/reads [EXCLUDE_REPS OUT_NAME [LABELS]]

LABELS (default db/marker_labels.tsv) selects the label table. With db/marker_labels_v2.tsv (corrected
CoxL forms, fix_coxL_v2.py; review round 1) group 2a [NiFe] is reported in its own column hydB_2a and is
no longer part of the hydB_other_nife control (R1-02). Use "-" for EXCLUDE_REPS to exclude nothing.

With EXCLUDE_REPS (one rep id or "rep|protein" id per line), reads whose best hit is an EQ protein of an
excluded rep, or an excluded protein, are discarded before normalisation (contamination and weak-call
sensitivity); output goes to OUT_NAME.
"""
import csv
import statistics
import sys
from collections import Counter
from pathlib import Path

HA = {"_NiFe__Group_1h": "hydB_1h", "_NiFe__Group_1l": "hydB_1l"}


def main(root: Path, exclude: set, out_name: str, labels: str = "db/marker_labels.tsv") -> None:
    v2 = "_v2" in Path(labels).name    # also the _v2_motifonly sensitivity labels
    ha = dict(HA, **({"_NiFe__Group_2a": "hydB_2a"} if v2 else {}))
    lab = {r["subject"]: (r["marker"], r["label"])
           for r in csv.DictReader(open(root / labels), delimiter="\t")}
    nt = {r["subject"]: 3 * int(r["length_aa"])
          for r in csv.DictReader(open(root / "db/subject_lengths.tsv"), delimiter="\t")}
    cols = ["sample", "id_floor", "reads_screened", "scg_genes_detected", "scg_median_rpk",
            "hydB_ha", "hydB_1l", "hydB_1h", "hydB_other_nife", "coxL_formI", "coxL_other",
            "reads_hydB_1l", "reads_hydB_1h", "reads_hydB_other_nife", "reads_coxL_formI",
            "reads_coxL_other", "reads_scg"] + (["hydB_2a", "reads_hydB_2a"] if v2 else [])
    out = open(root / out_name, "w")
    out.write("\t".join(cols) + "\n")
    n_done = 0
    for sd in sorted((root / "per_sample").iterdir()):
        if not (sd / "subject_counts.tsv").exists():
            continue
        meta = next(csv.DictReader(open(sd / "meta.tsv"), delimiter="\t"))
        rows = [line.rstrip("\n").split("\t") for line in open(sd / "subject_counts.tsv")]
        n_done += 1
        for k, floor in ((1, 50), (2, 60), (3, 70)):
            rpk, reads = Counter(), Counter()
            for r in rows:
                if exclude and r[0].startswith("EQ|") and (r[0].split("|")[2] in exclude or "|".join(r[0].split("|")[2:4]) in exclude):
                    continue
                m, label = lab[r[0]]
                c = int(r[k])
                if not c:
                    continue
                key = (m, label)
                if m == "hydB":
                    key = ("hydB", ha.get(label, "hydB_other_nife" if "NiFe" in label else "other"))
                elif m == "coxL":
                    key = ("coxL", "coxL_formI" if label == "form_I" else "coxL_other")
                rpk[key] += c / (nt[r[0]] / 1000)
                reads[key] += c
            scg = [v for (m, _), v in rpk.items() if m == "SCG"]
            med = statistics.median(scg) if len(scg) >= 10 else float("nan")
            ge = lambda key: rpk[key] / med if med == med and med > 0 else float("nan")
            ha_ge = ge(("hydB", "hydB_1l")) + ge(("hydB", "hydB_1h"))
            out.write("\t".join(str(x) for x in [
                meta["sample"], floor, meta["reads_screened"], len(scg), f"{med:.2f}",
                f"{ha_ge:.5f}", f"{ge(('hydB', 'hydB_1l')):.5f}", f"{ge(('hydB', 'hydB_1h')):.5f}",
                f"{ge(('hydB', 'hydB_other_nife')):.5f}", f"{ge(('coxL', 'coxL_formI')):.5f}",
                f"{ge(('coxL', 'coxL_other')):.5f}",
                reads[("hydB", "hydB_1l")], reads[("hydB", "hydB_1h")], reads[("hydB", "hydB_other_nife")],
                reads[("coxL", "coxL_formI")], reads[("coxL", "coxL_other")],
                sum(v for (m, _), v in reads.items() if m == "SCG")]
                + ([f"{ge(('hydB', 'hydB_2a')):.5f}", reads[("hydB", "hydB_2a")]] if v2 else [])) + "\n")
    out.close()
    print(f"aggregated {n_done} samples -> {root / out_name}")


if __name__ == "__main__":
    excl = set(open(sys.argv[2]).read().split()) if len(sys.argv) > 2 and sys.argv[2] != "-" else set()
    main(Path(sys.argv[1]), excl, sys.argv[3] if len(sys.argv) > 3 else "guild_genome_equivalents.tsv",
         sys.argv[4] if len(sys.argv) > 4 else "db/marker_labels.tsv")
