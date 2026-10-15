#!/usr/bin/env python3
"""Life on air (P3): corrected CoxL form calls (review round 1, finding R1-01; 12 October 2026).

The original typing (IBEX lifeonair_clean/scripts/lifeonair_2_type_guild.sh) had the active-site motifs
inverted: it called AYRGAGR form I and AYR[CG]SFR form II. Form I CoxL carries AYXCSFR and form II carries
AYRGAGR (King & Weber 2007, Nat Rev Microbiol 5:107); the three form-I references of Cordero et al. 2019
used here carry AYRCSFR, AYACSFR and AYRCSFR. Corrected rule, per protein:
  - [AS]Y.CS[FL]R present                 -> form_I (AYXCSFR and the variants SYRCSFR, SYACSFR and AYRCSLR
                                            reported by Imaura et al. 2026; review round 2)
  - AYRGAGR present (and no AY.CSFR)      -> form_II
  - neither motif (motif region missing)  -> form_I if > 55% identity and >= 90% query coverage to the
                                             Cordero form-I references, otherwise unresolved
Reads are relabelled without re-running DIAMOND: every EQ CoxL protein was already a separate subject in
the competitive database, so only the label attached to each subject changes.

    python3 lifeonair/code/fix_coxL_v2.py <coxL.faa from IBEX lifeonair_clean/queries>

Writes lifeonair/results/inputs/coxL_class_v2.tsv, lifeonair/results/inputs/guild_coxL_reps_v2.txt and
lifeonair/results/reads/db/marker_labels_v2.tsv (also separates group 2a hydrogenases, R1-02, by keeping
their HydDB label; r3_aggregate.py maps it to its own column).
"""
import csv
import re
import sys
from collections import Counter

IN, DB = "lifeonair/results/inputs", "lifeonair/results/reads/db"
FORM_I, FORM_II = re.compile(r"[AS]Y.CS[FL]R"), re.compile(r"AYRGAGR")


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


def main(faa):
    old = {r["query"]: r for r in csv.DictReader(open(f"{IN}/coxL_class.tsv"), delimiter="\t")}
    calls, n = {}, Counter()
    with open(f"{IN}/coxL_class_v2.tsv", "w") as o:
        o.write("query\tlength\tpident_formI\tqcov_formI\tmotif\tform_call\told_call\n")
        for name, seq in read_fasta(faa):
            r = old[name]
            pid, qc = float(r["pident_formI"]), float(r["qcov_formI"])
            if FORM_I.search(seq):
                mot, call = "AYxCSFR", "form_I"
            elif FORM_II.search(seq):
                mot, call = "AYRGAGR", "form_II"
            else:
                mot = "none"
                call = "form_I" if (pid > 55.0 and qc >= 90.0) else "unresolved"
            calls[name] = call
            n[(mot, call, r["form_I_call"])] += 1
            o.write(f"{name}\t{len(seq)}\t{pid:.1f}\t{qc:.1f}\t{mot}\t{call}\t{r['form_I_call']}\n")
    for k, v in sorted(n.items(), key=lambda x: -x[1]):
        print(*k, v, sep="\t")
    reps = sorted({q.split("|")[0] for q, c in calls.items() if c == "form_I"})
    open(f"{IN}/guild_coxL_reps_v2.txt", "w").write("\n".join(reps) + "\n")
    print(f"form_I proteins {sum(c == 'form_I' for c in calls.values())}; form_I reps {len(reps)}")
    rows = list(csv.DictReader(open(f"{DB}/marker_labels.tsv"), delimiter="\t"))
    with open(f"{DB}/marker_labels_v2.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            if r["marker"] == "coxL" and r["subject"].startswith("EQ|"):
                q = "|".join(r["subject"].split("|")[2:4])
                r["label"] = {"form_I": "form_I", "form_II": "form_II"}.get(calls[q], "coxL_unresolved")
            w.writerow(r)


if __name__ == "__main__":
    main(sys.argv[1])
