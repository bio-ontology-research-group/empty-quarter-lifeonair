#!/bin/bash
# Life on air (P3): held-out benchmark of the read classifier (review round 1, finding R1-20).
#
# How often is a 100-nt read from a hydrogenase or CoxL gene assigned to the right label when its source
# is absent from the competitive database? Two hold-out schemes, five folds each:
#   genome  hold out all marker proteins of 20% of the representatives (relatives stay in the database)
#   clu50   cluster all EQ marker proteins at 50% identity (MMseqs2, coverage 0.5) and hold out 20% of the
#           clusters, so that no database protein shares >= 50% identity with the read's source protein
# For each fold the database is rebuilt without the held-out EQ proteins (HydDB and Cordero references
# and the SCG set stay), 100-nt fragments are tiled every 30 nt over each held-out gene (alternating
# strand), and 50 random non-marker CDS per held-out genome are tiled the same way as negatives. Reads are
# searched exactly as in r2_reads_array.sh (DIAMOND blastx --id 50 --evalue 1e-5 --max-target-seqs 1,
# best hit with alignment >= 25 aa) and the assigned label is tabulated against the true label at the
# 50/60/70% identity floors. Labels are those of marker_labels_v2.tsv (corrected CoxL forms, group 2a).
#
#   sbatch b1_classifier_benchmark.sh      (workbench 06_atmospheric/lifeonair_review_r1)
#
#SBATCH --job-name=loa_bench
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1/logs/bench_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1/logs/bench_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
W=$R/06_atmospheric/lifeonair_review_r1
DB=$R/06_atmospheric/lifeonair_reads/db
B=$W/bench
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
DIAMOND=$CONDA/eggnog/bin/diamond
MMSEQS=$CONDA/mmseqs2/bin/mmseqs
T=${SLURM_CPUS_PER_TASK:-32}
source "$CONDA/../bin/activate" bakta   # any env with python3
mkdir -p "$B" "$W/logs" "$W/results"
LAB=$W/inputs/marker_labels_v2.tsv       # copied from lifeonair/results/reads/db/marker_labels_v2.tsv
[ -s "$LAB" ] || { echo "ERROR: $LAB missing"; exit 1; }

echo "[$(date -Iseconds)] cluster EQ marker proteins at 50% identity"
awk '/^>/{p = ($0 ~ /^>EQ\|/)} p' "$DB/markers.faa" > "$B/eq_markers.faa"
"$MMSEQS" easy-cluster "$B/eq_markers.faa" "$B/clu50" "$B/tmp_clu" --min-seq-id 0.5 -c 0.5 \
    --threads "$T" >/dev/null 2>&1
rm -rf "$B/tmp_clu"

echo "[$(date -Iseconds)] folds, gene nucleotide sequences and fragments"
python3 - "$R" "$B" "$LAB" <<'PY'
import csv, os, random, sys
R, B, LAB = sys.argv[1:4]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}

def resolve(rep):                              # same resolver as r1_build_db.sh
    sample = rep.split("__")[0]; m = meta.get(sample)
    if "__" not in rep or not m:
        return None
    tok, local = rep.split("__", 1)[1].split("_", 1)
    anno = "single_sample" if tok == "single" else "self_supervised"
    base = f"{TREE}/{m['trip']}/site{m['site']}/{sample}/wgs/annotation/{anno}/bakta"
    for lc in [local] + ([local[:-4]] if local.endswith("_sub") else []):
        p = f"{base}/{lc}/{lc}.ffn"
        if os.path.exists(p):
            return p
    return None

def read_fa(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

label = {r["subject"]: (r["marker"], r["label"]) for r in csv.DictReader(open(LAB), delimiter="\t")
         if r["subject"].startswith("EQ|")}
clu = {m: c for c, m in (l.rstrip("\n").split("\t") for l in open(f"{B}/clu50_cluster.tsv"))}
reps = sorted({s.split("|")[2] for s in label})
rng = random.Random(7)
rfold = {r: i % 5 for i, r in enumerate(rng.sample(reps, len(reps)))}
cl = sorted(set(clu.values()))
cfold = {c: i % 5 for i, c in enumerate(rng.sample(cl, len(cl)))}
with open(f"{B}/folds.tsv", "w") as o:
    o.write("subject\tmarker\tlabel\tcluster50\tfold_genome\tfold_clu50\n")
    for s, (m, l) in sorted(label.items()):
        o.write(f"{s}\t{m}\t{l}\t{clu[s]}\t{rfold[s.split('|')[2]]}\t{cfold[clu[s]]}\n")

def tile(seq, step=30, L=100):
    rc = str.maketrans("ACGTN", "TGCAN")
    for k, i in enumerate(range(0, max(len(seq) - L, 0) + 1, step)):
        f = seq[i:i + L]
        yield f if k % 2 == 0 else f.translate(rc)[::-1]

by_rep = {}
for s in label:
    by_rep.setdefault(s.split("|")[2], set()).add(s.split("|")[3])
nmiss = 0
with open(f"{B}/fragments.fa", "w") as o, open(f"{B}/fragments_truth.tsv", "w") as t:
    t.write("read\tsource\tmarker\tlabel\n")
    for rep, pids in sorted(by_rep.items()):
        p = resolve(rep)
        if not p:
            nmiss += 1; continue
        seqs = dict(read_fa(p))
        for pid in sorted(pids):
            subj = [k for k in (f"EQ|hydB|{rep}|{pid}", f"EQ|coxL|{rep}|{pid}") if k in label][0]
            g = seqs.get(pid)
            if not g:
                continue
            for j, f in enumerate(tile(g.upper())):
                rid = f"{subj}#{j}"
                o.write(f">{rid}\n{f}\n"); t.write(f"{rid}\t{subj}\t{label[subj][0]}\t{label[subj][1]}\n")
        others = [k for k in seqs if k not in pids and len(seqs[k]) >= 300]
        for pid in rng.sample(others, min(50, len(others))):
            for j, f in enumerate(tile(seqs[pid].upper(), step=150)):
                rid = f"NEG|{rep}|{pid}#{j}"
                o.write(f">{rid}\n{f}\n"); t.write(f"{rid}\tNEG|{rep}|{pid}\tnone\tnone\n")
print(f"reps {len(reps)}; unresolved {nmiss}; clusters {len(cl)}", flush=True)
PY

for SCHEME in genome clu50; do
  for F in 0 1 2 3 4; do
    D=$B/${SCHEME}_f$F; mkdir -p "$D"
    [ -s "$D/hits.tsv" ] && continue
    COL=$([ "$SCHEME" = genome ] && echo 5 || echo 6)
    awk -F'\t' -v c=$COL -v f=$F 'NR > 1 && $c == f {print $1}' "$B/folds.tsv" > "$D/heldout.txt"
    awk 'NR == FNR {h[">" $1]; next} /^>/{p = !(($1) in h)} p' "$D/heldout.txt" "$DB/loa_reads.faa" > "$D/db.faa"
    "$DIAMOND" makedb --in "$D/db.faa" -d "$D/db" -p "$T" --quiet
    # reads: fragments of the held-out marker genes plus negatives from the genomes of held-out proteins
    python3 - "$D/heldout.txt" "$B/fragments.fa" "$D/reads.fa" <<'PY'
import sys
held = {l.strip() for l in open(sys.argv[1])}
reps = {h.split("|")[2] for h in held}
keep = False
with open(sys.argv[3], "w") as o:
    for line in open(sys.argv[2]):
        if line.startswith(">"):
            rid = line[1:].strip(); src = rid.split("#")[0]
            keep = src in held or (src.startswith("NEG|") and src.split("|")[1] in reps)
        if keep:
            o.write(line)
PY
    "$DIAMOND" blastx -q "$D/reads.fa" -d "$D/db" -p "$T" --id 50 --evalue 1e-5 --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length evalue bitscore --quiet -o "$D/hits.tsv"
    rm -f "$D/db.faa" "$D/db.dmnd"
    echo "[$(date -Iseconds)] $SCHEME fold $F: $(grep -c '^>' "$D/reads.fa") reads, $(cut -f1 "$D/hits.tsv" | sort -u | wc -l) hit"
  done
done

echo "[$(date -Iseconds)] tabulate"
python3 - "$B" "$DB/marker_labels.tsv" "$LAB" "$W/results" <<'PY'
import csv, sys
from collections import Counter
B, LAB0, LAB, OUT = sys.argv[1:5]
lab = {r["subject"]: r["label"] for r in csv.DictReader(open(LAB0), delimiter="\t")}
lab.update({r["subject"]: r["label"] for r in csv.DictReader(open(LAB), delimiter="\t")})
def group(l):
    for k, g in (("Group_1l", "hydB_1l"), ("Group_1h", "hydB_1h"), ("Group_2a", "hydB_2a")):
        if l.endswith(k): return g
    if l == "form_I": return "coxL_formI"
    if l in ("form_II", "coxL_other", "coxL_unresolved"): return "coxL_other"
    if l in ("none",): return "none"
    if "NiFe" in l or "FeFe" in l or "Fe" in l: return "other_nife"
    return "SCG"
truth = {r["read"]: group(r["label"]) for r in csv.DictReader(open(f"{B}/fragments_truth.tsv"), delimiter="\t")}
rows = []
for scheme in ("genome", "clu50"):
    for floor in (50, 60, 70):
        c = Counter(); n = Counter()
        for f in range(5):
            D = f"{B}/{scheme}_f{f}"
            held = {l.strip() for l in open(f"{D}/heldout.txt")}
            reps = {h.split("|")[2] for h in held}
            best = {}
            for line in open(f"{D}/hits.tsv"):
                q, s, pid, aln = line.split("\t")[:4]
                if q not in best and int(aln) >= 25:
                    best[q] = (s, float(pid))
            for rid, t in truth.items():
                src = rid.split("#")[0]
                if not (src in held or (src.startswith("NEG|") and src.split("|")[1] in reps)):
                    continue
                n[t] += 1
                h = best.get(rid)
                a = group(lab[h[0]]) if h and h[1] >= floor else "unassigned"
                c[(t, a)] += 1
        for (t, a), k in sorted(c.items()):
            rows.append((scheme, floor, t, a, k, n[t], k / n[t]))
with open(f"{OUT}/classifier_benchmark_confusion.tsv", "w") as o:
    o.write("scheme\tid_floor\ttrue\tassigned\treads\treads_true\tfraction\n")
    for r in rows: o.write("\t".join(map(str, r[:-1])) + f"\t{r[-1]:.4f}\n")
# precision / recall per marker group (negatives and SCG counted as possible false positives)
with open(f"{OUT}/classifier_benchmark_summary.tsv", "w") as o:
    o.write("scheme\tid_floor\tgroup\trecall\tprecision\tassigned_reads\n")
    for scheme in ("genome", "clu50"):
        for floor in (50, 60, 70):
            sub = [r for r in rows if r[0] == scheme and r[1] == floor]
            for g in ("hydB_1l", "hydB_1h", "hydB_2a", "other_nife", "coxL_formI", "coxL_other"):
                tp = sum(r[4] for r in sub if r[2] == g and r[3] == g)
                pos = sum(r[4] for r in sub if r[2] == g)
                ass = sum(r[4] for r in sub if r[3] == g)
                o.write(f"{scheme}\t{floor}\t{g}\t{tp / pos if pos else float('nan'):.4f}\t"
                        f"{tp / ass if ass else float('nan'):.4f}\t{ass}\n")
print(open(f"{OUT}/classifier_benchmark_summary.tsv").read())
PY
sha256sum "$W/results/classifier_benchmark_"*.tsv "$B/folds.tsv" > "$W/results/classifier_benchmark.sha256"
echo "[$(date -Iseconds)] done"
