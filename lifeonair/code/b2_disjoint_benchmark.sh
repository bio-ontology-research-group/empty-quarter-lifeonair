#!/bin/bash
# Life on air (P3): sequence-disjoint benchmark of the read classifier (review round 2, finding R2-03).
#
# The round-1 benchmark (b1_classifier_benchmark.sh) held out clusters of EQ marker proteins only; HydDB and
# Cordero references stayed in the database. Here ALL marker proteins of the competitive database (EQ
# catalogue markers, HydDB references, Cordero form-I CoxL references) are grouped into components by
# single-linkage over DIAMOND all-vs-all blastp hits (identity >= T and alignment >= 50% of the shorter
# protein), for T = 50 and 70. Whole components go to five folds (greedy size balancing). For each fold every
# protein of the held-out components, EQ and external alike, is removed from the database (SCG set kept),
# 100-nt fragments are tiled every 30 nt over the held-out EQ genes and 50 random non-marker CDS per affected
# genome are tiled as negatives, and reads are searched exactly as in r2_reads_array.sh (DIAMOND blastx
# --id 50 --evalue 1e-5 --max-target-seqs 1; best hit, alignment >= 25 aa; floors 50/60/70%).
# Disjointness is verified per fold by blastp of the held-out EQ proteins against the reduced marker set,
# recording the best identity (alignment >= 50% of the shorter protein) and the best local identity over
# any alignment >= 33 aa (one read). Recall is also reported by nearest-retained-neighbour identity and by
# source phylum/order (GTDB R232).
#
#   sbatch b2_disjoint_benchmark.sh      (workbench 06_atmospheric/lifeonair_review_r2)
#
#SBATCH --job-name=loa_disj
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r2/logs/disj_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r2/logs/disj_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
W=$R/06_atmospheric/lifeonair_review_r2
DB=$R/06_atmospheric/lifeonair_reads/db
B=$W/disjoint
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
DIAMOND=$CONDA/eggnog/bin/diamond
T=${SLURM_CPUS_PER_TASK:-32}
source "$CONDA/../bin/activate" bakta   # any env with python3
mkdir -p "$B" "$W/logs" "$W/results"
LAB=$W/inputs/marker_labels_v2.tsv
SUBX=$W/inputs/sub_excluded_db_proteins.txt

echo "[$(date -Iseconds)] all-vs-all blastp of marker proteins"
cp "$DB/markers.faa" "$B/markers.faa"
"$DIAMOND" makedb --in "$B/markers.faa" -d "$B/markers" -p "$T" --quiet
"$DIAMOND" blastp -q "$B/markers.faa" -d "$B/markers" -p "$T" --id 40 --evalue 1e-5 --more-sensitive \
    --max-target-seqs 0 --outfmt 6 qseqid sseqid pident length qlen slen --quiet -o "$B/ava.tsv"

echo "[$(date -Iseconds)] components, folds, fragments"
python3 - "$R" "$B" "$LAB" "$SUBX" "$W/inputs" <<'PY'
import csv, os, random, sys
from collections import defaultdict
R, B, LAB, SUBX, INP = sys.argv[1:6]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}
subx = {l.strip() for l in open(SUBX) if l.strip()}

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

prots = [n for n, _ in read_fa(f"{B}/markers.faa")]
label = {r["subject"]: (r["marker"], r["label"]) for r in csv.DictReader(open(LAB), delimiter="\t")}
for T in (50, 70):
    par = {p: p for p in prots}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for line in open(f"{B}/ava.tsv"):
        q, s, pid, aln, ql, sl = line.split("\t")
        if q != s and float(pid) >= T and int(aln) >= 0.5 * min(int(ql), int(sl)):
            a, b = find(q), find(s)
            if a != b: par[a] = b
    comp = defaultdict(list)
    for p in prots: comp[find(p)].append(p)
    # greedy balancing on the number of EQ proteins per component
    order = sorted(comp, key=lambda c: (-sum(x.startswith("EQ|") for x in comp[c]), c))
    load = [0] * 5; fold = {}
    for c in order:
        f = min(range(5), key=lambda i: (load[i], i)); fold[c] = f
        load[f] += sum(x.startswith("EQ|") for x in comp[c])
    with open(f"{B}/comp{T}.tsv", "w") as o:
        o.write("subject\tmarker\tlabel\tcomponent\tcomponent_size\tfold\n")
        for c, mem in comp.items():
            for p in mem:
                m, l = label[p]
                o.write(f"{p}\t{m}\t{l}\t{c}\t{len(mem)}\t{fold[c]}\n")
    big = sorted(comp.values(), key=len, reverse=True)[:5]
    print(f"[comp{T}] components {len(comp)}; largest {[len(b) for b in big]}; EQ per fold {load}", flush=True)

def tile(seq, step=30, L=100):
    rc = str.maketrans("ACGTN", "TGCAN")
    for k, i in enumerate(range(0, max(len(seq) - L, 0) + 1, step)):
        f = seq[i:i + L]
        yield f if k % 2 == 0 else f.translate(rc)[::-1]

eq = [p for p in prots if p.startswith("EQ|") and "|".join(p.split("|")[2:4]) not in subx]
by_rep = defaultdict(set)
for s in eq: by_rep[s.split("|")[2]].add(s.split("|")[3])
rng = random.Random(7); nmiss = 0
with open(f"{B}/fragments.fa", "w") as o, open(f"{B}/fragments_truth.tsv", "w") as t, \
     open(f"{B}/eq_heldout_candidates.faa", "w") as pq:
    t.write("read\tsource\tmarker\tlabel\n")
    pseq = dict(read_fa(f"{B}/markers.faa"))
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
            pq.write(f">{subj}\n{pseq[subj]}\n")
            for j, f in enumerate(tile(g.upper())):
                rid = f"{subj}#{j}"
                o.write(f">{rid}\n{f}\n"); t.write(f"{rid}\t{subj}\t{label[subj][0]}\t{label[subj][1]}\n")
        others = [k for k in seqs if k not in pids and len(seqs[k]) >= 300]
        for pid in rng.sample(others, min(50, len(others))):
            for j, f in enumerate(tile(seqs[pid].upper(), step=150)):
                rid = f"NEG|{rep}|{pid}#{j}"
                o.write(f">{rid}\n{f}\n"); t.write(f"{rid}\tNEG|{rep}|{pid}\tnone\tnone\n")
print(f"EQ marker proteins (retained contigs) {len(eq)}; reps {len(by_rep)}; unresolved {nmiss}", flush=True)
PY

for SCHEME in comp50 comp70; do
  for F in 0 1 2 3 4; do
    D=$B/${SCHEME}_f$F; mkdir -p "$D"
    [ -s "$D/hits.tsv" ] && [ -s "$D/verify.tsv" ] && continue
    awk -F'\t' -v f=$F 'NR > 1 && $6 == f {print $1}' "$B/$SCHEME.tsv" > "$D/heldout.txt"
    awk 'NR == FNR {h[">" $1]; next} /^>/{p = !(($1) in h)} p' "$D/heldout.txt" "$DB/loa_reads.faa" > "$D/db.faa"
    "$DIAMOND" makedb --in "$D/db.faa" -d "$D/db" -p "$T" --quiet
    python3 - "$D/heldout.txt" "$B/fragments.fa" "$D/reads.fa" "$B/eq_heldout_candidates.faa" "$D/heldout_eq.faa" <<'PY'
import sys
held = {l.strip() for l in open(sys.argv[1])}
reps = {h.split("|")[2] for h in held if h.startswith("EQ|")}
keep = False
with open(sys.argv[3], "w") as o:
    for line in open(sys.argv[2]):
        if line.startswith(">"):
            rid = line[1:].strip(); src = rid.split("#")[0]
            keep = src in held or (src.startswith("NEG|") and src.split("|")[1] in reps)
        if keep:
            o.write(line)
keep = False
with open(sys.argv[5], "w") as o:
    for line in open(sys.argv[4]):
        if line.startswith(">"):
            keep = line[1:].strip() in held
        if keep:
            o.write(line)
PY
    "$DIAMOND" blastx -q "$D/reads.fa" -d "$D/db" -p "$T" --id 50 --evalue 1e-5 --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length evalue bitscore --quiet -o "$D/hits.tsv"
    # disjointness check: held-out EQ proteins against the retained marker proteins
    grep -A1 --no-group-separator -E '^>(EQ|REF)\|' "$D/db.faa" > "$D/db_markers.faa" || true
    "$DIAMOND" makedb --in "$D/db_markers.faa" -d "$D/dbm" -p "$T" --quiet
    "$DIAMOND" blastp -q "$D/heldout_eq.faa" -d "$D/dbm" -p "$T" --id 20 --evalue 10 --more-sensitive \
        --max-target-seqs 0 --outfmt 6 qseqid sseqid pident length qlen slen --quiet -o "$D/verify.tsv"
    rm -f "$D/db.faa" "$D/db.dmnd" "$D/db_markers.faa" "$D/dbm.dmnd"
    echo "[$(date -Iseconds)] $SCHEME fold $F: $(grep -c '^>' "$D/reads.fa") reads, $(cut -f1 "$D/hits.tsv" | sort -u | wc -l) hit"
  done
done

echo "[$(date -Iseconds)] tabulate"
python3 - "$B" "$LAB" "$W/results" "$W/inputs" <<'PY'
import csv, sys
from collections import Counter, defaultdict
B, LAB, OUT, INP = sys.argv[1:5]
lab = {r["subject"]: r["label"] for r in csv.DictReader(open(LAB), delimiter="\t")}
def group(l):
    for k, g in (("Group_1l", "hydB_1l"), ("Group_1h", "hydB_1h"), ("Group_2a", "hydB_2a")):
        if l.endswith(k): return g
    if l == "form_I": return "coxL_formI"
    if l in ("form_II", "coxL_other", "coxL_unresolved"): return "coxL_other"
    if l == "none": return "none"
    if "NiFe" in l: return "other_nife"
    if "Fe" in l or "HydDB" in l: return "other_hyd"
    return "SCG"
tax = {}
for f in ("r232_bac120.summary.tsv", "r232_ar53.summary.tsv"):
    for r in csv.DictReader(open(f"{INP}/{f}"), delimiter="\t"):
        t = dict(x.split("__", 1) for x in r["classification"].split(";") if "__" in x)
        tax[r["user_genome"]] = (t.get("p", "") or "unclassified", t.get("o", "") or "unclassified")
truth = {r["read"]: (r["source"], group(r["label"])) for r in csv.DictReader(open(f"{B}/fragments_truth.tsv"), delimiter="\t")}
GROUPS = ("hydB_1l", "hydB_1h", "hydB_2a", "other_nife", "coxL_formI", "coxL_other")
def nnbin(x):
    return "none" if x is None else ("<50" if x < 50 else "50-60" if x < 60 else "60-70" if x < 70 else "70-80" if x < 80 else ">=80")
conf, nn_rows, lin, ver = [], Counter(), Counter(), []
for scheme in ("comp50", "comp70"):
    for floor in (50, 60, 70):
        c = Counter(); n = Counter()
        for f in range(5):
            D = f"{B}/{scheme}_f{f}"
            held = {l.strip() for l in open(f"{D}/heldout.txt")}
            reps = {h.split("|")[2] for h in held if h.startswith("EQ|")}
            best = {}
            for line in open(f"{D}/hits.tsv"):
                q, s, pid, aln = line.split("\t")[:4]
                if q not in best and int(aln) >= 25:
                    best[q] = (s, float(pid))
            gl, loc = {}, {}
            for line in open(f"{D}/verify.tsv"):
                q, s, pid, aln, ql, sl = line.split("\t")
                pid, aln = float(pid), int(aln)
                if aln >= 0.5 * min(int(ql), int(sl)): gl[q] = max(gl.get(q, 0), pid)
                if aln >= 33: loc[q] = max(loc.get(q, 0), pid)
            if floor == 50:
                for q in {h for h in held if h.startswith("EQ|")}:
                    ver.append((scheme, f, q, group(lab[q]), gl.get(q), loc.get(q)))
            for rid, (src, t) in truth.items():
                if not (src in held or (src.startswith("NEG|") and src.split("|")[1] in reps)):
                    continue
                n[t] += 1
                h = best.get(rid)
                a = group(lab[h[0]]) if h and h[1] >= floor else "unassigned"
                c[(t, a)] += 1
                if t in GROUPS:
                    nn_rows[(scheme, floor, t, nnbin(gl.get(src)), "n")] += 1
                    nn_rows[(scheme, floor, t, nnbin(gl.get(src)), "tp")] += (a == t)
                    if t in ("hydB_1l", "coxL_formI"):
                        p, o = tax.get(src.split("|")[2], ("NA", "NA"))
                        lin[(scheme, floor, t, p, o, "n")] += 1
                        lin[(scheme, floor, t, p, o, "tp")] += (a == t)
        for (t, a), k in sorted(c.items()):
            conf.append((scheme, floor, t, a, k, n[t]))
with open(f"{OUT}/disjoint_benchmark_confusion.tsv", "w") as o:
    o.write("scheme\tid_floor\ttrue\tassigned\treads\treads_true\tfraction\n")
    for r in conf: o.write("\t".join(map(str, r)) + f"\t{r[4] / r[5]:.4f}\n")
with open(f"{OUT}/disjoint_benchmark_summary.tsv", "w") as o:
    o.write("scheme\tid_floor\tgroup\trecall\tprecision\tassigned_reads\ttrue_reads\tnegatives_assigned\n")
    for scheme in ("comp50", "comp70"):
        for floor in (50, 60, 70):
            sub = [r for r in conf if r[0] == scheme and r[1] == floor]
            for g in GROUPS:
                tp = sum(r[4] for r in sub if r[2] == g and r[3] == g)
                pos = sum(r[4] for r in sub if r[2] == g)
                ass = sum(r[4] for r in sub if r[3] == g)
                neg = sum(r[4] for r in sub if r[2] == "none" and r[3] == g)
                o.write(f"{scheme}\t{floor}\t{g}\t{tp / pos if pos else float('nan'):.4f}\t"
                        f"{tp / ass if ass else float('nan'):.4f}\t{ass}\t{pos}\t{neg}\n")
with open(f"{OUT}/disjoint_recall_by_nn_identity.tsv", "w") as o:
    o.write("scheme\tid_floor\tgroup\tnearest_retained_identity\treads\trecall\n")
    for k in sorted({k[:4] for k in nn_rows}):
        nn = nn_rows[k + ("n",)]
        o.write("\t".join(map(str, k)) + f"\t{nn}\t{nn_rows[k + ('tp',)] / nn:.4f}\n")
with open(f"{OUT}/disjoint_recall_by_lineage.tsv", "w") as o:
    o.write("scheme\tid_floor\tgroup\tphylum\torder\treads\trecall\n")
    for k in sorted({k[:5] for k in lin}):
        nn = lin[k + ("n",)]
        o.write("\t".join(map(str, k)) + f"\t{nn}\t{lin[k + ('tp',)] / nn:.4f}\n")
with open(f"{OUT}/disjoint_verification.tsv", "w") as o:
    o.write("scheme\tfold\tsubject\tgroup\tmax_identity_cov50\tmax_identity_local33aa\n")
    for r in ver: o.write("\t".join("NA" if x is None else str(x) for x in r) + "\n")
for scheme in ("comp50", "comp70"):
    v = [r for r in ver if r[0] == scheme]
    g = [r[4] for r in v if r[4] is not None]; l = [r[5] for r in v if r[5] is not None]
    print(f"[verify] {scheme}: held-out EQ {len(v)}; max cov50 identity {max(g) if g else None}; "
          f"n >= {scheme[4:]}% (cov50) {sum(x >= int(scheme[4:]) for x in g)}; max local33 {max(l) if l else None}; "
          f"n local33 >= 50% {sum(x >= 50 for x in l)}", flush=True)
print(open(f"{OUT}/disjoint_benchmark_summary.tsv").read())
PY
sha256sum "$W/results/disjoint_"*.tsv "$B/comp50.tsv" "$B/comp70.tsv" > "$W/results/disjoint_benchmark.sha256"
echo "[$(date -Iseconds)] done"
