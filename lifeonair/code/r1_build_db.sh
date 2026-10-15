#!/bin/bash
# Life on air (P3), read-based guild abundance -- stage R1: build the competitive protein database.
#
# Database = labelled proteins; every read is assigned to its single best hit, so related but
# non-atmospheric enzymes compete with the high-affinity ones:
#   hydB   EQ K06281 proteins of the 5,438 reps, labelled by HydDB class (lifeonair_clean typing)
#          + all HydDB reference hydrogenases, labelled by group
#   coxL   EQ K03520 proteins, labelled form_I / coxL_other (Cordero identity + AYRGAGR concordance)
#          + the three Cordero form-I references
#   SCG    14 universal single-copy ribosomal proteins from the same 5,438 reps (eggNOG KO),
#          clustered at 90% identity, used to normalise to genome equivalents
# Outputs under 06_atmospheric/lifeonair_reads/db/.
#
#   sbatch r1_build_db.sh
#
#SBATCH --job-name=loa_reads_db
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/db_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/db_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
L=$R/06_atmospheric/lifeonair_clean
REFS=$R/06_atmospheric/refs
A=$R/06_atmospheric/lifeonair_reads
DB=$A/db
mkdir -p "$DB" "$A/logs"
[ -s "$DB/DB_DONE" ] && { echo "db already built ($DB/DB_DONE); delete it to rebuild"; exit 0; }

CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
DIAMOND=$CONDA/eggnog/bin/diamond
MMSEQS=$CONDA/mmseqs2/bin/mmseqs
T=${SLURM_CPUS_PER_TASK:-16}
source "$CONDA/../bin/activate" bakta   # any env with python3

HYDDB=$(find "$REFS/HydDB" -name 'HydDB_all_hydrogenases.faa' | head -1)
[ -s "$HYDDB" ] || { echo "ERROR: HydDB fasta not found under $REFS/HydDB"; exit 1; }

echo "[$(date -Iseconds)] extract SCG ribosomal proteins from the 5,438 reps"
python3 - "$R" "$DB" <<'PY'
import os, csv, glob, sys
R, OUTD = sys.argv[1], sys.argv[2]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
REPS = f"{R}/05_refine/_drep/out/dereplicated_genomes"
SCG = {"K02886": "rplB", "K02906": "rplC", "K02926": "rplD", "K02931": "rplE", "K02933": "rplF",
       "K02867": "rplK", "K02874": "rplN", "K02878": "rplP", "K02982": "rpsC", "K02988": "rpsE",
       "K02994": "rpsH", "K02946": "rpsJ", "K02948": "rpsK", "K02965": "rpsS"}
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}

def resolve(rep, kind):                       # identical resolver to lifeonair_2_type_guild.sh
    sample = rep.split("__")[0]
    m = meta.get(sample)
    if "__" not in rep or not m:
        return None
    rest = rep.split("__", 1)[1]
    tok, local = rest.split("_", 1)
    anno = "single_sample" if tok == "single" else "self_supervised"
    base = f"{TREE}/{m['trip']}/site{m['site']}/{sample}/wgs/annotation/{anno}"
    cands = [local] + ([local[:-4]] if local.endswith("_sub") else [])
    for lc in cands:
        p = f"{base}/eggnog/{lc}/{lc}.emapper.annotations" if kind == "eggnog" \
            else f"{base}/bakta/{lc}/{lc}.faa"
        if os.path.exists(p):
            return p
    return None

def scg_proteins(ann_path):
    hits, hdr = {}, None
    with open(ann_path) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.rstrip("\n").lstrip("#").split("\t"); continue
            if line.startswith("#") or not hdr:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            ko = row.get("KEGG_ko", "-")
            if not ko or ko == "-":
                continue
            ks = [k.replace("ko:", "") for k in ko.split(",")]
            g = [SCG[k] for k in ks if k in SCG]
            if len(g) == 1:
                hits[row["query"]] = g[0]
    return hits

def read_faa(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None: yield name, "".join(seq)
                name = line[1:].split()[0]; seq = []
            else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

reps = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{REPS}/*.fa"))
out = open(f"{OUTD}/scg_raw.faa", "w")
n_rep = n_prot = miss = 0
for rep in reps:
    ap, fp = resolve(rep, "eggnog"), resolve(rep, "bakta")
    if not ap or not fp:
        miss += 1; continue
    hits = scg_proteins(ap)
    if not hits:
        continue
    n_rep += 1
    for pid, s in read_faa(fp):
        g = hits.get(pid)
        if g and len(s) >= 50:
            out.write(f">SCG|{g}|{rep}|{pid}\n{s}\n"); n_prot += 1
out.close()
print(f"[scg] reps={len(reps)} unresolved={miss} reps_with_scg={n_rep} proteins={n_prot}", flush=True)
PY

echo "[$(date -Iseconds)] cluster SCG at 90% identity"
"$MMSEQS" easy-cluster "$DB/scg_raw.faa" "$DB/scg_c90" "$DB/tmp_scg" --min-seq-id 0.9 -c 0.8 \
    --threads "$T" >/dev/null 2>&1
rm -rf "$DB/tmp_scg"

echo "[$(date -Iseconds)] label marker proteins"
python3 - "$L" "$HYDDB" "$REFS/coxL_formI_refs.faa" "$DB" <<'PY'
import sys, csv, re
L, HYDDB, CORD, OUTD = sys.argv[1:5]
def read_faa(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None: yield name, "".join(seq)
                name = line[1:].rstrip("\n"); seq = []
            else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)
safe = lambda s: re.sub(r"[^A-Za-z0-9_.-]", "_", s)

hyd_class = {}
for r in csv.DictReader(open(f"{L}/classify/atmo_classification.tsv"), delimiter="\t"):
    hyd_class[(r["rep_id"], r["protein_id"], r["marker"])] = r["assignment"]
cox_call = {}
for r in csv.DictReader(open(f"{L}/classify/coxL_class.tsv"), delimiter="\t"):
    cox_call[r["query"]] = r["form_I_call"]

out = open(f"{OUTD}/markers.faa", "w"); lab = open(f"{OUTD}/marker_labels.tsv", "w")
lab.write("subject\tmarker\tlabel\tsource\n")
n = {}
def emit(sid, seq, marker, label, src):
    out.write(f">{sid}\n{seq}\n"); lab.write(f"{sid}\t{marker}\t{label}\t{src}\n")
    n[(marker, label, src)] = n.get((marker, label, src), 0) + 1

for name, s in read_faa(f"{L}/queries/hydB.faa"):
    rep, pid = name.split()[0].split("|")
    a = hyd_class.get((rep, pid, "hydB"))
    if a is None:
        continue
    emit(f"EQ|hydB|{rep}|{pid}", s, "hydB", safe(a), "EQ")
for name, s in read_faa(HYDDB):
    # HydDB headers carry the group label; take the [NiFe]/[FeFe]/[Fe] group token
    m = re.search(r"(\[(?:NiFe|FeFe|Fe)\][ _]Group[ _][0-9A-Za-z]+)", name)
    label = safe(m.group(1).replace(" ", "_")) if m else "HydDB_unlabelled"
    emit(f"REF|hydB|{safe(name.split()[0])}", s, "hydB", label, "HydDB")
for name, s in read_faa(f"{L}/queries/coxL.faa"):
    q = name.split()[0]
    c = cox_call.get(q)
    if c is None:
        continue
    rep, pid = q.split("|")
    emit(f"EQ|coxL|{rep}|{pid}", s, "coxL", "form_I" if c == "form_I" else "coxL_other", "EQ")
for name, s in read_faa(CORD):
    emit(f"REF|coxL|{safe(name.split()[0])}", s, "coxL", "form_I", "Cordero2019")
out.close(); lab.close()
for k in sorted(n): print("[label]", *k, n[k], flush=True)
PY

python3 - "$DB" <<'PY'
import sys
D = sys.argv[1]
with open(f"{D}/scg_c90_rep_seq.fasta") as fh, open(f"{D}/marker_labels.tsv", "a") as lab:
    for line in fh:
        if line.startswith(">"):
            sid = line[1:].split()[0]
            lab.write(f"{sid}\tSCG\t{sid.split('|')[1]}\tEQ\n")
PY
cat "$DB/markers.faa" "$DB/scg_c90_rep_seq.fasta" > "$DB/loa_reads.faa"
"$DIAMOND" makedb --in "$DB/loa_reads.faa" -d "$DB/loa_reads" -p "$T"
python3 - "$DB/loa_reads.faa" > "$DB/subject_lengths.tsv" <<'PY'
import sys
n, L = None, 0
print("subject\tlength_aa")
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if n: print(f"{n}\t{L}")
        n, L = line[1:].split()[0], 0
    else: L += len(line.strip())
if n: print(f"{n}\t{L}")
PY
cut -f2,3 "$DB/marker_labels.tsv" | sort | uniq -c | sort -k2,2 -k1,1nr
sha256sum "$DB/loa_reads.faa" "$DB/marker_labels.tsv" "$DB/subject_lengths.tsv" > "$DB/db.sha256"
date -Iseconds > "$DB/DB_DONE"
echo "[$(date -Iseconds)] done"
