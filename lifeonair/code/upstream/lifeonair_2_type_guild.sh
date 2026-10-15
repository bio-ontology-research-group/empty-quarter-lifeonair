#!/bin/bash
# lifeonair step 2 -- TYPING + GUILD. Types hydB/coxL of the current dRep catalogue and builds the
# high-affinity trace-gas guild. Run after step 1 (census).
# Fuses scripts/ibex/82b_atmo_extract.sh (extract) + 82c_atmo_classify.sh (classify),
# repointed from the OLD 15,469 tables to the 5,438 basis:
#   reps        : 05_refine/_drep/out/dereplicated_genomes/*.fa                 (the 5,438)
#   annotations : reused per-MAG eggNOG   (SAME resolver as scripts/eco_probes/gather_rep_ko.py)
#   proteome    : companion Bakta .faa    (sibling of the eggnog dir: .../bakta/<local>/<local>.faa)
#   refs        : 06_atmospheric/refs/    (HydDB + 3 Cordero form-I refs — ALREADY BUILT, reused as-is)
# Nothing is re-annotated; only the phylogenetic typing (DIAMOND vs refs) is new compute.
# Outputs under 06_atmospheric/lifeonair_clean/{queries,classify}/ .
#
#   sbatch scripts/type_atmo_5438.sh
#
#SBATCH --job-name=atmo_type_5438
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_clean/logs/type_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_clean/logs/type_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
A=$R/06_atmospheric/lifeonair_clean
REFS=$R/06_atmospheric/refs
Q=$A/queries
OUT=$A/classify
mkdir -p "$Q" "$OUT" "$A/logs"
[ -s "$OUT/CLASSIFY_DONE" ] && { echo "typing already complete ($OUT/CLASSIFY_DONE) -- delete it to rerun"; exit 0; }

CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
DIAMOND=$CONDA/eggnog/bin/diamond
MMSEQS=$CONDA/mmseqs2/bin/mmseqs
MAFFT=$CONDA/busco/bin/mafft
FASTTREE=$CONDA/busco/bin/FastTree
T=${SLURM_CPUS_PER_TASK:-16}

[ -s "$REFS/REFS_DONE" ] || { echo "ERROR: refs not built (06_atmospheric/refs; run 82a)"; exit 1; }
[ -x "$DIAMOND" ] || { echo "ERROR: missing $DIAMOND"; exit 1; }

# ================================================================= STAGE B: extract (5,438 basis)
echo "[$(date -Iseconds)] Stage B: extract hydB(K06281)/coxL(K03520) proteins from the 5,438 reps"
source "$CONDA/../bin/activate" bakta   # any env with python3
python3 - "$R" "$Q" <<'PY'
import os, csv, glob, sys
R, OUTD = sys.argv[1], sys.argv[2]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
REPS = f"{R}/05_refine/_drep/out/dereplicated_genomes"
KO   = {"K06281": "hydB", "K03520": "coxL"}
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}

def resolve(rep, kind):                       # kind = "eggnog" (.emapper.annotations) | "bakta" (.faa)
    sample = rep.split("__")[0]
    m = meta.get(sample)
    if "__" not in rep or not m:
        return None
    rest = rep.split("__", 1)[1]              # single_SemiBin_0 | self_SemiBin_162[_sub]
    tok, local = rest.split("_", 1)           # tok=single/self ; local=SemiBin_0[...]
    anno = "single_sample" if tok == "single" else "self_supervised"
    base = f"{TREE}/{m['trip']}/site{m['site']}/{sample}/wgs/annotation/{anno}"
    cands = [local] + ([local[:-4]] if local.endswith("_sub") else [])   # strip trailing _sub fallback
    for lc in cands:
        p = f"{base}/eggnog/{lc}/{lc}.emapper.annotations" if kind == "eggnog" \
            else f"{base}/bakta/{lc}/{lc}.faa"
        if os.path.exists(p):
            return p
    return None

def marker_proteins(ann_path):                # protein_id -> marker for rows carrying K06281/K03520
    hits, hdr = {}, None
    with open(ann_path) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.rstrip("\n").lstrip("#").split("\t"); continue
            if line.startswith("#") or not hdr:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            ko = row.get("KEGG_ko", "-")
            if not ko or ko in ("-", ""):
                continue
            for k in ko.split(","):
                kk = k.replace("ko:", "")
                if kk in KO:
                    hits[row["query"]] = KO[kk]
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
print(f"[extract] {len(reps)} reps", flush=True)
outs = {"hydB": open(f"{OUTD}/hydB.faa", "w"), "coxL": open(f"{OUTD}/coxL.faa", "w")}
lens = open(f"{OUTD}/query_lengths.tsv", "w"); lens.write("rep\tprotein_id\tmarker\tlength\n")
mh   = open(f"{OUTD}/marker_hits.tsv", "w");   mh.write("rep\tprotein_id\tKO\tmarker\n")
counts = {"hydB": 0, "coxL": 0}; reps_with = {"hydB": set(), "coxL": set()}
no_ann = no_faa = 0
for rep in reps:
    ap = resolve(rep, "eggnog")
    if not ap:
        no_ann += 1; continue
    hits = marker_proteins(ap)
    if not hits:
        continue
    fp = resolve(rep, "bakta")
    if not fp:
        no_faa += 1; continue
    for pid, s in read_faa(fp):
        m = hits.get(pid)
        if m:
            outs[m].write(f">{rep}|{pid}\n{s}\n")
            lens.write(f"{rep}\t{pid}\t{m}\t{len(s)}\n")
            mh.write(f"{rep}\t{pid}\t{'K06281' if m=='hydB' else 'K03520'}\t{m}\n")
            counts[m] += 1; reps_with[m].add(rep)
for o in outs.values(): o.close()
lens.close(); mh.close()
print(f"[extract] no_ann={no_ann} no_faa={no_faa}", flush=True)
print(f"[extract] hydB: {counts['hydB']} proteins across {len(reps_with['hydB'])} reps  (census expected ~1102)", flush=True)
print(f"[extract] coxL: {counts['coxL']} proteins across {len(reps_with['coxL'])} reps  (census expected ~1818)", flush=True)
PY
date -Iseconds > "$Q/EXTRACT_DONE"
grep -c '^>' "$Q"/hydB.faa "$Q"/coxL.faa || true

# ================================================================= STAGE C: classify (reuse 82c logic)
for b in "$DIAMOND" "$MMSEQS" "$MAFFT" "$FASTTREE"; do [ -x "$b" ] || echo "WARN missing $b (tree step may skip)"; done

echo "[$(date -Iseconds)] hydB: DIAMOND vs HydDB"
"$DIAMOND" blastp -q "$Q/hydB.faa" -d "$REFS/hyddb" -p "$T" --ultra-sensitive --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore -o "$OUT/hydB_diamond.tsv" --quiet
awk -F'\t' 'BEGIN{OFS="\t"}{split($2,a,/\|/); grp=a[3]; if(grp=="")grp=$2; print $1, grp, $3, $4}' \
    "$OUT/hydB_diamond.tsv" > "$OUT/hydB_class.tsv"           # qseqid grp pident qcov
echo "  hydB classified: $(wc -l < "$OUT/hydB_class.tsv") proteins"
echo "  --- HydDB class distribution (top) ---"; cut -f2 "$OUT/hydB_class.tsv" | sort | uniq -c | sort -rn | head -25 || true

echo "[$(date -Iseconds)] coxL: DIAMOND vs form-I refs"
"$DIAMOND" blastp -q "$Q/coxL.faa" -d "$REFS/coxL_formI_refs" -p "$T" --ultra-sensitive --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length evalue bitscore -o "$OUT/coxL_diamond.tsv" --quiet
echo "  coxL with a form-I hit: $(cut -f1 "$OUT/coxL_diamond.tsv" | sort -u | wc -l)"

python3 - "$Q" "$OUT" <<'PY'
import sys, os, re
Q, OUT = sys.argv[1], sys.argv[2]
def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None: yield name, "".join(seq)
                name = line[1:].rstrip("\n"); seq = []
            else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)
cox = {}
dp = os.path.join(OUT, "coxL_diamond.tsv")
if os.path.exists(dp):
    with open(dp) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t"); q = p[0]
            if q in cox: continue
            cox[q] = {"pident": float(p[2]), "qcov": float(p[3])}
FORM_I  = re.compile(r"AYRGAGR")
FORM_II = re.compile(r"AYR[CG]SFR")
def motif_call(seq):
    if FORM_I.search(seq):  return "AYRGAGR(formI)"
    if FORM_II.search(seq): return "AYRxSFR(formII)"
    return "none"
with open(os.path.join(OUT, "coxL_class.tsv"), "w") as o:
    o.write("query\tlength\tpident_formI\tqcov_formI\tmotif\tform_I_call\n")
    n_formI = n_motifI = total = 0
    for name, seq in read_fasta(os.path.join(Q, "coxL.faa")):
        total += 1
        h = cox.get(name, {}); pid = h.get("pident", 0.0); qc = h.get("qcov", 0.0)
        mot = motif_call(seq); thr = (pid > 55.0 and qc >= 90.0)
        call = "form_I" if (thr or mot == "AYRGAGR(formI)") else ("form_II/other" if mot == "AYRxSFR(formII)" else "ambiguous")
        if call == "form_I": n_formI += 1
        if mot == "AYRGAGR(formI)": n_motifI += 1
        o.write(f"{name}\t{len(seq)}\t{pid:.1f}\t{qc:.1f}\t{mot}\t{call}\n")
    print(f"[coxL] {total} queries; form_I={n_formI}; AYRGAGR motif={n_motifI}", flush=True)
PY
echo "  --- coxL form call distribution ---"; tail -n +2 "$OUT/coxL_class.tsv" | cut -f6 | sort | uniq -c | sort -rn

echo "[$(date -Iseconds)] writing atmo_classification.tsv (per-protein)"
python3 - "$OUT" <<'PY'
import os, sys
OUT = sys.argv[1]; rows = []
with open(os.path.join(OUT, "hydB_class.tsv")) as fh:
    for line in fh:
        q, grp, pid, qc = line.rstrip("\n").split("\t")
        rep, prot = q.split("|", 1)
        ha = 1 if any(t in grp.lower() for t in ("1h", "1l")) and "nife" in grp.lower() else (
             1 if any(t in grp for t in ("Group 1h", "Group 1l")) else 0)
        rows.append((rep, prot, "hydB", grp, pid, qc, ha))
with open(os.path.join(OUT, "coxL_class.tsv")) as fh:
    next(fh)
    for line in fh:
        q, L, pid, qc, mot, call = line.rstrip("\n").split("\t")
        rep, prot = q.split("|", 1)
        ha = 1 if call == "form_I" else 0
        rows.append((rep, prot, "coxL", f"{call}|{mot}", pid, qc, ha))
with open(os.path.join(OUT, "atmo_classification.tsv"), "w") as o:
    o.write("rep_id\tprotein_id\tmarker\tassignment\tpident\tqcov\tis_high_affinity\n")
    for r in rows: o.write("\t".join(str(x) for x in r) + "\n")
print(f"[combined] {len(rows)} classified proteins -> atmo_classification.tsv", flush=True)
PY

# ---- guild = reps with >=1 high-affinity hydB (1h/1l) and/or >=1 form-I coxL (per-species, non-redundant)
echo "[$(date -Iseconds)] building high-affinity guild"
awk -F'\t' 'NR>1 && $7==1 && $3=="hydB"{print $1}' "$OUT/atmo_classification.tsv" | sort -u > "$OUT/guild_hydB_reps.txt"
awk -F'\t' 'NR>1 && $7==1 && $3=="coxL"{print $1}' "$OUT/atmo_classification.tsv" | sort -u > "$OUT/guild_coxL_reps.txt"
sort -u "$OUT/guild_hydB_reps.txt" "$OUT/guild_coxL_reps.txt" > "$OUT/guild_reps.txt"
comm -12 "$OUT/guild_hydB_reps.txt" "$OUT/guild_coxL_reps.txt" > "$OUT/guild_both_reps.txt"
echo "  high-affinity hydB reps : $(wc -l < "$OUT/guild_hydB_reps.txt")"
echo "  form-I coxL reps        : $(wc -l < "$OUT/guild_coxL_reps.txt")"
echo "  guild (union) reps      : $(wc -l < "$OUT/guild_reps.txt")"
echo "  both H2+CO reps         : $(wc -l < "$OUT/guild_both_reps.txt")"

# ---- dereplicate + tree per marker (figure inputs; identical to 82c, optional)
for m in coxL hydB; do
    [ -x "$MMSEQS" ] || { echo "  (skip tree: no mmseqs)"; break; }
    echo "[$(date -Iseconds)] $m: dereplicate + tree"
    "$MMSEQS" easy-cluster "$Q/$m.faa" "$OUT/${m}_clu" "$OUT/tmp_$m" --min-seq-id 0.7 -c 0.8 --threads "$T" >/dev/null 2>&1 \
        || { echo "  (mmseqs cluster skipped for $m)"; continue; }
    reps="$OUT/${m}_clu_rep_seq.fasta"; nrep=$(grep -c '^>' "$reps" 2>/dev/null || echo 0)
    echo "  $m reps: $nrep"
    if [ "$nrep" -gt 2 ] && [ "$nrep" -lt 6000 ] && [ -x "$MAFFT" ] && [ -x "$FASTTREE" ]; then
        # marker-matched references (82c added coxL refs to BOTH trees -- fixed here)
        case "$m" in
          coxL) MREF="$REFS/coxL_formI_refs.faa" ;;
          hydB) MREF=$(find "$REFS/HydDB" -name 'HydDB_all_hydrogenases.faa' | head -1) ;;
        esac
        if [ -s "$MREF" ]; then cat "$reps" "$MREF" > "$OUT/${m}_reps_plus_refs.faa"; else cp "$reps" "$OUT/${m}_reps_plus_refs.faa"; fi
        "$MAFFT" --auto --thread "$T" "$OUT/${m}_reps_plus_refs.faa" > "$OUT/${m}.aln" 2>/dev/null
        "$FASTTREE" "$OUT/${m}.aln" > "$OUT/${m}.nwk" 2>/dev/null && echo "  tree -> $OUT/${m}.nwk"
    else
        echo "  (skip tree for $m: nrep=$nrep out of range or aligner missing)"
    fi
    rm -rf "$OUT/tmp_$m"
done

echo "[$(date -Iseconds)] Stage C done"
date -Iseconds > "$OUT/CLASSIFY_DONE"
