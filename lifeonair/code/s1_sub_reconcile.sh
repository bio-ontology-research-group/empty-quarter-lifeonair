#!/bin/bash
# Life on air (P3): reconcile "_sub" representatives with their retained contigs (review round 1, R1-19).
#
# DAS_Tool refinement emitted some bins as contig subsets ("..._sub"). Their annotation was resolved to the
# parent bin (the resolver strips "_sub"), so typed markers, SCG proteins and the gRodon CDS set could include
# genes on contigs that the representative no longer contains. Bakta renames contigs, so retained contigs are
# matched by sequence (MD5 of the upper-case contig sequence) between the representative FASTA and the parent
# Bakta .fna. Outputs (workbench 06_atmospheric/lifeonair_review_r1/results/):
#   sub_protein_retention.tsv   rep, protein_id, contig, retained (all CDS of all _sub reps)
#   sub_excluded_db_proteins.txt  "rep|protein" ids of database subjects (markers and SCG) on removed contigs
#   rep_growth_sub_fixed.tsv    gRodon2 re-estimated for _sub reps from retained CDS only (same settings as
#                               g1_growth_rate.sh)
#
#   sbatch s1_sub_reconcile.sh
#
#SBATCH --job-name=loa_sub
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1/logs/sub_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_review_r1/logs/sub_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
W=$R/06_atmospheric/lifeonair_review_r1
DB=$R/06_atmospheric/lifeonair_reads/db
ENV=$R/06_atmospheric/lifeonair_reads/envs/grodon
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3
T=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$W/results" "$W/sub_ffn"
source "$CONDA/bin/activate" bakta

python3 - "$R" "$W" "$DB" <<'PY'
import csv, glob, hashlib, os, re, sys
R, W, DB = sys.argv[1:4]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
REPS = f"{R}/05_refine/_drep/out/dereplicated_genomes"
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}
RIBO = re.compile(r"^((30S|50S) ribosomal protein|(Large|Small) ribosomal subunit protein) ")

def parent(rep):
    sample = rep.split("__")[0]; m = meta[sample]
    tok, local = rep.split("__", 1)[1].split("_", 1)
    anno = "single_sample" if tok == "single" else "self_supervised"
    lc = local[:-4]
    return f"{TREE}/{m['trip']}/site{m['site']}/{sample}/wgs/annotation/{anno}/bakta/{lc}/{lc}"

def fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None: yield name, "".join(seq)
            name, seq = line[1:].rstrip("\n"), []
        else: seq.append(line.strip())
    if name is not None: yield name, "".join(seq)

md5 = lambda s: hashlib.md5(s.upper().encode()).hexdigest()
subs = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{REPS}/*_sub.fa"))
db_ids = set()
for line in open(f"{DB}/loa_reads.faa"):
    if line.startswith(">EQ|"):
        _, _, rep, pid = line[1:].split()[0].split("|"); db_ids.add((rep, pid))
    elif line.startswith(">SCG|"):
        f = line[1:].split()[0].split("|"); db_ids.add((f[2], f[3]))
nrep = nprot = nret = 0; excl = []
with open(f"{W}/results/sub_protein_retention.tsv", "w") as o:
    o.write("rep\tprotein_id\tcontig\tretained\n")
    for rep in subs:
        p = parent(rep)
        if not os.path.exists(p + ".gff3"):
            print("missing parent annotation", rep, flush=True); continue
        keep = {md5(s) for _, s in fasta(f"{REPS}/{rep}.fa")}
        kept_contigs = {n.split()[0] for n, s in fasta(p + ".fna") if md5(s) in keep}
        ret = set()
        for line in open(p + ".gff3"):
            if line.startswith("##FASTA"): break
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            m = re.search(r"ID=([^;]+)", f[8])
            if not m: continue
            pid = m.group(1); r = f[0] in kept_contigs
            o.write(f"{rep}\t{pid}\t{f[0]}\t{int(r)}\n"); nprot += 1
            if r: ret.add(pid); nret += 1
            elif (rep, pid) in db_ids: excl.append(f"{rep}|{pid}")
        nrep += 1
        cds = {l[1:].split()[0] for l in open(p + ".faa") if l.startswith(">")} & ret
        with open(f"{W}/sub_ffn/{rep}.ffn", "w") as g:
            k = False
            for line in open(p + ".ffn"):
                if line.startswith(">"):
                    tag, _, prod = line[1:].strip().partition(" ")
                    k = tag in cds
                    if k: g.write(f">{tag} {'HEG' if RIBO.match(prod) else 'x'}\n")
                elif k: g.write(line)
open(f"{W}/results/sub_excluded_db_proteins.txt", "w").write("\n".join(sorted(excl)) + "\n")
print(f"[sub] reps {nrep}/{len(subs)}; CDS {nprot}; retained {nret}; database subjects on removed contigs {len(excl)}", flush=True)
PY

conda activate base
"$ENV/bin/Rscript" "$R/06_atmospheric/lifeonair_reads/growth/run_grodon.R" "$W/sub_ffn" "$W/results/rep_growth_sub_fixed.tsv" "$T"
rm -rf "$W/sub_ffn"
sha256sum "$W/results/sub_"* "$W/results/rep_growth_sub_fixed.tsv" > "$W/results/sub_reconcile.sha256"
echo "[$(date -Iseconds)] done"
