#!/bin/bash
# Life on air (P3), analysis 6: figure trees for the typed large subunits.
#   hydB: EQ hydB proteins of the 5,438 reps clustered at 70% identity (lifeonair_clean) + every
#         HydDB [NiFe] Group 1 and Group 2 reference; mafft --auto, trimAl-free gap filter (>50%
#         gaps removed), FastTree -lg -gamma.
#   coxL: EQ coxL proteins >= 600 aa clustered at 70% identity + the three Cordero form-I references.
# Tip labels keep the source (EQ|rep|protein or REF|accession|organism|group).
#
#   sbatch t1_trees.sh
#
#SBATCH --job-name=loa_trees
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/trees_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/trees_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
L=$R/06_atmospheric/lifeonair_clean
REFS=$R/06_atmospheric/refs
OUT=$R/06_atmospheric/lifeonair_reads/trees
mkdir -p "$OUT"
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3/envs
MMSEQS=$CONDA/mmseqs2/bin/mmseqs
MAFFT=$CONDA/busco/bin/mafft
FASTTREE=$CONDA/busco/bin/FastTree
T=${SLURM_CPUS_PER_TASK:-16}
source "$CONDA/../bin/activate" bakta
HYDDB=$(find "$REFS/HydDB" -name 'HydDB_all_hydrogenases.faa' | head -1)

python3 - "$L" "$HYDDB" "$REFS/coxL_formI_refs.faa" "$OUT" <<'PY'
import sys
L, HYDDB, CORD, OUT = sys.argv[1:5]
def read_faa(p):   # sanitised: selenocysteine U -> X, stop/gap characters removed
    n, s = None, []
    clean = lambda x: "".join(x).upper().replace("U", "X").replace("*", "").replace("-", "")
    for line in open(p):
        if line.startswith(">"):
            if n: yield n, clean(s)
            n, s = line[1:].strip(), []
        else: s.append(line.strip())
    if n: yield n, clean(s)
with open(f"{OUT}/hydB_eq.faa", "w") as o:
    for n, s in read_faa(f"{L}/queries/hydB.faa"):
        if len(s) >= 400: o.write(f">EQ|{n.split()[0]}\n{s}\n")
with open(f"{OUT}/hydB_refs.faa", "w") as o:
    for n, s in read_faa(HYDDB):
        if "[NiFe]_Group_1" in n or "[NiFe]_Group_2" in n:
            o.write(">REF|" + n.split()[0].replace("[", "").replace("]", "") + f"\n{s}\n")
with open(f"{OUT}/coxL_eq.faa", "w") as o:
    for n, s in read_faa(f"{L}/queries/coxL.faa"):
        if len(s) >= 600: o.write(f">EQ|{n.split()[0]}\n{s}\n")
with open(f"{OUT}/coxL_refs.faa", "w") as o:
    for n, s in read_faa(CORD):
        o.write(">REF|" + n.split()[0] + f"|formI\n{s}\n")
PY

gapfilter() {  # drop alignment columns with > 50% gaps
python3 - "$1" "$2" <<'PY'
import sys
names, seqs = [], []
for line in open(sys.argv[1]):
    if line.startswith(">"): names.append(line[1:].strip()); seqs.append([])
    else: seqs[-1].append(line.strip())
seqs = ["".join(s) for s in seqs]; n = len(seqs)
keep = [i for i in range(len(seqs[0])) if sum(s[i] == "-" for s in seqs) <= 0.5 * n]
with open(sys.argv[2], "w") as o:
    for nm, s in zip(names, seqs): o.write(f">{nm}\n{''.join(s[i] for i in keep)}\n")
print(f"[gapfilter] {len(seqs[0])} -> {len(keep)} columns, {n} sequences")
PY
}

for m in hydB coxL; do
  "$MMSEQS" easy-cluster "$OUT/${m}_eq.faa" "$OUT/${m}_c70" "$OUT/tmp_$m" --min-seq-id 0.7 -c 0.8 \
      --threads "$T" >/dev/null 2>&1
  rm -rf "$OUT/tmp_$m"
  cat "$OUT/${m}_c70_rep_seq.fasta" "$OUT/${m}_refs.faa" > "$OUT/${m}_tree_input.faa"
  echo "[$m] $(grep -c '>' "$OUT/${m}_tree_input.faa") sequences"
  "$MAFFT" --auto --thread "$T" "$OUT/${m}_tree_input.faa" > "$OUT/${m}.aln" 2> "$OUT/${m}.mafft.log" \
      || { tail -20 "$OUT/${m}.mafft.log"; exit 1; }
  gapfilter "$OUT/${m}.aln" "$OUT/${m}.trim.aln"
  "$FASTTREE" -lg -gamma -quiet "$OUT/${m}.trim.aln" > "$OUT/${m}.nwk"
done
date -Iseconds > "$OUT/TREES_DONE"
