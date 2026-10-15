#!/bin/bash
# Life on air (P3), mechanism proxy for H1: predicted maximum growth rate (copiotrophy) per representative.
#   Li et al. 2023 (ISME J) tie the loss of trace-gas oxidisers under shrubs to organic carbon favouring
#   copiotrophs. Soil organic carbon was never measured here, so we use the same life-history proxy:
#   gRodon2 codon-usage-bias minimum doubling time (Weissman et al. 2021), mode "partial" (MAGs),
#   highly expressed set = Bakta CDS annotated as ribosomal proteins under either naming scheme
#   ("50S ribosomal protein L2" or "Large ribosomal subunit protein uL2"); the first run (53541016) matched only the
#   old names and left 2,824 reps without an estimate.
#   Input: the 5,438 ANI95 representatives, Bakta .ffn/.faa resolved as in c1_operon_context.sh.
#   Output: growth/rep_growth.tsv (rep, n_cds, n_ribo, d_hours, lower, upper, cub_hegs).
#
#   sbatch g1_growth_rate.sh
#
#SBATCH --job-name=loa_growth
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/growth_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/growth_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$R/06_atmospheric/lifeonair_reads/growth
ENV=$R/06_atmospheric/lifeonair_reads/envs/grodon
mkdir -p "$OUT/ffn" "$(dirname "$ENV")"
CONDA=/ibex/scratch/projects/c2014/rob/miniconda3
source "$CONDA/bin/activate" base
T=${SLURM_CPUS_PER_TASK:-16}

if [[ ! -x "$ENV/bin/Rscript" ]]; then
  conda create -y -p "$ENV" -c conda-forge -c bioconda r-base=4.3 bioconductor-biostrings \
      bioconductor-cordon r-matrixstats r-remotes r-data.table >/dev/null
  "$ENV/bin/Rscript" -e 'remotes::install_github("jlw-ecoevo/gRodon2", upgrade = "never", quiet = TRUE)'
fi
"$ENV/bin/Rscript" -e 'cat("gRodon", as.character(packageVersion("gRodon")), packageDescription("gRodon")$RemoteSha, "\n")' | tee "$OUT/grodon_version.txt"

# 1. per rep: CDS-only ffn (ids present in .faa), ribosomal-protein flag in the header
conda activate bakta
python3 - "$R" "$OUT" <<'PY'
import csv, glob, os, re, sys
R, OUT = sys.argv[1:3]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
REPS = f"{R}/05_refine/_drep/out/dereplicated_genomes"
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}
def resolve(rep):                             # identical resolver to c1_operon_context.sh
    sample = rep.split("__")[0]; m = meta.get(sample)
    if "__" not in rep or not m: return None
    tok, local = rep.split("__", 1)[1].split("_", 1)
    anno = "single_sample" if tok == "single" else "self_supervised"
    base = f"{TREE}/{m['trip']}/site{m['site']}/{sample}/wgs/annotation/{anno}"
    for lc in [local] + ([local[:-4]] if local.endswith("_sub") else []):
        p = f"{base}/bakta/{lc}/{lc}"
        if os.path.exists(p + ".ffn") and os.path.exists(p + ".faa"): return p
    return None
RIBO = re.compile(r"^((30S|50S) ribosomal protein|(Large|Small) ribosomal subunit protein) ")   # old and new UniProt names
reps = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{REPS}/*.fa"))
miss = 0
for rep in reps:
    p = resolve(rep)
    if not p: miss += 1; continue
    cds = {l[1:].split()[0] for l in open(p + ".faa") if l.startswith(">")}
    with open(f"{OUT}/ffn/{rep}.ffn", "w") as o:
        keep = False
        for line in open(p + ".ffn"):
            if line.startswith(">"):
                tag, _, prod = line[1:].strip().partition(" ")
                keep = tag in cds
                if keep: o.write(f">{tag} {'HEG' if RIBO.match(prod) else 'x'}\n")
            elif keep: o.write(line)
print(f"[ffn] reps={len(reps)} unresolved={miss}", flush=True)
PY

# 2. gRodon per rep, parallel
cat > "$OUT/run_grodon.R" <<'RS'
suppressMessages({library(gRodon); library(Biostrings); library(parallel)})
args <- commandArgs(TRUE); files <- list.files(args[1], "\\.ffn$", full.names = TRUE)
one <- function(f) {
  s <- readDNAStringSet(f); s <- s[width(s) %% 3 == 0 & width(s) >= 240]
  heg <- grepl(" HEG$", names(s))
  rep <- sub("\\.ffn$", "", basename(f))
  if (sum(heg) < 10) return(data.frame(rep, n_cds = length(s), n_ribo = sum(heg), d_hours = NA, lower = NA, upper = NA, cub_hegs = NA))
  g <- tryCatch(predictGrowth(s, heg, mode = "partial"), error = function(e) NULL)
  if (is.null(g)) return(data.frame(rep, n_cds = length(s), n_ribo = sum(heg), d_hours = NA, lower = NA, upper = NA, cub_hegs = NA))
  data.frame(rep, n_cds = length(s), n_ribo = sum(heg), d_hours = g$d, lower = g$LowerCI, upper = g$UpperCI, cub_hegs = g$CUBHE)
}
res <- do.call(rbind, mclapply(files, one, mc.cores = as.integer(args[3])))
write.table(res, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat("[grodon] reps", nrow(res), "with estimate", sum(!is.na(res$d_hours)), "\n")
RS
"$ENV/bin/Rscript" "$OUT/run_grodon.R" "$OUT/ffn" "$OUT/rep_growth.tsv" "$T"
rm -rf "$OUT/ffn"
sha256sum "$OUT/rep_growth.tsv" > "$OUT/rep_growth.sha256"
date -Iseconds > "$OUT/GROWTH_DONE"
