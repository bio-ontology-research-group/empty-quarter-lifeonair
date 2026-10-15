#!/bin/bash
# Life on air (P3), analysis 5: genomic context of the typed hydB and coxL genes in the 5,438 reps.
# For every typed hydB (all HydDB classes) and coxL (form I and other) protein, takes the +/-10
# gene neighbourhood on the same contig from the representative's Bakta GFF3, annotates each
# neighbour with Bakta gene/product and eggNOG KO/Preferred_name, and scores operon features:
#   hydB: small subunit (K06282 hydA/hupS/hhyS), maturation hyp/hup (K04651-K04656, K03605 hupD/hybD,
#         K03192 ureJ/hupE), plus Bakta product text matches
#   coxL: coxM (K03519), coxS (K03518), accessory coxD/E/F/G (product text)
# Writes one row per marker protein plus the neighbour table. Small; runs as a single job.
#
#   sbatch c1_operon_context.sh
#
#SBATCH --job-name=loa_operon
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/operon_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_reads/logs/operon_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=batch

set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics
L=$R/06_atmospheric/lifeonair_clean
OUT=$R/06_atmospheric/lifeonair_reads/operon
mkdir -p "$OUT"
source /ibex/scratch/projects/c2014/rob/miniconda3/bin/activate bakta

python3 - "$R" "$L" "$OUT" <<'PY'
import csv, os, re, sys
R, L, OUT = sys.argv[1:4]
TREE = "/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites"
meta = {r["sample"]: r for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}
W = 10

def resolve(rep, kind):                       # identical resolver to lifeonair_2_type_guild.sh
    sample = rep.split("__")[0]
    m = meta.get(sample)
    if "__" not in rep or not m:
        return None
    rest = rep.split("__", 1)[1]
    tok, local = rest.split("_", 1)
    anno = "single_sample" if tok == "single" else "self_supervised"
    base = f"{TREE}/{m['trip']}/site{m['site']}/{sample}/wgs/annotation/{anno}"
    for lc in [local] + ([local[:-4]] if local.endswith("_sub") else []):
        p = {"eggnog": f"{base}/eggnog/{lc}/{lc}.emapper.annotations",
             "gff": f"{base}/bakta/{lc}/{lc}.gff3"}[kind]
        if os.path.exists(p):
            return p
    return None

# typed marker proteins
targets = {}
for r in csv.DictReader(open(f"{L}/classify/atmo_classification.tsv"), delimiter="\t"):
    if r["marker"] == "hydB":
        targets.setdefault(r["rep_id"], {})[r["protein_id"]] = ("hydB", r["assignment"])
for r in csv.DictReader(open(f"{L}/classify/coxL_class.tsv"), delimiter="\t"):
    rep, pid = r["query"].split("|")
    targets.setdefault(rep, {})[pid] = ("coxL", "form_I" if r["form_I_call"] == "form_I" else "coxL_other")

HYD_SMALL = {"K06282", "K05927", "K00437"}
HYD_MAT = {"K04651", "K04652", "K04653", "K04654", "K04655", "K04656", "K03605", "K03192", "K03618", "K03619"}
HYD_MAT_RE = re.compile(r"\bhyp[A-F]\b|\bhup[A-Z]\b|hydrogenase (maturation|expression|formation)|"
                        r"nickel incorporation|hydrogenase.{0,20}(endopeptidase|protease)", re.I)
HYD_SMALL_RE = re.compile(r"hydrogenase.{0,30}small (subunit|chain)", re.I)
COX = {"K03519": "coxM", "K03518": "coxS"}
COX_ACC_RE = re.compile(r"\bcox[DEFG]\b|CoxD|CoxE|CoxF|CoxG|carbon monoxide dehydrogenase.{0,40}(accessory|D|E|F|G)\b", re.I)
COX_SUB_RE = re.compile(r"carbon[- ]monoxide dehydrogenase.{0,30}(small|medium)", re.I)

def gff_genes(path):
    contigs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            pid = a.get("locus_tag") or a.get("ID", "")
            contigs.setdefault(p[0], []).append((int(p[3]), p[6], pid, a.get("gene", ""),
                                                 a.get("product", "").replace("%2C", ",")))
    for c in contigs.values():
        c.sort()
    return contigs

def eggnog(path):
    ann, hdr = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                hdr = line.rstrip("\n").lstrip("#").split("\t"); continue
            if line.startswith("#") or not hdr:
                continue
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            kos = {k.replace("ko:", "") for k in row.get("KEGG_ko", "-").split(",") if k not in ("-", "")}
            ann[row["query"]] = (kos, row.get("Preferred_name", "-"))
    return ann

per = open(f"{OUT}/marker_context.tsv", "w")
per.write("rep\tprotein_id\tmarker\tclass\tcontig\tcontig_genes\tposition\tgenes_to_contig_end\t"
          "small_subunit_adjacent\tn_maturation\tcoxM\tcoxS\tn_cox_accessory\tn_hypothetical\tgene_order\n")
nb = open(f"{OUT}/neighbours.tsv", "w")
nb.write("rep\tprotein_id\toffset\tneighbour\tstrand\tbakta_gene\tbakta_product\tkos\teggnog_name\n")
miss = n = 0
for rep, prots in sorted(targets.items()):
    g, e = resolve(rep, "gff"), resolve(rep, "eggnog")
    if not g:
        miss += 1; continue
    contigs = gff_genes(g); ann = eggnog(e) if e else {}
    loc = {pid: (c, i) for c, genes in contigs.items() for i, x in enumerate(genes) for pid in [x[2]]}
    for pid, (marker, cls) in prots.items():
        if pid not in loc:
            continue
        c, i = loc[pid]; genes = contigs[c]
        lo, hi = max(0, i - W), min(len(genes), i + W + 1)
        small = mat = cm = cs = acc = hyp = 0; order = []
        for j in range(lo, hi):
            _, strand, q, gene, prod = genes[j]
            kos, name = ann.get(q, (set(), "-"))
            off = j - i
            nb.write(f"{rep}\t{pid}\t{off}\t{q}\t{strand}\t{gene}\t{prod}\t{','.join(sorted(kos))}\t{name}\n")
            if off == 0:
                order.append(f"[{marker}]"); continue
            lab = name if name not in ("-", "") else (gene or ("hyp" if prod == "hypothetical protein" else "?"))
            order.append(lab)
            if prod == "hypothetical protein" and not kos:
                hyp += 1
            if marker == "hydB":
                if abs(off) <= 2 and (kos & HYD_SMALL or HYD_SMALL_RE.search(prod)):
                    small = 1
                if kos & HYD_MAT or HYD_MAT_RE.search(prod) or HYD_MAT_RE.search(gene):
                    mat += 1
            else:
                if "K03519" in kos or (COX_SUB_RE.search(prod) and "medium" in prod.lower()):
                    cm = 1
                if "K03518" in kos or (COX_SUB_RE.search(prod) and "small" in prod.lower()):
                    cs = 1
                if COX_ACC_RE.search(prod) or COX_ACC_RE.search(gene) or COX_ACC_RE.search(name):
                    acc += 1
        per.write(f"{rep}\t{pid}\t{marker}\t{cls}\t{c}\t{len(genes)}\t{i}\t{min(i, len(genes) - 1 - i)}\t"
                  f"{small}\t{mat}\t{cm}\t{cs}\t{acc}\t{hyp}\t{' '.join(order)}\n")
        n += 1
per.close(); nb.close()
print(f"[operon] marker proteins with context: {n}; reps without GFF: {miss}", flush=True)
PY
date -Iseconds > "$OUT/OPERON_DONE"
