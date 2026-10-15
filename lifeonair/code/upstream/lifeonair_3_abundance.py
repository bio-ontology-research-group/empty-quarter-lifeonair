#!/usr/bin/env python3
"""§2.2 abundance — MAG-direct (CoverM) weighting of the high-affinity trace-gas guild
on the CURRENT 5,438 catalogue. Fully local (no amplicon data needed).

Reads:
  guild reps          : classify/guild_reps.txt  (+ guild_hydB_reps.txt / guild_coxL_reps.txt)
  per-sample CoverM   : rebuild_20260729/picrust2_current_20260803/coverm_current/
                        profiles_full_recompute/<sample>/abund.tsv
                          col1 = Genome (5,438 rep IDs + an 'unmapped' row)
                          col2 = "reps.mmi/<sample>_R1.fastq.gz Relative Abundance (%)"  (per-sample label; use position)
                          values are % summing to 100 INCLUDING unmapped
  canonical 551       : metadata/sample_metadata.tsv  (col1 = sample)

For each canonical sample that has a profile (expected 548): guild abundance is reported two ways:
  whole  = Σ rel_abund(guild reps)                              (fraction of the WHOLE community, incl. unmapped)
  mapped = whole / (100 - unmapped) * 100                       (fraction of the MAG-mapped community)
Reports median / IQR / mean across samples for hydB-only, coxL-only, both, and the union guild.

Writes (under lifeonair_clean/results/):
  guild_abundance_by_sample.tsv , guild_abundance_summary.tsv
"""
import os, csv, glob, statistics, sys

R    = "/ibex/scratch/projects/c2014/rob/metagenomics"
A    = f"{R}/06_atmospheric/lifeonair_clean"
PROF = f"{R}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"
CLS  = f"{A}/classify"
RES  = os.environ.get("RES_OUT", f"{A}/results")   # v551 rerun: RES_OUT=.../results/v551 (2026-09-24)
os.makedirs(RES, exist_ok=True)

def load_set(p):
    if not os.path.exists(p):
        return set()
    with open(p) as fh:
        return {l.strip() for l in fh if l.strip()}

guild      = load_set(f"{CLS}/guild_reps.txt")
guild_hyd  = load_set(f"{CLS}/guild_hydB_reps.txt")
guild_cox  = load_set(f"{CLS}/guild_coxL_reps.txt")
guild_both = guild_hyd & guild_cox
if not guild:
    sys.exit(f"ERROR: no guild reps at {CLS}/guild_reps.txt — run type_atmo_5438.sh first")

canonical = {r["sample"] for r in csv.DictReader(open(f"{R}/metadata/sample_metadata.tsv"), delimiter="\t")}
have      = {os.path.basename(d) for d in glob.glob(f"{PROF}/*") if os.path.isdir(d)}
samples   = sorted(canonical & have)
missing   = sorted(canonical - have)
print(f"[abund] guild union={len(guild)} (hydB={len(guild_hyd)} coxL={len(guild_cox)} both={len(guild_both)})")
print(f"[abund] canonical={len(canonical)} profiled={len(have)} usable={len(samples)} missing_profile={len(missing)}: {missing}")
# (degenerate-profile exclusion happens below, in the sample loop, and is reported after it)

def read_profile(sample):
    """return (dict rep->pct, unmapped_pct)"""
    p = f"{PROF}/{sample}/abund.tsv"
    ab, unmapped = {}, 0.0
    with open(p) as fh:
        next(fh)                                   # header (col2 label is per-sample)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2:
                continue
            g = c[0]
            try:
                v = float(c[1])
            except ValueError:
                continue
            if g == "unmapped":
                unmapped = v
            else:
                ab[g] = v
    return ab, unmapped

def frac(ab, unmapped, gset):
    whole = sum(ab.get(g, 0.0) for g in gset)       # % of whole community (incl. unmapped)
    mapped_tot = 100.0 - unmapped
    mapped = (whole / mapped_tot * 100.0) if mapped_tot > 0 else 0.0
    return whole, mapped

rows = []
degenerate = []
for s in samples:
    ab, un = read_profile(s)
    if un >= 100.0:   # 0 reads mapped (e.g. F5Dr1, 1,232 total reads / 0 mapped) -> every genome row is NaN
        degenerate.append(s)
        continue
    row = {"sample": s, "unmapped_pct": un, "n_guild_present": sum(1 for g in guild if ab.get(g, 0) > 0)}
    for tag, gs in (("union", guild), ("hydB", guild_hyd), ("coxL", guild_cox), ("both", guild_both)):
        w, m = frac(ab, un, gs)
        row[f"{tag}_whole"] = w
        row[f"{tag}_mapped"] = m
    rows.append(row)

cols = ["sample", "unmapped_pct", "n_guild_present",
        "union_whole", "union_mapped", "hydB_whole", "hydB_mapped",
        "coxL_whole", "coxL_mapped", "both_whole", "both_mapped"]
print(f"[abund] degenerate profiles excluded (0 mapped reads): {degenerate}; n_used={len(rows)}")
with open(f"{RES}/guild_abundance_by_sample.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

def summ(vals):
    vals = sorted(vals)
    q = statistics.quantiles(vals, n=4) if len(vals) >= 2 else [vals[0]] * 3
    return statistics.median(vals), q[0], q[2], statistics.mean(vals), min(vals), max(vals)

with open(f"{RES}/guild_abundance_summary.tsv", "w") as o:
    o.write("guild\tdenominator\tmedian\tq1\tq3\tmean\tmin\tmax\tn_samples\n")
    for tag in ("union", "hydB", "coxL", "both"):
        for den in ("whole", "mapped"):
            md, q1, q3, mn, lo, hi = summ([r[f"{tag}_{den}"] for r in rows])
            o.write(f"{tag}\t{den}\t{md:.3f}\t{q1:.3f}\t{q3:.3f}\t{mn:.3f}\t{lo:.3f}\t{hi:.3f}\t{len(rows)}\n")
            print(f"[abund] {tag:5s} {den:6s}: median={md:.2f}%  IQR {q1:.2f}-{q3:.2f}  mean={mn:.2f}%  (n={len(rows)})")

med_un = statistics.median([r["unmapped_pct"] for r in rows])
print(f"[abund] median unmapped = {med_un:.1f}%  -> mapped-to-MAGs community ~{100-med_un:.1f}%")
print(f"[abund] wrote {RES}/guild_abundance_by_sample.tsv + guild_abundance_summary.tsv")
