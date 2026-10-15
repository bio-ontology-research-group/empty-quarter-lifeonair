#!/usr/bin/env python3
"""Life on air (P3): in-silico mock communities for the genome-equivalent estimator (review round 2, R2-05).

    python3 b3_mock_sim.py select  WORKDIR N_VISITS SEED
    python3 b3_mock_sim.py simulate WORKDIR COMMUNITY N_READS SEED
    python3 b3_mock_sim.py novel_db WORKDIR COMMUNITY
    python3 b3_mock_sim.py compare WORKDIR

select    picks N_VISITS site visits with both a bulk (Deep) and a root-adjacent (PlantRoot) sample in the
          CoverM table (stratified over expeditions) and writes one community per sample: the sample's CoverM
          relative abundances over the 5,438 representatives (unmapped reads dropped), WORKDIR/communities.tsv
          and WORKDIR/comm/<community>.abund.tsv.
simulate  draws N_READS reads from the representative genome FASTA files, genome g receiving reads in
          proportion to its relative abundance (CoverM relative abundance is a share of reads, i.e. proportional
          to cells x genome length), uniform positions on contigs >= 150 nt, random strand, 150 nt with
          Poisson(0.75) substitutions (0.5% per base), then truncated to the first 100 nt exactly as
          r2_reads_array.sh does (seqtk trimfq -L 100). The true copies per genome of each group is the
          cell-weighted mean marker copy number, sum_g c_g x m_g / sum_g c_g with c_g = RA_g / L_g and m_g the
          number of typed marker proteins of g on retained contigs (marker_labels_v2, minus the _sub proteins on
          removed contigs).
novel_db  writes the marker+SCG database without any EQ marker protein or SCG cluster representative that
          comes from a genome of this community (realistic novelty: reads must be matched by relatives).
compare   joins the r3_aggregate.py output with the truth table.
"""
import csv
import gzip
import math
import os
import random
import sys
from collections import Counter, defaultdict

R = "/ibex/scratch/projects/c2014/rob/metagenomics"
REPS = f"{R}/05_refine/_drep/out/dereplicated_genomes"
COVERM = f"{R}/06_atmospheric/lifeonair_reads/coverm_rep_abundance_long.tsv.gz"
DB = f"{R}/06_atmospheric/lifeonair_reads/db"
GROUPS = ("hydB_1l", "hydB_1h", "hydB_2a", "hydB_other_nife", "coxL_formI", "coxL_other")


def group(marker, label):                    # identical mapping to r3_aggregate.py with v2 labels
    if marker == "hydB":
        for k, g in (("_NiFe__Group_1h", "hydB_1h"), ("_NiFe__Group_1l", "hydB_1l"), ("_NiFe__Group_2a", "hydB_2a")):
            if label == k:
                return g
        return "hydB_other_nife" if "NiFe" in label else "hydB_other"
    if marker == "coxL":
        return "coxL_formI" if label == "form_I" else "coxL_other"
    return marker


def read_fa(path):
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


def copies(work):
    inp = f"{work}/../inputs"
    subx = {l.strip() for l in open(f"{inp}/sub_excluded_db_proteins.txt") if l.strip()}
    m = defaultdict(Counter)
    for r in csv.DictReader(open(f"{inp}/marker_labels_v2.tsv"), delimiter="\t"):
        s = r["subject"]
        if not s.startswith("EQ|"):
            continue
        _, _, rep, pid = s.split("|")
        if f"{rep}|{pid}" in subx:
            continue
        m[rep][group(r["marker"], r["label"])] += 1
    return m


def select(work, nvis, seed):
    samples = {r["sample_id"]: r for r in csv.DictReader(open(f"{work}/../inputs/samples.tsv"), delimiter="\t")}
    ab = defaultdict(dict)
    with gzip.open(COVERM, "rt") as fh:
        for line in fh:
            s, g, v = line.rstrip("\n").split("\t")
            if g != "unmapped" and float(v) > 0:
                ab[s][g] = float(v)
    vis = defaultdict(dict)
    for s, r in samples.items():
        if s in ab and r["compartment"] in ("Deep", "PlantRoot"):
            vis[(r["site"], r["trip"])].setdefault(r["compartment"], []).append(s)
    ok = sorted(k for k, v in vis.items() if "Deep" in v and "PlantRoot" in v)
    rng = random.Random(seed)
    by_trip = defaultdict(list)
    for k in ok:
        by_trip[k[1]].append(k)
    chosen = []
    trips = sorted(by_trip)
    while len(chosen) < nvis and any(by_trip.values()):
        for t in trips:
            if by_trip[t] and len(chosen) < nvis:
                chosen.append(by_trip[t].pop(rng.randrange(len(by_trip[t]))))
    os.makedirs(f"{work}/comm", exist_ok=True)
    with open(f"{work}/communities.tsv", "w") as o:
        o.write("community\tsample\tsite\ttrip\tcompartment\tn_genomes\tmapped_pct\n")
        for site, trip in chosen:
            for comp in ("Deep", "PlantRoot"):
                s = sorted(vis[(site, trip)][comp])[0]
                c = f"m_{s}"
                with open(f"{work}/comm/{c}.abund.tsv", "w") as a:
                    for g, v in sorted(ab[s].items()):
                        a.write(f"{g}\t{v}\n")
                o.write(f"{c}\t{s}\t{site}\t{trip}\t{comp}\t{len(ab[s])}\t{sum(ab[s].values()):.2f}\n")
    print(f"selected {len(chosen)} visits from {len(ok)}")


POIS = [math.exp(-0.75) * 0.75 ** k / math.factorial(k) for k in range(8)]
POIS_CUM = [sum(POIS[:k + 1]) for k in range(8)]
RC = str.maketrans("ACGTN", "TGCAN")


def simulate(work, comm, nreads, seed):
    rng = random.Random(seed)
    ab = [l.rstrip("\n").split("\t") for l in open(f"{work}/comm/{comm}.abund.tsv")]
    tot = sum(float(v) for _, v in ab)
    m = copies(work)
    num = Counter(); den = 0.0; carr = 0.0; nwritten = 0
    out = open(f"{work}/comm/{comm}.reads.fa", "w")
    for g, v in ab:
        n = int(round(nreads * float(v) / tot))
        contigs = [s.upper() for _, s in read_fa(f"{REPS}/{g}.fa")]
        L = sum(len(s) for s in contigs)
        c = float(v) / L
        den += c
        for k in GROUPS:
            num[k] += c * m[g][k]
        carr += c * (m[g]["hydB_1l"] > 0)
        usable = [s for s in contigs if len(s) >= 150]
        if not usable or n == 0:
            continue
        cw, acc = [], 0
        for s in usable:
            acc += len(s) - 149; cw.append(acc)
        for i, s in enumerate(rng.choices(usable, cum_weights=cw, k=n)):
            p = rng.randrange(len(s) - 149)
            r = s[p:p + 150]
            if rng.random() < 0.5:
                r = r.translate(RC)[::-1]
            u = rng.random(); ne = next(k for k, q in enumerate(POIS_CUM) if u <= q or k == 7)
            if ne:
                r = list(r)
                for j in rng.sample(range(150), ne):
                    r[j] = rng.choice([b for b in "ACGT" if b != r[j]])
                r = "".join(r)
            out.write(f">{comm}_{nwritten}\n{r[:100]}\n"); nwritten += 1
    out.close()
    with open(f"{work}/comm/{comm}.truth.tsv", "w") as o:
        o.write("community\treads\t" + "\t".join(f"true_{k}" for k in GROUPS) + "\ttrue_carrier_frac_1l\n")
        o.write(f"{comm}\t{nwritten}\t" + "\t".join(f"{num[k] / den:.5f}" for k in GROUPS) + f"\t{carr / den:.5f}\n")
    print(f"{comm}: {nwritten} reads; true 1l {num['hydB_1l'] / den:.4f}")


def novel_db(work, comm):
    reps = {l.split("\t")[0] for l in open(f"{work}/comm/{comm}.abund.tsv")}
    keep = False; n = drop = 0
    with open(f"{work}/comm/{comm}.novel.faa", "w") as o:
        for line in open(f"{DB}/loa_reads.faa"):
            if line.startswith(">"):
                f = line[1:].split()[0].split("|")
                src = f[2] if f[0] in ("EQ", "SCG") else None
                keep = src not in reps
                n += 1; drop += not keep
            if keep:
                o.write(line)
    print(f"{comm}: novel DB drops {drop} of {n} proteins")


def compare(work):
    comms = {r["community"]: r for r in csv.DictReader(open(f"{work}/communities.tsv"), delimiter="\t")}
    truth = {}
    for c in comms:
        truth[c] = next(csv.DictReader(open(f"{work}/comm/{c}.truth.tsv"), delimiter="\t"))
    rows = []
    for variant in ("full", "novel"):
        f = f"{work}/agg_{variant}/mock_ge.tsv"
        if not os.path.exists(f):
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            c = r["sample"]
            meta = comms[c]
            for k in GROUPS:
                rows.append([variant, r["id_floor"], c, meta["sample"], meta["site"], meta["trip"], meta["compartment"],
                             k, truth[c][f"true_{k}"], r[k], r["scg_genes_detected"], r["reads_screened"]])
    with open(f"{work}/mock_ge_vs_truth.tsv", "w") as o:
        o.write("db\tid_floor\tcommunity\tsample\tsite\ttrip\tcompartment\tgroup\ttrue_ge\test_ge\tscg_genes\treads\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {len(rows)} rows")


if __name__ == "__main__":
    cmd, work = sys.argv[1], sys.argv[2]
    if cmd == "select":
        select(work, int(sys.argv[3]), int(sys.argv[4]))
    elif cmd == "simulate":
        simulate(work, sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    elif cmd == "novel_db":
        novel_db(work, sys.argv[3])
    elif cmd == "compare":
        compare(work)
