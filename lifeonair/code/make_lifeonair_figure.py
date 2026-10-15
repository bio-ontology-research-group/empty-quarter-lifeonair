#!/usr/bin/env python3
"""Section 2.2 figure — the trace-gas guild: who carries it, and how much of the community it is.

DESIGN DECISION (2026-08-23, session 09-display). Two panels on the 5,438 catalogue / 548 profiled samples:

  (a) MEMBERSHIP: species carrying the high-affinity [NiFe]-hydrogenase (HydDB 1h/1l), form-I CODH,
      either, and both — counts and % of the catalogue, direct-labelled. The 1l-over-1h dominance
      (1,124 vs 445 proteins) is annotated once; it is a description of the pool, not an aridity claim.
  (b) ABUNDANCE: per-sample share of the guild (union) as ECDFs under both denominators, whole community
      and mapped community. Two denominators are two readings of one guild, so one hue with the mapped
      reading de-emphasised in grey; medians (18.4% / 43.9%) are marked and direct-labelled.

Variant written to figures/options/: _opt_tree adds the 269-tip hydrogenase-marker tree coloured by
HydDB call (high-affinity 1l/1h in blue, circles vs triangles; other groups grey) — evidence the calls
sit in coherent clades, not scattered. The tree holds catalogue proteins only (no reference tips), so it
shows classification structure, not reference placement.

Colour: guild/focal #4477aa; grey #999999 for the alternative denominator and non-high-affinity groups.
Reads  M/06_atmospheric/lifeonair_clean/classify/{atmo_classification.tsv,hydB_class.tsv,hydB.nwk}
       and results/guild_abundance_by_sample.tsv.
Writes lifeonair/figures/fig_lifeonair_guild.{pdf,png} and lifeonair/figures/options/fig_lifeonair_opt_tree.pdf.
Membership counts and abundance medians are asserted against the registered §2.2 numbers
(968 / 1,170 / 1,860 / 278; 1l 1,124 vs 1h 445 proteins; medians 18.4% whole, 43.9% mapped).
"""
import os, csv, statistics
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
FIG = os.path.join(REPO, "lifeonair", "figures")
OPT = os.path.join(FIG, "options")
WB = "/ibex/scratch/projects/c2014/rob/metagenomics/06_atmospheric/lifeonair_clean"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
import sys; sys.path.insert(0, os.path.join(REPO, "code", "figures"))
from palette import C_DARK, INK, MUTED, GRID, GREY
N_CAT = 5438

# Membership comes from the run's shipped guild lists (the registered basis), not a re-derivation:
# atmo_classification.tsv carries a handful of extra high-affinity rows relative to the final lists.
def reps_of(fn):
    return set(x.strip() for x in open(os.path.join(WB, "classify", fn)) if x.strip())
hyd, cox = reps_of("guild_hydB_reps.txt"), reps_of("guild_coxL_reps.txt")
either, both = reps_of("guild_reps.txt"), reps_of("guild_both_reps.txt")
n_1l = n_1h = 0
for r in csv.DictReader(open(os.path.join(WB, "classify", "atmo_classification.tsv")), delimiter="\t"):
    if r["is_high_affinity"] == "1" and r["marker"] == "hydB":
        if r["assignment"].endswith("_1l"): n_1l += 1
        if r["assignment"].endswith("_1h"): n_1h += 1
assert (len(hyd), len(cox), len(either), len(both)) == (968, 1170, 1860, 278), \
    (len(hyd), len(cox), len(either), len(both))
assert (n_1l, n_1h) == (1124, 445), (n_1l, n_1h)

ab = list(csv.DictReader(open(os.path.join(WB, "results", "guild_abundance_by_sample.tsv")), delimiter="\t"))
kw = [k for k in ab[0] if "union" in k and "whole" in k][0]
km = [k for k in ab[0] if "union" in k and "mapped" in k][0]
whole = sorted(float(r[kw]) for r in ab)
mapped = sorted(float(r[km]) for r in ab)
assert len(whole) == 548 and round(statistics.median(whole), 1) == 18.4 \
    and round(statistics.median(mapped), 1) == 43.9, (len(whole), statistics.median(whole), statistics.median(mapped))
print(f"asserted: guild {len(hyd)}/{len(cox)}/{len(either)}/{len(both)}; 1l/1h {n_1l}/{n_1h}; "
      f"medians {statistics.median(whole):.1f}/{statistics.median(mapped):.1f} (n=548)")

def panel_membership(a):
    rows = [("high-affinity [NiFe]-hydrogenase (1h/1l)", len(hyd)),
            ("form-I CO dehydrogenase", len(cox)),
            ("either marker", len(either)),
            ("both markers", len(both))]
    ys = range(len(rows) - 1, -1, -1)
    a.barh(list(ys), [n for _, n in rows], height=0.62, color=C_DARK, lw=0)
    for y, (lab, n) in zip(ys, rows):
        a.text(n + 25, y, f"{n:,}  ({100*n/N_CAT:.0f}%)", va="center", fontsize=8.5, color=INK)
    a.set_yticks(list(ys)); a.set_yticklabels([r[0] for r in rows], fontsize=8.5)
    a.set_xlabel(f"species representatives (of {N_CAT:,})")
    a.set_xlim(0, 2400)
    a.grid(axis="x", color=GRID, lw=0.6); a.set_axisbelow(True)
    a.set_title("a  A third of the catalogue carries the machinery", loc="left", fontsize=10)
    a.text(0.97, 0.20, f"hydrogenase pool: group 1l {n_1l:,} proteins,\ngroup 1h {n_1h} "
           "(arid-soil form predominates)", transform=a.transAxes, fontsize=8, color=MUTED, ha="right")

def panel_abundance(b):
    n = len(whole)
    for xs, col, lab in ((whole, C_DARK, "of the whole community"),
                         (mapped, GREY, "of the mapped community")):
        b.step(xs, [i / n for i in range(1, n + 1)], where="post", color=col, lw=1.8)
        med = statistics.median(xs)
        b.scatter([med], [0.5], s=40, color=col, zorder=5, edgecolor="white", lw=0.8)
        if col == C_DARK:
            b.text(0.35, 0.13, f"median {med:.1f}%\n{lab}", transform=b.transAxes,
                   fontsize=8, color=col, ha="right")
        else:
            b.annotate(f"median {med:.1f}%\n{lab}", (med, 0.5), xytext=(10, 8),
                       textcoords="offset points", fontsize=8, color=MUTED)
    b.set_xlabel("guild share per sample (%)")
    b.set_ylabel("cumulative fraction of 548 samples")
    b.set_xlim(0, 90); b.set_ylim(0, 1.02)
    b.grid(color=GRID, lw=0.6); b.set_axisbelow(True)
    b.set_title("b  A pervasive baseline, not a niche", loc="left", fontsize=10)

def read_tree_panel(c):
    from Bio import Phylo
    groups = {}
    for row in csv.reader(open(os.path.join(WB, "classify", "hydB_class.tsv")), delimiter="\t"):
        if len(row) >= 2:
            groups[row[0]] = row[1]
    tree = Phylo.read(os.path.join(WB, "classify", "hydB.nwk"), "newick")
    tree.ladderize()
    tips = tree.get_terminals()
    depth = tree.depths()
    ymap = {t: i for i, t in enumerate(tips)}
    # branch lines
    def draw(cl, x0):
        x1 = depth[cl]
        y = ymap[cl] if cl.is_terminal() else statistics.mean(
            (min(ymap[t] for t in cl.get_terminals()), max(ymap[t] for t in cl.get_terminals())))
        c.plot([x0, x1], [y, y], color="#bbbbbb", lw=0.5, zorder=1)
        if not cl.is_terminal():
            kids = cl.clades
            yy = [ymap[k] if k.is_terminal() else statistics.mean(
                (min(ymap[t] for t in k.get_terminals()), max(ymap[t] for t in k.get_terminals())))
                for k in kids]
            c.plot([x1, x1], [min(yy), max(yy)], color="#bbbbbb", lw=0.5, zorder=1)
            for k in kids:
                draw(k, x1)
        return y
    import sys
    sys.setrecursionlimit(10000)
    draw(tree.root, 0.0)
    xmax = max(depth.values())
    n_query = sum(1 for t in tips if t.name and "|" in t.name)
    for t in tips:
        if t.name and "|" not in t.name:      # HydDB reference sequence retained in the alignment
            c.scatter([depth[t]], [ymap[t]], s=16, marker="x", color=INK, lw=1.0, zorder=4)
            continue
        g = groups.get(t.name, "")
        hi = g.endswith("_1l") or g.endswith("_1h")
        c.scatter([depth[t]], [ymap[t]], s=10 if hi else 6,
                  marker="o" if g.endswith("_1l") else ("^" if g.endswith("_1h") else "o"),
                  color=C_DARK if hi else GREY, alpha=0.9 if hi else 0.5, lw=0, zorder=3)
    c.scatter([], [], s=14, marker="o", color=C_DARK, label="group 1l (high-affinity)")
    c.scatter([], [], s=14, marker="^", color=C_DARK, label="group 1h (high-affinity)")
    c.scatter([], [], s=10, marker="o", color=GREY, label="other [NiFe] groups")
    c.scatter([], [], s=16, marker="x", color=INK, lw=1.0, label="HydDB reference")
    c.legend(loc="lower right", frameon=False, fontsize=8)
    c.set_xlim(-0.02 * xmax, 1.05 * xmax)
    c.set_xticks([]); c.set_yticks([])
    for s in ("left", "bottom"):
        c.spines[s].set_visible(False)
    c.set_title(f"c  High-affinity calls form coherent clades ({n_query} cluster representatives)",
                loc="left", fontsize=10)

os.makedirs(FIG, exist_ok=True); os.makedirs(OPT, exist_ok=True)
fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 4.6))
panel_membership(a); panel_abundance(b)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_lifeonair_guild.{ext}"), bbox_inches="tight", dpi=200)
plt.close(fig)

fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.6, 4.6))
panel_membership(a); panel_abundance(b); read_tree_panel(c)
fig.tight_layout()
fig.savefig(os.path.join(OPT, "fig_lifeonair_opt_tree.pdf"), bbox_inches="tight")
plt.close(fig)
print("wrote lifeonair/figures/fig_lifeonair_guild.{pdf,png} + 1 option")
