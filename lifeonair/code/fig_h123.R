#!/usr/bin/env Rscript
# Life on air (P3): Figure 1 (H1-H3) and supplementary CoxL tree.
#   a  hydB large-subunit tree: EQ cluster representatives (70% id) + HydDB [NiFe] group 1/2 refs
#   b  PlantRoot vs Deep (same 5-10 cm depth): per site x trip log2 fold change per marker (H1)
#   c  Surface vs Deep: per site x trip log2 fold change per marker (H2)
#   d  Deep site-median group 1l genome equivalents vs ERA5 mean annual precipitation (H3)
# Colours: categorical slots 1 (blue, group 1l) and 2 (orange, group 1h), validated with the
# dataviz palette checker (light surface); other markers in neutral grey.
#
#   Rscript lifeonair/code/fig_h123.R

suppressMessages(library(ape))
IN <- "lifeonair/results"; OUT <- "lifeonair/figures"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.delim(f, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
C1L <- "#2a78d6"; C1H <- "#eb6834"; GREY <- "#9a9a94"; INK <- "#2b2b28"; MUTED <- "#6b6b66"

cls <- rd(file.path(IN, "inputs/atmo_classification.tsv"))
cls$tip <- paste("EQ", cls$rep_id, cls$protein_id, sep = "|")
tip_group <- function(t) {
  ifelse(startsWith(t, "EQ|"), sub("^\\[NiFe\\]_Group_", "", cls$assignment[match(t, cls$tip)]),
         sub(".*NiFe_Group_", "", t))
}
tip_col <- function(g) ifelse(g == "1l", C1L, ifelse(g == "1h", C1H, GREY))

c1 <- rd(file.path(IN, "h123/h1_cell_log2fc_plantroot_vs_deep.tsv"))
ge <- rd(file.path(IN, "reads/guild_genome_equivalents.tsv")); ge <- ge[ge$id_floor == 60 & ge$scg_genes_detected >= 10, ]
meta <- rd(file.path(IN, "inputs/samples.tsv")); d <- merge(meta, ge, by.x = "sample_id", by.y = "sample")
site_env <- rd(file.path(IN, "inputs/env_site.tsv"))
cells <- function(a, b, vars) {
  k <- paste(d$site, d$trip); out <- NULL
  for (z in unique(k)) { x <- d[k == z & d$compartment == a, ]; y <- d[k == z & d$compartment == b, ]
    if (nrow(x) && nrow(y)) out <- rbind(out, sapply(vars, function(v) log2((median(x[[v]]) + .01) / (median(y[[v]]) + .01)))) }
  as.data.frame(out)
}
MK <- c("hydB_1l", "hydB_1h", "coxL_formI", "hydB_other_nife", "coxL_other")
LAB <- c("[NiFe] 1l", "[NiFe] 1h", "CoxL form I", "other [NiFe]", "other CoxL")
COL <- c(C1L, C1H, INK, GREY, GREY)
c2 <- cells("Surface", "Deep", MK)

strip <- function(df, title, n) {
  set.seed(1)
  plot(NA, xlim = c(0.5, 5.5), ylim = c(-7, 4), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  abline(h = seq(-6, 4, 2), col = "#ececea", lwd = 1); abline(h = 0, col = MUTED, lwd = 1)
  axis(2, at = seq(-6, 4, 2), las = 1, col = NA, col.ticks = NA, col.axis = MUTED, cex.axis = 0.8)
  mtext("log2 fold change", side = 2, line = 2.2, cex = 0.75, col = MUTED)
  for (i in seq_along(MK)) {
    y <- df[[MK[i]]]
    points(i + runif(length(y), -0.22, 0.22), y, pch = 16, cex = 0.55, col = adjustcolor(COL[i], 0.45))
    segments(i - 0.3, median(y), i + 0.3, median(y), lwd = 2.5, col = INK)
  }
  text(1:5, -7, LAB, cex = 0.72, col = INK, srt = 30, adj = c(1, 0.5), xpd = NA)
  title(main = title, adj = 0, cex.main = 0.9, font.main = 1, col.main = INK)
  mtext(sprintf("%d site x trip pairs; bar = median", n), side = 3, adj = 0, line = -0.2, cex = 0.62, col = MUTED)
}

draw <- function(dev_open) {
  dev_open()
  layout(matrix(c(1, 2, 1, 3, 1, 4), nrow = 3, byrow = TRUE), widths = c(1.25, 1))
  par(mar = c(1, 1, 2, 1), family = "sans")
  tr <- read.tree(file.path(IN, "trees/hydB.nwk"))
  # root on the clade that best separates the HydDB group 2 references (F1 over all clades and complements)
  n <- length(tr$tip.label); g2 <- grepl("Group_2", tr$tip.label); best <- c(0, NA)
  for (nd in (n + 2):(n + tr$Nnode)) { inn <- tr$tip.label %in% extract.clade(tr, nd)$tip.label
    for (s in list(inn, !inn)) { f1 <- 2 * sum(s & g2) / (sum(s) + sum(g2)); if (f1 > best[1]) best <- c(f1, nd) } }
  tr <- ladderize(root(tr, node = best[2], resolve.root = TRUE))
  g <- tip_group(tr$tip.label); eq <- startsWith(tr$tip.label, "EQ|")
  plot(tr, type = "fan", show.tip.label = FALSE, edge.color = "#c9c9c4", edge.width = 0.5, no.margin = FALSE)
  tiplabels(pch = ifelse(eq, 16, 1), col = tip_col(g), cex = ifelse(eq, 0.45, 0.35))
  title(main = "a  [NiFe] large subunits, groups 1 and 2, rooted on group 2", adj = 0, cex.main = 0.9, font.main = 1, col.main = INK)
  legend("bottomleft", bty = "n", cex = 0.7, text.col = INK,
         legend = c(sprintf("group 1l (EQ %d)", sum(eq & g == "1l")), sprintf("group 1h (EQ %d)", sum(eq & g == "1h")),
                    "other groups", "filled = Empty Quarter, open = HydDB"),
         pch = c(16, 16, 16, NA), col = c(C1L, C1H, GREY, NA))
  par(mar = c(4.2, 3.6, 2.2, 0.8))
  strip(c1[, MK], "b  root-adjacent vs bulk soil, both 5-10 cm", nrow(c1))
  strip(c2, "c  surface (0 cm) vs bulk soil (5-10 cm)", nrow(c2))
  s <- aggregate(d$hydB_1l[d$compartment == "Deep"], list(site = d$site[d$compartment == "Deep"]), median)
  s <- merge(s, site_env, by = "site")
  rho <- cor.test(s$x, s$map_era5_mm, method = "spearman", exact = FALSE)
  plot(s$map_era5_mm, s$x, log = "x", pch = 16, col = adjustcolor(C1L, 0.8), cex = 0.9, bty = "n",
       xlab = "", ylab = "", las = 1, col.axis = MUTED, cex.axis = 0.8, fg = MUTED)
  mtext("ERA5 mean annual precipitation 1995-2024 (mm, log)", side = 1, line = 2.3, cex = 0.65, col = MUTED)
  mtext("group 1l copies per genome", side = 2, line = 2.4, cex = 0.7, col = MUTED)
  title(main = "d  bulk soil, site medians", adj = 0, cex.main = 0.9, font.main = 1, col.main = INK)
  mtext(sprintf("Spearman rho %.2f, p = %.1e, n = %d sites", rho$estimate, rho$p.value, nrow(s)),
        side = 3, adj = 0, line = -0.2, cex = 0.62, col = MUTED)
  dev.off()
}
draw(function() pdf(file.path(OUT, "fig_h123.pdf"), width = 10, height = 8.5))
draw(function() png(file.path(OUT, "fig_h123.png"), width = 2000, height = 1700, res = 200))

# supplementary: CoxL tree coloured by form-I call
cx <- rd(file.path(IN, "inputs/coxL_class.tsv"))
tr <- ladderize(read.tree(file.path(IN, "trees/coxL.nwk")))
call <- ifelse(startsWith(tr$tip.label, "REF|"), "form_I", cx$form_I_call[match(sub("^EQ\\|", "", tr$tip.label), cx$query)])
pdf(file.path(OUT, "figS_coxL_tree.pdf"), width = 7, height = 7)
par(mar = c(1, 1, 2, 1))
plot(tr, type = "fan", show.tip.label = FALSE, edge.color = "#c9c9c4", edge.width = 0.5)
tiplabels(pch = ifelse(startsWith(tr$tip.label, "REF|"), 1, 16), col = ifelse(call == "form_I", C1L, GREY), cex = 0.45)
legend("bottomleft", bty = "n", cex = 0.75, legend = c(sprintf("form I (EQ %d)", sum(call == "form_I" & startsWith(tr$tip.label, "EQ|"), na.rm = TRUE)),
       "other / ambiguous", "open = Cordero form-I references"), pch = c(16, 16, NA), col = c(C1L, GREY, NA), text.col = INK)
dev.off()
cat("figures written to", OUT, "\n")
