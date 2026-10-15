#!/usr/bin/env Rscript
# Life on air (P3): figures after review round 1 (supersedes fig_h123.R and fig_mechanism.R).
#   Figure 1  a [NiFe] large-subunit tree (groups 1l, 1h, 2a); b root-adjacent vs bulk and c surface vs bulk
#             log2 ratios per site visit; d excess change beyond the controls with site-bootstrap 95% CIs
#   Figure 2  a group 1l GE against extract DNA concentration; b genus depletion beside roots against the
#             fraction of the genus carrying group 1l; c root effect on group 1l under successive adjustments;
#             d Spearman correlation with ERA5 mean annual precipitation per position, with 6-block bootstrap CIs
#   Figure S1 CoxL tree with corrected form calls (R1-01)
#   Figure S2 within-site change against antecedent precipitation (trip-adjusted), scaled to the within-site
#             interquartile range, with the +/-0.5 log2 equivalence band
# Inputs: lifeonair/results/h123_v2, h1_mechanism_v2, reads/guild_genome_equivalents_v2.tsv, inputs/.
#   Rscript lifeonair/code/fig_v2.R
suppressMessages(library(ape))
IN <- "lifeonair/results"; OUT <- "lifeonair/figures"
rd <- function(f) read.delim(f, quote = "", stringsAsFactors = FALSE, check.names = FALSE)
C1L <- "#2a78d6"; C1H <- "#eb6834"; C2A <- "#1f9e74"; CFI <- "#8a5cc2"; GREY <- "#9a9a94"; INK <- "#2b2b28"; MUTED <- "#6b6b66"
POSCOL <- c(Surface = "#d9a441", Deep = "#6b6b66", PlantRoot = "#1f9e74")
MK <- c("hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI", "hydB_other_nife", "coxL_other")
LAB <- c("[NiFe] 1l", "[NiFe] 1h", "[NiFe] 2a", "CoxL form I", "other [NiFe]", "other CoxL")
COL <- c(C1L, C1H, C2A, CFI, GREY, GREY)
ttl <- function(s) title(main = s, adj = 0, cex.main = 1.05, font.main = 2, col.main = INK)
ORPHANS <- readLines(file.path(IN, "reads/orphan_sub_reps.txt"))

cls <- rd(file.path(IN, "inputs/atmo_classification.tsv"))
cls$tip <- paste("EQ", cls$rep_id, cls$protein_id, sep = "|")
grp <- function(t) ifelse(startsWith(t, "EQ|"), sub("^\\[NiFe\\]_Group_", "", cls$assignment[match(t, cls$tip)]),
                          sub(".*NiFe_Group_", "", t))
gcol <- function(g) ifelse(g == "1l", C1L, ifelse(g == "1h", C1H, ifelse(g == "2a", C2A, GREY)))

c1 <- rd(file.path(IN, "h123_v2/h1_visit_log2fc_root_vs_bulk.tsv"))
c2 <- rd(file.path(IN, "h123_v2/h1_visit_log2fc_surface_vs_bulk.tsv"))
ct <- rd(file.path(IN, "h123_v2/contrasts_by_position.tsv"))

strip <- function(df, title, n) {
  set.seed(1)
  plot(NA, xlim = c(0.5, 6.5), ylim = c(-7, 4), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  abline(h = seq(-6, 4, 2), col = "#ececea"); abline(h = 0, col = MUTED)
  axis(2, at = seq(-6, 4, 2), las = 1, col = NA, col.ticks = NA, col.axis = INK, cex.axis = 0.95)
  mtext(expression(log[2] ~ "ratio of genome equivalents"), side = 2, line = 2.3, cex = 0.8, col = INK)
  for (i in seq_along(MK)) {
    y <- df[[MK[i]]]
    points(i + runif(length(y), -0.22, 0.22), y, pch = 16, cex = 0.6, col = adjustcolor(COL[i], 0.45))
    segments(i - 0.3, median(y), i + 0.3, median(y), lwd = 3, col = INK)
  }
  text(1:6, -7.3, LAB, cex = 0.9, col = INK, srt = 35, adj = c(1, 0.5), xpd = NA)
  ttl(title); mtext(sprintf("%d site visits; bar = median", n), side = 3, adj = 0, line = -0.3, cex = 0.9, col = MUTED)
}

fig1 <- function(dev_open) {
  dev_open()
  layout(matrix(c(1, 2, 1, 3, 4, 4), nrow = 3, byrow = TRUE), widths = c(1.15, 1), heights = c(1, 1, 1.2))
  par(mar = c(1, 1, 2.2, 1), family = "sans", cex = 1)
  tr <- read.tree(file.path(IN, "trees/hydB.nwk"))
  n <- length(tr$tip.label); g2 <- grepl("Group_2", tr$tip.label); best <- c(0, NA)
  for (nd in (n + 2):(n + tr$Nnode)) { inn <- tr$tip.label %in% extract.clade(tr, nd)$tip.label
    for (s in list(inn, !inn)) { f1 <- 2 * sum(s & g2) / (sum(s) + sum(g2)); if (f1 > best[1]) best <- c(f1, nd) } }
  tr <- ladderize(root(tr, node = best[2], resolve.root = TRUE, edgelabel = TRUE))   # support values belong to edges
  g <- grp(tr$tip.label); eq <- startsWith(tr$tip.label, "EQ|")
  plot(tr, type = "fan", show.tip.label = FALSE, edge.color = "#c9c9c4", edge.width = 0.5)
  tiplabels(pch = ifelse(eq, 16, 1), col = gcol(g), cex = ifelse(eq, 0.5, 0.4))
  # FastTree SH-like support of the smallest clade holding all tips of each focal group (R2-26)
  sup <- suppressWarnings(as.numeric(tr$node.label))
  for (gg in c("1l", "1h", "2a")) { nd <- getMRCA(tr, tr$tip.label[g == gg])
    nodelabels(sprintf("%s %.2f", gg, sup[nd - n]), node = nd, frame = "none", cex = 0.8, col = gcol(gg), font = 2) }
  ttl("a  [NiFe] large-subunit tree")
  legend("bottomleft", bty = "n", cex = 0.9, text.col = INK,
         legend = c(sprintf("group 1l (EQ clusters %d)", sum(eq & g == "1l")), sprintf("group 1h (EQ %d)", sum(eq & g == "1h")),
                    sprintf("group 2a (EQ %d)", sum(eq & g == "2a")), "other groups", "filled = Rub' al-Khali, open = HydDB"),
         pch = c(16, 16, 16, 16, NA), col = c(C1L, C1H, C2A, GREY, NA))
  mtext("node labels: group and SH-like support", side = 1, adj = 0, line = -0.2, cex = 0.9, col = MUTED)
  par(mar = c(4.1, 3.5, 2.4, 0.5))
  strip(c1, "b  root-adjacent / bulk", nrow(c1))
  strip(c2, "c  surface / bulk", nrow(c2))
  # d excess vs controls, three contrasts
  par(mar = c(3.6, 11, 2.4, 1))
  ex <- ct[grepl("^excess_", ct$quantity), ]
  ex$m <- sub("^excess_", "", ex$quantity)
  ord <- expand.grid(m = c("hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI"), c = c("root/bulk", "surface/bulk", "root/surface"), stringsAsFactors = FALSE)
  ex <- ex[match(paste(ord$c, ord$m), paste(ex$contrast, ex$m)), ]
  y <- rev(seq_len(nrow(ex)) + rep(0:2, each = 4) * 0.8)
  plot(NA, xlim = c(-3, 1), ylim = range(y) + c(-0.5, 0.5), yaxt = "n", xlab = "", ylab = "", bty = "n", las = 1, cex.axis = 0.95)
  abline(v = 0, col = MUTED); abline(v = seq(-3, 1, 0.5), col = "#ececea")
  mc <- setNames(COL[1:4], MK[1:4])
  segments(ex$ci_lo, y, ex$ci_hi, y, lwd = 2.5, col = mc[ex$m]); points(ex$median_log2, y, pch = 16, cex = 1.3, col = mc[ex$m])
  text(par("usr")[1] - 0.05, y, paste(ex$contrast, LAB[match(ex$m, MK)]), adj = 1, cex = 0.9, xpd = NA, col = INK)
  mtext(expression("excess " * log[2] * " ratio (site-bootstrap 95% CI)"), side = 1, line = 2.4, cex = 0.8, col = INK)
  ttl("d  excess change beyond controls")
  dev.off()
}
# sized so that text is >= 7 pt when the figure is set at the text width (R2-26)
fig1(function() pdf(file.path(OUT, "fig1_guilds.pdf"), width = 6.37, height = 7.1, pointsize = 9))
fig1(function() png(file.path(OUT, "fig1_guilds.png"), width = 1593, height = 1775, res = 250, pointsize = 9))

# ---------------------------------------------------------------------------------------------- Figure 2
ge <- rd(file.path(IN, "reads/guild_genome_equivalents_v2.tsv")); ge <- ge[ge$id_floor == 60 & ge$scg_genes_detected >= 10, ]
meta <- rd(file.path(IN, "inputs/samples.tsv")); meta <- meta[!is.na(meta$site) & meta$site %in% 1:60, ]
d <- merge(meta, ge, by.x = "sample_id", by.y = "sample")
led <- rd(file.path(IN, "inputs/sample_ledger_dna.tsv")); led$dna <- suppressWarnings(as.numeric(led$dna_concentration))
nconc <- tapply(led$dna, led$sample_id, function(x) length(unique(round(x[is.finite(x)], 6))))
led <- led[led$sample_id %in% names(nconc)[nconc == 1] & is.finite(led$dna) & led$dna > 0, ]
d$dna <- led$dna[match(d$sample_id, led$sample_id)]
a1 <- rd(file.path(IN, "h123_v2/A1_dna_scaling.tsv")); ra <- rd(file.path(IN, "h123_v2/A1_root_effect_dna_adjusted.tsv"))
gen <- rd(file.path(IN, "h1_mechanism_v2/m2d_genus_root_depletion.tsv"))
m2d <- rd(file.path(IN, "h1_mechanism_v2/m2d_genus_models.tsv"))
m2c <- rd(file.path(IN, "h1_mechanism_v2/m2c_excess_loss_vs_copiotrophy_shift.tsv"))
h3 <- rd(file.path(IN, "h123_v2/H3_spatial_by_marker_position.tsv")); h3 <- h3[h3$precip == "map_era5_mm", ]

fig2 <- function(dev_open) {
  dev_open()
  par(mfrow = c(2, 2), mar = c(4.6, 4.8, 2.8, 1), family = "sans"); par(cex = 1)
  x <- d[is.finite(d$dna), ]
  plot(log2(x$dna), log2(x$hydB_1l + 0.01), pch = 16, cex = 0.8, col = adjustcolor(POSCOL[x$compartment], 0.7), bty = "n",
       xlab = expression(log[2] ~ "extract DNA concentration (ng/" * mu * "L)"), ylab = expression(log[2] * "(group 1l GE + 0.01)"),
       las = 1, cex.lab = 1.05, cex.axis = 0.95)
  b <- a1[a1$response == "hydB_1l", ]
  legend("bottomleft", bty = "n", cex = 0.9, pch = 16, col = POSCOL, legend = c("surface", "bulk", "root-adjacent"))
  mtext(sprintf("slope %.2f (SE %.2f); n = %d", b$slope_ldna, b$se, b$n),
        side = 3, adj = 0, line = -0.3, cex = 0.9, col = MUTED)
  ttl("a  group 1l and extract DNA")
  g <- gen[is.finite(gen$log2fc_root) & is.finite(gen$frac_1l), ]
  plot(g$frac_1l + runif(nrow(g), -0.015, 0.015), g$log2fc_root, pch = 16, cex = 0.4 + 6 * sqrt(g$deep_mean_pct / max(g$deep_mean_pct)),
       col = adjustcolor(C1L, 0.45), bty = "n", las = 1, xlab = "genus fraction carrying group 1l",
       ylab = expression(log[2] ~ "ratio, root-adjacent vs bulk"), cex.lab = 1.05, cex.axis = 0.95)
  abline(h = 0, col = MUTED)
  k <- m2d[m2d$model == "growth + carriage + phylum" & m2d$term == "frac_1l", ]
  mtext(sprintf("%d genera; coefficient %.2f (SE %.2f)", k$n_genera, k$Estimate, k[["Std..Error"]]),
        side = 3, adj = 0, line = -0.3, cex = 0.9, col = MUTED)
  ttl("b  genus carriage and depletion")
  # c root effect under adjustments
  ex <- ct[ct$contrast == "root/bulk" & ct$quantity == "excess_hydB_1l", ]
  raw <- ct[ct$contrast == "root/bulk" & ct$quantity == "hydB_1l", ]
  iz <- m2c[grepl("lmer slope", m2c$test) & grepl("excess", m2c$response), ]
  bt <- as.numeric(strsplit(iz$intercept_boot95, " to ")[[1]])
  est <- data.frame(lab = c(sprintf("visit median ratio\n%d visits; bootstrap", raw$n_visits),
                            sprintf("unadjusted model\n%d samples; Wald", ra$n[1]),
                            sprintf("DNA-adjusted model\n%d samples; Wald", ra$n[2]),
                            sprintf("excess vs controls\n%d visits; bootstrap", ex$n_visits),
                            sprintf("zero-shift excess\n%d visits; bootstrap", iz$n)),
                    m = c(raw$median_log2, ra$estimate[1], ra$estimate[2], ex$median_log2, iz$intercept_at_no_shift),
                    lo = c(raw$ci_lo, ra$estimate[1] - 1.96 * ra$se[1], ra$estimate[2] - 1.96 * ra$se[2], ex$ci_lo, bt[1]),
                    hi = c(raw$ci_hi, ra$estimate[1] + 1.96 * ra$se[1], ra$estimate[2] + 1.96 * ra$se[2], ex$ci_hi, bt[2]))
  par(mar = c(4.6, 9, 2.8, 0.5)); y <- rev(seq_len(nrow(est)))
  plot(NA, xlim = c(-3.2, 0.3), ylim = c(0.5, nrow(est) + 0.5), yaxt = "n", bty = "n", las = 1, xlab = expression("group 1l root effect (" * log[2] * ")"), ylab = "", cex.lab = 1.05)
  abline(v = 0, col = MUTED); segments(est$lo, y, est$hi, y, lwd = 2.5, col = C1L); points(est$m, y, pch = 16, cex = 1.3, col = C1L)
  axis(2, at = y, labels = est$lab, las = 1, tick = FALSE, cex.axis = 0.85, padj = 0.5)
  ttl("c  root effects (95% CI)")
  # d aridity
  par(mar = c(4.6, 8, 2.8, 0.5))
  o <- expand.grid(m = c("hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI"), c = c("Deep", "Surface", "PlantRoot"), stringsAsFactors = FALSE)
  h <- h3[match(paste(o$c, o$m), paste(h3$compartment, h3$marker)), ]
  ci <- do.call(rbind, lapply(strsplit(h$block6_ci, " to "), as.numeric))
  y <- rev(seq_len(nrow(h)) + rep(0:2, each = 4) * 0.8)
  plot(NA, xlim = c(-0.9, 0.5), ylim = range(y) + c(-0.5, 0.5), yaxt = "n", bty = "n", las = 1, xlab = "Spearman rho (6-block 95% CI)", ylab = "", cex.lab = 1.05)
  abline(v = 0, col = MUTED); mc <- setNames(COL[1:4], MK[1:4])
  segments(ci[, 1], y, ci[, 2], y, lwd = 2.5, col = mc[h$marker]); points(h$rho, y, pch = 16, cex = 1.3, col = mc[h$marker])
  pos <- c(Deep = "bulk", Surface = "surface", PlantRoot = "root-adj.")
  axis(2, at = y, labels = paste(pos[h$compartment], LAB[match(h$marker, MK)]), las = 1, tick = FALSE, cex.axis = 0.9)
  ttl("d  ERA5 precipitation association")
  dev.off()
}
fig2(function() pdf(file.path(OUT, "fig2_scaling.pdf"), width = 6.37, height = 6.4, pointsize = 9))
fig2(function() png(file.path(OUT, "fig2_scaling.png"), width = 1593, height = 1600, res = 250, pointsize = 9))

# ---------------------------------------------------------------------------------------------- Figure S1
cx <- rd(file.path(IN, "inputs/coxL_class_v2.tsv"))
tr <- ladderize(read.tree(file.path(IN, "trees/coxL.nwk")))
ref <- startsWith(tr$tip.label, "REF|")
call <- ifelse(ref, "form_I", cx$form_call[match(sub("^EQ\\|", "", tr$tip.label), cx$query)])
pdf(file.path(OUT, "figS1_coxL_tree.pdf"), width = 7.5, height = 7.5)
par(mar = c(1, 1, 2, 1))
plot(tr, type = "fan", show.tip.label = FALSE, edge.color = "#c9c9c4", edge.width = 0.5)
tiplabels(pch = ifelse(ref, 1, 16), col = ifelse(call == "form_I", CFI, ifelse(call == "form_II", C1H, GREY)), cex = 0.5)
legend("bottomleft", bty = "n", cex = 0.85, text.col = INK, pch = c(16, 16, 16, NA), col = c(CFI, C1H, GREY, NA),
       legend = c(sprintf("form I, AYXCSFR (EQ %d)", sum(call == "form_I" & !ref, na.rm = TRUE)), sprintf("form II, AYRGAGR (EQ %d)", sum(call == "form_II" & !ref, na.rm = TRUE)),
                  "unresolved", "open = Cordero et al. form-I references"))
dev.off()

# ---------------------------------------------------------------------------------------------- Figure S2
tm <- rd(file.path(IN, "h123_v2/H3_temporal_within_between.tsv")); tm <- tm[tm$exposure %in% c("P90_era5", "P365_era5"), ]
pdf(file.path(OUT, "figS2_temporal.pdf"), width = 6.37, height = 6.5, pointsize = 10)
par(mar = c(4.6, 12, 2.5, 1))
y <- rev(seq_len(nrow(tm)))
plot(NA, xlim = range(c(tm$change_across_IQR_lo, tm$change_across_IQR_hi, -0.6, 0.6)), ylim = c(0.5, nrow(tm) + 0.5), yaxt = "n", bty = "n", las = 1,
     xlab = expression("change across exposure IQR (" * log[2] * ")"), ylab = "")
rect(-0.5, 0, 0.5, nrow(tm) + 1, col = "#efefec", border = NA); abline(v = 0, col = MUTED)
mc <- setNames(COL[1:4], MK[1:4]); mid <- tm$within_slope * tm$within_IQR
segments(tm$change_across_IQR_lo, y, tm$change_across_IQR_hi, y, lwd = 2.2, col = mc[tm$marker]); points(mid, y, pch = 16, col = mc[tm$marker])
pos <- c(Deep = "bulk", Surface = "surface", PlantRoot = "root-adj.")
axis(2, at = y, labels = paste(pos[tm$compartment], LAB[match(tm$marker, MK)], sub("_era5", "", tm$exposure)), las = 1, tick = FALSE, cex.axis = 0.9)
title(main = "Expedition-adjusted change; band = exploratory +/-0.5", adj = 0, cex.main = 0.95, font.main = 1)
dev.off()
cat("figures written to", OUT, "\n")
