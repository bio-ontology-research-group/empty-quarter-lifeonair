#!/usr/bin/env Rscript
# Life on air (P3): Figure 2, proxies for the mechanism of the root-adjacent group 1l loss.
#   a  per site x trip pair: 1l excess loss (root-adjacent vs bulk, beyond controls) against NDVI
#   b  per pair: shift of the community-weighted log doubling time (root-adjacent minus bulk)
#   c  per pair: 1l excess loss against that shift, with the mixed-model line
#   d  per genus: root-adjacent log2 fold change against the fraction of genomes carrying group 1l
# Inputs from h1_mechanism.R. Colours as in fig_h123.R (group 1l = categorical slot 1).
#
#   Rscript lifeonair/code/fig_mechanism.R   (from the repository root)

IN <- "lifeonair/results/h1_mechanism"; OUT <- "lifeonair/figures"
rd <- function(f) read.delim(file.path(IN, f), quote = "", check.names = FALSE)
C1L <- "#2a78d6"; GREY <- "#9a9a94"; INK <- "#2b2b28"; MUTED <- "#6b6b66"; GRID <- "#ececea"

pp <- rd("m_pairs_plantroot_vs_deep.tsv")
m1 <- rd("m1_root_effect_vs_ndvi.tsv"); m2b <- rd("m2b_community_doubling_time_shift.tsv")
m2c <- rd("m2c_excess_loss_vs_copiotrophy_shift.tsv"); gt <- rd("m2d_genus_root_depletion.tsv")
gm <- rd("m2d_genus_models.tsv")

frame <- function(x, y, xlab, ylab, title, sub, ...) {
  plot(x, y, type = "n", bty = "n", xlab = "", ylab = "", las = 1, col.axis = MUTED, cex.axis = 0.8, fg = MUTED, ...)
  abline(h = 0, col = MUTED, lwd = 1)
  mtext(xlab, side = 1, line = 2.3, cex = 0.65, col = MUTED)
  mtext(ylab, side = 2, line = 2.4, cex = 0.65, col = MUTED)
  title(main = title, adj = 0, cex.main = 0.9, font.main = 1, col.main = INK)
  mtext(sub, side = 3, adj = 0, line = -0.2, cex = 0.62, col = MUTED)
}

draw <- function(dev_open) {
  dev_open()
  par(mfrow = c(2, 2), mar = c(4, 4, 2.6, 1), family = "sans")
  # a
  s <- aggregate(cbind(excess_1l, ndvi_W60d) ~ site, pp, median)
  frame(pp$ndvi_W60d, pp$excess_1l, "NDVI, 60 days before sampling", "group 1l excess loss (log2)",
        "a  root contrast against vegetation greenness",
        sprintf("site medians: Spearman rho %.2f, p = %.3f, n = %d sites", m1$estimate[2], m1$p_or_se[2], m1$n[2]))
  points(pp$ndvi_W60d, pp$excess_1l, pch = 16, cex = 0.55, col = adjustcolor(GREY, 0.6))
  points(s$ndvi_W60d, s$excess_1l, pch = 16, cex = 0.9, col = adjustcolor(C1L, 0.85))
  legend("bottomright", bty = "n", cex = 0.7, text.col = INK, pch = 16, col = c(GREY, C1L),
         legend = c("site x trip pair", "site median"))
  # b
  set.seed(1)
  frame(c(0.5, 1.5), range(pp$d_cwm_logd), "", "shift in log doubling time",
        "b  community growth rate, root-adjacent minus bulk",
        sprintf("%d pairs; %.0f%% faster beside roots; site Wilcoxon p = %.1e", m2b$n_pairs[1],
                100 * m2b$frac_pairs_faster[1], m2b$site_wilcox_p[1]), xaxt = "n", xlim = c(0.5, 1.5))
  points(1 + runif(nrow(pp), -0.2, 0.2), pp$d_cwm_logd, pch = 16, cex = 0.6, col = adjustcolor(INK, 0.45))
  segments(0.7, median(pp$d_cwm_logd), 1.3, median(pp$d_cwm_logd), lwd = 2.5, col = INK)
  text(1.32, median(pp$d_cwm_logd), "median", adj = 0, cex = 0.65, col = MUTED)
  text(0.55, min(pp$d_cwm_logd), "faster growers\nbeside roots", adj = c(0, 0), cex = 0.65, col = MUTED)
  # c
  b0 <- m2c$intercept_at_no_shift[2]; b1 <- m2c$estimate[2]
  frame(pp$d_cwm_logd, pp$excess_1l, "shift in log doubling time (root-adjacent minus bulk)",
        "group 1l excess loss (log2)", "c  1l loss against the growth-rate shift",
        sprintf("Spearman rho %.2f, p = %.3f; model value at zero shift %.2f", m2c$estimate[1], m2c$p_or_se[1], b0))
  abline(v = 0, col = GRID, lwd = 1)
  points(pp$d_cwm_logd, pp$excess_1l, pch = 16, cex = 0.6, col = adjustcolor(C1L, 0.6))
  xr <- range(pp$d_cwm_logd); lines(xr, b0 + b1 * xr, lwd = 2, col = INK)
  points(0, b0, pch = 21, bg = "white", col = INK, cex = 1.1, lwd = 1.5)
  # d
  f <- gm[gm$model == "growth + carriage + phylum", ]
  frame(gt$frac_1l, gt$log2fc_root, "fraction of genus representatives carrying group 1l",
        "genus log2 fold change, root-adjacent vs bulk", "d  genus depletion beside roots",
        sprintf("%d genera; 1l carriage %.2f (p = %.0e) given growth rate and phylum", nrow(gt),
                f$Estimate[f$term == "frac_1l"], f$`Pr...t..`[f$term == "frac_1l"]))
  sz <- 0.4 + 1.6 * sqrt(gt$deep_mean_pct / max(gt$deep_mean_pct))
  points(jitter(gt$frac_1l, amount = 0.01), gt$log2fc_root, pch = 16, cex = sz, col = adjustcolor(C1L, 0.45))
  legend("topright", bty = "n", cex = 0.7, text.col = INK, legend = "symbol area: mean share in bulk soil")
  dev.off()
}
draw(function() pdf(file.path(OUT, "fig_mechanism.pdf"), width = 9, height = 7.5))
draw(function() png(file.path(OUT, "fig_mechanism.png"), width = 1800, height = 1500, res = 200))
cat("figure written to", OUT, "\n")
