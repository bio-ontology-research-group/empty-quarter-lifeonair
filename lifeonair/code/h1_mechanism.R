#!/usr/bin/env Rscript
# Life on air (P3): mechanism proxies for H1 (root-adjacent loss of group 1l), 10 Oct 2026.
# Li et al. 2023 (ISME J) attribute the loss of trace-gas oxidisers under desert shrubs to organic carbon
# favouring copiotrophs. Soil organic carbon was never measured in this campaign (P1 has XRF, pH and
# weather only), so two proxies are tested:
#   M1 plant input: NDVI (60-day window, site x trip) against 1l per compartment and against the root effect.
#   M2 copiotrophy: gRodon2 predicted minimum doubling time per representative (g1_growth_rate.sh).
#      a  are 1l carriers slower growers than non-carriers within order?
#      b  is the abundance-weighted community doubling time shorter beside roots (paired site x trip)?
#      c  does the per-pair 1l excess loss scale with that shift?
#      d  genus level: does carrying 1l predict root depletion beyond growth rate?
# Capacity, not activity. Run from the repository root:
#   micromamba run -n eq-lifeonair Rscript lifeonair/code/h1_mechanism.R

suppressMessages(library(lme4))
IN <- "lifeonair/results/inputs"; OUT <- "lifeonair/results/h1_mechanism"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f, ...) read.delim(f, check.names = FALSE, quote = "", ...)
wt <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
sp <- function(x, y) { ok <- is.finite(x) & is.finite(y); r <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  c(rho = unname(r$estimate), p = r$p.value, n = sum(ok)) }
log_lines <- character(); say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

meta <- rd(file.path(IN, "samples.tsv")); meta <- meta[meta$trip != "" & meta$compartment != "", ]
ge <- rd("lifeonair/results/reads/guild_genome_equivalents.tsv"); ge <- ge[ge$id_floor == 60 & ge$scg_genes_detected >= 10, ]
env <- rd(file.path(IN, "env_site_trip.tsv")); site_env <- rd(file.path(IN, "env_site.tsv"))
d <- merge(meta, ge, by.x = "sample_id", by.y = "sample")
d <- merge(d, env, by = c("site", "trip")); d <- merge(d, site_env[, c("site", "map_era5_mm")], by = "site")
c1 <- rd("lifeonair/results/h123/h1_cell_log2fc_plantroot_vs_deep.tsv")

# ---- M1 NDVI --------------------------------------------------------------------------------------
say("M1 NDVI (60-day window) vs group 1l genome equivalents; site NDVI vs MAP rho %.2f",
    sp(aggregate(ndvi_W60d ~ site, d, median)$ndvi_W60d, aggregate(map_era5_mm ~ site, d, median)$map_era5_mm)[1])
m1 <- do.call(rbind, lapply(c("Deep", "PlantRoot", "Surface"), function(cp) {
  x <- d[d$compartment == cp, ]
  m <- lmer(log2(hydB_1l + .01) ~ scale(ndvi_W60d) + scale(log(map_era5_mm)) + (1 | site) + (1 | trip), x)
  cf <- summary(m)$coefficients
  data.frame(compartment = cp, n = nrow(x), ndvi_per_sd = cf[2, 1], ndvi_se = cf[2, 2],
             logmap_per_sd = cf[3, 1], logmap_se = cf[3, 2]) }))
wt(m1, "m1_ndvi_by_compartment.tsv"); print(m1, digits = 3)
cp <- merge(c1, env, by = c("site", "trip"))
mm <- summary(lmer(excess_1l ~ scale(ndvi_W60d) + (1 | site) + (1 | trip), cp))$coefficients
s_cp <- aggregate(cbind(excess_1l, ndvi_W60d) ~ site, cp, median)
m1b <- data.frame(test = c("pair excess_1l ~ NDVI (spearman)", "site-median excess_1l ~ NDVI (spearman)", "pair excess_1l ~ NDVI (lmer, per SD)"),
                  estimate = c(sp(cp$excess_1l, cp$ndvi_W60d)[1], sp(s_cp$excess_1l, s_cp$ndvi_W60d)[1], mm[2, 1]),
                  p_or_se = c(sp(cp$excess_1l, cp$ndvi_W60d)[2], sp(s_cp$excess_1l, s_cp$ndvi_W60d)[2], mm[2, 2]),
                  n = c(nrow(cp), nrow(s_cp), nrow(cp)))
wt(m1b, "m1_root_effect_vs_ndvi.tsv"); print(m1b, digits = 3)

# ---- M2 growth rate -------------------------------------------------------------------------------
gr_file <- file.path(IN, "rep_growth.tsv")
if (!file.exists(gr_file)) { writeLines(log_lines, file.path(OUT, "SUMMARY.txt")); stop("rep_growth.tsv not yet available") }
gr <- rd(gr_file); gr <- gr[is.finite(gr$d_hours) & gr$d_hours > 0, ]
cls <- rd(file.path(IN, "atmo_classification.tsv")); cx <- rd(file.path(IN, "coxL_class.tsv"))
r1l <- unique(cls$rep_id[cls$marker == "hydB" & grepl("Group_1l$", cls$assignment)])
r1h <- unique(cls$rep_id[cls$marker == "hydB" & grepl("Group_1h$", cls$assignment)])
rfI <- unique(sub("\\|.*", "", cx$query[cx$form_I_call == "form_I"]))
tax <- rbind(rd(file.path(IN, "r232_bac120.summary.tsv"))[, 1:2], rd(file.path(IN, "r232_ar53.summary.tsv"))[, 1:2])
rk <- function(cl, p) ifelse(grepl(paste0(p, "__[^;]"), cl), sub(paste0(".*", p, "__([^;]*).*"), "\\1", cl), "unassigned")
gr$phylum <- rk(tax$classification, "p")[match(gr$rep, tax$user_genome)]
gr$order <- paste(gr$phylum, rk(tax$classification, "o")[match(gr$rep, tax$user_genome)], sep = ";")
gr$genus <- paste(gr$order, rk(tax$classification, "g")[match(gr$rep, tax$user_genome)], sep = ";")
gr$c1l <- gr$rep %in% r1l; gr$c1h <- gr$rep %in% r1h; gr$cfI <- gr$rep %in% rfI
gr$logd <- log(gr$d_hours)
say("M2 gRodon estimates for %d representatives (median doubling %.1f h)", nrow(gr), median(gr$d_hours))

# a carriers vs non-carriers, within order
ma <- lmer(logd ~ c1l + c1h + cfI + (1 | order), gr)
fa <- summary(ma)$coefficients
m2a <- data.frame(term = rownames(fa), fa, row.names = NULL,
                  median_d_carriers = c(NA, median(gr$d_hours[gr$c1l]), median(gr$d_hours[gr$c1h]), median(gr$d_hours[gr$cfI])),
                  median_d_noncarriers = c(NA, median(gr$d_hours[!gr$c1l]), median(gr$d_hours[!gr$c1h]), median(gr$d_hours[!gr$cfI])))
wt(m2a, "m2a_doubling_time_by_carrier_status.tsv"); print(m2a, digits = 3)
say("M2a 1l carriers: log doubling time %+.3f within order (SE %.3f); medians %.1f vs %.1f h",
    fa[2, 1], fa[2, 2], m2a$median_d_carriers[2], m2a$median_d_noncarriers[2])

# b community-weighted mean log doubling time per sample (CoverM, reps with an estimate)
ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE, col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund) & ab$abund > 0, ]
ab$logd <- gr$logd[match(ab$rep, gr$rep)]
cov_frac <- tapply((!is.na(ab$logd)) * ab$abund, ab$sample, sum) / tapply(ab$abund, ab$sample, sum)
ab <- ab[!is.na(ab$logd), ]
cwm <- data.frame(sample_id = names(cov_frac), cwm_logd = as.numeric(tapply(ab$logd * ab$abund, ab$sample, sum)[names(cov_frac)] /
                                                     tapply(ab$abund, ab$sample, sum)[names(cov_frac)]),
                  covered = as.numeric(cov_frac))
d2 <- merge(d, cwm, by = "sample_id")
say("M2b CWM coverage of mapped abundance: median %.2f", median(d2$covered))
cells <- function(x, a, b, v) { k <- paste(x$site, x$trip); out <- NULL
  for (z in unique(k)) { p <- x[k == z & x$compartment == a, ]; q <- x[k == z & x$compartment == b, ]
    if (nrow(p) && nrow(q)) out <- rbind(out, data.frame(site = p$site[1], trip = p$trip[1], diff = median(p[[v]]) - median(q[[v]]))) }
  out }
db <- cells(d2, "PlantRoot", "Deep", "cwm_logd"); ds <- cells(d2, "Surface", "Deep", "cwm_logd")
sb <- aggregate(diff ~ site, db, median); ss <- aggregate(diff ~ site, ds, median)
m2b <- data.frame(contrast = c("PlantRoot - Deep", "Surface - Deep"), n_pairs = c(nrow(db), nrow(ds)),
                  median_diff_logd = c(median(db$diff), median(ds$diff)),
                  frac_pairs_faster = c(mean(db$diff < 0), mean(ds$diff < 0)),
                  site_wilcox_p = c(wilcox.test(sb$diff)$p.value, wilcox.test(ss$diff)$p.value))
wt(m2b, "m2b_community_doubling_time_shift.tsv"); print(m2b, digits = 3)

# c per-pair 1l excess loss vs the copiotrophy shift
pc <- merge(c1, setNames(db, c("site", "trip", "d_cwm_logd")), by = c("site", "trip"))
mc <- summary(lmer(excess_1l ~ d_cwm_logd + (1 | site) + (1 | trip), pc))$coefficients
mc2 <- summary(lmer(hydB_1l ~ d_cwm_logd + (1 | site) + (1 | trip), pc))$coefficients
m2c <- data.frame(response = c("excess_1l", "excess_1l", "hydB_1l log2FC"),
                  test = c("spearman", "lmer slope per log-unit", "lmer slope per log-unit"),
                  estimate = c(sp(pc$excess_1l, pc$d_cwm_logd)[1], mc[2, 1], mc2[2, 1]),
                  p_or_se = c(sp(pc$excess_1l, pc$d_cwm_logd)[2], mc[2, 2], mc2[2, 2]), n = nrow(pc),
                  intercept_at_no_shift = c(NA, mc[1, 1], mc2[1, 1]))
wt(m2c, "m2c_excess_loss_vs_copiotrophy_shift.tsv"); print(m2c, digits = 3)
wt(merge(pc[, c("site", "trip", "hydB_1l", "excess_1l", "d_cwm_logd")], env[, c("site", "trip", "ndvi_W60d")],
         by = c("site", "trip"), all.x = TRUE), "m_pairs_plantroot_vs_deep.tsv")   # per-pair table for the figure

# d genus level: root depletion ~ growth rate + 1l carriage
mapped <- tapply(ab$abund, ab$sample, sum)
ab$genus <- gr$genus[match(ab$rep, gr$rep)]
gm <- xtabs(abund ~ sample + genus, ab); gm <- sweep(as.matrix(unclass(gm)), 1, rowSums(gm), "/") * 100
gm <- as.data.frame(gm); gm$sample_id <- rownames(gm); gm <- merge(meta, gm, by = "sample_id")
gen <- setdiff(names(gm), names(meta))
deep_mean <- colMeans(gm[gm$compartment == "Deep", gen]); gen <- gen[deep_mean >= 0.05]
k <- paste(gm$site, gm$trip); keys <- intersect(k[gm$compartment == "PlantRoot"], k[gm$compartment == "Deep"])
lfc <- sapply(gen, function(g) median(sapply(keys, function(z)
  log2((median(gm[[g]][k == z & gm$compartment == "PlantRoot"]) + .01) / (median(gm[[g]][k == z & gm$compartment == "Deep"]) + .01)))))
gt <- data.frame(genus = gen, deep_mean_pct = deep_mean[gen], log2fc_root = lfc,
                 median_logd = tapply(gr$logd, gr$genus, median)[gen],
                 frac_1l = tapply(gr$c1l, gr$genus, mean)[gen], frac_1h = tapply(gr$c1h, gr$genus, mean)[gen],
                 frac_formI = tapply(gr$cfI, gr$genus, mean)[gen], n_reps = as.numeric(table(gr$genus)[gen]))
gt$phylum <- sub(";.*", "", gt$genus)
wt(gt, "m2d_genus_root_depletion.tsv")
fit0 <- lm(log2fc_root ~ median_logd, gt, weights = sqrt(deep_mean_pct))
fit1 <- lm(log2fc_root ~ median_logd + frac_1l + frac_1h + frac_formI, gt, weights = sqrt(deep_mean_pct))
fit2 <- lm(log2fc_root ~ median_logd + frac_1l + frac_1h + frac_formI + phylum, gt, weights = sqrt(deep_mean_pct))
co <- function(f, nm) { s <- summary(f)$coefficients; s <- s[intersect(c("median_logd", "frac_1l", "frac_1h", "frac_formI"), rownames(s)), , drop = FALSE]
  data.frame(model = nm, term = rownames(s), s, r2 = summary(f)$r.squared, n_genera = nrow(gt), row.names = NULL) }
m2d <- rbind(co(fit0, "growth only"), co(fit1, "growth + carriage"), co(fit2, "growth + carriage + phylum"))
wt(m2d, "m2d_genus_models.tsv"); print(m2d, digits = 3)

say("M2b root-adjacent CWM log doubling time shift %+.3f (%.0f%% of %d pairs faster; site p %.2g)",
    m2b$median_diff_logd[1], 100 * m2b$frac_pairs_faster[1], m2b$n_pairs[1], m2b$site_wilcox_p[1])
say("M2c excess_1l vs shift: spearman rho %.2f (p %.2g); lmer slope %.2f (SE %.2f), intercept at no shift %.2f",
    m2c$estimate[1], m2c$p_or_se[1], m2c$estimate[2], m2c$p_or_se[2], m2c$intercept_at_no_shift[2])
f2 <- m2d[m2d$model == "growth + carriage + phylum", ]
say("M2d genus (n=%d): growth slope %.2f (p %.2g); frac_1l %.2f (p %.2g) with phylum",
    nrow(gt), f2$Estimate[f2$term == "median_logd"], f2$Pr...t..[f2$term == "median_logd"],
    f2$Estimate[f2$term == "frac_1l"], f2$Pr...t..[f2$term == "frac_1l"])
writeLines(log_lines, file.path(OUT, "SUMMARY.txt"))
