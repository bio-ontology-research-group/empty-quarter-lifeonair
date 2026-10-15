#!/usr/bin/env Rscript
# Life on air (P3): registered analyses for the three hypotheses (decided 9 October 2026).
#
#   H1  Root-adjacent soil selectively loses group 1l H2 scavengers. PlantRoot and Deep are both
#       5-10 cm (data-paper methods), so PlantRoot/Deep is a same-depth root contrast. The 1l
#       depletion must exceed the depletion of control markers (non-high-affinity [NiFe] large
#       subunits and non-form-I CoxL) sharing the same SCG denominator.
#   H2  H2 and CO scavenging capacity are decoupled across microhabitats: high-affinity hydB varies
#       with compartment, form-I CoxL does not. Primary test uses the raw read ratio
#       (1h + 1l reads) / (form-I coxL reads), which has no shared normaliser.
#   H3  Long-term aridity (ERA5 1995-2024 and NASA POWER 2001-2020 mean annual precipitation) sets the
#       spatial 1l baseline, while within-site, between-trip change does not track antecedent
#       precipitation (space-for-time substitution fails).
#
# Primary measure: read-based genome equivalents (GE) at 60% identity (r3_aggregate.py). GE are
# relative capacity (copies per genome equivalent), not cell fractions or activity.
# Sensitivity variants: identity 50/70; reads to contamination-flagged reps removed; reads to group 1l
# proteins typed below 50% identity or 80% coverage removed (strict_1l_calls); trip 5
# excluded (high 1l median, possible batch); trip 4 excluded (101-nt libraries, truncated to 100 nt).
#
#   Rscript lifeonair/code/h123_analysis.R   (repository root; env lifeonair/code/environment.yml)

suppressMessages({ library(lme4); library(vegan) })
set.seed(20261009)
IN  <- "lifeonair/results/inputs"; RD <- "lifeonair/results/reads"
OUT <- "lifeonair/results/h123"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd  <- function(f, ...) read.delim(f, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE,
                                   check.names = FALSE, ...)
EPS <- 0.01
NB  <- 2000
MARKERS  <- c("hydB_1l", "hydB_1h", "coxL_formI", "hydB_other_nife", "coxL_other")
CONTROLS <- c("hydB_other_nife", "coxL_other")

meta <- rd(file.path(IN, "samples.tsv"))
meta <- meta[meta$trip != "" & meta$compartment != "", ]
site_env <- rd(file.path(IN, "env_site.tsv"))
site_env <- merge(site_env, rd(file.path(IN, "nasa_power_precip_climatology.tsv"))[, c("site", "map_power_mm")], by = "site")
st_env <- rd(file.path(IN, "env_site_trip.tsv"))

load_ge <- function(file, floor, drop_trips = character()) {
  g <- rd(file.path(RD, file)); g <- g[g$id_floor == floor & g$scg_genes_detected >= 10, ]
  d <- merge(meta, g, by.x = "sample_id", by.y = "sample")
  d[!d$trip %in% drop_trips, ]
}
variants <- list(
  primary_id60        = list(file = "guild_genome_equivalents.tsv", floor = 60, drop = character()),
  id50                = list(file = "guild_genome_equivalents.tsv", floor = 50, drop = character()),
  id70                = list(file = "guild_genome_equivalents.tsv", floor = 70, drop = character()),
  no_flagged_reps     = list(file = "guild_genome_equivalents_noflagged.tsv", floor = 60, drop = character()),
  strict_1l_calls     = list(file = "guild_genome_equivalents_strict1l.tsv", floor = 60, drop = character()),
  without_trip5       = list(file = "guild_genome_equivalents.tsv", floor = 60, drop = "trip5"),
  without_trip4       = list(file = "guild_genome_equivalents.tsv", floor = 60, drop = "trip4"))

# cell = site x trip; per-cell compartment medians of each marker
cells <- function(d, a, b, vars) {
  key <- paste(d$site, d$trip, sep = "|")
  out <- NULL
  for (k in unique(key)) {
    x <- d[key == k & d$compartment == a, ]; y <- d[key == k & d$compartment == b, ]
    if (!nrow(x) || !nrow(y)) next
    row <- data.frame(site = x$site[1], trip = x$trip[1])
    for (v in vars) row[[v]] <- log2((median(x[[v]]) + EPS) / (median(y[[v]]) + EPS))
    out <- rbind(out, row)
  }
  out
}
site_boot <- function(df, f) {          # cluster bootstrap over sites
  s <- unique(df$site)
  replicate(NB, { b <- sample(s, replace = TRUE); f(do.call(rbind, lapply(b, function(z) df[df$site == z, ]))) })
}
ci <- function(x) sprintf("[%.2f, %.2f]", quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
site_wilcox <- function(df, v) { s <- tapply(df[[v]], df$site, mean); wilcox.test(s)$p.value }

# ================================================================== H1
h1_core <- function(d) {
  c1 <- cells(d, "PlantRoot", "Deep", MARKERS)
  c1$excess_1l <- c1$hydB_1l - rowMeans(c1[, CONTROLS])
  c1$excess_1h <- c1$hydB_1h - rowMeans(c1[, CONTROLS])
  c1$excess_formI <- c1$coxL_formI - rowMeans(c1[, CONTROLS])
  c1
}
c1 <- h1_core(load_ge(variants$primary_id60$file, 60))
write.table(c1, file.path(OUT, "h1_cell_log2fc_plantroot_vs_deep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
h1_tab <- do.call(rbind, lapply(c(MARKERS, "excess_1l", "excess_1h", "excess_formI"), function(v) {
  b <- site_boot(c1, function(x) median(x[[v]]))
  data.frame(quantity = v, n_cells = nrow(c1), n_sites = length(unique(c1$site)), median_log2fc = median(c1[[v]]),
             boot95 = ci(b), frac_cells_below0 = mean(c1[[v]] < 0), site_level_wilcoxon_p = site_wilcox(c1, v))
}))
h1_trip <- aggregate(cbind(hydB_1l, excess_1l) ~ trip, data = c1, FUN = median)
h1_trip$n_cells <- as.numeric(table(c1$trip)[h1_trip$trip])

# marker x compartment mixed model (long format, PlantRoot + Deep)
d0 <- load_ge(variants$primary_id60$file, 60)
dl <- do.call(rbind, lapply(MARKERS, function(v) {
  x <- d0[d0$compartment %in% c("PlantRoot", "Deep"), c("sample_id", "site", "trip", "compartment")]
  x$marker <- v; x$y <- log2(d0[[v]][d0$compartment %in% c("PlantRoot", "Deep")] + EPS); x }))
dl$marker <- relevel(factor(dl$marker), ref = "hydB_other_nife")
dl$compartment <- relevel(factor(dl$compartment), ref = "Deep")
m_full <- lmer(y ~ marker * compartment + trip + (1 | site) + (1 | site:trip) + (1 | sample_id), data = dl, REML = FALSE)
m_red  <- lmer(y ~ marker + compartment + trip + (1 | site) + (1 | site:trip) + (1 | sample_id), data = dl, REML = FALSE)
a <- anova(m_red, m_full)
fx <- summary(m_full)$coefficients
h1_mm <- data.frame(term = rownames(fx)[grepl(":", rownames(fx))], fx[grepl(":", rownames(fx)), , drop = FALSE], row.names = NULL)
h1_mm$interaction_LRT_chisq <- a$Chisq[2]; h1_mm$interaction_LRT_df <- a$Df[2]; h1_mm$interaction_LRT_p <- a$`Pr(>Chisq)`[2]
write.table(h1_mm, file.path(OUT, "h1_marker_by_compartment_mixed_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# lineages carrying the 1l signal (CoverM, PlantRoot vs Deep at order level)
ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE, col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund) & ab$abund > 0, ]
cls <- rd(file.path(IN, "atmo_classification.tsv"))
reps_1l <- unique(cls$rep_id[cls$marker == "hydB" & grepl("Group_1l$", cls$assignment)])
tax <- rbind(rd(file.path(IN, "r232_bac120.summary.tsv"))[, 1:2], rd(file.path(IN, "r232_ar53.summary.tsv"))[, 1:2])
rk  <- function(cl, p) ifelse(grepl(paste0(p, "__[^;]"), cl), sub(paste0(".*", p, "__([^;]*).*"), "\\1", cl), "unassigned")
tax$phylum <- rk(tax$classification, "p"); tax$order <- rk(tax$classification, "o"); tax$genus <- rk(tax$classification, "g")
a1 <- ab[ab$rep %in% reps_1l & ab$sample %in% d0$sample_id, ]
a1$order <- paste(tax$phylum, tax$order, sep = ";")[match(a1$rep, tax$user_genome)]
mapped <- tapply(ab$abund, ab$sample, sum)
a1$pct <- 100 * a1$abund / as.numeric(mapped[a1$sample])
om <- xtabs(pct ~ sample + order, data = a1); om <- as.data.frame.matrix(unclass(om))
om$sample_id <- rownames(om); om <- merge(meta, om, by = "sample_id")
ords <- setdiff(names(om), names(meta))
top <- names(sort(colSums(om[om$compartment == "Deep", ords]), decreasing = TRUE))[1:10]
c_ord <- cells(om, "PlantRoot", "Deep", top)
h1_lin <- data.frame(order = top, deep_mean_pct_mapped = colMeans(om[om$compartment == "Deep", top]),
                     median_log2fc_plantroot_vs_deep = apply(c_ord[, top], 2, median),
                     frac_cells_below0 = apply(c_ord[, top], 2, function(v) mean(v < 0)), row.names = NULL)
write.table(h1_lin, file.path(OUT, "h1_lineages_1l_carriers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ================================================================== H2
d0$lr_H2_CO <- log2((d0$reads_hydB_1l + d0$reads_hydB_1h + 0.5) / (d0$reads_coxL_formI + 0.5))
d0$lr_1l_CO <- log2((d0$reads_hydB_1l + 0.5) / (d0$reads_coxL_formI + 0.5))
h2_rows <- NULL
for (pr in list(c("PlantRoot", "Deep"), c("Surface", "Deep"), c("PlantRoot", "Surface"))) {
  cc <- cells(d0, pr[1], pr[2], c("hydB_1l", "hydB_1h", "coxL_formI"))
  # raw read ratio difference (log2 of ratio-of-ratios), no shared SCG normaliser
  key <- paste(d0$site, d0$trip)
  rr <- do.call(rbind, lapply(split(d0, key), function(x) {
    a <- x[x$compartment == pr[1], ]; b <- x[x$compartment == pr[2], ]
    if (!nrow(a) || !nrow(b)) return(NULL)
    data.frame(site = a$site[1], trip = a$trip[1], d_lr_H2_CO = median(a$lr_H2_CO) - median(b$lr_H2_CO))
  }))
  for (v in c("hydB_1l", "hydB_1h", "coxL_formI")) {
    b <- site_boot(cc, function(x) median(x[[v]]))
    h2_rows <- rbind(h2_rows, data.frame(contrast = paste(pr, collapse = "/"), quantity = paste0("log2FC_GE_", v),
      n_cells = nrow(cc), median = median(cc[[v]]), boot95 = ci(b), site_level_wilcoxon_p = site_wilcox(cc, v)))
  }
  b <- site_boot(rr, function(x) median(x$d_lr_H2_CO))
  h2_rows <- rbind(h2_rows, data.frame(contrast = paste(pr, collapse = "/"), quantity = "delta_log2_raw_read_ratio_H2HA_to_formI_CO",
    n_cells = nrow(rr), median = median(rr$d_lr_H2_CO), boot95 = ci(b), site_level_wilcoxon_p = site_wilcox(rr, "d_lr_H2_CO")))
}
# per-marker compartment effect size in the same mixed-model frame
h2_var <- do.call(rbind, lapply(c("hydB_1l", "hydB_1h", "coxL_formI", CONTROLS), function(v) {
  x <- d0; x$y <- log2(x[[v]] + EPS)
  f1 <- lmer(y ~ compartment + trip + (1 | site), data = x, REML = FALSE)
  f0 <- lmer(y ~ trip + (1 | site), data = x, REML = FALSE)
  vf <- var(as.vector(model.matrix(f1)[, grep("^compartment", colnames(model.matrix(f1)))] %*%
                      fixef(f1)[grep("^compartment", names(fixef(f1)))]))
  vt <- vf + var(as.vector(model.matrix(f1)[, grep("^trip", colnames(model.matrix(f1)))] %*%
                           fixef(f1)[grep("^trip", names(fixef(f1)))])) + sum(as.data.frame(VarCorr(f1))$vcov)
  data.frame(marker = v, compartment_LRT_p = anova(f0, f1)$`Pr(>Chisq)`[2], compartment_var_share = vf / vt)
}))
# carrier redundancy: effective number of genera (exp Shannon) carrying each capacity per sample
hyd_reps <- unique(cls$rep_id[cls$marker == "hydB" & grepl("Group_1[hl]$", cls$assignment)])
cox_reps <- readLines(file.path(IN, "guild_coxL_reps.txt"))
effn <- function(reps) {
  x <- ab[ab$rep %in% reps & ab$sample %in% d0$sample_id, ]
  x$genus <- tax$genus[match(x$rep, tax$user_genome)]
  g <- tapply(x$abund, list(x$sample, x$genus), sum); g[is.na(g)] <- 0
  setNames(exp(diversity(g, "shannon")), rownames(g))
}
eh <- effn(hyd_reps); ec <- effn(cox_reps)
red <- merge(data.frame(sample_id = names(eh), eff_genera_H2 = eh), data.frame(sample_id = names(ec), eff_genera_CO = ec))
red <- merge(meta[, c("sample_id", "compartment")], red)
h2_red <- aggregate(cbind(eff_genera_H2, eff_genera_CO) ~ compartment, data = red, FUN = median)
write.table(h2_rows, file.path(OUT, "h2_decoupling_contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(h2_var, file.path(OUT, "h2_compartment_effect_by_marker.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(h2_red, file.path(OUT, "h2_carrier_redundancy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ================================================================== H3
prank <- function(x, y, Z) {   # partial Spearman: correlate rank residuals after regressing on Z ranks
  rx <- rank(x); ry <- rank(y); Zr <- apply(as.matrix(Z), 2, rank)
  ex <- resid(lm(rx ~ Zr)); ey <- resid(lm(ry ~ Zr)); ct <- cor.test(ex, ey)
  c(rho = unname(ct$estimate), p = ct$p.value)
}
site_env$block <- kmeans(scale(site_env[, c("lat", "lon")]), centers = 6, nstart = 50)$cluster
# community-composition axis (Deep samples, site means; PCoA axis 1 of Bray-Curtis)
deep_ids <- meta$sample_id[meta$compartment == "Deep"]
cm <- xtabs(abund ~ sample + rep, data = ab[ab$sample %in% deep_ids, ]); cm <- unclass(cm); cm <- cm / rowSums(cm)
cm_site <- rowsum(cm, meta$site[match(rownames(cm), meta$sample_id)]) /
           as.numeric(table(meta$site[match(rownames(cm), meta$sample_id)])[rownames(rowsum(cm, meta$site[match(rownames(cm), meta$sample_id)]))])
pc1 <- cmdscale(vegdist(cm_site, "bray"), k = 1)[, 1]
site_env$comp_pc1_deep <- pc1[as.character(site_env$site)]

h3_rows <- NULL; h3_lobo <- NULL
for (comp in c("Deep", "Surface", "PlantRoot")) for (v in c("hydB_1l", "hydB_1h", "coxL_formI")) {
  x <- d0[d0$compartment == comp, ]
  s <- aggregate(x[[v]], list(site = x$site), median); names(s)[2] <- "ge"
  s <- merge(s, site_env, by = "site")
  for (pv in c("map_era5_mm", "map_power_mm")) {
    r0 <- cor.test(s$ge, s[[pv]], method = "spearman", exact = FALSE)
    pg <- prank(s$ge, s[[pv]], s[, c("lat", "lon", "elevation_m")])
    sc <- s[!is.na(s$comp_pc1_deep), ]
    pc <- prank(sc$ge, sc[[pv]], sc[, c("comp_pc1_deep"), drop = FALSE])
    bb <- replicate(NB, { b <- sample(unique(s$block), replace = TRUE)
                          z <- do.call(rbind, lapply(b, function(k) s[s$block == k, ]))
                          suppressWarnings(cor(z$ge, z[[pv]], method = "spearman")) })
    lobo <- sapply(unique(s$block), function(k) { z <- s[s$block != k, ]; cor(z$ge, z[[pv]], method = "spearman") })
    h3_rows <- rbind(h3_rows, data.frame(compartment = comp, marker = v, precip = pv, n_sites = nrow(s),
      rho = unname(r0$estimate), p = r0$p.value, block_boot95 = ci(bb),
      leave_one_block_out_range = sprintf("%.2f to %.2f", min(lobo), max(lobo)),
      partial_rho_latlonelev = pg[["rho"]], partial_p_latlonelev = pg[["p"]],
      partial_rho_compPC1 = pc[["rho"]], partial_p_compPC1 = pc[["p"]]))
  }
}
# mixed model with log MAP (ERA5) across compartments
d3 <- merge(d0, site_env[, c("site", "map_era5_mm", "map_power_mm")], by = "site")
h3_mm <- do.call(rbind, lapply(c("hydB_1l", "hydB_1h", "coxL_formI"), function(v) {
  d3$y <- log2(d3[[v]] + EPS); d3$lmap <- as.numeric(scale(log(d3$map_era5_mm)))
  f1 <- lmer(y ~ compartment + trip + lmap + (1 | site), data = d3, REML = FALSE)
  f0 <- lmer(y ~ compartment + trip + (1 | site), data = d3, REML = FALSE)
  data.frame(marker = v, beta_per_sd_logMAP = fixef(f1)[["lmap"]], LRT_p = anova(f0, f1)$`Pr(>Chisq)`[2])
}))
# temporal: within-site, between-trip change vs antecedent precipitation (ERA5 record ends 2024)
h3_t <- NULL
for (comp in c("Deep", "Surface", "PlantRoot")) for (v in c("hydB_1l", "hydB_1h", "coxL_formI")) {
  x <- d0[d0$compartment == comp, ]
  cc <- aggregate(log2(x[[v]] + EPS), list(site = x$site, trip = x$trip), median); names(cc)[3] <- "y"
  cc <- merge(cc, st_env, by = c("site", "trip"))
  for (pv in c("P90_era5", "P365_era5")) {
    z <- cc[!is.na(cc[[pv]]), ]; z$lp <- log(z[[pv]] + 1)
    z <- z[z$site %in% names(which(table(z$site) >= 2)), ]
    z$y_d  <- z$y  - ave(z$y,  z$site); z$lp_d <- z$lp - ave(z$lp, z$site)
    tr <- cor.test(z$y_d, z$lp_d, method = "spearman", exact = FALSE)
    slope_t <- coef(lm(y_d ~ lp_d, data = z))[2]
    bt <- replicate(NB, { b <- sample(unique(z$site), replace = TRUE)
                          w <- do.call(rbind, lapply(b, function(k) z[z$site == k, ])); coef(lm(y_d ~ lp_d, data = w))[2] })
    sm <- aggregate(cbind(y, lp) ~ site, data = z, FUN = mean)
    sm <- merge(sm, site_env[, c("site", "map_era5_mm")]); slope_s <- coef(lm(y ~ log(map_era5_mm), data = sm))[2]
    h3_t <- rbind(h3_t, data.frame(compartment = comp, marker = v, antecedent = pv, n_cells = nrow(z),
      n_sites = length(unique(z$site)), within_site_rho = unname(tr$estimate), within_site_p = tr$p.value,
      temporal_slope = slope_t, temporal_slope_boot95 = ci(bt),
      min_detectable_slope_80pct = 2.8 * sd(bt), spatial_slope_logMAP = slope_s))
  }
}
# temporal, adjusted for trip (antecedent precipitation varies mostly between trips, so a
# within-site slope can be aliased with trip / batch): y ~ log(P365 + 1) + trip + (1 | site)
h3_tfe <- NULL
for (comp in c("Deep", "Surface", "PlantRoot")) for (v in c("hydB_1l", "hydB_1h", "coxL_formI")) {
  x <- merge(d0[d0$compartment == comp, ], st_env, by = c("site", "trip"))
  x <- x[!is.na(x$P365_era5), ]; x$y <- log2(x[[v]] + EPS); x$lp <- log(x$P365_era5 + 1)
  f1 <- lmer(y ~ lp + trip + (1 | site), data = x, REML = FALSE); f0 <- lmer(y ~ trip + (1 | site), data = x, REML = FALSE)
  g1 <- lmer(y ~ lp + (1 | site), data = x, REML = FALSE);        g0 <- lmer(y ~ (1 | site), data = x, REML = FALSE)
  h3_tfe <- rbind(h3_tfe, data.frame(compartment = comp, marker = v, n_samples = nrow(x),
    beta_logP365_no_trip = fixef(g1)[["lp"]], p_no_trip = anova(g0, g1)$`Pr(>Chisq)`[2],
    beta_logP365_trip_adjusted = fixef(f1)[["lp"]], p_trip_adjusted = anova(f0, f1)$`Pr(>Chisq)`[2]))
}
write.table(h3_tfe, file.path(OUT, "h3_temporal_trip_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(h3_rows, file.path(OUT, "h3_spatial_aridity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(h3_mm, file.path(OUT, "h3_mixed_model_logMAP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(h3_t, file.path(OUT, "h3_temporal_space_for_time.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(site_env, file.path(OUT, "h3_site_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ================================================================== sensitivity
sens <- do.call(rbind, lapply(names(variants), function(nm) {
  v <- variants[[nm]]; d <- load_ge(v$file, v$floor, v$drop)
  c1v <- h1_core(d)
  s <- aggregate(d$hydB_1l[d$compartment == "Deep"], list(site = d$site[d$compartment == "Deep"]), median)
  s <- merge(s, site_env, by = "site")
  dd <- d; dd$lr <- log2((dd$reads_hydB_1l + dd$reads_hydB_1h + 0.5) / (dd$reads_coxL_formI + 0.5))
  rr <- cells(transform(dd, lr = lr + 0), "PlantRoot", "Deep", "lr")
  data.frame(variant = nm, n_samples = nrow(d),
             H1_1l_log2fc = median(c1v$hydB_1l), H1_excess_1l = median(c1v$excess_1l),
             H1_site_p = site_wilcox(c1v, "excess_1l"),
             H2_formI_log2fc = median(c1v$coxL_formI),
             H3_deep_1l_rho_ERA5 = cor(s$x, s$map_era5_mm, method = "spearman"),
             H3_deep_1l_rho_POWER = cor(s$x, s$map_power_mm, method = "spearman"))
}))
write.table(sens, file.path(OUT, "sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ================================================================== summary
sink(file.path(OUT, "SUMMARY.txt"))
cat("Life on air (P3) H1-H3, primary measure: read-based GE at 60% identity\n\n== H1 PlantRoot vs Deep (same 5-10 cm depth), per site x trip ==\n")
print(h1_tab, row.names = FALSE, digits = 3)
cat("\nper trip:\n"); print(h1_trip, row.names = FALSE, digits = 3)
cat(sprintf("\nmarker x compartment interaction LRT: chisq %.1f, df %d, p %.2g\n", a$Chisq[2], a$Df[2], a$`Pr(>Chisq)`[2]))
print(h1_mm[, 1:3], row.names = FALSE, digits = 3)
cat("\n1l-carrier orders (CoverM, % of mapped):\n"); print(h1_lin, row.names = FALSE, digits = 3)
cat("\n== H2 decoupling ==\n"); print(h2_rows, row.names = FALSE, digits = 3)
cat("\n"); print(h2_var, row.names = FALSE, digits = 3)
cat("\neffective number of carrier genera (median per sample):\n"); print(h2_red, row.names = FALSE, digits = 3)
cat("\n== H3 spatial aridity (site medians) ==\n"); print(h3_rows, row.names = FALSE, digits = 3)
cat("\n"); print(h3_mm, row.names = FALSE, digits = 3)
cat("\n== H3 temporal (within-site, between-trip) ==\n"); print(h3_t, row.names = FALSE, digits = 3)
cat("\ntrip-adjusted temporal models:\n"); print(h3_tfe, row.names = FALSE, digits = 3)
cat("\n== sensitivity ==\n"); print(sens, row.names = FALSE, digits = 3)
cat("\n"); print(sessionInfo()$otherPkgs |> names())
sink()
cat(readLines(file.path(OUT, "SUMMARY.txt")), sep = "\n")
