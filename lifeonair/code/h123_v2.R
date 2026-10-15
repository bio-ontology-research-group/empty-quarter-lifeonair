#!/usr/bin/env Rscript
# Life on air (P3), revised analyses after review round 1 (12 October 2026). Supersedes h123_analysis.R,
# which is kept for provenance (outputs in results/h123/).
#
# Corrections (finding IDs from lifeonair/reviews/round1/findings.tsv):
#   R1-01  corrected CoxL form calls (fix_coxL_v2.py; guild_genome_equivalents_v2.tsv)
#   R1-02  group 2a [NiFe] reported separately, not part of the hydrogenase control
#   R1-19  reads to proteins of two refinement-dropped _sub representatives excluded
#   R1-06  canonical transect only (sites 1-60); 62 and 64 are special-interest locations
#   R1-03  pure ERA5 1995-2025 (fetch_era5_2025.py), so antecedent rain covers expedition 5
#   R1-11  temporal test as an explicit within/between-site decomposition with an equivalence bound
#   R1-12  site-visit random effects wherever samples share a visit
#   R1-13  residual spatial autocorrelation and block-scale sensitivity
#   R1-15/16 PERMANOVA with restricted permutations; synchrony within position
#   R1-18  Benjamini-Hochberg within declared families of secondary tests
#   R1-21  pseudocount sensitivity for paired contrasts
#   R1-31  marginal R2 (Nakagawa) instead of the ad hoc variance share
# New analyses (round-1 consultation, lifeonair/reviews/round1_consult/):
#   A1  scaling of each marker with DNA yield (biomass proxy)
#   A2  lineage turnover vs within-genus carrier fraction (root contrast and aridity)
#   A3  root vs surface filters across markers; surface contrast vs temperature and soil moisture
#   A4  carbon-fixation (CBB) co-carriage by hydrogenase group
#   A6  carrier representation vs copies per carrier
#   A7  geochemistry (XRF) as an alternative to aridity
# Round-2 corrections (lifeonair/reviews/round2/findings.tsv):
#   R2-01  topography had two records for site 19; site-level covariates are averaged per site and every
#          site-level join is checked for unique keys
#   R2-09  contamination sensitivity uses the union of genus-based and batch-frequency flags
#   R2-12  aridity decomposition zero-fills absent genus x visit combinations before averaging
#   R2-14  DNA sensitivities keep kit adjustment and visit random effects
#   R2-15/16/17  temporal models: NASA POWER as an independent exposure, leave-one-expedition-out,
#          ERA5-grid-cell bootstrap, pseudocount sensitivity
#   R2-18  PERMANOVA on visit-level profiles with visits exchangeable only within site
#   R2-19/20/21/28  bootstrap of the turnover share, clipping sensitivity, CBB completeness sensitivity,
#          surface temperature residual support
#
# GE (genome equivalents) are relative capacity: marker copies per ribosomal-protein genome equivalent.
#
#   micromamba run -n eq-lifeonair Rscript lifeonair/code/h123_v2.R   (repository root)

suppressMessages({ library(lme4); library(vegan); library(ape); library(permute) })
set.seed(20261012)
IN  <- "lifeonair/results/inputs"; RD <- "lifeonair/results/reads"
OUT <- "lifeonair/results/h123_v2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd  <- function(f, ...) read.delim(f, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE,
                                   check.names = FALSE, ...)
wr  <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)
EPS <- 0.01; NB <- 2000
FOCAL    <- c("hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI")
CONTROLS <- c("hydB_other_nife", "coxL_other")
MARKERS  <- c(FOCAL, CONTROLS)
ORPHANS  <- readLines(file.path(RD, "orphan_sub_reps.txt"))

meta <- rd(file.path(IN, "samples.tsv"))
meta <- meta[meta$trip != "" & meta$compartment != "" & meta$site %in% 1:60, ]
site_env <- rd(file.path(IN, "env_site.tsv"))
site_env <- merge(site_env, rd(file.path(IN, "nasa_power_precip_climatology.tsv"))[, c("site", "map_power_mm")], by = "site")
st_env <- rd(file.path(IN, "env_site_trip.tsv"))
tmp <- rd(file.path(IN, "per_trip_site_temperature.tsv")); tmp$trip <- paste0("trip", tmp$trip)
st_env <- merge(st_env, tmp, by = c("site", "trip"), all.x = TRUE)

load_ge <- function(file = "guild_genome_equivalents_v2.tsv", floor = 60, drop_trips = character()) {
  g <- rd(file.path(RD, file)); g <- g[g$id_floor == floor & g$scg_genes_detected >= 10, ]
  d <- merge(meta, g, by.x = "sample_id", by.y = "sample")
  d[!d$trip %in% drop_trips, ]
}
d0 <- load_ge()
d0$visit <- paste(d0$site, d0$trip, sep = "|")

# DNA yield (extract concentration), biomass proxy; conflicting duplicate records and non-positive values dropped
led <- rd(file.path(IN, "sample_ledger_dna.tsv"))
led$dna <- suppressWarnings(as.numeric(led$dna_concentration))
nconc <- tapply(led$dna, led$sample_id, function(x) length(unique(round(x[is.finite(x)], 6))))
led <- led[led$sample_id %in% names(nconc)[nconc == 1] & is.finite(led$dna) & led$dna > 0, ]
led <- led[!duplicated(led$sample_id), c("sample_id", "dna", "dna_kit")]
d0 <- merge(d0, led, by = "sample_id", all.x = TRUE)

# ------------------------------------------------------------------ helpers
cells <- function(d, a, b, vars, eps = EPS) {          # visit-level log2 ratio of replicate medians
  key <- paste(d$site, d$trip, sep = "|"); out <- NULL
  for (k in unique(key)) {
    x <- d[key == k & d$compartment == a, ]; y <- d[key == k & d$compartment == b, ]
    if (!nrow(x) || !nrow(y)) next
    row <- data.frame(site = x$site[1], trip = x$trip[1])
    for (v in vars) row[[v]] <- log2((median(x[[v]]) + eps) / (median(y[[v]]) + eps))
    out <- rbind(out, row)
  }
  out
}
add_excess <- function(c1) {
  ctl <- rowMeans(c1[, CONTROLS])
  for (v in FOCAL) c1[[paste0("excess_", v)]] <- c1[[v]] - ctl
  c1
}
site_boot <- function(df, f) {           # cluster bootstrap over sites (row indices, fast)
  idx <- split(seq_len(nrow(df)), df$site); s <- names(idx)
  replicate(NB, f(df[unlist(idx[sample(s, replace = TRUE)], use.names = FALSE), , drop = FALSE]))
}
q95 <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
site_wilcox <- function(df, v) { s <- tapply(df[[v]], df$site, mean); wilcox.test(s, exact = FALSE)$p.value }
contrast_table <- function(cc, vars, label) do.call(rbind, lapply(vars, function(v) {
  b <- q95(site_boot(cc, function(x) median(x[[v]])))
  data.frame(contrast = label, quantity = v, n_visits = nrow(cc), n_sites = length(unique(cc$site)),
             median_log2 = median(cc[[v]]), ci_lo = b[1], ci_hi = b[2], frac_below0 = mean(cc[[v]] < 0),
             p_site_wilcoxon = site_wilcox(cc, v), row.names = NULL)
}))
r2m <- function(f) {    # Nakagawa marginal R2 for a Gaussian lmer
  vf <- var(as.vector(model.matrix(f) %*% fixef(f))); vr <- sum(as.data.frame(VarCorr(f))$vcov)
  vf / (vf + vr)
}

# ================================================================== H1 root contrast (relative representation)
c1 <- add_excess(cells(d0, "PlantRoot", "Deep", MARKERS))
wr(c1, "h1_visit_log2fc_root_vs_bulk.tsv")
h1 <- contrast_table(c1, c(MARKERS, paste0("excess_", FOCAL)), "root/bulk")
cs <- add_excess(cells(d0, "Surface", "Deep", MARKERS))
wr(cs, "h1_visit_log2fc_surface_vs_bulk.tsv")
hs <- contrast_table(cs, c(MARKERS, paste0("excess_", FOCAL)), "surface/bulk")
crs <- add_excess(cells(d0, "PlantRoot", "Surface", MARKERS))
hrs <- contrast_table(crs, c(MARKERS, paste0("excess_", FOCAL)), "root/surface")
contr <- rbind(h1, hs, hrs)
contr$q_BH_within_contrast <- ave(contr$p_site_wilcoxon, contr$contrast, FUN = function(p) p.adjust(p, "BH"))
wr(contr, "contrasts_by_position.tsv")

# pseudocount sensitivity (no paired visit median is zero for 1l; ratios are finite without a pseudocount)
pc_sens <- do.call(rbind, lapply(c(0, 0.001, 0.01, 0.1), function(e) {
  x <- add_excess(cells(d0, "PlantRoot", "Deep", MARKERS, eps = e))
  fm <- function(v) { z <- x[[v]][is.finite(x[[v]])]; median(z) }
  data.frame(pseudocount = e, n_visits = nrow(x), n_finite_1l = sum(is.finite(x$hydB_1l)),
             n_finite_excess_1l = sum(is.finite(x$excess_hydB_1l)), root_1l = fm("hydB_1l"),
             root_excess_1l = fm("excess_hydB_1l"), root_1h = fm("hydB_1h"), root_2a = fm("hydB_2a"),
             root_formI = fm("coxL_formI"))
}))
wr(pc_sens, "h1_pseudocount_sensitivity.tsv")

h1_trip <- do.call(rbind, lapply(split(c1, c1$trip), function(x)
  data.frame(trip = x$trip[1], n_visits = nrow(x), median_1l = median(x$hydB_1l), median_excess_1l = median(x$excess_hydB_1l))))
wr(h1_trip, "h1_by_trip.tsv")

# marker x position joint model (root and bulk; all six markers; sample effect for within-sample correlation)
long <- function(d, comps) do.call(rbind, lapply(MARKERS, function(v) {
  x <- d[d$compartment %in% comps, c("sample_id", "site", "trip", "visit", "compartment")]
  x$marker <- v; x$y <- log2(d[[v]][d$compartment %in% comps] + EPS); x }))
joint <- function(comps, ref_comp = "Deep") {
  dl <- long(d0, comps)
  dl$marker <- relevel(factor(dl$marker), ref = "hydB_other_nife")
  dl$compartment <- relevel(factor(dl$compartment), ref = ref_comp)
  f1 <- lmer(y ~ marker * compartment + trip + (1 | site) + (1 | visit) + (1 | sample_id), data = dl, REML = FALSE)
  f0 <- lmer(y ~ marker + compartment + trip + (1 | site) + (1 | visit) + (1 | sample_id), data = dl, REML = FALSE)
  a <- anova(f0, f1); fx <- summary(f1)$coefficients; k <- grepl(":", rownames(fx))
  data.frame(positions = paste(comps, collapse = "+"), term = rownames(fx)[k], fx[k, 1:2, drop = FALSE],
             LRT_chisq = a$Chisq[2], LRT_df = a$Df[2], LRT_p = a$`Pr(>Chisq)`[2], row.names = NULL)
}
jm <- rbind(joint(c("PlantRoot", "Deep")), joint(c("Surface", "Deep")))
wr(jm, "joint_marker_by_position_models.tsv")

# per-marker position effect: marginal R2 and LRT, with visit random effect (R1-12, R1-31)
pos_eff <- do.call(rbind, lapply(MARKERS, function(v) {
  x <- d0; x$y <- log2(x[[v]] + EPS)
  f1 <- lmer(y ~ compartment + trip + (1 | site) + (1 | visit), data = x, REML = FALSE)
  f0 <- lmer(y ~ trip + (1 | site) + (1 | visit), data = x, REML = FALSE)
  data.frame(marker = v, position_LRT_p = anova(f0, f1)$`Pr(>Chisq)`[2],
             marginal_R2_full = r2m(f1), marginal_R2_trip_only = r2m(f0), delta_marginal_R2_position = r2m(f1) - r2m(f0))
}))
wr(pos_eff, "position_effect_by_marker.tsv")

# ================================================================== A1 scaling with DNA yield
dd <- d0[is.finite(d0$dna), ]
dd$ldna <- log2(dd$dna)
dd$ldna_visit <- ave(dd$ldna, dd$visit); dd$ldna_within <- dd$ldna - dd$ldna_visit
kit_ok <- names(which(tapply(dd$dna_kit, dd$visit, function(k) length(unique(k)) == 1)))
a1 <- NULL
for (v in c(MARKERS, "r_1l_1h", "r_1l_2a", "r_1l_formI")) {
  y <- switch(v, r_1l_1h = log2(dd$hydB_1l / dd$hydB_1h), r_1l_2a = log2(dd$hydB_1l / dd$hydB_2a),
              r_1l_formI = log2(dd$hydB_1l / dd$coxL_formI), log2(dd[[v]] + EPS))
  x <- dd; x$y <- y; x <- x[is.finite(x$y), ]
  f <- lmer(y ~ ldna + compartment + trip + dna_kit + (1 | site) + (1 | visit), data = x, REML = FALSE)
  g <- lmer(y ~ ldna_within + ldna_visit + compartment + trip + dna_kit + (1 | site) + (1 | visit), data = x, REML = FALSE)
  xs <- x[x$visit %in% kit_ok, ]
  h <- lmer(y ~ ldna + compartment + trip + dna_kit + (1 | site) + (1 | visit), data = xs, REML = FALSE)
  xb <- x[x$compartment == "Deep", ]
  b <- lmer(y ~ ldna + trip + dna_kit + (1 | site) + (1 | visit), data = xb, REML = FALSE)
  cf <- function(m, t) summary(m)$coefficients[t, 1:2]
  a1 <- rbind(a1, data.frame(response = v, n = nrow(x), slope_ldna = cf(f, "ldna")[1], se = cf(f, "ldna")[2],
    slope_within_visit = cf(g, "ldna_within")[1], se_within = cf(g, "ldna_within")[2],
    slope_between_visit = cf(g, "ldna_visit")[1], se_between = cf(g, "ldna_visit")[2],
    slope_same_kit_visits = cf(h, "ldna")[1], se_same_kit = cf(h, "ldna")[2], n_same_kit = nrow(xs),
    slope_bulk_only = cf(b, "ldna")[1], se_bulk_only = cf(b, "ldna")[2],
    # errors-in-variables: classical error with reliability lambda attenuates slope by lambda
    slope_eiv_rel0.8 = cf(f, "ldna")[1] / 0.8, slope_eiv_rel0.6 = cf(f, "ldna")[1] / 0.6))
}
a1$implied_alpha <- 1 + a1$slope_ldna        # carriers ~ biomass^alpha (for single markers)
wr(a1, "A1_dna_scaling.tsv")
# root effect on 1l before and after DNA adjustment, same complete-case cohort (root + bulk)
rb <- dd[dd$compartment %in% c("PlantRoot", "Deep"), ]; rb$y <- log2(rb$hydB_1l + EPS)
m0 <- lmer(y ~ compartment + trip + dna_kit + (1 | site) + (1 | visit), data = rb, REML = FALSE)
m1 <- lmer(y ~ compartment + ldna + trip + dna_kit + (1 | site) + (1 | visit), data = rb, REML = FALSE)
a1b <- data.frame(model = c("root_effect_unadjusted", "root_effect_DNA_adjusted"),
                  n = nrow(rb), estimate = c(fixef(m0)["compartmentPlantRoot"], fixef(m1)["compartmentPlantRoot"]),
                  se = c(summary(m0)$coefficients["compartmentPlantRoot", 2], summary(m1)$coefficients["compartmentPlantRoot", 2]))
# per-position DNA ratio and DNA-scaled (per-gram proxy) contrasts; product built from GE x DNA (descriptive only)
dvis <- cells(transform(dd, dna_ = dna), "PlantRoot", "Deep", "dna_", eps = 0)
a1b <- rbind(a1b, data.frame(model = "root_bulk_log2_DNA_ratio_median", n = nrow(dvis),
                             estimate = median(dvis$dna_), se = NA))
wr(a1b, "A1_root_effect_dna_adjusted.tsv")
# aridity and DNA: site medians in bulk soil
bs <- aggregate(cbind(ldna, l1 = log2(hydB_1l + EPS)) ~ site, dd[dd$compartment == "Deep", ], median)
bs <- merge(bs, site_env, by = "site")
a1c <- rbind(
  data.frame(test = "bulk site log DNA vs log MAP (ERA5)", n_sites = nrow(bs),
             rho = cor(bs$ldna, log(bs$map_era5_mm), method = "spearman"),
             p = cor.test(bs$ldna, log(bs$map_era5_mm), method = "spearman", exact = FALSE)$p.value),
  data.frame(test = "bulk site 1l vs log DNA", n_sites = nrow(bs),
             rho = cor(bs$l1, bs$ldna, method = "spearman"),
             p = cor.test(bs$l1, bs$ldna, method = "spearman", exact = FALSE)$p.value))
x <- dd[dd$compartment == "Deep", ]; x$y <- log2(x$hydB_1l + EPS); x$lmap <- as.numeric(scale(log(site_env$map_era5_mm[match(x$site, site_env$site)])))
f0 <- lmer(y ~ lmap + trip + dna_kit + (1 | site) + (1 | visit), data = x, REML = FALSE)
f1 <- lmer(y ~ lmap + ldna + trip + dna_kit + (1 | site) + (1 | visit), data = x, REML = FALSE)
a1c <- rbind(a1c, data.frame(test = c("bulk 1l per SD log MAP, DNA cohort", "bulk 1l per SD log MAP, + log DNA"),
                             n_sites = length(unique(x$site)), rho = c(fixef(f0)["lmap"], fixef(f1)["lmap"]),
                             p = c(summary(f0)$coefficients["lmap", 2], summary(f1)$coefficients["lmap", 2])))
names(a1c)[3:4] <- c("estimate", "p_or_se")
wr(a1c, "A1_aridity_and_dna.tsv")

# ================================================================== CoverM-based lineage analyses (A2, A4, A6)
ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE, col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund) & ab$abund > 0 & ab$sample %in% d0$sample_id & !ab$rep %in% ORPHANS, ]
mapped <- tapply(ab$abund, ab$sample, sum); ab$share <- ab$abund / as.numeric(mapped[ab$sample])
cls <- rd(file.path(IN, "atmo_classification_v2.tsv")); cls <- cls[!cls$rep_id %in% ORPHANS, ]
grp_reps <- function(g) unique(cls$rep_id[cls$marker == "hydB" & grepl(paste0("Group_", g, "$"), cls$assignment)])
R1L <- grp_reps("1l"); R1H <- grp_reps("1h"); R2A <- grp_reps("2a")
RFI <- setdiff(readLines(file.path(IN, "guild_coxL_reps_v2.txt")), ORPHANS)
tax <- rbind(rd(file.path(IN, "r232_bac120.summary.tsv"))[, 1:2], rd(file.path(IN, "r232_ar53.summary.tsv"))[, 1:2])
rk  <- function(cl, p) ifelse(grepl(paste0(p, "__[^;]"), cl), sub(paste0(".*", p, "__([^;]*).*"), "\\1", cl), NA)
tax$phylum <- rk(tax$classification, "p"); tax$order <- rk(tax$classification, "o"); tax$genus <- rk(tax$classification, "g")
ab$genus <- tax$genus[match(ab$rep, tax$user_genome)]; ab$genus[is.na(ab$genus)] <- paste0("unclassified_", ab$rep[is.na(ab$genus)])
ab$c1l <- ab$rep %in% R1L
sm <- merge(ab, meta[, c("sample_id", "site", "trip", "compartment")], by.x = "sample", by.y = "sample_id")
sm$visit <- paste(sm$site, sm$trip, sep = "|")

# A2 (a) decomposition of the root/bulk change in 1l-carrier share into genus-abundance and within-genus terms
gv <- aggregate(cbind(share, s1l = share * c1l) ~ visit + compartment + genus, data = sm, FUN = sum)
# per visit and position: mean over replicate samples of genus share and genus 1l share
nrep <- aggregate(sample ~ visit + compartment, data = unique(sm[, c("sample", "visit", "compartment")]), FUN = length)
gv <- merge(gv, nrep, by = c("visit", "compartment")); gv$w <- gv$share / gv$sample; gv$w1 <- gv$s1l / gv$sample
dec <- NULL
for (v in unique(gv$visit)) {
  r <- gv[gv$visit == v & gv$compartment == "PlantRoot", ]; b <- gv[gv$visit == v & gv$compartment == "Deep", ]
  if (!nrow(r) || !nrow(b)) next
  g <- merge(r[, c("genus", "w", "w1")], b[, c("genus", "w", "w1")], by = "genus", all = TRUE, suffixes = c("_r", "_b"))
  g[is.na(g)] <- 0
  q_r <- ifelse(g$w_r > 0, g$w1_r / g$w_r, NA); q_b <- ifelse(g$w_b > 0, g$w1_b / g$w_b, NA)
  q_r[is.na(q_r)] <- q_b[is.na(q_r)]; q_b[is.na(q_b)] <- q_r[is.na(q_b)]; q_r[is.na(q_r)] <- 0; q_b[is.na(q_b)] <- 0
  tot_r <- sum(g$w1_r); tot_b <- sum(g$w1_b)
  # Bennet decomposition: delta = sum(dw * qbar) + sum(wbar * dq)
  dec <- rbind(dec, data.frame(visit = v, share1l_root = tot_r, share1l_bulk = tot_b,
    delta = tot_r - tot_b, genus_turnover = sum((g$w_r - g$w_b) * (q_r + q_b) / 2),
    within_genus = sum((g$w_r + g$w_b) / 2 * (q_r - q_b))))
}
dec$site <- sub("\\|.*", "", dec$visit)
a2a <- data.frame(n_visits = nrow(dec), median_share1l_bulk = median(dec$share1l_bulk), median_share1l_root = median(dec$share1l_root),
                  median_delta = median(dec$delta), median_genus_turnover = median(dec$genus_turnover),
                  median_within_genus = median(dec$within_genus),
                  frac_delta_from_turnover = sum(dec$genus_turnover) / sum(dec$delta),
                  frac_turnover_site_boot95 = { b <- q95(site_boot(dec, function(x) sum(x$genus_turnover) / sum(x$delta)))
                                                sprintf("%.3f to %.3f", b[1], b[2]) },
                  p_within_genus_site = { s <- tapply(dec$within_genus, dec$site, mean); wilcox.test(s, exact = FALSE)$p.value })
wr(dec, "A2_root_bulk_decomposition_by_visit.tsv"); wr(a2a, "A2_root_bulk_decomposition_summary.tsv")

# A2 (b) aridity: across bulk sites, genus-weight term vs within-genus term for the 1l-carrier share
bulk <- gv[gv$compartment == "Deep", ]; bulk$site <- as.integer(sub("\\|.*", "", bulk$visit))
# mean over ALL bulk visits of a site: genera absent from a visit contribute zero (R2-12)
nvis <- tapply(bulk$visit, bulk$site, function(v) length(unique(v)))
bsg <- aggregate(cbind(w, w1) ~ site + genus, data = bulk, FUN = sum)
bsg$w <- bsg$w / as.numeric(nvis[as.character(bsg$site)]); bsg$w1 <- bsg$w1 / as.numeric(nvis[as.character(bsg$site)])
stopifnot(all(abs(tapply(bsg$w, bsg$site, sum) - 1) < 1e-6))
qbar <- aggregate(cbind(w, w1) ~ genus, data = bsg, FUN = mean); qbar$q <- ifelse(qbar$w > 0, qbar$w1 / qbar$w, 0)
bsg$qbar <- qbar$q[match(bsg$genus, qbar$genus)]
site_dec <- aggregate(cbind(obs = w1, turnover_only = w * qbar) ~ site, data = bsg, FUN = sum)
site_dec <- merge(site_dec, site_env[, c("site", "map_era5_mm")], by = "site")
a2b <- data.frame(n_sites = nrow(site_dec),
  rho_observed_vs_logMAP = cor(site_dec$obs, log(site_dec$map_era5_mm), method = "spearman"),
  rho_turnover_only_vs_logMAP = cor(site_dec$turnover_only, log(site_dec$map_era5_mm), method = "spearman"),
  rho_within_residual_vs_logMAP = cor(site_dec$obs - site_dec$turnover_only, log(site_dec$map_era5_mm), method = "spearman"),
  p_within_residual = cor.test(site_dec$obs - site_dec$turnover_only, log(site_dec$map_era5_mm), method = "spearman", exact = FALSE)$p.value)
# within mixed genera: carrier fraction vs log MAP (genera with carriers and non-carriers in the catalogue)
gen_reps <- split(tax$user_genome, tax$genus)
mixed <- names(Filter(function(r) { r <- intersect(r, unique(ab$rep)); any(r %in% R1L) && any(!r %in% R1L) }, gen_reps))
mg <- bsg[bsg$genus %in% mixed & bsg$w > 0, ]; mg$q <- mg$w1 / mg$w
mg$lmap <- as.numeric(scale(log(site_env$map_era5_mm[match(mg$site, site_env$site)])))
fm <- lmer(qlogis(pmin(pmax(q, 0.005), 0.995)) ~ lmap + (1 | genus) + (1 | site), data = mg, weights = sqrt(w), REML = FALSE)
a2b$n_mixed_genera <- length(unique(mg$genus)); a2b$within_mixed_genera_logit_slope_per_SD_logMAP <- fixef(fm)["lmap"]
a2b$within_mixed_genera_se <- summary(fm)$coefficients["lmap", 2]
a2b$mixed_genera_share_of_bulk_1l <- sum(bsg$w1[bsg$genus %in% mixed]) / sum(bsg$w1)
# clipping sensitivity of the weighted Gaussian model on clipped logits (R2-20)
for (cl in c(0.001, 0.02)) {
  fc <- lmer(qlogis(pmin(pmax(q, cl), 1 - cl)) ~ lmap + (1 | genus) + (1 | site), data = mg, weights = sqrt(w), REML = FALSE)
  a2b[[paste0("slope_clip", cl)]] <- fixef(fc)["lmap"]; a2b[[paste0("se_clip", cl)]] <- summary(fc)$coefficients["lmap", 2]
}
wr(site_dec, "A2_aridity_site_decomposition.tsv"); wr(a2b, "A2_aridity_decomposition_summary.tsv")

# A6 carrier representation vs copies per carrier (1l)
ncopy <- table(cls$rep_id[cls$marker == "hydB" & grepl("Group_1l$", cls$assignment)])
sm$copies <- as.numeric(ncopy[sm$rep]); sm$copies[is.na(sm$copies)] <- 0
cp <- aggregate(cbind(F = share * c1l, M = share * copies) ~ sample + visit + compartment + site + trip, data = sm, FUN = sum)
cp$C <- ifelse(cp$F > 0, cp$M / cp$F, NA)
cpc <- cells(cp[is.finite(cp$C), ], "PlantRoot", "Deep", c("F", "M", "C"), eps = 1e-6)
cps <- cells(cp[is.finite(cp$C), ], "Surface", "Deep", c("F", "M", "C"), eps = 1e-6)
wr(rbind(contrast_table(cpc, c("F", "M", "C"), "root/bulk"), contrast_table(cps, c("F", "M", "C"), "surface/bulk")),
   "A6_copy_number_vs_carrier_representation.tsv")

# A4 CBB co-carriage (RuBisCO large + small subunit + phosphoribulokinase) by hydrogenase group
kt <- rd(file.path(IN, "rep_cbb_ko.tsv")); kt <- kt[!kt$rep %in% ORPHANS, ]
kt$cbb <- as.integer(kt$K01601 == 1 & kt$K01602 == 1 & kt$K00855 == 1)
kt$group <- ifelse(kt$rep %in% R1L & kt$rep %in% R1H, "1l+1h", ifelse(kt$rep %in% R1L, "1l", ifelse(kt$rep %in% R1H, "1h",
            ifelse(kt$rep %in% R2A, "2a_only", "no_high_affinity_or_2a"))))
kt$order <- tax$order[match(kt$rep, tax$user_genome)]; kt$phylum <- tax$phylum[match(kt$rep, tax$user_genome)]
gs <- if (file.exists(file.path(IN, "rep_genome_size.tsv"))) rd(file.path(IN, "rep_genome_size.tsv")) else NULL
a4 <- aggregate(cbind(n = 1, cbb) ~ group, data = kt, FUN = sum); a4$frac_cbb <- a4$cbb / a4$n
k2 <- kt[kt$group %in% c("1l", "1h") & !is.na(kt$order), ]; k2$is1l <- as.integer(k2$group == "1l")
gm <- glmer(cbb ~ is1l + (1 | order), data = k2, family = binomial)
a4$within_order_log_odds_1l_vs_1h <- fixef(gm)["is1l"]; a4$se <- summary(gm)$coefficients["is1l", 2]
a4$p <- summary(gm)$coefficients["is1l", 4]
if (!is.null(gs)) {
  k2$comp <- gs$completeness[match(k2$rep, gs$genome_id)]
  gm2 <- glmer(cbb ~ is1l + scale(comp) + (1 | order), data = k2, family = binomial)
  a4$completeness_adjusted_log_odds <- fixef(gm2)["is1l"]; a4$completeness_adjusted_p <- summary(gm2)$coefficients["is1l", 4]
  # high-completeness genomes only (>= 90%), where absence of a CBB gene is less likely to be missed (R2-21)
  k3 <- k2[is.finite(k2$comp) & k2$comp >= 90, ]
  gm3 <- glmer(cbb ~ is1l + (1 | order), data = k3, family = binomial)
  a4$hc90_n_1l <- sum(k3$is1l == 1); a4$hc90_n_1h <- sum(k3$is1l == 0)
  a4$hc90_frac_1l <- mean(k3$cbb[k3$is1l == 1]); a4$hc90_frac_1h <- mean(k3$cbb[k3$is1l == 0])
  a4$hc90_log_odds <- fixef(gm3)["is1l"]; a4$hc90_se <- summary(gm3)$coefficients["is1l", 2]; a4$hc90_p <- summary(gm3)$coefficients["is1l", 4]
}
wr(a4, "A4_cbb_cocarriage_by_group.tsv")
# abundance of 1l carriers with and without CBB: root contrast and aridity
cbb_reps <- kt$rep[kt$cbb == 1]
sm$sub <- ifelse(sm$c1l & sm$rep %in% cbb_reps, "l1_cbb", ifelse(sm$c1l, "l1_nocbb", "other"))
sp <- aggregate(share ~ sample + visit + compartment + site + trip + sub, data = sm[sm$sub != "other", ], FUN = sum)
sw <- reshape(sp, idvar = c("sample", "visit", "compartment", "site", "trip"), timevar = "sub", direction = "wide")
sw[is.na(sw)] <- 0; names(sw) <- sub("share\\.", "", names(sw))
sw$frac_cbb <- sw$l1_cbb / (sw$l1_cbb + sw$l1_nocbb)
a4c <- rbind(contrast_table(cells(sw, "PlantRoot", "Deep", c("l1_cbb", "l1_nocbb"), eps = 1e-5), c("l1_cbb", "l1_nocbb"), "root/bulk"),
             contrast_table(cells(sw, "Surface", "Deep", c("l1_cbb", "l1_nocbb"), eps = 1e-5), c("l1_cbb", "l1_nocbb"), "surface/bulk"))
swb <- sw[sw$compartment == "Deep", ]; sb <- aggregate(cbind(frac_cbb, l1_cbb, l1_nocbb) ~ site, swb, median)
sb <- merge(sb, site_env, by = "site")
a4c <- rbind(a4c, data.frame(contrast = "bulk sites vs log MAP (Spearman rho in median_log2)", quantity = c("frac_cbb", "l1_cbb", "l1_nocbb"),
  n_visits = NA, n_sites = nrow(sb), median_log2 = sapply(c("frac_cbb", "l1_cbb", "l1_nocbb"), function(v) cor(sb[[v]], log(sb$map_era5_mm), method = "spearman")),
  ci_lo = NA, ci_hi = NA, frac_below0 = NA,
  p_site_wilcoxon = sapply(c("frac_cbb", "l1_cbb", "l1_nocbb"), function(v) cor.test(sb[[v]], log(sb$map_era5_mm), method = "spearman", exact = FALSE)$p.value)))
wr(a4c, "A4_cbb_subsets_position_and_aridity.tsv")

# ================================================================== A3 surface filter vs exposure
sx <- merge(cs, st_env, by = c("site", "trip"))
a3 <- do.call(rbind, lapply(FOCAL, function(v) do.call(rbind, lapply(c("T_d30", "smap_sm_W30d", "P90_era5"), function(e) {
  x <- sx[is.finite(sx[[e]]), ]; x$y <- x[[paste0("excess_", v)]]; x$e <- as.numeric(scale(if (e == "P90_era5") log1p(x[[e]]) else x[[e]]))
  f1 <- lmer(y ~ e + trip + (1 | site), data = x, REML = FALSE); f0 <- lmer(y ~ trip + (1 | site), data = x, REML = FALSE)
  data.frame(marker = v, exposure = e, n_visits = nrow(x), beta_per_SD = fixef(f1)["e"], se = summary(f1)$coefficients["e", 2],
             LRT_p = anova(f0, f1)$`Pr(>Chisq)`[2])
}))))
a3$q_BH <- p.adjust(a3$LRT_p, "BH")
wr(a3, "A3_surface_excess_vs_exposure.tsv")
# residual support of the temperature exposure (R2-28): variance explained by expedition, residual SD, and the
# 1h estimate leaving out one expedition at a time
xt <- sx[is.finite(sx$T_d30), ]
a3t <- data.frame(n_visits = nrow(xt), sd_T = sd(xt$T_d30), sd_T_resid_trip = sd(resid(lm(T_d30 ~ trip, xt))),
                  r2_T_by_trip = summary(lm(T_d30 ~ trip, xt))$r.squared)
for (tp in sort(unique(xt$trip))) {
  x <- xt[xt$trip != tp, ]; x$e <- as.numeric(scale(x$T_d30)); x$y <- x$excess_hydB_1h
  f <- lmer(y ~ e + trip + (1 | site), data = x, REML = FALSE)
  a3t[[paste0("beta_1h_without_", tp)]] <- fixef(f)["e"]
}
wr(a3t, "A3_temperature_residual_support.tsv")

# ================================================================== H3 spatial (bulk site medians)
prank <- function(x, y, Z) {
  rx <- rank(x); ry <- rank(y); Zr <- apply(as.matrix(Z), 2, rank)
  ex <- resid(lm(rx ~ Zr)); ey <- resid(lm(ry ~ Zr)); ct <- cor.test(ex, ey)
  c(rho = unname(ct$estimate), p = ct$p.value)
}
xrf <- rd(file.path(IN, "xrf_site_compartment_panel.tsv")); xrf <- xrf[xrf$compartment == "deep", ]
topo <- rd(file.path(IN, "topography.tsv"))
topo <- aggregate(cbind(twi, slope_deg) ~ site, data = topo, FUN = mean)   # site 19 has two records (R2-01)
site_env <- merge(site_env, xrf[, c("site", "Ca", "Cl", "S", "Fe", "Al", "Si", "K", "Na")], by = "site", all.x = TRUE)
site_env <- merge(site_env, topo[, c("site", "twi", "slope_deg")], by = "site", all.x = TRUE)
# geochemical axes: PCA of centred log-ratio XRF (major elements, zeros replaced by half the minimum positive)
el <- c("Al", "Ca", "Cl", "Fe", "K", "Mg", "Na", "P", "S", "Si", "Ti", "Mn", "Sr", "Zr")
xm <- as.matrix(xrf[, el]); for (j in seq_len(ncol(xm))) { z <- xm[, j]; z[z <= 0] <- min(z[z > 0]) / 2; xm[, j] <- z }
clr <- log(xm) - rowMeans(log(xm)); pca <- prcomp(clr, scale. = TRUE)
geo <- data.frame(site = xrf$site, geoPC1 = pca$x[, 1], geoPC2 = pca$x[, 2])
site_env <- merge(site_env, geo, by = "site", all.x = TRUE)
stopifnot(!anyDuplicated(site_env$site))
# community-composition axis (bulk site means of CoverM profiles)
deep_ids <- meta$sample_id[meta$compartment == "Deep"]
cm <- xtabs(abund ~ sample + rep, data = ab[ab$sample %in% deep_ids, ]); cm <- unclass(cm); cm <- cm / rowSums(cm)
sid <- meta$site[match(rownames(cm), meta$sample_id)]
cm_site <- rowsum(cm, sid) / as.numeric(table(sid)[rownames(rowsum(cm, sid))])
site_env$comp_pc1_deep <- cmdscale(vegdist(cm_site, "bray"), k = 1)[, 1][as.character(site_env$site)]

blocks <- function(k) { set.seed(k); kmeans(scale(site_env[, c("lat", "lon")]), centers = k, nstart = 50)$cluster }
for (k in c(4, 6, 8, 10)) site_env[[paste0("block", k)]] <- blocks(k)
h3 <- NULL
for (comp in c("Deep", "Surface", "PlantRoot")) for (v in FOCAL) {
  x <- d0[d0$compartment == comp, ]
  s <- aggregate(x[[v]], list(site = x$site), median); names(s)[2] <- "ge"
  s <- merge(s, site_env, by = "site")
  for (pv in c("map_era5_mm", "map_power_mm", "map_bestmatch_mm")) {
    r0 <- cor.test(s$ge, s[[pv]], method = "spearman", exact = FALSE)
    pg <- prank(s$ge, s[[pv]], s[, c("lat", "lon", "elevation_m")])
    sc <- s[is.finite(s$comp_pc1_deep), ]; pc <- prank(sc$ge, sc[[pv]], sc[, "comp_pc1_deep", drop = FALSE])
    sx2 <- s[is.finite(s$geoPC1) & is.finite(s$Ca), ]
    pgeo <- prank(sx2$ge, sx2[[pv]], sx2[, c("geoPC1", "geoPC2")]); psalt <- prank(sx2$ge, sx2[[pv]], sx2[, c("Ca", "Cl", "S")])
    bb <- sapply(c(4, 6, 8, 10), function(k) { bk <- s[[paste0("block", k)]]
      q95(replicate(NB, { b <- sample(unique(bk), replace = TRUE)
        z <- do.call(rbind, lapply(b, function(j) s[bk == j, ])); suppressWarnings(cor(z$ge, z[[pv]], method = "spearman")) })) })
    # residual spatial autocorrelation of rank regression (Moran's I, inverse distance)
    rr <- resid(lm(rank(s$ge) ~ rank(s[[pv]])))
    dm <- as.matrix(dist(cbind(s$lon * cos(mean(s$lat) * pi / 180), s$lat))); w <- 1 / dm; w[!is.finite(w)] <- 0; diag(w) <- 0
    mi <- Moran.I(rr, w)
    h3 <- rbind(h3, data.frame(compartment = comp, marker = v, precip = pv, n_sites = nrow(s),
      rho = unname(r0$estimate), p = r0$p.value,
      block4_ci = sprintf("%.2f to %.2f", bb[1, 1], bb[2, 1]), block6_ci = sprintf("%.2f to %.2f", bb[1, 2], bb[2, 2]),
      block8_ci = sprintf("%.2f to %.2f", bb[1, 3], bb[2, 3]), block10_ci = sprintf("%.2f to %.2f", bb[1, 4], bb[2, 4]),
      moran_I_resid = mi$observed, moran_p = mi$p.value,
      partial_rho_latlonelev = pg[["rho"]], partial_p_latlonelev = pg[["p"]],
      partial_rho_geoPC12 = pgeo[["rho"]], partial_p_geoPC12 = pgeo[["p"]],
      partial_rho_CaClS = psalt[["rho"]], partial_p_CaClS = psalt[["p"]],
      partial_rho_compPC1 = pc[["rho"]], partial_p_compPC1 = pc[["p"]]))
  }
}
h3$q_BH_family <- ave(h3$p, h3$precip, FUN = function(p) p.adjust(p, "BH"))
wr(h3, "H3_spatial_by_marker_position.tsv")
x3 <- merge(d0, site_env[, c("site", "map_era5_mm")], by = "site")
h3mm <- do.call(rbind, lapply(FOCAL, function(v) {
  x3$y <- log2(x3[[v]] + EPS); x3$lmap <- as.numeric(scale(log(x3$map_era5_mm)))
  f1 <- lmer(y ~ compartment + trip + lmap + (1 | site) + (1 | visit), data = x3, REML = FALSE)
  f0 <- lmer(y ~ compartment + trip + (1 | site) + (1 | visit), data = x3, REML = FALSE)
  data.frame(marker = v, beta_per_SD_logMAP = fixef(f1)[["lmap"]], se = summary(f1)$coefficients["lmap", 2],
             LRT_p = anova(f0, f1)$`Pr(>Chisq)`[2])
}))
wr(h3mm, "H3_mixed_model_logMAP.tsv")

# ================================================================== H3 temporal: within/between decomposition (R1-03, R1-04, R1-11)
EQB <- 0.5   # equivalence bound: |change| < 0.5 log2 across the within-site IQR of the exposure
tcells <- function(comp, v, pv, eps = EPS) {
  x <- d0[d0$compartment == comp, ]
  cc <- aggregate(log2(x[[v]] + eps), list(site = x$site, trip = x$trip), median); names(cc)[3] <- "y"
  cc <- merge(cc, st_env, by = c("site", "trip")); cc <- cc[is.finite(cc[[pv]]), ]
  cc$lp <- log1p(cc[[pv]])
  cc <- cc[cc$site %in% names(which(table(cc$site) >= 2)), ]
  cc$lp_between <- ave(cc$lp, cc$site); cc$lp_within <- cc$lp - cc$lp_between
  cc$cell <- site_env$era5_cell[match(cc$site, site_env$site)]
  cc
}
wfit <- function(cc) tryCatch(fixef(lmer(y ~ lp_within + lp_between + trip + (1 | site), data = cc, REML = FALSE))[["lp_within"]],
                              error = function(e) NA)
cboot <- function(cc, unit = "site", nb = NB) q95(replicate(nb, {   # cluster bootstrap over sites or ERA5 cells
  s <- sample(unique(cc[[unit]]), replace = TRUE)
  wfit(do.call(rbind, lapply(seq_along(s), function(i) { z <- cc[cc[[unit]] == s[i], ]; z$site <- paste0(z$site, "_", i); z })))
}))
h3t <- NULL
for (comp in c("Deep", "Surface", "PlantRoot")) for (v in FOCAL)
  for (pv in c("P365_era5", "P90_era5", "P365_bestmatch", "P365_power", "P90_power")) {
  cc <- tcells(comp, v, pv)
  f <- lmer(y ~ lp_within + lp_between + trip + (1 | site), data = cc, REML = FALSE)
  f0 <- lmer(y ~ lp_between + trip + (1 | site), data = cc, REML = FALSE)
  b <- summary(f)$coefficients["lp_within", 1:2]
  iqr <- IQR(cc$lp_within)
  bt <- cboot(cc)
  h3t <- rbind(h3t, data.frame(compartment = comp, marker = v, exposure = pv, n_visits = nrow(cc), n_sites = length(unique(cc$site)),
    trips = paste(sort(unique(cc$trip)), collapse = ","), date_range = paste(range(cc$trip_date), collapse = " to "),
    within_slope = b[1], within_se = b[2], within_LRT_p = anova(f0, f)$`Pr(>Chisq)`[2],
    within_boot_lo = bt[1], within_boot_hi = bt[2], within_IQR = iqr,
    change_across_IQR_lo = bt[1] * iqr, change_across_IQR_hi = bt[2] * iqr,
    equivalent_within_0.5 = (abs(bt[1] * iqr) < EQB & abs(bt[2] * iqr) < EQB),
    between_slope = fixef(f)[["lp_between"]]))
}
h3t$q_BH_family <- ave(h3t$within_LRT_p, h3t$exposure, FUN = function(p) p.adjust(p, "BH"))
wr(h3t, "H3_temporal_within_between.tsv")
# robustness of the within-site estimates (R2-15/16/17): ERA5-cell bootstrap, leave-one-expedition-out,
# pseudocount 0.001 and 0.1 (change across the IQR, 95% site-bootstrap interval), 500 resamples each
NB2 <- 500; h3r <- NULL
for (comp in c("Deep", "Surface", "PlantRoot")) for (v in FOCAL) for (pv in c("P90_era5", "P365_era5", "P90_power")) {
  cc <- tcells(comp, v, pv); iqr <- IQR(cc$lp_within)
  lo <- sapply(sort(unique(cc$trip)), function(tp) { z <- cc[cc$trip != tp, ]
    z <- z[z$site %in% names(which(table(z$site) >= 2)), ]; z$lp_between <- ave(z$lp, z$site); z$lp_within <- z$lp - z$lp_between
    wfit(z) * iqr })
  gb <- cboot(cc, "cell", NB2) * iqr
  p1 <- tcells(comp, v, pv, 0.001); p2 <- tcells(comp, v, pv, 0.1)
  b1 <- cboot(p1, "site", NB2) * IQR(p1$lp_within); b2 <- cboot(p2, "site", NB2) * IQR(p2$lp_within)
  h3r <- rbind(h3r, data.frame(compartment = comp, marker = v, exposure = pv, n_visits = nrow(cc), n_sites = length(unique(cc$site)),
    n_era5_cells = length(unique(cc$cell)), change_IQR = wfit(cc) * iqr,
    cellboot_lo = gb[1], cellboot_hi = gb[2], loo_min = min(lo, na.rm = TRUE), loo_max = max(lo, na.rm = TRUE),
    loo_detail = paste(sprintf("%s:%.2f", names(lo), lo), collapse = ";"),
    pc0.001_lo = b1[1], pc0.001_hi = b1[2], pc0.1_lo = b2[1], pc0.1_hi = b2[2],
    frac_visits_at_or_below_0.01 = mean(cc$y <= log2(2 * EPS))))
}
wr(h3r, "H3_temporal_robustness.tsv")
# agreement of the two products on identical visits (within-site deviations of log1p antecedent precipitation)
ag <- do.call(rbind, lapply(c(90, 365), function(w) {
  z <- st_env[is.finite(st_env[[paste0("P", w, "_era5")]]) & is.finite(st_env[[paste0("P", w, "_power")]]), ]
  z <- z[z$site %in% names(which(table(z$site) >= 2)), ]
  e <- log1p(z[[paste0("P", w, "_era5")]]); p <- log1p(z[[paste0("P", w, "_power")]])
  de <- e - ave(e, z$site); dp <- p - ave(p, z$site)
  data.frame(window_days = w, n_visits = nrow(z), rho_levels = cor(e, p, method = "spearman"),
             rho_within_site = cor(de, dp, method = "spearman"), frac_same_sign_within = mean(sign(de) == sign(dp)))
}))
wr(ag, "H3_temporal_product_agreement.tsv")

# ================================================================== stability (R1-15, R1-16)
pm_out <- NULL
# Units are visit-level profiles (mean relative abundance over replicate samples of a site, expedition and
# position), so replicates never exchange separately (R2-18). Expedition is tested with whole visits permuted
# only within site (and within position when positions are pooled); dispersion uses the same design.
for (comp in c("all", "Deep", "Surface", "PlantRoot")) {
  ids <- if (comp == "all") d0$sample_id else d0$sample_id[d0$compartment == comp]
  m <- xtabs(abund ~ sample + rep, data = ab[ab$sample %in% ids, ]); m <- unclass(m); m <- m / rowSums(m)
  md <- d0[match(rownames(m), d0$sample_id), ]
  cellk <- paste(md$site, md$trip, md$compartment, sep = "|")
  mv <- rowsum(m, cellk) / as.numeric(table(cellk)[rownames(rowsum(m, cellk))])
  mdv <- md[match(rownames(mv), cellk), c("site", "trip", "compartment")]
  keep <- ave(seq_len(nrow(mdv)), paste(mdv$site, mdv$compartment), FUN = length) >= 2   # sites with repeat visits
  mv <- mv[keep, ]; mdv <- mdv[keep, ]
  bc <- vegdist(mv, "bray")
  blk <- factor(paste(mdv$site, mdv$compartment))
  ctrl_trip <- how(nperm = 999, blocks = blk)
  f <- if (comp == "all") bc ~ factor(trip) + factor(site) + compartment else bc ~ factor(trip) + factor(site)
  a_all <- adonis2(f, data = mdv, by = "margin", permutations = 999)
  a_trip <- adonis2(f, data = mdv, by = "margin", permutations = ctrl_trip)
  disp <- permutest(betadisper(bc, mdv$trip), permutations = ctrl_trip)
  pm_out <- rbind(pm_out, data.frame(subset = comp, n_visit_profiles = nrow(mv), n_sites = length(unique(mdv$site)),
    term = rownames(a_all), Df = a_all$Df, R2 = a_all$R2,
    R2_per_df = a_all$R2 / a_all$Df, p_free = a_all$`Pr(>F)`,
    p_trip_within_site = ifelse(rownames(a_all) == "factor(trip)", a_trip$`Pr(>F)`[rownames(a_trip) == "factor(trip)"], NA),
    trip_dispersion_p_within_site = disp$tab$`Pr(>F)`[1], row.names = NULL))
}
wr(pm_out, "stability_permanova_restricted.tsv")
# synchrony of 1l-carrier genera across trips within bulk soil only, sites with >= 4 trips
car <- sm[sm$c1l & sm$compartment == "Deep", ]
gt <- aggregate(share ~ site + trip + genus, data = car, FUN = mean)
syn <- do.call(rbind, lapply(split(gt, gt$site), function(x) {
  w <- xtabs(share ~ trip + genus, data = x); w <- unclass(w)
  if (nrow(w) < 4 || ncol(w) < 2) return(NULL)
  phi <- function(w) var(rowSums(w)) / sum(apply(w, 2, sd))^2
  null <- replicate(199, phi(apply(w, 2, sample)))
  data.frame(site = x$site[1], n_trips = nrow(w), n_genera = ncol(w), phi = phi(w), phi_null_median = median(null))
}))
wr(syn, "stability_synchrony_bulk_1l_genera.tsv")

# ================================================================== sensitivity of headline estimates
variants <- list(primary = list(file = "guild_genome_equivalents_v2.tsv", floor = 60, drop = character()),
                 id50 = list(file = "guild_genome_equivalents_v2.tsv", floor = 50, drop = character()),
                 id70 = list(file = "guild_genome_equivalents_v2.tsv", floor = 70, drop = character()),
                 without_trip5 = list(file = "guild_genome_equivalents_v2.tsv", floor = 60, drop = "trip5"),
                 without_trip4 = list(file = "guild_genome_equivalents_v2.tsv", floor = 60, drop = "trip4"))
for (f in c("guild_genome_equivalents_v2_noflagged.tsv", "guild_genome_equivalents_v2_strict1l.tsv",
            "guild_genome_equivalents_v2_motifonly.tsv", "guild_genome_equivalents_v2_coxLcorr.tsv"))
  if (file.exists(file.path(RD, f))) variants[[sub("guild_genome_equivalents_v2_", "", sub(".tsv", "", f))]] <- list(file = f, floor = 60, drop = character())
sens <- do.call(rbind, lapply(names(variants), function(nm) {
  v <- variants[[nm]]; d <- load_ge(v$file, v$floor, v$drop)
  r <- add_excess(cells(d, "PlantRoot", "Deep", MARKERS)); s <- add_excess(cells(d, "Surface", "Deep", MARKERS))
  b <- aggregate(d$hydB_1l[d$compartment == "Deep"], list(site = d$site[d$compartment == "Deep"]), median)
  b <- merge(b, site_env, by = "site")
  data.frame(variant = nm, n_samples = nrow(d), root_1l = median(r$hydB_1l), root_excess_1l = median(r$excess_hydB_1l),
             root_excess_1l_p = site_wilcox(r, "excess_hydB_1l"), root_excess_1h = median(r$excess_hydB_1h),
             root_excess_formI = median(r$excess_coxL_formI),
             root_excess_1l_vs_coxLother = median(r$hydB_1l - r$coxL_other),     # R2: other-[NiFe] is mostly off-target
             root_excess_1l_vs_nife = median(r$hydB_1l - r$hydB_other_nife),
             surface_excess_1l = median(s$excess_hydB_1l),
             surface_excess_formI = median(s$excess_coxL_formI),
             rho_bulk_1l_ERA5 = cor(b$x, b$map_era5_mm, method = "spearman"),
             rho_bulk_1l_POWER = cor(b$x, b$map_power_mm, method = "spearman"))
}))
wr(sens, "sensitivity.tsv")

sink(file.path(OUT, "SUMMARY.txt"))
cat("Life on air (P3), revised analyses (h123_v2.R), sites 1-60, GE v2 at 60% identity\n")
cat(sprintf("samples %d (surface %d, bulk %d, root-adjacent %d) at %d sites; with DNA yield %d\n", nrow(d0),
    sum(d0$compartment == "Surface"), sum(d0$compartment == "Deep"), sum(d0$compartment == "PlantRoot"),
    length(unique(d0$site)), sum(is.finite(d0$dna))))
cat("\nbulk medians (GE):\n"); print(round(sapply(MARKERS, function(v) median(d0[d0$compartment == "Deep", v])), 3))
for (nm in c("contr", "pc_sens", "h1_trip", "jm", "pos_eff", "a1", "a1b", "a1c", "a2a", "a2b", "a4", "a4c", "a3",
             "h3", "h3mm", "h3t", "h3r", "ag", "a3t", "pm_out", "sens")) { cat("\n==", nm, "==\n"); print(get(nm), row.names = FALSE, digits = 3) }
cat(sprintf("\nsynchrony (bulk, 1l genera): median phi %.3f vs null %.3f at %d sites\n", median(syn$phi), median(syn$phi_null_median), nrow(syn)))
sink()
cat(readLines(file.path(OUT, "SUMMARY.txt")), sep = "\n")
