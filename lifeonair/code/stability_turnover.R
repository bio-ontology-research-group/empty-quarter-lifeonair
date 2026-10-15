#!/usr/bin/env Rscript
# Life on air (P3), analyses 3 and 4: stability of guild capacity versus turnover of composition.
#
# Measures (per sample):
#   reads   read-based genome equivalents (r3_aggregate.py), identity floor 60%:
#           hydB_ha (group 1h + 1l), coxL_formI
#   mags    CoverM MAG-direct guild abundance, % of mapped community, with and without the
#           contamination-flagged representatives (decontam_screen.R)
# Analysis 3 (stability): lme4 mixed models  y ~ trip + compartment + (1 | site), trip and
#   compartment tested by likelihood ratio (ML fits); marginal variance shares; repeated on the
#   37 sites sampled in >= 4 trips.
# Analysis 4 (turnover): PERMANOVA (vegan::adonis2, Bray-Curtis on Hellinger-free relative
#   abundances, terms trip + site + compartment, by = "margin") for the whole mapped community and
#   for guild carriers only; guild synchrony per site (Loreau & de Mazancourt 2008 phi) across trips
#   at genus level, with a within-site trip-label permutation null.
#
#   Rscript lifeonair/code/stability_turnover.R   (from the repository root)

suppressMessages({ library(lme4); library(vegan) })
set.seed(20261009)
IN  <- "lifeonair/results/inputs"
OUT <- "lifeonair/results/stability"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f, ...) read.delim(f, sep = "\t", quote = "", comment.char = "",
                                  stringsAsFactors = FALSE, check.names = FALSE, ...)

meta <- rd(file.path(IN, "samples.tsv"))
meta <- meta[meta$trip != "" & meta$compartment != "", ]
meta$trip <- factor(meta$trip); meta$compartment <- factor(meta$compartment)
meta$site <- factor(meta$site)
trips_per_site <- tapply(meta$trip, meta$site, function(x) length(unique(x)))
long_sites <- names(trips_per_site)[trips_per_site >= 4]

# ---- measures ------------------------------------------------------------------------------------
measures <- list()
mg <- rd("lifeonair/results/decontam/guild_abundance_with_without_flagged.tsv")
for (v in c("union_mapped", "union_mapped_clean", "hydB_mapped_clean", "coxL_mapped_clean"))
  measures[[paste0("mags_", v)]] <- setNames(mg[[v]], mg$sample)
ge_f <- "lifeonair/results/reads/guild_genome_equivalents.tsv"
if (file.exists(ge_f)) {
  ge <- rd(ge_f); ge <- ge[ge$id_floor == 60, ]
  for (v in c("hydB_ha", "hydB_1l", "hydB_1h", "coxL_formI"))
    measures[[paste0("reads_", v)]] <- setNames(ge[[v]], ge$sample)
}

# ---- analysis 3: stability ----------------------------------------------------------------------
lrt <- function(full, reduced) { a <- anova(reduced, full); c(chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2]) }
stab <- list()
for (nm in names(measures)) for (subset in c("all", "sites_ge4_trips")) {
  d <- meta; d$y <- measures[[nm]][d$sample_id]
  d <- d[is.finite(d$y), ]
  if (subset != "all") d <- d[d$site %in% long_sites, ]
  d <- droplevels(d)
  if (nrow(d) < 50) next
  f0 <- lmer(y ~ trip + compartment + (1 | site), data = d, REML = FALSE)
  ft <- lmer(y ~ compartment + (1 | site), data = d, REML = FALSE)
  fc <- lmer(y ~ trip + (1 | site), data = d, REML = FALSE)
  vc <- as.data.frame(VarCorr(f0)); vfix <- var(as.vector(model.matrix(f0) %*% fixef(f0)))
  vtrip <- var(as.vector(model.matrix(f0)[, grep("^trip", colnames(model.matrix(f0))), drop = FALSE] %*%
                         fixef(f0)[grep("^trip", names(fixef(f0)))]))
  tot <- vfix + sum(vc$vcov)
  trip_med <- tapply(d$y, d$trip, median)
  stab[[length(stab) + 1]] <- data.frame(
    measure = nm, subset = subset, n_samples = nrow(d), n_sites = nlevels(d$site),
    median = median(d$y), q25 = quantile(d$y, .25), q75 = quantile(d$y, .75),
    trip_LRT_chisq = lrt(f0, ft)[1], trip_df = lrt(f0, ft)[2], trip_p = lrt(f0, ft)[3],
    compartment_LRT_chisq = lrt(f0, fc)[1], compartment_p = lrt(f0, fc)[3],
    var_share_trip = vtrip / tot, var_share_site = vc$vcov[vc$grp == "site"] / tot,
    var_share_residual = vc$vcov[vc$grp == "Residual"] / tot,
    trip_medians = paste(sprintf("%s=%.3f", names(trip_med), trip_med), collapse = ";"),
    row.names = NULL)
}
stab <- do.call(rbind, stab)
write.table(stab, file.path(OUT, "stability_mixed_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- analysis 4: turnover --------------------------------------------------------------------------
ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE,
         col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund) & ab$abund > 0 & ab$sample %in% meta$sample_id, ]
flags <- rd("lifeonair/results/decontam/rep_contamination_flags.tsv")
ab <- ab[!ab$rep %in% flags$rep[flags$flag_any == "TRUE" | flags$flag_any == TRUE], ]
guild <- readLines(file.path(IN, "guild_reps.txt"))
mat <- xtabs(abund ~ sample + rep, data = ab, sparse = FALSE)
mat <- unclass(mat); mat <- mat[rowSums(mat) > 0, ]
m <- meta[match(rownames(mat), meta$sample_id), ]

perm <- function(X, label) {
  X <- X[rowSums(X) > 0, , drop = FALSE]; mm <- m[match(rownames(X), m$sample_id), ]
  X <- X / rowSums(X)
  a <- adonis2(vegdist(X, "bray") ~ trip + site + compartment, data = mm, by = "margin", permutations = 999)
  data.frame(community = label, term = rownames(a)[1:3], R2 = a$R2[1:3], F = a$F[1:3], p = a$`Pr(>F)`[1:3],
             n_samples = nrow(X), n_reps = sum(colSums(X) > 0))
}
pm_f <- file.path(OUT, "permanova_turnover.tsv")
if (file.exists(pm_f) && Sys.getenv("LOA_RERUN_PERMANOVA") == "") {
  pm <- rd(pm_f)          # cached: 3 x 999 permutations take ~1.5 h; set LOA_RERUN_PERMANOVA=1 to redo
} else {
  pm <- rbind(perm(mat, "whole_mapped_community"),
              perm(mat[, colnames(mat) %in% guild, drop = FALSE], "guild_carriers"),
              perm(mat[, !colnames(mat) %in% guild, drop = FALSE], "non_guild"))
  write.table(pm, pm_f, sep = "\t", quote = FALSE, row.names = FALSE)
}

# synchrony of guild-carrier genera across trips within long-sampled sites.
# Per sample, genus abundance as % of the mapped community; per site x trip, the mean over the
# replicate samples of that trip (all compartments); phi = var(total) / (sum_i sd_i)^2 across trips.
mapped <- tapply(ab$abund, ab$sample, sum)
gl <- ab[ab$rep %in% guild, ]
gen <- flags$genus[match(gl$rep, flags$rep)]; gen[is.na(gen) | gen == ""] <- "unassigned"
gl$genus <- gen
gl$pct <- 100 * gl$abund / as.numeric(mapped[gl$sample])
sg <- xtabs(pct ~ sample + genus, data = gl); sg <- unclass(sg)
sg_meta <- meta[match(rownames(sg), meta$sample_id), ]
phi_of <- function(P) var(rowSums(P)) / sum(apply(P, 2, sd))^2
syn <- do.call(rbind, lapply(long_sites, function(s) {
  idx <- which(sg_meta$site == s)
  if (!length(idx)) return(NULL)
  tr <- droplevels(sg_meta$trip[idx])
  if (nlevels(tr) < 4) return(NULL)
  P <- rowsum(sg[idx, , drop = FALSE], as.character(tr)) / as.numeric(table(as.character(tr)))
  P <- P[, colSums(P) > 0, drop = FALSE]
  tot <- rowSums(P)
  phi <- phi_of(P)
  null <- replicate(999, phi_of(apply(P, 2, sample)))
  cvg <- apply(P, 2, function(v) sd(v) / mean(v))
  data.frame(site = s, n_trips = nrow(P), n_samples = length(idx), n_genera = ncol(P),
             guild_total_mean = mean(tot), cv_guild_total = sd(tot) / mean(tot),
             median_cv_genus = median(cvg), phi = phi, phi_null_median = median(null),
             p_phi_lower = mean(null <= phi))
}))
write.table(syn, file.path(OUT, "guild_synchrony_by_site.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(OUT, "SUMMARY.txt"))
cat("== stability (mixed models; LRT for trip and compartment) ==\n")
print(stab[, c("measure", "subset", "n_samples", "n_sites", "median", "trip_p", "compartment_p",
               "var_share_trip", "var_share_site", "var_share_residual")], row.names = FALSE, digits = 3)
cat("\n== turnover (PERMANOVA, marginal terms) ==\n"); print(pm, row.names = FALSE, digits = 3)
cat(sprintf("\n== guild synchrony across trips, %d sites with >= 4 trips ==\n", nrow(syn)))
cat(sprintf("median phi %.3f (null median %.3f); sites with phi below null at p<0.05: %d\n",
            median(syn$phi), median(syn$phi_null_median), sum(syn$p_phi_lower < 0.05)))
cat(sprintf("median CV guild total %.3f vs median mean-CV of carrier genera %.3f\n",
            median(syn$cv_guild_total), median(syn$median_cv_genus)))
sink()
cat(readLines(file.path(OUT, "SUMMARY.txt")), sep = "\n")
