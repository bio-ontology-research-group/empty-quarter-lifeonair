#!/usr/bin/env Rscript
# Life on air (P3): genomic traits of atmospheric-guild carriers (review round 1, consultation A5).
# Do group 1l carriers differ from 1h carriers (and non-carriers) in estimated genome size, GC content or
# gRodon doubling time, within order? Genome size is completeness-adjusted (length / completeness, CheckM2),
# from IBEX 06_atmospheric/lifeonair_review_r1 job 53611804 (results/inputs/rep_genome_size.tsv).
#   Rscript lifeonair/code/a5_genome_traits.R   -> lifeonair/results/a5_genome_traits/
suppressPackageStartupMessages(library(lme4))
IN <- "lifeonair/results/inputs"; OUT <- "lifeonair/results/a5_genome_traits"; dir.create(OUT, showWarnings = FALSE)
rd <- function(f, ...) read.delim(f, check.names = FALSE, quote = "", ...)
ORPHANS <- readLines("lifeonair/results/reads/orphan_sub_reps.txt")
gs <- rd(file.path(IN, "rep_genome_size.tsv")); gs <- gs[!gs$genome_id %in% ORPHANS, ]
gr <- rd(file.path(IN, "rep_growth_v2.tsv"))
gs$logd <- log(gr$d_hours[match(gs$genome_id, gr$rep)])
gs$est_mb <- gs$length_bp / (gs$completeness / 100) / 1e6
cls <- rd(file.path(IN, "atmo_classification_v2.tsv")); cls <- cls[!cls$rep_id %in% ORPHANS, ]
g <- function(x) unique(cls$rep_id[cls$marker == "hydB" & grepl(paste0("Group_", x, "$"), cls$assignment)])
gs$c1l <- gs$genome_id %in% g("1l"); gs$c1h <- gs$genome_id %in% g("1h"); gs$c2a <- gs$genome_id %in% g("2a")
gs$cfI <- gs$genome_id %in% setdiff(readLines(file.path(IN, "guild_coxL_reps_v2.txt")), ORPHANS)
tax <- rbind(rd(file.path(IN, "r232_bac120.summary.tsv"))[, 1:2], rd(file.path(IN, "r232_ar53.summary.tsv"))[, 1:2])
rk <- function(cl, p) ifelse(grepl(paste0(p, "__[^;]"), cl), sub(paste0(".*", p, "__([^;]*).*"), "\\1", cl), "unassigned")
cl <- tax$classification[match(gs$genome_id, tax$user_genome)]
gs$order <- paste(rk(cl, "p"), rk(cl, "o"), sep = ";")
gs$group <- ifelse(gs$c1l & gs$c1h, "1l+1h", ifelse(gs$c1l, "1l_only", ifelse(gs$c1h, "1h_only", "neither")))
desc <- do.call(rbind, lapply(split(gs, gs$group), function(d) data.frame(group = d$group[1], n = nrow(d),
  median_est_Mb = median(d$est_mb), median_gc = median(d$gc), median_doubling_h = median(exp(d$logd), na.rm = TRUE),
  median_completeness = median(d$completeness))))
fit <- function(resp) {
  m <- lmer(as.formula(paste(resp, "~ c1l + c1h + c2a + cfI + (1 | order)")), gs)
  co <- summary(m)$coefficients
  data.frame(response = resp, term = rownames(co)[-1], estimate = co[-1, 1], se = co[-1, 2],
             p_wald = 2 * pnorm(-abs(co[-1, 3])))
}
gs$log_est_mb <- log(gs$est_mb); gs$gc100 <- 100 * gs$gc
mod <- do.call(rbind, lapply(c("log_est_mb", "gc100", "logd"), fit))
mod$q_BH <- p.adjust(mod$p_wald, "BH")
write.table(desc, file.path(OUT, "trait_medians_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mod, file.path(OUT, "trait_models_within_order.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(desc); print(mod, digits = 3)
