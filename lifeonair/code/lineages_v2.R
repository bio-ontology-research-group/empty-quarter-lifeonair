#!/usr/bin/env Rscript
# Life on air (P3): orders of group 1l carriers beside roots (Table S3), revised after review round 1
# (sites 1-60, carriers restricted to retained contigs). Share of the mapped community held by the group 1l
# carriers of each order per sample (CoverM), median over replicates per site visit, log2 ratio of
# root-adjacent to bulk (pseudocount 0.001 percentage points), for the ten orders with the highest bulk share.
#   Rscript lifeonair/code/lineages_v2.R   -> lifeonair/results/h123_v2/h1_lineages_1l_carriers.tsv
IN <- "lifeonair/results/inputs"
rd <- function(f, ...) read.delim(f, quote = "", comment.char = "", check.names = FALSE, ...)
ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE, col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund), ]
tot <- tapply(ab$abund, ab$sample, sum)
cls <- rd(file.path(IN, "atmo_classification_v2.tsv"))
r1l <- unique(cls$rep_id[cls$marker == "hydB" & grepl("Group_1l$", cls$assignment)])
tax <- rbind(rd(file.path(IN, "r232_bac120.summary.tsv"))[, 1:2], rd(file.path(IN, "r232_ar53.summary.tsv"))[, 1:2])
rk <- function(cl, p) ifelse(grepl(paste0(p, "__[^;]"), cl), sub(paste0(".*", p, "__([^;]*).*"), "\\1", cl), "unassigned")
ab <- ab[ab$rep %in% r1l, ]
cl <- tax$classification[match(ab$rep, tax$user_genome)]
ab$order <- paste(rk(cl, "p"), rk(cl, "o"), sep = ";")
meta <- rd(file.path(IN, "samples.tsv")); meta <- meta[!is.na(meta$site) & meta$site %in% 1:60, ]
ge <- rd("lifeonair/results/reads/guild_genome_equivalents_v2.tsv"); ge <- ge[ge$id_floor == 60 & ge$scg_genes_detected >= 10, ]
meta <- meta[meta$sample_id %in% ge$sample & meta$sample_id %in% names(tot), ]
sh <- tapply(ab$abund, list(ab$sample, ab$order), sum); sh[is.na(sh)] <- 0
sh <- 100 * sh[match(meta$sample_id, rownames(sh)), , drop = FALSE] / as.numeric(tot[meta$sample_id]); sh[is.na(sh)] <- 0
deep <- meta$compartment == "Deep"
top <- names(sort(colMeans(sh[deep, ]), decreasing = TRUE))[1:10]
vis <- paste(meta$site, meta$trip)
out <- do.call(rbind, lapply(top, function(o) {
  md <- tapply(sh[, o], list(vis, meta$compartment), median)
  ok <- is.finite(md[, "Deep"]) & is.finite(md[, "PlantRoot"])
  l <- log2((md[ok, "PlantRoot"] + 0.001) / (md[ok, "Deep"] + 0.001))
  data.frame(order = o, deep_mean_pct_mapped = mean(sh[deep, o]), n_visits = sum(ok),
             median_log2fc_plantroot_vs_deep = median(l), frac_cells_below0 = mean(l < 0))
}))
write.table(out, "lifeonair/results/h123_v2/h1_lineages_1l_carriers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out)
