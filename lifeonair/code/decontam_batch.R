#!/usr/bin/env Rscript
# Life on air (P3): batch-stratified decontam frequency screen (review round 1, finding R1-23).
# decontam_screen.R pooled all expeditions and extraction kits, so its frequency test compared samples
# whose DNA yield differs by batch. Here the two models (abundance independent of DNA concentration vs
# proportional to 1/concentration) are compared after centring within expedition x extraction kit, with a
# pooled F statistic over strata (df = n - number of strata). Flags at score < 0.1 with >= 10 samples.
#   Rscript lifeonair/code/decontam_batch.R   -> lifeonair/results/decontam/rep_frequency_by_batch.tsv
IN <- "lifeonair/results/inputs"; OUT <- "lifeonair/results/decontam"
rd <- function(f, ...) read.delim(f, quote = "", comment.char = "", check.names = FALSE, ...)
ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE, col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund) & ab$abund > 0, ]
led <- rd(file.path(IN, "sample_ledger_dna.tsv")); led$conc <- suppressWarnings(as.numeric(led$dna_concentration))
led <- led[is.finite(led$conc) & led$conc > 0, ]
nu <- tapply(led$conc, led$sample_id, function(x) length(unique(round(x, 6))))
led <- led[!led$sample_id %in% names(nu)[nu > 1] & !duplicated(led$sample_id), ]
led$batch <- paste(tolower(led$trip), led$dna_kit)
ab <- ab[ab$sample %in% led$sample_id, ]
ab$batch <- led$batch[match(ab$sample, led$sample_id)]; ab$x <- log(led$conc[match(ab$sample, led$sample_id)])
cat("samples", length(unique(ab$sample)), "batches:\n"); print(table(led$batch[led$sample_id %in% ab$sample]))
res <- do.call(rbind, lapply(split(ab, ab$rep), function(d) {
  n <- nrow(d); k <- length(unique(d$batch))
  if (n < 10 || n - k < 5) return(data.frame(rep = d$rep[1], n_present = n, n_batches = k, score = NA))
  y <- log(d$abund); r1 <- y + d$x
  cy <- y - ave(y, d$batch); cr <- r1 - ave(r1, d$batch)
  data.frame(rep = d$rep[1], n_present = n, n_batches = k, score = pf(sum(cr^2) / sum(cy^2), n - k, n - k))
}))
res$flag <- !is.na(res$score) & res$score < 0.1
write.table(res, file.path(OUT, "rep_frequency_by_batch.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
g <- unique(c(readLines(file.path(IN, "guild_hydB_reps.txt")), readLines(file.path(IN, "guild_coxL_reps_v2.txt"))))
cat(sprintf("testable %d; flagged %d (guild carriers %d)\n", sum(!is.na(res$score)), sum(res$flag), sum(res$flag & res$rep %in% g)))
writeLines(res$rep[res$flag], file.path(OUT, "flagged_reps_by_batch.txt"))
