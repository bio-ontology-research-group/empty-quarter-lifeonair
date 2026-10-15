#!/usr/bin/env Rscript
# Life on air (P3), analysis 1: MAG-level contamination screen for the 5,438 ANI95 representatives.
#
# Two independent flags, reported separately and as a union:
#   blank_genus  the representative's GTDB R232 genus (suffix stripped) matches the genus of a
#                Trip-5 16S feature called candidate_contaminant from 17 extraction blanks
#                (data-paper evidence/control-audit/trip5_primary_contaminant_calls.tsv)
#   frequency    decontam frequency method (Davis et al. 2018) re-implemented: among samples with a
#                recorded extraction DNA concentration where the representative is present, compare
#                log(abundance) ~ 1 against log(abundance) ~ offset(-log(conc));
#                score = pf(SS_contam / SS_null, n - 1, n - 1); flagged when score < 0.1, n >= 10.
# Guild abundance (union, hydB, coxL) is then recomputed per sample with and without flagged reps.
#
#   Rscript lifeonair/code/decontam_screen.R   (from the repository root)

suppressMessages({ library(stats); library(utils) })
IN  <- "lifeonair/results/inputs"
OUT <- "lifeonair/results/decontam"
CTL <- path.expand("~/Public/software/empty-quarter-paper-repos/data-paper/evidence/control-audit/trip5_primary_contaminant_calls.tsv")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(f, ...) read.delim(f, sep = "\t", quote = "", comment.char = "",
                                  stringsAsFactors = FALSE, check.names = FALSE, ...)

# ---- inputs -------------------------------------------------------------------------------------
tax <- rbind(rd(file.path(IN, "r232_bac120.summary.tsv"))[, c("user_genome", "classification")],
             rd(file.path(IN, "r232_ar53.summary.tsv"))[, c("user_genome", "classification")])
rank_of <- function(cl, p) { x <- sub(paste0(".*", p, "__([^;]*).*"), "\\1", cl); ifelse(grepl(paste0(p, "__"), cl), x, "") }
tax$genus  <- sub("_[A-Z]+$", "", rank_of(tax$classification, "g"))
tax$phylum <- rank_of(tax$classification, "p")

ab <- rd(gzfile(file.path(IN, "coverm_rep_abundance_long.tsv.gz")), header = FALSE,
         col.names = c("sample", "rep", "abund"))
ab <- ab[ab$rep != "unmapped" & is.finite(ab$abund) & ab$abund > 0, ]  # F5Dr1 profile is all-NaN
samples <- rd(file.path(IN, "samples.tsv"))
reps <- sort(unique(ab$rep))
stopifnot(length(reps) <= 5438)
tax <- tax[tax$user_genome %in% reps, ]

guild <- list(union = readLines(file.path(IN, "guild_reps.txt")),
              hydB  = readLines(file.path(IN, "guild_hydB_reps.txt")),
              coxL  = readLines(file.path(IN, "guild_coxL_reps.txt")))

# ---- flag 1: genera of blank-enriched features --------------------------------------------------
ctl <- rd(CTL)
ctl <- ctl[ctl$call == "candidate_contaminant", ]
ctl_genus <- vapply(strsplit(ctl$taxonomy, ";"), function(x) if (length(x) >= 6) x[6] else "", "")
ctl_genus <- sort(unique(ctl_genus[!ctl_genus %in% c("", "NA") & !grepl("^uncultured", ctl_genus)]))
writeLines(ctl_genus, file.path(OUT, "blank_contaminant_genera.txt"))
flag_blank <- tax$user_genome[tax$genus %in% ctl_genus]

# ---- flag 2: decontam frequency method ----------------------------------------------------------
led <- rd(file.path(IN, "sample_ledger_dna.tsv"))
led$dna_concentration <- suppressWarnings(as.numeric(led$dna_concentration))
led <- led[is.finite(led$dna_concentration) & led$dna_concentration > 0, ]
conc_n <- tapply(led$dna_concentration, led$sample_id, function(x) length(unique(round(x, 6))))
conflict <- names(conc_n)[conc_n > 1]
led <- led[!led$sample_id %in% conflict & !duplicated(led$sample_id), ]
conc <- setNames(led$dna_concentration, led$sample_id)
conc <- conc[names(conc) %in% samples$sample_id]

abc <- ab[ab$sample %in% names(conc), ]
freq <- do.call(rbind, lapply(split(abc, abc$rep), function(d) {
  n <- nrow(d)
  if (n < 10) return(data.frame(rep = d$rep[1], n_present = n, score = NA_real_))
  y <- log(d$abund); x <- log(conc[d$sample])
  ss0 <- sum((y - mean(y))^2)                 # non-contaminant: abundance independent of input
  r1  <- y + x; ss1 <- sum((r1 - mean(r1))^2) # contaminant: abundance proportional to 1 / input
  data.frame(rep = d$rep[1], n_present = n, score = pf(ss1 / ss0, n - 1, n - 1))
}))
flag_freq <- freq$rep[!is.na(freq$score) & freq$score < 0.1]

flags <- data.frame(rep = reps)
flags <- merge(flags, tax[, c("user_genome", "phylum", "genus")], by.x = "rep", by.y = "user_genome", all.x = TRUE)
flags <- merge(flags, freq, by = "rep", all.x = TRUE)
flags$flag_blank_genus <- flags$rep %in% flag_blank
flags$flag_frequency   <- flags$rep %in% flag_freq
flags$flag_any         <- flags$flag_blank_genus | flags$flag_frequency
flags$in_guild         <- flags$rep %in% guild$union
write.table(flags, file.path(OUT, "rep_contamination_flags.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- guild abundance with and without flagged reps ----------------------------------------------
tot <- tapply(ab$abund, ab$sample, sum)                        # mapped % of reads (whole = 100)
gab <- function(set, excl = character()) {
  d <- ab[ab$rep %in% setdiff(set, excl), ]
  v <- tapply(d$abund, d$sample, sum); v[setdiff(names(tot), names(v))] <- 0; v[names(tot)]
}
res <- data.frame(sample = names(tot), mapped_pct = as.numeric(tot))
for (g in names(guild)) {
  res[[paste0(g, "_whole")]]          <- as.numeric(gab(guild[[g]]))
  res[[paste0(g, "_whole_clean")]]    <- as.numeric(gab(guild[[g]], flags$rep[flags$flag_any]))
}
mapped_clean <- as.numeric(tot - gab(flags$rep[flags$flag_any]))
for (g in names(guild)) {
  res[[paste0(g, "_mapped")]]       <- 100 * res[[paste0(g, "_whole")]] / res$mapped_pct
  res[[paste0(g, "_mapped_clean")]] <- 100 * res[[paste0(g, "_whole_clean")]] / mapped_clean
}
res <- res[res$mapped_pct > 0, ]
write.table(res, file.path(OUT, "guild_abundance_with_without_flagged.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- summary -------------------------------------------------------------------------------------
q <- function(x) { x <- x[is.finite(x)]; sprintf("%.2f (IQR %.2f-%.2f; n=%d)", median(x), quantile(x, .25), quantile(x, .75), length(x)) }
sink(file.path(OUT, "SUMMARY.txt"))
cat("Life on air: contamination screen of the ANI95 representatives\n")
cat(sprintf("representatives with CoverM abundance: %d\n", length(reps)))
cat(sprintf("blank-enriched candidate_contaminant features: %d; distinct genera: %d\n", nrow(ctl), length(ctl_genus)))
cat(sprintf("samples with a single unambiguous DNA concentration: %d (conflicting values dropped: %d)\n", length(conc), length(conflict)))
cat(sprintf("reps testable by frequency (present in >= 10 such samples): %d\n", sum(!is.na(flags$score))))
cat(sprintf("flagged: blank genus %d, frequency %d, either %d (of which in guild %d of %d)\n",
            sum(flags$flag_blank_genus), sum(flags$flag_frequency), sum(flags$flag_any),
            sum(flags$flag_any & flags$in_guild), length(guild$union)))
cat(sprintf("samples: %d\n", nrow(res)))
for (g in names(guild)) {
  cat(sprintf("%-5s whole-community %%: all %s | clean %s\n", g, q(res[[paste0(g, "_whole")]]), q(res[[paste0(g, "_whole_clean")]])))
  cat(sprintf("%-5s mapped-community %%: all %s | clean %s\n", g, q(res[[paste0(g, "_mapped")]]), q(res[[paste0(g, "_mapped_clean")]])))
}
cat("\nflagged guild members by genus:\n")
fg <- flags[flags$flag_any & flags$in_guild, ]
print(sort(table(paste(fg$phylum, fg$genus)), decreasing = TRUE))
sink()
cat(readLines(file.path(OUT, "SUMMARY.txt")), sep = "\n")
