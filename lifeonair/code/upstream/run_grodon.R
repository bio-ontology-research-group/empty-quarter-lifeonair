suppressMessages({library(gRodon); library(Biostrings); library(parallel)})
args <- commandArgs(TRUE); files <- list.files(args[1], "\\.ffn$", full.names = TRUE)
one <- function(f) {
  s <- readDNAStringSet(f); s <- s[width(s) %% 3 == 0 & width(s) >= 240]
  heg <- grepl(" HEG$", names(s))
  rep <- sub("\\.ffn$", "", basename(f))
  if (sum(heg) < 10) return(data.frame(rep, n_cds = length(s), n_ribo = sum(heg), d_hours = NA, lower = NA, upper = NA, cub_hegs = NA))
  g <- tryCatch(predictGrowth(s, heg, mode = "partial"), error = function(e) NULL)
  if (is.null(g)) return(data.frame(rep, n_cds = length(s), n_ribo = sum(heg), d_hours = NA, lower = NA, upper = NA, cub_hegs = NA))
  data.frame(rep, n_cds = length(s), n_ribo = sum(heg), d_hours = g$d, lower = g$LowerCI, upper = g$UpperCI, cub_hegs = g$CUBHE)
}
res <- do.call(rbind, mclapply(files, one, mc.cores = as.integer(args[3])))
write.table(res, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat("[grodon] reps", nrow(res), "with estimate", sum(!is.na(res$d_hours)), "\n")
