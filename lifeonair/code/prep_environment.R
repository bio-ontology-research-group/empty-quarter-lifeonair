#!/usr/bin/env Rscript
# Life on air (P3): compact environmental inputs for H1-H3, frozen as TSV with checksums.
# Sources (amplicon repository cache; sha256 recorded in lifeonair/results/inputs/env_sources.sha256):
#   climate_historical_1995_2024.parquet  daily T and P per site (Open-Meteo archive default "best_match",
#                                         i.e. ERA5 blended with ECMWF IFS from 2017; kept as map_bestmatch_mm)
#   inputs/era5_daily_1995_2025.tsv       pure ERA5 daily P 1995-01-01..2025-10-31 (fetch_era5_2025.py; review
#                                         round 1, R1-03): primary source of MAP (1995-2024) and antecedent sums
#   metadata_with_rainfall.parquet        per amplicon sample: trip_date (sampling date)
#   ndvi_by_sample.parquet, smap_by_sample.parquet   MODIS NDVI (60 d) and SMAP soil moisture
#   srtm_topo.parquet                     site lat, lon, elevation
# Metagenome sample ids equal the amplicon ids without the "eNNNN_" prefix.
# Outputs: inputs/env_site.tsv (site climate normals), inputs/env_site_trip.tsv (sampling date,
# antecedent precipitation, NDVI/SMAP per site x trip).
#
#   Rscript lifeonair/code/prep_environment.R

suppressMessages(library(arrow))
C   <- path.expand("~/Documents/papers/empty-quarter-amplicon/repository/cache")
IN  <- "lifeonair/results/inputs"
rd  <- function(f) read.delim(f, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)

clim <- as.data.frame(read_parquet(file.path(C, "climate_historical_1995_2024.parquet")))
clim$date <- as.Date(clim$date)
clim$year <- as.integer(format(clim$date, "%Y"))
ann <- aggregate(P ~ site + year, data = clim, FUN = sum)
site <- data.frame(site = sort(unique(clim$site)))
site$map_bestmatch_mm <- tapply(ann$P, ann$site, mean)[as.character(site$site)]
e5 <- rd(file.path(IN, "era5_daily_1995_2025.tsv")); e5$date <- as.Date(e5$date)
e5$year <- as.integer(format(e5$date, "%Y"))
ann5 <- aggregate(P ~ site + year, data = e5[e5$year <= 2024, ], FUN = sum)
site$map_era5_mm  <- tapply(ann5$P, ann5$site, mean)[as.character(site$site)]
site$mat_era5_c   <- tapply(clim$T, clim$site, mean)[as.character(site$site)]
site$rainy_days_per_yr <- tapply(clim$P >= 1, clim$site, sum)[as.character(site$site)] / length(unique(clim$year))
topo <- as.data.frame(read_parquet(file.path(C, "srtm_topo.parquet")))
site <- merge(site, topo, by = "site", all.x = TRUE)
# ERA5 grid cell shared by sites (identical daily series), for grid-aware resampling (review round 2, R2-16)
sig <- tapply(e5$P, e5$site, function(p) paste(round(p, 2), collapse = ","))
site$era5_cell <- as.integer(factor(sig))[match(site$site, as.integer(names(sig)))]

md <- as.data.frame(read_parquet(file.path(C, "metadata_with_rainfall.parquet")))
md$sample_id <- sub("^e[0-9]+_", "", md$sample)
md$trip_date <- as.Date(md$trip_date)
nd <- as.data.frame(read_parquet(file.path(C, "ndvi_by_sample.parquet")))
sm <- as.data.frame(read_parquet(file.path(C, "smap_by_sample.parquet")))
md <- merge(md, nd, by = "sample", all.x = TRUE); md <- merge(md, sm, by = "sample", all.x = TRUE)

# antecedent precipitation over 90 and 365 days before sampling (NA past the record end); pure ERA5 is
# primary, the best_match record (ends 2024) is kept for sensitivity
ante <- function(s, d, days, cl = e5) {
  x <- cl[cl$site == s & cl$date < d & cl$date >= d - days, ]
  if (d - 1 > max(cl$date)) return(NA_real_)
  sum(x$P)
}
# per site x trip: sampling date (earliest amplicon date in that cell) and site-level remote sensing
st <- aggregate(trip_date ~ site + trip, data = md, FUN = min)
st$trip_date <- as.Date(st$trip_date, origin = "1970-01-01")
rs <- aggregate(cbind(ndvi_W60d, evi_W60d, smap_sm_W30d) ~ site + trip, data = md, FUN = median, na.action = na.pass)
st <- merge(st, rs, by = c("site", "trip"), all.x = TRUE)
# antecedent ERA5 precipitation over 90 and 365 days before sampling (NA past the record end)
st$P90_era5  <- mapply(ante, st$site, st$trip_date, 90)
st$P365_era5 <- mapply(ante, st$site, st$trip_date, 365)
st$P365_bestmatch <- mapply(function(a, b) ante(a, b, 365, clim), st$site, st$trip_date)
# independent product for the temporal models (review round 2, R2-15): NASA POWER PRECTOTCORR (MERRA-2,
# GPCP-corrected), same visits and windows (fetch_power_daily.py)
pw <- rd(file.path(IN, "power_daily_2021_2025.tsv")); pw$date <- as.Date(pw$date); pw <- pw[is.finite(pw$P), ]
st$P90_power  <- mapply(function(a, b) ante(a, b, 90, pw), st$site, st$trip_date)
st$P365_power <- mapply(function(a, b) ante(a, b, 365, pw), st$site, st$trip_date)
st$trip <- paste0("trip", st$trip)

write.table(site, file.path(IN, "env_site.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(st, file.path(IN, "env_site_trip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
src <- c(file.path(IN, "era5_daily_1995_2025.tsv"), file.path(IN, "power_daily_2021_2025.tsv"), file.path(C, c("climate_historical_1995_2024.parquet", "metadata_with_rainfall.parquet",
                      "ndvi_by_sample.parquet", "smap_by_sample.parquet", "srtm_topo.parquet")))
writeLines(vapply(src, function(f) system2("sha256sum", f, stdout = TRUE), ""), file.path(IN, "env_sources.sha256"))
mg <- rd(file.path(IN, "samples.tsv"))
cov <- sum(paste(mg$site, mg$trip) %in% paste(st$site, st$trip))
cat(sprintf("sites %d; site x trip cells %d (P365 available %d); metagenome samples covered %d of %d\n",
            nrow(site), nrow(st), sum(!is.na(st$P365_era5)), cov, nrow(mg)))
