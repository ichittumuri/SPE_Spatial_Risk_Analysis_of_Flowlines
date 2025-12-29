# =========================
# Script 2: split raw by spills IDs, clean pool, impute, outlier-check, recombine
# =========================
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(FNN)
library(tibble)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# -------------------------
# 1) Load
# -------------------------
spills      <- st_read("spills_w_flowline_attributes.geojson")
matched_raw <- st_read("flowlines_matched_raw.geojson")
unmatched_raw <- st_read("flowlines_unmatched_raw.geojson")  # saved from Script 1

stopifnot("unique_id" %in% names(spills))
stopifnot(all(c("unique_id","operator_name","operator_number","flowline_id","location_id",
                "max_operating_pressure","diameter_in") %in% names(matched_raw)))

# -------------------------
# 2) Split matched_raw by spills IDs
# -------------------------
ids_kept <- unique(spills$unique_id[!is.na(spills$unique_id)])

matched_kept_raw    <- matched_raw %>% filter(unique_id %in% ids_kept)
matched_notkept_raw <- matched_raw %>% filter(!(unique_id %in% ids_kept))

# Save splits for traceability
st_write(matched_kept_raw,    "flowlines_matched_kept_raw.geojson",    delete_dsn = TRUE)
st_write(matched_notkept_raw, "flowlines_matched_notkept_raw.geojson", delete_dsn = TRUE)

cat("Matched RAW kept (by spills IDs):   ", nrow(matched_kept_raw), "\n")
cat("Matched RAW not-kept:               ", nrow(matched_notkept_raw), "\n")

# -------------------------
# 3) Clean the pool = notkept_raw + unmatched_raw
#     Drop rows with NA or zero in required key columns
# -------------------------
required_nonzero <- c("max_operating_pressure","diameter_in")
required_nonNA   <- c("operator_name","operator_number","flowline_id","location_id",
                      "max_operating_pressure","diameter_in")

pool_raw <- bind_rows(matched_notkept_raw, unmatched_raw)

# zero->NA for numeric required_nonzero, then drop rows with NA in required_nonNA
drop_bad <- function(gdf) {
  gdf2 <- gdf
  for (z in intersect(required_nonzero, names(gdf2))) {
    zvec <- suppressWarnings(as.numeric(gdf2[[z]]))
    gdf2[[z]] <- zvec
    gdf2[[z]][gdf2[[z]] == 0] <- NA_real_
  }
  gdf2 %>% filter(complete.cases(st_drop_geometry(select(., all_of(required_nonNA)))))
}

pool_clean <- drop_bad(pool_raw)

cat("Pool size (raw):                      ", nrow(pool_raw), "\n")
cat("Pool size after key-col NA/zero drop: ", nrow(pool_clean), "\n")

# -------------------------
# 4) Combine kept_raw + cleaned pool (pre-impute)
# -------------------------
combined_preimp <- bind_rows(matched_kept_raw, pool_clean)

# -------------------------
# 5) Impute remaining numerics + outlier filter (on imputed entries)
# -------------------------
# Helpers
scale_with <- function(M, center, scale) sweep(sweep(M, 2, center, `-`), 2, scale, `/`)
knn_impute_col <- function(df_num, target, k = 5, backup = c("median","mean"), ignore_zero_neighbors = TRUE) {
  backup <- match.arg(backup)
  feats <- setdiff(colnames(df_num), target)
  feat_complete <- complete.cases(df_num[, feats, drop = FALSE])
  obs <- which(!is.na(df_num[[target]]) & feat_complete)
  mis <- which( is.na(df_num[[target]]) & feat_complete)
  if (!length(mis) || !length(obs)) {
    if (anyNA(df_num[[target]])) {
      fill <- if (backup == "median") median(df_num[[target]], na.rm = TRUE) else mean(df_num[[target]], na.rm = TRUE)
      df_num[[target]][is.na(df_num[[target]])] <- fill
    }
    return(df_num[[target]])
  }
  X_obs <- as.matrix(df_num[obs, feats, drop = FALSE])
  X_mis <- as.matrix(df_num[mis, feats, drop = FALSE])
  ctr <- colMeans(X_obs)
  sds <- apply(X_obs, 2, sd); sds[!is.finite(sds) | sds == 0] <- 1
  X_obs_s <- scale_with(X_obs, ctr, sds)
  X_mis_s <- scale_with(X_mis, ctr, sds)
  k_use <- min(k, nrow(X_obs_s))
  nn <- FNN::get.knnx(data = X_obs_s, query = X_mis_s, k = k_use)$nn.index
  y_obs <- df_num[[target]][obs]
  y_nonzero <- y_obs[y_obs != 0]
  fallback <- if (ignore_zero_neighbors && length(y_nonzero)) median(y_nonzero, na.rm = TRUE) else median(y_obs, na.rm = TRUE)
  imputed_vals <- apply(nn, 1, function(idx) {
    vals <- y_obs[idx]
    if (ignore_zero_neighbors) vals <- vals[vals != 0]
    if (!length(vals) || all(!is.finite(vals))) return(fallback)
    mean(vals, na.rm = TRUE)
  })
  out <- df_num[[target]]
  out[mis] <- imputed_vals
  out[is.na(out)] <- fallback
  out
}

# Pick columns to impute (don’t touch IDs; allow physics cols to be imputed where present)
num_candidates <- intersect(c("max_operating_pressure","diameter_in","length_ft","line_age_yr"), names(combined_preimp))

# Encode cats to aid imputation if present
work <- combined_preimp
if ("material" %in% names(work)) work$material_encoded <- as.integer(factor(tidyr::replace_na(work$material, "")))
if ("fluid"    %in% names(work)) work$fluid_encoded    <- as.integer(factor(tidyr::replace_na(work$fluid, "")))

df_num <- work %>% st_drop_geometry()
keep_feats <- c(num_candidates, intersect(c("material_encoded","fluid_encoded"), names(df_num)))
df_num <- df_num %>% select(all_of(keep_feats))

# ensure numeric
for (nm in names(df_num)) df_num[[nm]] <- suppressWarnings(as.numeric(df_num[[nm]]))

df_pre <- df_num
for (tg in intersect(num_candidates, names(df_num))) {
  df_num[[tg]] <- knn_impute_col(df_num, target = tg, k = 5, backup = "median", ignore_zero_neighbors = TRUE)
}

# Outlier filter only where we imputed (per-column p1/p99 or |z|>3)
rows_to_drop <- integer(0)
for (tg in intersect(num_candidates, names(df_num))) {
  imp_idx <- which(is.na(df_pre[[tg]]) & !is.na(df_num[[tg]]))
  if (!length(imp_idx)) next
  obs <- df_num[[tg]][!is.na(df_pre[[tg]])]
  if (!length(obs)) next
  lo <- quantile(obs, 0.01, na.rm = TRUE); hi <- quantile(obs, 0.99, na.rm = TRUE)
  mu <- mean(obs, na.rm = TRUE); sg <- sd(obs, na.rm = TRUE); if (!is.finite(sg) || sg == 0) sg <- NA_real_
  val <- df_num[[tg]][imp_idx]
  bad <- (val < lo | val > hi) | if (is.na(sg)) FALSE else abs((val - mu)/sg) > 3
  rows_to_drop <- union(rows_to_drop, imp_idx[bad])
}

if (length(rows_to_drop)) {
  combined_preimp <- combined_preimp[-rows_to_drop, , drop = FALSE]
  df_num          <- df_num[-rows_to_drop, , drop = FALSE]
}

# Write imputed columns back
for (tg in intersect(num_candidates, names(df_num))) {
  combined_preimp[[tg]] <- df_num[[tg]]
}

final_combined <- combined_preimp

# -------------------------
# Ensure all spill IDs are preserved
# -------------------------
ids_final <- unique(final_combined$unique_id)
missing_from_final <- setdiff(ids_kept, ids_final)

if (length(missing_from_final)) {
  cat("\n[INFO] Adding back", length(missing_from_final),
      "missing spill IDs from matched_kept_raw before saving...\n")
  add_back <- matched_kept_raw %>%
    filter(unique_id %in% missing_from_final)
  final_combined <- bind_rows(final_combined, add_back)
} else {
  cat("\n[INFO] All spill IDs are already present in final_combined.\n")
}

# -------------------------
# Save + Checks
# -------------------------
st_write(final_combined, "combined_flowlines_final.geojson", delete_dsn = TRUE)

cat("\n--- Final Checks ---\n")
cat("Spills unique_id count:           ", length(ids_kept), "\n")
cat("Matched_kept_raw unique_id count: ", matched_kept_raw %>% st_drop_geometry() %>% pull(unique_id) %>% unique() %>% length(), "\n")
cat("Final combined rows:              ", nrow(final_combined), "\n")

