# =========================
# 0) Deps & Setup
# =========================
library(sf)
library(dplyr)
library(stringr)
library(tidyr)
library(FNN)
library(purrr)
library(tibble)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =========================
# 1) Load & normalize names
# =========================
flowlines <- st_read("clean_combined_flowlines.geojson") %>%
  mutate(operator_name = str_trim(str_to_upper(operator_name)))

spills <- st_read("clean_spills.geojson") %>%
  mutate(operator_name = str_trim(str_to_upper(operator_name)))

# =========================
# 2) Split matched / unmatched by operator_name
# =========================
common_ops <- intersect(unique(flowlines$operator_name), unique(spills$operator_name))
flowlines_matched   <- filter(flowlines, operator_name %in% common_ops)
flowlines_unmatched <- filter(flowlines, !operator_name %in% common_ops)

# Save the *raw* matched subset (pre-impute/QC)
st_write(flowlines_matched, "flowlines_matched_raw.geojson", delete_dsn = TRUE)
st_write(flowlines_unmatched, "flowlines_unmatched_raw.geojson", delete_dsn = TRUE)

# Keep only already-complete rows for unmatched (ignore geometry col)
flowlines_unmatched_clean <- flowlines_unmatched %>%
  filter(complete.cases(st_drop_geometry(.)))

# =========================
# 3) kNN imputation helper
# =========================
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

# =========================
# 4) Prepare matched subset + impute numerics
# =========================
cols_for_impute <- intersect(
  c("max_operating_pressure","diameter_in","length_ft","line_age_yr","material","fluid"),
  names(flowlines_matched)
)

zero_as_na  <- intersect(c("max_operating_pressure","diameter_in","length_ft"), cols_for_impute)
num_targets <- intersect(c("max_operating_pressure","diameter_in","length_ft","line_age_yr"), cols_for_impute)

# Build numeric frame (encode cats), drop geometry
df <- flowlines_matched[, cols_for_impute, drop = FALSE]

if ("material" %in% names(df)) df$material_encoded <- as.integer(factor(tidyr::replace_na(df$material, "")))
if ("fluid"    %in% names(df)) df$fluid_encoded    <- as.integer(factor(tidyr::replace_na(df$fluid, "")))
df <- df %>% select(-any_of(c("material","fluid"))) %>% st_drop_geometry()

# Coerce to numeric safely
df[] <- lapply(df, function(x) if (is.factor(x)) as.integer(x) else suppressWarnings(as.numeric(x)))

# Turn exact zeros to NA for selected cols (so they can be imputed)
for (zcol in zero_as_na) {
  if (zcol %in% names(df)) df[[zcol]][df[[zcol]] == 0] <- NA_real_
}

# Keep a copy to know which entries were imputed
df_pre <- df

# kNN impute selected numeric targets
for (tg in num_targets) {
  df[[tg]] <- knn_impute_col(df, target = tg, k = 5, backup = "median", ignore_zero_neighbors = TRUE)
}

# =========================
# 5) QC: drop rows with extreme imputed values (per-column)
# =========================
imp_mask <- setNames(lapply(num_targets, function(c) is.na(df_pre[[c]])), num_targets)

flags <- map(num_targets, function(c) {
  x <- df[[c]]
  imp <- imp_mask[[c]]
  obs <- x[!imp]
  if (!length(obs)) return(tibble(column=c, row_idx=integer(), value=numeric(), out_of_band=logical(), z_bad=logical()))
  
  lo <- quantile(obs, 0.01, na.rm=TRUE)
  hi <- quantile(obs, 0.99, na.rm=TRUE)
  mu <- mean(obs, na.rm=TRUE)
  sg <- sd(obs, na.rm=TRUE); if (!is.finite(sg) || sg == 0) sg <- NA_real_
  
  idx <- which(imp)
  if (!length(idx)) return(tibble(column=c, row_idx=integer(), value=numeric(), out_of_band=logical(), z_bad=logical()))
  
  val <- x[idx]
  out_of_band <- val < lo | val > hi
  z_bad <- if (is.na(sg)) rep(FALSE, length(val)) else abs((val - mu)/sg) > 3
  
  tibble(column=c, row_idx=idx, value=val, out_of_band, z_bad)
}) %>% bind_rows() %>% filter(out_of_band | z_bad)

rows_to_drop <- unique(flags$row_idx)

if (length(rows_to_drop)) {
  df <- df[-rows_to_drop, , drop = FALSE]
  flowlines_matched <- flowlines_matched[-rows_to_drop, , drop = FALSE]
}

# Write imputed values back; drop any residual NAs on matched
flowlines_matched[num_targets] <- df[, num_targets, drop = FALSE]
flowlines_matched_clean <- flowlines_matched %>%
  filter(complete.cases(st_drop_geometry(.)))

# =========================
# 6) Combine & Save
# =========================
combined_flowlines <- bind_rows(flowlines_unmatched_clean, flowlines_matched_clean)

st_write(flowlines_matched_clean,   "flowlines_matched_clean.geojson",   delete_dsn = TRUE)
st_write(flowlines_unmatched_clean, "flowlines_unmatched_clean.geojson", delete_dsn = TRUE)
st_write(combined_flowlines, "interpolated_clean_flowlines.geojson", delete_dsn = TRUE)

# --- Print counts ---
cat("Unmatched kept (NA-free):", nrow(flowlines_unmatched_clean), "\n")
cat("Matched kept (cleaned):  ", nrow(flowlines_matched_clean), "\n")
cat("Combined total:          ", nrow(combined_flowlines), "\n")

 