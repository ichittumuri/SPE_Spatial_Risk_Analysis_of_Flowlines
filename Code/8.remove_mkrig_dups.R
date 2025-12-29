# =========================
# Script 3: clean points (flowlines + spills), drop dupes, drop overlaps (keep spills), combine & save
# =========================
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
# (purrr optional) library(purrr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# -------------------------
# 0) Filepaths & knobs
# -------------------------
FLOWLINES_FP         <- "flowline_points_50m_dedup.geojson"
SPILLS_FP            <- "spills_w_flowline_attributes.geojson"

FLOWLINES_CLEAN_FP   <- "flowlines_dedup_nonNA.geojson"
SPILLS_CLEAN_FP      <- "spills_dedup_nonNA.geojson"

FLOWLINES_NO_OVER_FP <- "flowlines_no_overlap.geojson"
COMBINED_RAW_FP      <- "combined_flowlines_spills_raw.geojson"
COMBINED_FP          <- "combined_flowlines_spills.geojson"   # final

# columns to require non-NA (adjust as needed)
cols_to_check <- c("flowline_action","fluid","material")

round_dp <- 6  # coordinate rounding precision for "R-level" dup & overlap checks

# -------------------------
# helpers
# -------------------------
# collapse rows of a numeric matrix to strings (mKrig-style keys)
cat.matrix <- function(M) apply(M, 1, function(r) paste(r, collapse = " "))

# force valid POINT geoms; drop non-points/empties
force_points <- function(gdf, label = "layer") {
  if (!inherits(gdf, "sf")) stop(label, ": not an sf object")
  gdf <- st_make_valid(gdf)
  # extract POINTS from collections; if still not POINT, try cast
  gt <- unique(as.character(st_geometry_type(gdf)))
  if (!all(gt %in% c("POINT"))) {
    pts <- trySuppressWarnings(st_collection_extract(gdf, "POINT"))
    if (nrow(pts) == 0) {
      # some drivers store points as MULTIPOINT; attempt st_cast
      pts <- trySuppressWarnings(st_cast(gdf, "POINT", warn = FALSE))
    }
    gdf <- pts
  }
  # drop empties
  empties <- st_is_empty(gdf)
  if (any(empties)) gdf <- gdf[!empties, , drop = FALSE]
  if (nrow(gdf) == 0) stop(label, ": no POINT geometries after cleaning")
  gdf
}

trySuppressWarnings <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

# safe coordinate extractor (X,Y)
get_xy <- function(gdf) {
  cm <- st_coordinates(gdf)
  cm[, c("X","Y"), drop = FALSE]
}

# Audit + clean one layer, then save
clean_and_save <- function(gdf, out_fp, label, cols_to_check, round_dp = 6) {
  message("\n========== ", toupper(label), " ==========")
  gdf <- force_points(gdf, label)
  coords <- get_xy(gdf)
  
  # R-level dups (rounded)
  coord_rounded <- round(coords, round_dp)
  dupe_rows <- duplicated(as.data.frame(coord_rounded))
  cat("R-level duplicates (rounded):", sum(dupe_rows), "\n")
  
  # mKrig-style dups
  coord_strings <- cat.matrix(coords)
  dups_mkrig <- duplicated(coord_strings)
  cat("mKrig-level duplicates:", sum(dups_mkrig), "\n")
  
  # Drop mKrig duplicates
  gdf_clean <- gdf[!dups_mkrig, , drop = FALSE]
  cat("Rows after mKrig-dup removal:", nrow(gdf_clean), "\n")
  
  # NA audit on original (just info)
  na_counts <- colSums(is.na(st_drop_geometry(gdf)))
  cat("\n--- NA counts per column (original) ---\n")
  print(na_counts)
  
  # Drop NAs on selected columns (only if columns exist)
  cols_use <- intersect(cols_to_check, names(gdf_clean))
  if (length(cols_use)) {
    rows_before <- nrow(gdf_clean)
    gdf_clean <- gdf_clean[
      complete.cases(st_drop_geometry(gdf_clean)[, cols_use, drop = FALSE]),
    ]
    cat("Rows before NA drop:", rows_before, " | after:", nrow(gdf_clean), "\n")
  } else {
    cat("(No overlap between cols_to_check and columns in ", label, "; skipping NA drop)\n", sep = "")
  }
  
  st_write(gdf_clean, out_fp, delete_dsn = TRUE, quiet = TRUE)
  cat("Saved:", out_fp, "\n")
  
  gdf_clean
}

# -------------------------
# 1) Load
# -------------------------
flowlines_raw <- st_read(FLOWLINES_FP, quiet = TRUE)
spills_raw    <- st_read(SPILLS_FP,    quiet = TRUE)

# -------------------------
# 2) Clean each layer & save
# -------------------------
flowlines_clean <- clean_and_save(flowlines_raw, FLOWLINES_CLEAN_FP, "flowlines", cols_to_check, round_dp)
spills_clean    <- clean_and_save(spills_raw,    SPILLS_CLEAN_FP,    "spills",    cols_to_check, round_dp)

# -------------------------
# 3) Remove overlaps (keep spills)
#     Define overlaps by rounded (X,Y) key at round_dp
# -------------------------
fl_xy <- get_xy(flowlines_clean)
sp_xy <- get_xy(spills_clean)

fl_key <- apply(round(fl_xy, round_dp), 1, paste, collapse = ",")
sp_key <- apply(round(sp_xy, round_dp), 1, paste, collapse = ",")

common_keys <- intersect(fl_key, sp_key)

cat("\n=== OVERLAP AUDIT (rounded ", round_dp, " d.p.) ===\n", sep = "")
cat("Flowlines BEFORE:", nrow(flowlines_clean), "\n")
cat("Spills:          ", nrow(spills_clean),    "\n")
cat("Overlapping keys:", length(common_keys),   "\n")

flowlines_no_overlap <- flowlines_clean[!(fl_key %in% common_keys), , drop = FALSE]
cat("Flowlines AFTER removing overlaps:", nrow(flowlines_no_overlap), "\n")

st_write(flowlines_no_overlap, FLOWLINES_NO_OVER_FP, delete_dsn = TRUE, quiet = TRUE)
cat("Saved:", FLOWLINES_NO_OVER_FP, "\n")

# -------------------------
# 4) Combine (union columns)
# -------------------------
# unify CRS
if (st_crs(flowlines_no_overlap) != st_crs(spills_clean)) {
  spills_clean <- st_transform(spills_clean, st_crs(flowlines_no_overlap))
}

# strip geometry
fl_df <- st_drop_geometry(flowlines_no_overlap)
sp_df <- st_drop_geometry(spills_clean)

fl_g  <- st_geometry(flowlines_no_overlap)
sp_g  <- st_geometry(spills_clean)

combined_df <- dplyr::bind_rows(fl_df, sp_df)  # unions columns
combined_raw <- st_sf(combined_df,
                      geometry = c(fl_g, sp_g),
                      crs = st_crs(flowlines_no_overlap))

st_write(combined_raw, COMBINED_RAW_FP, delete_dsn = TRUE, quiet = TRUE)
cat("Rows combined (raw):", nrow(combined_raw), "\n")
cat("Saved:", COMBINED_RAW_FP, "\n")

# -------------------------
# 5) Final combined duplicate check + save final
# -------------------------
cat("\n========== COMBINED (FINAL AUDIT) ==========\n")
combined <- force_points(combined_raw, "combined")
coords_c <- get_xy(combined)

# R-level dups (rounded)
coord_rounded_c <- round(coords_c, round_dp)
dupe_rows_c <- duplicated(as.data.frame(coord_rounded_c))
cat("R-level duplicates (rounded):", sum(dupe_rows_c), "\n")

# mKrig-style dups
coord_strings_c <- cat.matrix(coords_c)
dups_mkrig_c <- duplicated(coord_strings_c)
cat("mKrig-level duplicates:", sum(dups_mkrig_c), "\n")

combined_clean <- combined[!dups_mkrig_c, , drop = FALSE]
cat("Final rows after mKrig-dup removal:", nrow(combined_clean), "\n")

st_write(combined_clean, COMBINED_FP, delete_dsn = TRUE, quiet = TRUE)
cat("Saved FINAL:", COMBINED_FP, "\n")

