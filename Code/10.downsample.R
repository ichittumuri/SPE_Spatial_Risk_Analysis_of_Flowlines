# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R
# example https://simonbrewer.github.io/geog5160/GEOG_5160_6160_lab04.html

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
# 1) Setup --------------------------------------------------------------------
library(sf)
library(dplyr)
library(ggplot2)
library(glmnet)
library(spdep)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================

usable_sf <- st_read("final_dataset_subset.geojson")

usable_sf[] <- lapply(usable_sf, function(col) {
  if (is.factor(col) || is.character(col)) {
    col <- gsub(" ", "_", col)      # replace spaces in the values
    col <- factor(col)              # ensure factors stay factors
  }
  col
})

# =============================================================================
# 3) Prep model data + quick class balance
# =============================================================================
coords_mat        <- st_coordinates(usable_sf)[, c("X","Y")]
usable_sf$lon     <- coords_mat[, 1]
usable_sf$lat     <- coords_mat[, 2]

df <- usable_sf |>
  st_drop_geometry() |>
  select(-unique_id) |>
  filter(
    lon >= -105.85,
    lon <= -104.79 + 0.2,
    lat >=  39.42,
    lat <=  40.68
  )

# # =============================================================================
# # 4) Use full dataset (no downsampling)
# # =============================================================================
# df_balanced <- df
# 
# cat("\n--- FULL DATASET (no downsampling) ---\n")
# df_balanced %>%
#   count(risk) %>%
#   mutate(prop = n / sum(n)) %>%
#   print()
# 
# write.csv(df_balanced, "df_balanced.csv", row.names = FALSE)

# =============================================================================
# 4) Downsample 0-risk points to ~5% spills, (stratified sampling or weighted sample???)
# =============================================================================
set.seed(123)

spills <- df %>% filter(risk == 1)
non_spills <- df %>% filter(risk == 0)

n_spills <- nrow(spills)
target_total <- ceiling(n_spills / 0.05)  # total rows for ~5% spills
n_non_spills_needed <- target_total - n_spills

non_spills_sample <- non_spills %>%
  sample_n(size = min(n_non_spills_needed, nrow(non_spills)))

df_balanced <- bind_rows(spills, non_spills_sample)

cat("\n--- NEW RATIO (~5% spills) ---\n")
df_balanced %>%
  count(risk) %>%
  mutate(prop = n / sum(n)) %>%
  print()

write.csv(df_balanced, "df_balanced.csv", row.names = FALSE)
