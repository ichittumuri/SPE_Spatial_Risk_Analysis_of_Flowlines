# 15.combine_predictions_wide.R
# =============================================================================
# Combine prediction tables into WIDE format:
# lon, lat, risk, pred_GLM, pred_SpatialOnly, pred_SpatialPlusZ, pred_Weighted
# =============================================================================

library(dplyr)
library(readr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# Load CSVs -----------------------------------------------------------
pred_glm       <- read_csv("predictions_glm_base.csv")
pred_spatial   <- read_csv("predictions_spatial_only.csv")
pred_spatial_Z <- read_csv("predictions_spatial_plus_z.csv")
pred_weighted  <- read_csv("predictions_spatial_plus_z_weighted.csv")

# Build WIDE dataset --------------------------------------------------
all_preds <- pred_glm %>%
  select(lon, lat, risk) %>%        # keep shared columns
  mutate(
    pred_GLM              = pred_glm$pred,
    pred_SpatialOnly      = pred_spatial$pred,
    pred_SpatialPlusZ     = pred_spatial_Z$pred,
    pred_SpatialPlusZ_wt  = pred_weighted$pred
  )

# Save ---------------------------------------------------------------
write_csv(all_preds, "all_predictions.csv")

cat("\nSaved prediction table to: all_predictions.csv\n")
