# ============================================================
# 1) Setup
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# ============================================================
# 2) Load combined predictions
# ============================================================
all_preds <- read_csv("all_predictions.csv")

# ============================================================
# 3) Helper: deviance residuals from y and p
# ============================================================
deviance_residuals <- function(y, p) {
  eps <- 1e-15                     # avoid log(0)
  p <- pmin(pmax(p, eps), 1 - eps)
  ifelse(y == 1,
         sqrt(2 * abs(log(p))),
         -sqrt(2 * abs(log(1 - p))))
}

# ============================================================
# 4) Compute deviance residuals + linear predictor for GLM
# ============================================================
all_preds <- all_preds %>%
  mutate(
    devres = deviance_residuals(risk, pred_SpatialPlusZ_wt), # pred_GLM, pred_SpatialOnly, pred_SpatialPlusZ, 
    # pred_SpatialPlusZ_wt
    eta    = log(pred_SpatialPlusZ_wt / (1 - pred_SpatialPlusZ_wt))   # = qlogis(pred_GLM)
  )

# ============================================================
# 5) Plot: deviance residuals vs linear predictor (GLM)
# ============================================================
ggplot(all_preds, aes(x = eta, y = devres)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess") +
  labs(
    x = "Linear predictor",
    y = "Deviance residuals",
    title = "Residuals Plot"
  )

# ============================================================
# 6) Calibration Curve
# ============================================================
all_preds %>%
  ggplot(aes(x = pred_SpatialPlusZ_wt, y = risk)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess") +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash") +
  coord_cartesian(ylim = c(-0.5, 1)) +   # ← extend downward
  labs(
    title = "Calibration Curve",
    x = "Predicted Probability",
    y = "Observed Risk"
  )


# ============================================================
# 7) Brier score
# ============================================================

brier_score <- function(y, p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  mean((y - p)^2)
}

brier_result <- with(all_preds, brier_score(risk, pred_SpatialPlusZ_wt))
print(brier_result)

# ============================================================
# 8) Log-likelihood
# ============================================================

loglik_binom <- function(y, p) {
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  
  sum(y * log(p) + (1 - y) * log(1 - p))
}

loglik_result <- with(all_preds, loglik_binom(risk, pred_SpatialPlusZ_wt))
print(loglik_result)


# ============================================================
# 9) Threshold sweep for SpatialPlusZ_wt
# ============================================================

# sequence of thresholds to try
thresh <- seq(0.01, 0.5, 0.01)

Sensitivity <- numeric(length(thresh))
Specificity <- numeric(length(thresh))
PPV <- numeric(length(thresh))
NPV <- numeric(length(thresh))
F1 <- numeric(length(thresh))

for (j in seq_along(thresh)) {
  # classify as 1 (spill) if predicted prob > threshold
  pp <- ifelse(all_preds$pred_SpatialPlusZ_wt > thresh[j], 1, 0)
  
  # confusion matrix: rows = true (risk), cols = predicted (pp)
  xx <- xtabs(~ risk + pp, all_preds)
  
  # assuming levels: row1 = 0, row2 = 1; col1 = 0, col2 = 1
  TN <- xx[1, 1]
  FP <- xx[1, 2]
  FN <- xx[2, 1]
  TP <- xx[2, 2]
  
  Specificity[j] <- TN / (TN + FP)          # true negative rate
  Sensitivity[j] <- TP / (TP + FN)          # true positive rate (recall)
  PPV[j]        <- TP / (TP + FP)           # precision
  NPV[j]        <- TN / (TN + FN)
  F1[j]         <- TP / (TP + 0.5 * (FP + FN))
}

full_df <- data.frame(
  thresh = thresh,
  sens   = Sensitivity,
  spec   = Specificity,
  PPV    = PPV,
  NPV    = NPV,
  F1     = F1
)

full_df$best <- NA
full_df$best[which.max(F1)] <- "Best"
best_thresh <- thresh[which.max(F1)]
print(best_thresh)

# ============================================================
# Plots
# ============================================================

# Sensitivity & Specificity vs threshold
p1 <- ggplot(full_df, aes(x = thresh)) +
  geom_line(aes(y = sens, linetype = "Sensitivity")) +
  geom_line(aes(y = spec, linetype = "Specificity")) +
  geom_vline(xintercept = best_thresh, linetype = "longdash") +
  theme(legend.position = "bottom") +
  labs(
    title = "Sensitivity and Specificity Across Classification Thresholds",
    x = "Classification Threshold",
    y = "Rate (Sensitivity / Specificity)"
  )

print(p1)

# ROC-style curve (Sensitivity vs 1 - Specificity)
p2 <- ggplot(full_df, aes(x = 1 - spec, y = sens)) +
  geom_line() +
  geom_point(aes(color = best)) +
  geom_abline(slope = 1, intercept = 0, linetype = "longdash") +
  theme(legend.position = "bottom") +
  labs(
    title = "ROC Curve",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  )

print(p2)
# grid.arrange(p1, p2, ncol = 2)
