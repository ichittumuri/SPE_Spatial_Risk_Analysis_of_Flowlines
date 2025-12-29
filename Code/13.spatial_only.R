# New cell= command, option, I or alt, windows symbol, I 
# Run = windows symbol, enter 
# Run rscript = Command + Option + R

# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(fields)
library(caret)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")
source("OGlogisticSmoother.R")
source("OGLambdaGridSearch.R")
source("make_split.R")

# =============================================================================
# 2) Load Data & Train/Test Split
# =============================================================================
df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

split_obj  <- make_train_test_split(df_balanced)
train_data <- split_obj$train
test_data  <- split_obj$test

# =============================================================================
# 6) Smoother fit + lambda search
# =============================================================================
train_coord_model <- as.matrix(train_data[, c("lon","lat")])
train_risk_model  <- train_data$risk

profile    <- LambdaGridSearch(train_coord_model, train_risk_model)
bestLambda <- profile$bestLambda

train_MLEFit <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda = bestLambda
)

test_coord_model <- as.matrix(test_data[, c("lon","lat")])
test_risk_model  <- test_data$risk

# =============================================================================
# 7) Predict from fitted model
# =============================================================================
nu_all_mle   <- predict(train_MLEFit, test_coord_model)
nu_all_mle   <- as.numeric(nu_all_mle)               # drop 1-col matrix -> vector
prob_all_mle <- plogis(nu_all_mle)                   # same as exp(nu)/(1+exp(nu))
prob_all_mle <- as.numeric(prob_all_mle)             # ensure plain numeric

test_data <- test_data %>%
  mutate(predicted_prob = prob_all_mle)

# =============================================================================
# 8) Total log-likelihood
# =============================================================================
nu_hat <- predict(train_MLEFit, test_coord_model)
p_hat  <- exp(nu_hat) / (1 + exp(nu_hat))  # same as plogis(nu_hat)
risk_model <- test_data$risk

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

loglik <- sum(loglik_i)

cat("Total log-likelihood =", round(loglik, 2), "\n")

# =============================================================================
# 9) Observed Risk plot (categorical 0/1)
# =============================================================================
ggplot() +
  geom_point(
    data = test_data %>% filter(risk == 0),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  geom_point(
    data = test_data %>% filter(risk == 1),
    aes(x = lon, y = lat, color = factor(risk)),
    size = 1.7, alpha = 0.8
  ) +
  scale_color_manual(
    values = c("0" = "#42a5f5", "1" = "#e53935"),  # blue vs red
    name   = "Observed Risk",
    labels = c("0" = "No Spill", "1" = "Spill")
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Observed Spill Risk",
    x = "Longitude",
    y = "Latitude"
  )

# =============================================================================
# 10) Predicted Risk plot (Turbo continuous scale)
# =============================================================================
# Predicted Map with Turbo scale + log-likelihood in title
pred_title <- paste0(
  "Predicted Spill Risk – Spatial Only\n",
  "Log-likelihood = ", round(loglik, 2)
)

ggplot(test_data, aes(x = lon, y = lat)) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = pred_title,
    x     = "Longitude",
    y     = "Latitude"
  )

# =============================================================================
# Save standardized prediction table (Spatial Only)
# =============================================================================
spatial_only_preds <- data.frame(
  model = "SpatialOnly",
  idx   = seq_len(nrow(test_data)),
  lon   = as.numeric(test_data$lon),
  lat   = as.numeric(test_data$lat),
  risk  = as.integer(test_data$risk),
  pred  = as.numeric(prob_all_mle)
)

write_csv(spatial_only_preds, "predictions_spatial_only.csv")
