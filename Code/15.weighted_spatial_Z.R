# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(caret)
library(fields)
library(readr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================
source("WeightsZlogisticSmoother.R")   # logisticSmoother with weights arg
source("WeightsZLambdaGridSearch.R")   # LambdaGridSearch with weights + doPlot arg
source("make_split.R")
source("SigmaGridSearch.R")   # <- sigma grid search function

df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

split_obj  <- make_train_test_split(df_balanced)
train_data <- split_obj$train
test_data  <- split_obj$test

# =============================================================================
# 3) Build X and Z matrices FOR EACH SPLIT
# =============================================================================
to.factor <- c("status","flowline_action","location_type","fluid","material")

## Training
train_data <- train_data %>% mutate(across(all_of(to.factor), as.factor))

X_train <- model.matrix(
  ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + elevation + line_age_yr,
  data = train_data
)[, -1]

Z_train <- X_train[, c("fluidOther", "elevation")]
Z_df_train <- as.data.frame(Z_train)
Z_df_train$risk <- train_data$risk

train_risk_model  <- train_data$risk
train_coord_model <- as.matrix(train_data[, c("lon","lat")])

## Testing
test_data <- test_data %>% mutate(across(all_of(to.factor), as.factor))

X_test <- model.matrix(
  ~ status + flowline_action + location_type + fluid + material +
    diameter_in + length_ft + elevation + line_age_yr,
  data = test_data
)[, -1]

Z_test <- X_test[, c("fluidOther", "elevation")]
Z_df_test <- as.data.frame(Z_test)
Z_df_test$risk <- test_data$risk

test_risk_model  <- test_data$risk
test_coord_model <- as.matrix(test_data[, c("lon","lat")])

# =============================================================================
# 4) Run sigma grid search (e.g., 5 to 100 by 5) for Spatial + Z + weights
# =============================================================================
sigma_grid <- seq(5, 100, by = 5)

sigma_fit     <- SigmaGridSearch_spatialZ(sigma_grid)
results_sigma <- sigma_fit$results_df
best_sigma    <- sigma_fit$best_sigma

print(results_sigma)
cat("Using best sigma =", best_sigma, "\n")

# Place label slightly below the max log-likelihood
label_y <- max(results_sigma$loglik) - 0.25

ggplot(results_sigma, aes(x = sigma, y = loglik)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = best_sigma,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = best_sigma,
    y = label_y,
    label = paste("best σ =", best_sigma),
    hjust = -0.1,
    vjust = 0.5,
    size = 3
  ) +
  theme_minimal() +
  labs(
    title = "Sigma grid search (Spatial + Z + distance-weights)",
    x     = expression(sigma),
    y     = "Test log-likelihood"
  )

# =============================================================================
# 5) Final model: use best_sigma, refit Spatial + Z + weights
# =============================================================================
sigma_dist <- best_sigma

d_tr <- train_data$match_distance_m

w_train <- ifelse(
  train_data$risk == 1 & !is.na(d_tr),
  exp(-d_tr / sigma_dist),
  1
)

summary(w_train)

# GLM for fixed Z-effects (weighted)
glm_sel <- glm(
  risk ~ . -1,
  data    = Z_df_train,
  family  = binomial(link = "logit"),
  weights = w_train
)

betaHat <- coef(glm_sel)

# Lambda search + final fit with Z, betaHat, and weights
# (Here doPlot = TRUE is fine: you get ONE λ-profile plot for the chosen σ)
profile <- LambdaGridSearch(
  coord_model = train_coord_model,
  risk_model  = train_risk_model,
  Z           = Z_train,
  betaHat     = betaHat,
  weights     = w_train,
  doPlot      = TRUE
)

bestLambda <- profile$bestLambda

MLEFit <- logisticSmoother(
  s       = train_coord_model,
  y       = train_risk_model,
  lambda  = bestLambda,
  Z       = Z_train,
  betaHat = betaHat,
  weights = w_train
)

# =============================================================================
# 6) Predict from fitted model (Spatial + Z + weights)
# =============================================================================
# spatial part at TEST locations
f_hat_test <- as.numeric(predict(MLEFit, test_coord_model))

# fixed-effects (GLM) part at TEST locations
eta_glm_test <- as.vector(Z_test %*% betaHat)

# full linear predictor: eta = Zβ + f(s)
nu_all_mle <- eta_glm_test + f_hat_test

# probabilities
prob_all_mle <- plogis(nu_all_mle)

test_data <- test_data %>%
  mutate(predicted_prob = prob_all_mle)

# =============================================================================
# 7) Total test log-likelihood for final model
# =============================================================================
eps   <- 1e-15
p_hat <- pmin(pmax(prob_all_mle, eps), 1 - eps)
risk_model <- test_data$risk

loglik_i <- risk_model * log(p_hat) +
  (1 - risk_model) * log(1 - p_hat)

loglik_weighted_spatialZ <- sum(loglik_i)

cat("Total log-likelihood (Spatial + Z + distance-weights, σ =", sigma_dist, ") =",
    round(loglik_weighted_spatialZ, 2), "\n")

# =============================================================================
# 8) Observed Risk plot (categorical 0/1)
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
    values = c("0" = "#42a5f5", "1" = "#e53935"),
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
# 9) Predicted Risk plot (Turbo continuous scale)
# =============================================================================
pred_title <- paste0(
  "Predicted Spill Risk – Spatial + Z + distance-weights\n",
  "σ = ", sigma_dist,
  ", Log-likelihood = ", round(loglik_weighted_spatialZ, 2)
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
# Save standardized prediction table (Spatial + Z + weights)
# =============================================================================
spatial_Z_weighted_preds <- data.frame(
  model = "SpatialPlusZ_weighted",
  idx   = seq_len(nrow(test_data)),
  lon   = as.numeric(test_data$lon),
  lat   = as.numeric(test_data$lat),
  risk  = as.integer(test_data$risk),
  pred  = as.numeric(prob_all_mle)
)

write_csv(spatial_Z_weighted_preds, "predictions_spatial_plus_z_weighted.csv")
