# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(caret)
library(fields)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

# =============================================================================
# 2) Sources & Data
# =============================================================================
source("ZlogisticSmoother.R")       
source("ZLambdaGridSearch.R")
source("make_split.R")

df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

split_obj  <- make_train_test_split(df_balanced)
train_data <- split_obj$train
test_data  <- split_obj$test

# =============================================================================
# 4) Build X and Z matrices FOR EACH SPLIT
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

train_risk_model   <- train_data$risk
train_coord_model  <- as.matrix(train_data[, c("lon","lat")])

glm_sel <- glm(risk ~ . -1, data = Z_df_train, family = binomial(link = "logit"))
betaHat <- coef(glm_sel)

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

test_risk_model   <- test_data$risk
test_coord_model  <- as.matrix(test_data[, c("lon","lat")])

# =============================================================================
# 7) Lambda search + final fit with Z and betaHat
# =============================================================================
profile <- LambdaGridSearch(
  train_coord_model,
  train_risk_model,
  Z       = Z_train,
  betaHat = betaHat
)

bestLambda <- profile$bestLambda

MLEFit <- logisticSmoother(
  train_coord_model,
  train_risk_model,
  lambda  = bestLambda,
  Z       = Z_train,
  betaHat = betaHat
)

# =============================================================================
# 7) Predict from fitted model (with Z)
# =============================================================================

# OLD CODE, incorrect predictions 
# nu_all_mle <- predict(MLEFit, test_coord_model, Z = Z_test)
# prob_all_mle <- plogis(nu_all_mle)

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
# 8) Total log-likelihood
# =============================================================================
nu_hat <- nu_all_mle
p_hat  <- plogis(nu_hat)
eps <- 1e-15
p_hat <- pmin(pmax(p_hat, eps), 1 - eps)

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
    values = c("0" = "#42a5f5", "1" = "#e53935"),  # blue for no spill, red for spill
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
pred_title <- paste0(
  "Predicted Spill Risk – Spatial + Covarites \n",
  "Log-likelihood = ", round(loglik, 2)
)

ggplot(test_data, aes(x = lon, y = lat)) +
  geom_point(aes(color = predicted_prob), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish,
    begin  = 0.12,
    end    = 0.9
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = pred_title,
    x     = "Longitude",
    y     = "Latitude"
  )

# =============================================================================
# Save standardized prediction table (Spatial + Z)
# =============================================================================
spatial_Z_preds <- data.frame(
  model = "SpatialPlusZ",
  idx   = seq_len(nrow(test_data)),
  lon   = as.numeric(test_data$lon),
  lat   = as.numeric(test_data$lat),
  risk  = as.integer(test_data$risk),
  pred  = as.numeric(prob_all_mle)
)

write_csv(spatial_Z_preds, "predictions_spatial_plus_z.csv")

