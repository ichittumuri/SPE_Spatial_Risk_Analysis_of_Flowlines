# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(dplyr)
library(caret)
library(readr)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")
source("make_split.R")

# =============================================================================
# 2) Load Data & Train/Test Split
# =============================================================================
df_balanced <- read.csv("df_balanced.csv", stringsAsFactors = FALSE)

split_obj  <- make_train_test_split(df_balanced)
train_data <- split_obj$train
test_data  <- split_obj$test

# =============================================================================
# 2b) Save train/test split for reuse in other scripts
# =============================================================================
readr::write_csv(train_data, "train_data_split.csv")
readr::write_csv(test_data, "test_data_split.csv")

# =============================================================================
# 3) One-hot encode categorical predictors
# =============================================================================
dummies <- dummyVars(risk ~ ., data = train_data, fullRank = TRUE)

X_train <- predict(dummies, newdata = train_data) %>% as.data.frame()
X_test  <- predict(dummies, newdata = test_data) %>% as.data.frame()

y_train <- train_data$risk
y_test  <- test_data$risk

# Use the two selected variables
Z_train <- X_train[, c("fluidOther", "elevation")]
Z_test  <- X_test[,  c("fluidOther", "elevation")]

train_sel <- cbind.data.frame(risk = y_train, Z_train)
test_sel  <- cbind.data.frame(risk = y_test,  Z_test)

# =============================================================================
# 4) Fit Logistic Regression (GLM)
# =============================================================================
glm_base <- glm(risk ~ fluidOther + elevation,
                data = train_sel,
                family = binomial())

# =============================================================================
# 5) Predict on TEST + Compute Log-likelihood
# =============================================================================
p_hat <- predict(glm_base, newdata = test_sel, type = "response")
# Safety clamp
eps <- 1e-15
p_hat <- pmin(pmax(p_hat, eps), 1 - eps)

y_true <- test_sel$risk

loglik_i <- y_true * log(p_hat) + (1 - y_true) * log(1 - p_hat)
loglik   <- sum(loglik_i)

cat("\nTest Log-likelihood:", round(loglik, 4), "\n")

# Build test_out (required)
test_out <- test_sel %>%
  mutate(pred_prob = p_hat) %>%
  bind_cols(test_data %>% dplyr::select(lon, lat))

# Standardized output for CSV
glm_preds <- test_out %>%
  transmute(
    model = "GLM",
    idx = row_number(),
    lon, lat,
    risk = as.integer(risk),
    pred = pred_prob
  )

# =============================================================================
# 6) Visualization (Turbo color scheme)
# =============================================================================
ggplot() +
  geom_point(data = glm_preds %>% filter(risk == 0),
             aes(x = lon, y = lat, color = factor(risk)),
             size = 1.7, alpha = 0.8) +
  geom_point(data = glm_preds %>% filter(risk == 1),
             aes(x = lon, y = lat, color = factor(risk)),
             size = 1.7, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#42a5f5", "1" = "#e53935"),  # blue vs red
                     name = "Observed Risk",
                     labels = c("No Spill", "Spill")) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Observed Spill Risk (TEST set)",
       x = "Longitude", y = "Latitude")

# Predicted Map with Turbo scale + log-likelihood in title
pred_title <- paste0(
  "Predicted Spill Risk – GLM (TEST set)\n",
  "Log-likelihood = ", round(loglik, 2)
)

ggplot(glm_preds %>% arrange(pred),
       aes(x = lon, y = lat)) +
  geom_point(aes(color = pred), size = 2) +
  scale_color_viridis_c(
    option = "turbo",
    name   = "Predicted Risk",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    oob    = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  labs(title = pred_title,
       x = "Longitude", y = "Latitude")

# =============================================================================
# 6b) Save standardized prediction table for GLM
# =============================================================================
glm_preds <- test_out %>%
  transmute(
    model = "GLM",
    idx   = row_number(),
    lon, lat,
    risk = as.integer(risk),
    pred = pred_prob
  )

write_csv(glm_preds, "predictions_glm_base.csv")




