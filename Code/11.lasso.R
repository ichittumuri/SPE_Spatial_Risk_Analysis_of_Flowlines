# =============================================================================
# 1) Setup & Libraries
# =============================================================================
library(dplyr)
library(ggplot2)
library(glmnet)
library(spdep)   # for Moran's I later

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")
source("make_split.R")

# =============================================================================
# 2) Load balanced data
# =============================================================================
df_balanced <- read.csv("df_balanced.csv")

# make sure risk is numeric 0/1
df_balanced$risk <- as.integer(df_balanced$risk)

# quick check
df_balanced %>%
  count(risk) %>%
  mutate(prop = n / sum(n)) %>%
  print()

# =============================================================================
# 3) Design matrix for LASSO
# =============================================================================
form <- risk ~ status + flowline_action + location_type + fluid + material +
  diameter_in + length_ft + max_operating_pressure + elevation + line_age_yr

y <- df_balanced$risk
X <- model.matrix(form, data = df_balanced)[, -1]  # drop intercept

# =============================================================================
# 4) CV-LASSO
# =============================================================================
set.seed(123)
cv_lasso <- cv.glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha  = 1,
  nfolds = 10
)

lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se

cat("\nSelected lambdas:\n",
    "  lambda.min =", signif(lambda_min, 3), "\n",
    "  lambda.1se =", signif(lambda_1se, 3), "\n")

# =============================================================================
# 5) Extract nonzero features at lambda.1se
# =============================================================================
lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
nz_idx      <- which(as.numeric(lasso_coefs) != 0)
sel_vars    <- setdiff(rownames(lasso_coefs)[nz_idx], "(Intercept)")

cat("\nNonzero features at lambda.1se:\n")
print(sel_vars)
