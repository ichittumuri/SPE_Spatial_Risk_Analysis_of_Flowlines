SigmaGridSearch_spatialZ <- function(sigma_grid) {
  
  val_loglik <- function(y, p) {
    eps <- 1e-15
    p <- pmin(pmax(p, eps), 1 - eps)
    sum(y * log(p) + (1 - y) * log(1 - p))
  }
  
  loglik_results <- numeric(length(sigma_grid))
  
  for (i in seq_along(sigma_grid)) {
    sigma_dist <- sigma_grid[i]
    cat("\n--- Sigma =", sigma_dist, "---\n")
    
    d_tr <- train_data$match_distance_m
    w_train <- ifelse(
      train_data$risk == 1 & !is.na(d_tr),
      exp(-d_tr / sigma_dist),
      1
    )
    
    glm_sel <- glm(
      risk ~ . -1,
      data    = Z_df_train,
      family  = binomial(link = "logit"),
      weights = w_train
    )
    betaHat <- coef(glm_sel)
    
    # TURN OFF PLOTTING HERE
    profile <- LambdaGridSearch(
      coord_model = train_coord_model,
      risk_model  = train_risk_model,
      Z           = Z_train,
      betaHat     = betaHat,
      weights     = w_train,
      doPlot      = FALSE      # <--- important
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
    
    f_hat_test   <- as.numeric(predict(MLEFit, test_coord_model))
    eta_glm_test <- as.vector(Z_test %*% betaHat)
    nu_all       <- eta_glm_test + f_hat_test
    p_all        <- plogis(nu_all)
    
    loglik_results[i] <- val_loglik(test_data$risk, p_all)
    cat("  → Test log-likelihood =", round(loglik_results[i], 2), "\n")
  }
  
  results_df <- data.frame(
    sigma  = sigma_grid,
    loglik = loglik_results
  )
  
  best_idx   <- which.max(loglik_results)
  best_sigma <- sigma_grid[best_idx]
  
  cat("\n====================================\n")
  cat("BEST σ =", best_sigma,
      "with log-likelihood =", round(loglik_results[best_idx], 2), "\n")
  cat("====================================\n")
  
  invisible(list(
    results_df = results_df,
    best_sigma = best_sigma
  ))
}
