logisticSmoother <- function(
    s,
    y,
    lambda,
    Z       = NULL,
    betaHat = NULL,
    nuOld   = NULL,
    weights = NULL
) {
  # IRLS + spatial smoothing on coords s + optional covariates Z
  # 'weights' are observation weights v_i (e.g. distance-based)
  
  n <- length(y)
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  
  # ---------------------------------------------------------------------------
  # Initial linear predictor
  # ---------------------------------------------------------------------------
  if (is.null(nuOld)) {
    if (!is.null(Z) && !is.null(betaHat)) {
      # linear predictor from GLM
      nuOld <- as.vector(Z %*% betaHat)
    } else {
      pStart  <- mean(y)
      nuStart <- log(pStart / (1 - pStart))
      nuOld   <- rep(nuStart, n)
    }
  }
  
  # ---------------------------------------------------------------------------
  # IRLS iterations
  # ---------------------------------------------------------------------------
  for (k in 1:20) {
    # Current probabilities
    pOld <- exp(nuOld) / (1 + exp(nuOld))
    
    # IRLS weights
    W_core <- c(pOld * (1 - pOld))  # variance part
    W      <- W_core * weights      # combine with obs weights v_i
    
    # Working response
    z <- nuOld + (y - pOld) / W_core
    
    # Spatial smoothing step on residual part
    if (!is.null(Z)) {
      z_spatial <- z - as.vector(Z %*% betaHat)
    } else {
      z_spatial <- z
    }
    
    # 'fields::spatialProcess' with weights
    tempObj <- spatialProcess(
      s,
      z_spatial,
      cov.function = "Tps.cov",
      weights      = W,
      lambda       = lambda
    )
    
    f_hat <- tempObj$fitted.values
    
    # Updated linear predictor
    if (!is.null(Z)) {
      nuNew <- as.vector(Z %*% betaHat) + f_hat
    } else {
      nuNew <- f_hat
    }
    
    # Convergence check
    testConv <- mean(abs(nuNew - nuOld))
    if (testConv < 1e-5) break
    
    nuOld <- nuNew
  }
  
  cat("[IRLS Iterations]:", k, "\n")
  
  # tempObj has fitted.values + summary["lnProfileLike.FULL"]
  return(tempObj)
}

