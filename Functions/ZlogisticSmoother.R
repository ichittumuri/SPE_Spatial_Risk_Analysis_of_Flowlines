logisticSmoother <- function(s, y, lambda, Z = NULL, betaHat = NULL, nuOld = NULL) {
  # IRLS + smoothing on spatial coords + fixed covariates Z
  if (is.null(nuOld)) {
    if (!is.null(Z) && !is.null(betaHat)) {
      nuOld <- as.vector(Z %*% betaHat)  # linear predictor from GLM
    } else {
      # no Z matrix given
      pStart <- mean(y)
      nuStart <- log(pStart / (1 - pStart))
      nuOld <- rep(nuStart, length(y))
    }
  }
  
  for (k in 1:20) {
    pOld <- exp(nuOld) / (1 + exp(nuOld))
    W <- c(pOld * (1 - pOld))
    z <- nuOld + (1 / W) * (y - pOld)
    
    # spatial smoothing step (on residual part)
    if (!is.null(Z)) {
      z_spatial <- z - as.vector(Z %*% betaHat)
    } else {
      z_spatial <- z
    }
    
    tempObj <- spatialProcess(
      s, z_spatial,
      cov.function = "Tps.cov",
      weights = W,
      lambda = lambda
    )
    
    f_hat <- tempObj$fitted.values
    
    if (!is.null(Z)) {
      nuNew <- as.vector(Z %*% betaHat) + f_hat
    } else {
      nuNew <- f_hat
    }
    
    testConv <- mean(abs(nuNew - nuOld))
    if (testConv < 1e-5) break
    
    nuOld <- nuNew
  }
  
  cat("[IRLS Iterations]:", k, "\n")
  return(tempObj)
}

