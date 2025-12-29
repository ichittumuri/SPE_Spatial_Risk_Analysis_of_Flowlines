LambdaGridSearch <- function(coord_model,
                             risk_model,
                             Z           = NULL,
                             betaHat     = NULL,
                             M           = 10,
                             exp_range   = c(-2, 2),
                             fine_range  = c(-2, 3),
                             fine_length = 250,
                             weights     = NULL,
                             doPlot      = TRUE) {   # <- NEW ARG
  
  lambdaGrid <- 10^seq(exp_range[1], exp_range[2], length.out = M)
  logLike    <- numeric(M)
  look2      <- NULL
  
  for (k in M:1) {
    
    if (k < M) {
      if (!is.null(Z) && !is.null(betaHat))
        nuOld <- look2$fitted.values + as.vector(Z %*% betaHat)
      else
        nuOld <- look2$fitted.values
    } else nuOld <- NULL
    
    look2 <- logisticSmoother(
      s       = coord_model,
      y       = risk_model,
      lambda  = lambdaGrid[k],
      Z       = Z,
      betaHat = betaHat,
      nuOld   = nuOld,
      weights = weights
    )
    
    logLike[k] <- look2$summary["lnProfileLike.FULL"]
  }
  
  # ---- ONLY PLOT IF doPlot = TRUE ----
  if (doPlot) {
    plot(
      log10(lambdaGrid), logLike, type = "b",
      xlab = "log10(λ)", ylab = "Profile log likelihood",
      ylim = c(min(logLike) - 2, max(logLike) + 1)
    )
    
    lGrid <- seq(fine_range[1], fine_range[2], length.out = fine_length)
    profile_spline <- splint(log10(lambdaGrid), logLike, lGrid)
    lines(lGrid, profile_spline, col = "blue")
    
    logLambdaHat <- lGrid[which.max(profile_spline)]
    abline(v = logLambdaHat, col = "blue", lty = 2)
    
  } else {
    # still compute the fine grid + best lambda, but NO PLOT
    lGrid <- seq(fine_range[1], fine_range[2], length.out = fine_length)
    profile_spline <- splint(log10(lambdaGrid), logLike, lGrid)
    logLambdaHat <- lGrid[which.max(profile_spline)]
  }
  
  bestLambda <- 10^logLambdaHat
  
  invisible(list(
    lambdaGrid      = lambdaGrid,
    logLike         = logLike,
    lGrid           = lGrid,
    profile_spline  = profile_spline,
    logLambdaHat    = logLambdaHat,
    bestLambda      = bestLambda
  ))
}
