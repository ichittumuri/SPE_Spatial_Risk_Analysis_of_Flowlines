LambdaGridSearch <- function(coord_model,
                             risk_model,
                             Z = NULL,
                             betaHat = NULL,
                             M = 10,
                             exp_range = c(-2, 2),
                             fine_range = c(-2, 3),
                             fine_length = 250) {
  lambdaGrid <- 10^seq(exp_range[1], exp_range[2], length.out = M)
  logLike <- numeric(M)
  look2 <- NULL
  
  for (k in M:1) {
    nuOld <- if (k < M) look2$fitted.values + if (!is.null(Z)) as.vector(Z %*% betaHat) else 0 else NULL
    
    look2 <- logisticSmoother(
      coord_model,
      risk_model,
      lambda = lambdaGrid[k],
      Z = Z,
      betaHat = betaHat,
      nuOld = nuOld
    )
    
    logLike[k] <- look2$summary["lnProfileLike.FULL"]
    cat("λ =", format(lambdaGrid[k], digits = 3),
        "→ logLik =", round(logLike[k], 3), "\n")
  }
  
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
  cat("Estimated log10(λ) =", round(logLambdaHat, 3), "\n")
  
  bestLambda <- 10^logLambdaHat
  cat("Ideal λ =", signif(bestLambda, 3), "\n")
  
  invisible(list(
    lambdaGrid = lambdaGrid,
    logLike = logLike,
    lGrid = lGrid,
    profile_spline = profile_spline,
    logLambdaHat = logLambdaHat,
    bestLambda = bestLambda
  ))
}
