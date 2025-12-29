LambdaGridSearch <- function(coord_model,
                             risk_model,
                             M          = 10,
                             exp_range  = c(-2,  2),
                             fine_range = c(-2,  3),
                             fine_length = 250,
                             ...) { # added "..."
  
  # 1) set up coarse grid
  lambdaGrid <- 10^seq(exp_range[1], exp_range[2], length.out = M)
  logLike    <- numeric(M)
  look2      <- NULL
  
  # --- grid-search over lambda on the TRAIN set to find MLE of λ ---
  for (k in M:1) {
    # warm-start
    if (k < M) {
      nuOld <- look2$fitted.values
    } else {
      nuOld <- NULL
    }
    
    # positional args: coord, response, then named lambda & nuOld
    look2 <- logisticSmoother(
      coord_model,
      risk_model,
      lambda = lambdaGrid[k],
      nuOld  = nuOld
    )
    
    logLike[k] <- look2$summary["lnProfileLike.FULL"]
    cat("λ =", format(lambdaGrid[k], digits = 3),
        "→ logLik =", round(logLike[k], 3), "\n")
  }
  
  # 2) plot coarse grid
  plot(
    log10(lambdaGrid), logLike, type = "b",
    xlab = "log10(λ)", ylab = "Profile log likelihood"
  )
  
  # 3) spline-interpolate on a finer grid
  lGrid          <- seq(fine_range[1], fine_range[2], length.out = fine_length)
  profile_spline <- splint(log10(lambdaGrid), logLike, lGrid)
  lines(lGrid, profile_spline, col = "blue")
  
  # 4) find the maximizer
  logLambdaHat <- lGrid[which.max(profile_spline)]
  abline(v = logLambdaHat, col = "blue", lty = 2)
  cat("Estimated log10(λ) =", round(logLambdaHat, 3), "\n")
  
  # 5) compute bestLambda
  bestLambda <- 10^logLambdaHat
  cat("Ideal λ =", signif(bestLambda, 3), "\n")
  
  # 6) return everything useful (invisibly)
  invisible(list(
    lambdaGrid     = lambdaGrid,
    logLike        = logLike,
    lGrid          = lGrid,
    profile_spline = profile_spline,
    logLambdaHat   = logLambdaHat,
    bestLambda     = bestLambda
  ))
}







# this is the code for the function:
#   # --- grid‐search over lambda on the TRAIN set to find MLE of λ ---
#   M <- 10
# lambdaGrid <- 10^seq(-2, 2, length.out = M)
# logLike    <- numeric(M)
# 
# for (k in M:1) {
#   # warm-start
#   if (k < M) {
#     nuOld <- look2$fitted.values
#   } else {
#     nuOld <- NULL
#   }
#   
#   # positional args: coord, response, then named lambda & nuOld
#   look2 <- logisticSmoother(
#     coord_model,
#     risk_model,
#     lambda = lambdaGrid[k],
#     nuOld  = nuOld
#   )
#   
#   logLike[k] <- look2$summary["lnProfileLike.FULL"]
#   cat("λ =", format(lambdaGrid[k], digits = 3),
#       "→ logLik =", round(logLike[k], 3), "\n")
# }
# 
# # plot coarse grid
# plot(
#   log10(lambdaGrid), logLike, type = "b",
#   xlab = "log10(λ)", ylab = "Profile log likelihood"
# )
# 
# # spline‐interpolate on a finer grid
# lGrid          <- seq(-2, 3, length.out = 250)
# profile_spline <- splint(log10(lambdaGrid), logLike, lGrid)
# lines(lGrid, profile_spline, col = "blue")
# 
# # find the maximizer
# logLambdaHat <- lGrid[which.max(profile_spline)]
# abline(v = logLambdaHat, col = "blue", lty = 2)
# cat("Estimated log10(λ) =", round(logLambdaHat, 3), "\n")
# 
# # re-fit at the best λ
# bestLambda <- 10^logLambdaHat
# cat("Ideal λ =", signif(bestLambda, 3), "\n")
# MLEFit     <- logisticSmoother(
#   coord_model,
#   risk_model,
#   lambda = bestLambda
# )