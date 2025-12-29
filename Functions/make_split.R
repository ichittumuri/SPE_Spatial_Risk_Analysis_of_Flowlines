make_train_test_split <- function(df_balanced, p_train = 0.70, seed = 42) {
  set.seed(seed)
  
  dat_fac <- df_balanced %>%
    mutate(risk = factor(risk, levels = c(0, 1)))
  
  train_idx <- createDataPartition(dat_fac$risk, p = p_train, list = FALSE)
  
  train_data <- df_balanced[train_idx, ]
  test_data  <- df_balanced[-train_idx, ]
  
  # quick prevalence check
  cat("\nPrevalence check (shared split):\n")
  cat("Overall p(Spill) =", mean(as.integer(as.character(dat_fac$risk)) == 1), "\n")
  cat("Train   p(Spill) =", mean(train_data$risk == 1), "\n")
  cat("Test    p(Spill) =", mean(test_data$risk == 1), "\n")
  
  list(
    train    = train_data,
    test     = test_data,
    dat_all  = dat_fac,
    idx      = train_idx
  )
}