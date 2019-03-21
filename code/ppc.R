## Posterior predictive checks

# Precompute residuals
compute_residuals <- function(X, y, yhat, linear_combs, iterations, burn) {
  
  # Effective iterations
  its <- (iterations - burn)
  
  # Precompute residuals
  resids <- array(0L, dim = c(nrow(X), its, 2))
  
  # For each desired sample, calculate
  #browser()
  for(i in 1:its) {
    
    # Calculate simulated resid
    resids[, i, 1] <- yhat[[i]] - linear_combs[,i]
    
    # Calculate observed residual
    resids[, i, 2] <- y - linear_combs[,i]
    
  }
  
  # Return
  return(resids)
  
}

# Heteroskedasticity test
test_for_homoskedasticity <- function(y, X) {
  
  # Regress residuals on X
  fit <- lm("y ~ .", as.data.frame(cbind(y, X[,-1])))
  
  # Return
  return(summary(fit)$r.squared)
  
}

# Model fit
rmse <- function(residuals) {
  
  sqrt(mean(residuals^2))
  
}
