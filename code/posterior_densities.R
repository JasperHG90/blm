# Posterior density functions

# Posterior mean helper function
posterior_mu <- function(X, xj, y, w, sigma, mu0, tau0) {
  
  # Numerator
  num1 <- sum(xj * (y - X %*% w)) / sigma
  num2 <- mu0 / tau0
  
  # Denominator
  den1 <- t(xj) %*% xj / sigma 
  den2 <- 1 / tau0
  
  # Return
  return( ( num1 + num2 ) / ( den1 + den2 ) )
  
}

# Posterior tau helper function
posterior_tau <- function(xj, sigma, tau0) {
  
  return( 1 / ( (t(xj) %*% xj/sigma) + (1/tau0) ) )
  
}

# Wrapper for the helper functions
posterior_coef <- function(X, y, w, j, prior, sigma) {
  
  # Rows and columns
  n <- nrow(X)
  m <- ncol(X)
  
  # Extract vectors
  xj <- matrix(X[,(j)])
  wj <- w[(j),]
  
  # Remove the jth column from X and row from w
  w <- matrix(w[-(j),])
  X <- matrix(X[,-(j)], nrow=n, ncol=m-1)
  
  # Calculate mu
  mu1 <- posterior_mu(X, xj, y, w, sigma, prior$mu, prior$sd)
  
  # Calculate tau
  tau1 <- posterior_tau(xj, sigma, prior$sd)
  
  # Return
  # (in reality, create an object)
  return(list(
    "mu_posterior" = mu1,
    "tau_posterior" = tau1
  ))
  
}

# Calculate posterior sigma
posterior_sigma <- function(X, y, w, priors) {
  
  # Get number of data points
  n <- nrow(X)
  
  # Calculate alpha1
  alpha1 <- (n / 2) + priors$alpha
  
  # Calculate beta1
  beta1 <- (sum((y - X %*% w)^2) / 2) + priors$beta
  
  # Return
  return(
    list(
      "alpha_posterior" = alpha1,
      "beta_posterior" = beta1
    )
  )
  
}
