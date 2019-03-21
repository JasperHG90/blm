# Gibbs sampler

# Initiate initial values
initialize_chain_values <- function(priors) {
  
  # For each prior, draw initial values and construct a weight matrix
  w <- matrix(0L, nrow = length(priors)-1, ncol=1)
  
  # To grid
  for(i in seq_along(priors)) {
    if(names(priors)[i] == "sigma") {
      sigma <- draw_value(priors[[i]])
    } else {
      w[i,1] <- draw_value(priors[[i]])
    }
  }
  
  # Return
  return(
    list(
      "w" = w,
      "sigma" = sigma
    )
  )
  
}

# One iteration of the sampler
gibbs_one_iteration <- function(X, y, w, sigma, priors) {
  
  # Update the parameters stepwise
  for(j in seq_along(priors)) {
    
    # Retrieve prior
    prior <- priors[[j]]
    
    # If sigma, do ...
    if(names(priors)[j] == "sigma") {
      
      # Get posterior for sigma
      posterior_sigma_values <- posterior_sigma(X, y, w, priors[["sigma"]])
      
      # Draw
      sigma <- 1 / rgamma(1, 
                          posterior_sigma_values$alpha_posterior, 
                          posterior_sigma_values$beta_posterior)
      
    } else {
        
      # Get posterior
      posterior_values <- posterior_coef(X, y, w, j, prior, sigma)
      
      # Draw from posterior
      w[j,1] <- rnorm(1, posterior_values$mu_posterior, posterior_values$tau_posterior)
      
    }
    
  }
  
  # Return values
  return(
    list(
      "w" = w,
      "sigma" = sigma
    )
  )
  
}

# Sample posterior distribution in R
#
# @param X Design matrix X
# @param y Outcome vector (dependent variable) y
# @param w coefficient weights
# @param initial_values initial values for the MCMC chain
# @param iterations number of iterations user desires to run the algorithm
# @param priors prior distribution parameters. Must be a list containing (#variables + 2) of class 'prior'
# @param burn number burn-in samples

gibbs_sampler <- function(X, y, initial_values, iterations, priors, burn) {
  
  # Set up result matrix
  results <- matrix(0L, nrow=iterations, ncol = length(priors))
  
  # Unroll initial values
  w <- initial_values$w
  sigma <- initial_values$sigma
  
  # For each pass in the iteration
  for(iteration in seq_along(1:iterations)) {
    
    # One step of the gibbs sampler
    params <- gibbs_one_iteration(X, y, w, sigma, priors)
    
    # Unroll new values for coefs and sigma
    w <- params$w
    sigma <- params$sigma
    
    # Append results to matrix
    results[iteration,1:(ncol(results)-1)] <- w[,1]
    results[iteration,ncol(results)] <- sigma
    
  }
  
  # Remove burn-in samples
  results <- results[-1:-burn,]
  
  # Return results
  return(results)
  
}

# Burn-in period diagnostics
# Run a linear regression on squared coefficient draws from posterior
burnin_diagnostic <- function(posterior) {
  
  out_post <- lapply(seq_along(posterior), function(chain_it) {
    
    # Get current chain
    chain <- posterior[[chain_it]]
    
    # For each column, compute linear coef
    out <- lapply(seq_along(1:ncol(chain)), function(x) {
      
      df <- data.frame(
        "y" = chain[,x]^2,
        "index"= (1:nrow(chain)) 
      )
      
      # Linear reg
      linr <- lm("y ~ index", data=df)
      
      # Get coef
      linr$coefficients["index"]
      
    })
    
    # Name out
    names(out) <- colnames(chain)
    
    # Bind
    out <- do.call(cbind.data.frame, out)
    row.names(out) <- paste0("chain ", chain_it)
    
    # Return
    return(out)
    
  })
  
  # Bind
  diagnostics <- do.call(rbind.data.frame, out_post)
  
  # Return
  return(diagnostics)
  
}
