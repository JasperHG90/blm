# Evaluation functions

# Calculate MAP values & SE over all chains
MAP <- function(posterior) {
  
  # Join the data row-wise
  D <- do.call(rbind.data.frame, posterior)
  
  # Get posterior MAP values
  MAP_R <- apply(D, 2, mean)
  SE_R <- apply(D, 2, sd)
  
  # Return
  return(
    list(
      "MAP" = MAP_R,
      "SE" = SE_R
    )
  )
  
}

# 95% credible interval
CI <- function(posterior) {
  
  # Bind data row-wise
  D <- do.call(rbind.data.frame, posterior)
  
  # Calculate 95% CI
  post_credint <- apply(D, 2, function(x) quantile(x, c(0.025, 0.975)))
  
  # Return
  post_credint
  
}

# Gelman-rubin statistic
#https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
# Compare the within / between variances
GR <- function(posterior, iterations) {
  
  ## For each chain, and each statistic, calculate the between and within
  ## variance
  M <- length(posterior)
  N <- iterations
  m <- ncol(posterior$chain_1)
  
  ## Open up results matrix
  RM <- matrix(0L, ncol = m, 
               nrow = 1,
               dimnames = list(
                 c(""),
                 colnames(posterior$chain_1)
               ))
  
  ## For each parameter
  for(j in seq_along(1:m)) {
  
      # Get chain means
      cm <- matrix(unname(
        sapply(posterior, function(x) mean(x[,j]))
        ), ncol=1)
      
      # Get chain vars
      cv <- matrix(unname(
        sapply(posterior, function(x) var(x[,j]))
      ), ncol=1)
      
      # Subtract parameter mean over all chains
      cmm <- cm - mean(cm[,1])
      
      # Between variance
      B <- (N / (M - 1)) * t(cmm) %*% cmm
      
      # Within variance
      W <- (1/M) * sum(cv)
      
      # Pooled variance
      V <- (((N-1) / N) * W) + (((M+1) / (M*N)) * B)
      
      # GR stat
      RM[1, j] <- V / W
      
  }
  
  # Return 
  return(RM)
    
}

# Calculate mode
calc_mode <- function(x) {
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

# Model DIC
# See: http://kylehardman.com/BlogPosts/View/6
DIC <- function(X, y, posterior) {
  
  ## Join posterior
  pb <- do.call(rbind.data.frame, posterior)
  
  ## MAP values
  MAP <- apply(pb, 2, mean)
  
  ## Coef separate from sigma
  sigma <- unname(MAP["sigma"])
  coefs <- matrix(MAP[-length(MAP)], ncol=1)
  
  ## Two parts to DIC ==> (1) sum of log of likelihood P(y|theta_MAP)
  
  # Mu for each point in the likelihood
  mu <- X %*% coefs
  # Draw from normal
  LL <- sum(log(dnorm(y, mean=mu, sd=sqrt(sigma))))
  
  ## Part two: effective number of parameters p_dic
  
  ## Multiply X by the posterior coefficients from the sample
  # Result: n x k matrix (n==examples, k==# gibbs samples across ALL chains)
  lincom <- X %*% t(pb[,-ncol(pb)])
  
  # For each gibbs sample, calculate the log likelihood
  P <- lapply(seq_along(1:ncol(lincom)), function(x) {
    
    # For each gibbs sample, draw sample 
    sum(log(dnorm(y, mean=lincom[,x], sd=sqrt(pb[x,ncol(pb)]))))
    
  })
  
  ## Average of all samples
  P <- mean(unlist(P))
  
  ## Effective Params
  P_eff <- 2*(LL-P)
  
  ## Return model fit
  return(
  do.call(cbind.data.frame,
          list(
            'DIC' = -2*LL + 2*P_eff,
            'LL' = LL,
            "Eff. P" = P_eff
          ))
  )
  
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
