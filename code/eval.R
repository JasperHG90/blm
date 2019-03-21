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


