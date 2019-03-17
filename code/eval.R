# Evaluation functions

# Calculate MAP values & SE over all chains
MAP <- function(posterior) {
  
  # MAP_R 
  MAP_R <- matrix(0L, ncol = ncol(posterior$chain_1), nrow = 1)
  SE_R <- matrix(0L, ncol = ncol(posterior$chain_1), nrow = 1)
  
  # For each chain, calculate stats
  for(i in seq_along(posterior)) {
    
    # MAP values
    MAP_R <- MAP_R + apply(posterior[[i]], 2, median)
    SE_R <- SE_R + apply(posterior[[i]], 2, sd)
  
  }
  
  # Scale by number of chains
  MAP_R <- MAP_R / length(posterior)
  SE_R <- SE_R / length(posterior)
  
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
  
  # Open results matrix
  post_credint <- matrix(0L, ncol = ncol(posterior$chain_1), nrow = 2,
                         dimnames = list(
                           c("2.5%","97.5%"),
                           colnames(posterior$chain_1)
                         ))
  
  # Calculate 95% credible interval for each chain
  for( i in seq_along(posterior) ) {
    post_credint <- post_credint + apply(posterior[[i]] , 2 , function(x) quantile(x, c(0.025, 0.975)))
  } 
  
  # Scale by number of chains
  return(
    post_credint / length(posterior)
  )
  
}
