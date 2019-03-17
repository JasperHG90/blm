## Functions

# Misc ----

# Generate a dataset with normal variables and outcome for testing purposes
#
# @param n number of examples
# @param j number of variables
# @param binary number of columns with 0/1 outcomes
# @seed seed for pseudo-random number generator
#
# @return data frame with n rows and j columns. 

generate_dataset <- function(n = 2000, j = 5, binary = 1, seed=NULL,
                             heteroskedastic = FALSE, center = FALSE,
                             ...) {
  
  opts <- list(...)
  if("degree" %in% names(opts)) {
    degree <- opts$degree
    if(!is.numeric(degree)) {
      stop("'degree' must be numeric")
    } else if(degree <= 0) {
      stop("'degree' cannot be less than or equal to 0")
    } else {
      # Continue
      degree <- degree
    }
  } else{
    degree <- 1
  }
  
  # Empty matrix
  X <- matrix(0L, ncol=j, nrow=n)
  # Determine which variables will be binary
  binary_j <- order(runif(j))[1:binary]
  
  # Populate matrix
  if(!is.null(seed)) {
    set.seed(seed)
  }
  for(i in 1:j) {
    if(!i %in% binary_j) {
      # Mean and var
      m <- runif(1, 1, 10)
      v <- runif(1, 1, 3) 
      X[,i] <- rnorm(n, m, v)
    } else {
      X[,i] <- rbinom(n, 1, 0.5)
    }
  }

  # Coefficients
  if(!is.null(seed)) {
    set.seed(seed)
  }
  # Intercept + j coefs
  coef <- c(rnorm(1, 0, 20), rnorm(j, 0, 10))
  
  # Make the example heteroskedastic
  if(heteroskedastic) {
    sigmaSq <- rep(1e-10, n)
    for(i in setdiff(1:j, binary_j)) {
      X[,i] <- X[,i] + seq(1,mean(X[,i]) + degree*sd(X[,i]), length.out = n)
      # Set sigma
      pw <- runif(1, 1, 1.6)
      sigmaSq <- sigmaSq + (abs(X[,i])^pw)
    }
    sigmaSq
  } else {
    # Sigma squared
    if(!is.null(seed)) {
      set.seed(seed)
    }
    sigmaSq <- runif(1, 1, 10)
  }
  # Center if desired
  if(center) {
    X[,setdiff(1:j, binary_j)] <- apply(X[,setdiff(1:j, binary_j)], 
                                        2, function(x) (x - mean(x)))
  }
  
  # Add intercept
  X <- cbind(rep(1, n), X)
  
  # Create y
  y <- rnorm(n, 
             mean = X %*% matrix(coef), 
             sd = sqrt(sigmaSq)) 
  
  # List of real values (for later)
  real <- list(
    coef = coef,
    sigma = sigmaSq
  )
  
  # Return
  return(
    list(
      "X" = X,
      "y" = y,
      "parameters" = real,
      "seed" = ifelse(is.null(seed), NA, seed)
    )
  )
  
}


