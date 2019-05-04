# Helper functions

# Calculate mode
calc_mode <- function(x) {
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

# Helper functions for sampling ----

# Initiate initial values for a Gibbs chain
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

# Helper function that calls the Julia MC sampler
mc_sampler <- function(X, y, initial_values, iterations, thinning, priors, samplers) {

  # Unroll initial values
  w <- initial_values$w
  # If intercept-only model, then w is a scalar. Coerce to array
  # JuliaCall screws this up in conversion!!!
  if(is.vector(w)) {
    w <- matrix(w, ncol=1, nrow=1)

  }

  sigma <- initial_values$sigma

  # Call MCMC sampler in Julia
  r <- .blm$julia$eval("MCMC_sampler")(X, as.numeric(y), w, sigma, as.integer(iterations),
                                       as.integer(thinning), unname(priors), unname(samplers))

  # Burn
  return(r)

}

# Helper functions for evaluation -----

# Helper function for autocorrelation
autocor <- function(x, n=10) {

  # Results
  res <- rep(0, n)

  # Lag for each n and calculate correlation
  for(i in 1:n) {
    res[i] <- cor(x, c(rep(NA, i), x[1:(length(x)-i)]),
                  use="complete.obs")
  }

  # Return
  return(
    data.frame(
      "lag" = 1:n,
      "correlation" = res
    )
  )

}

# Posterior predictive checks in julia
ppc_julia <- function(X, y, posterior_samples) {

  # Call Julia function for posterior predictive checks
  return(
    .blm$julia$eval("ppc_draws")(X, as.numeric(y), posterior_samples)
  )

}

# R-squared calculation in julia
bayes_R2 <- function(X, y, posterior_samples) {

  # Call Julia function for posterior predictive checks
  return(
    .blm$julia$eval("bayes_R2")(X, as.numeric(y), posterior_samples)
  )

}

# Model DIC
# See: http://kylehardman.com/BlogPosts/View/6
DIC <- function(X, y, posterior_samples) {

  ## Map values
  map <- MAP(posterior_samples)$MAP

  ## Bind data
  pb <- bind(posterior_samples)
  ## To matrix (expected by julia)
  pb <- as.matrix(pb)

  ## Coef separate from sigma
  sigma <- unname(map["sigma"])
  coefs <- matrix(map[-length(map)], ncol=1)

  ## Two parts to DIC ==> (1) sum of log of likelihood P(y|theta_MAP)

  # Call Julia implementation
  # NOTE: sigma here is the standard deviation. We use posterior values. These have been squared-rooted
  #  after sampling.
  r <- .blm$julia$eval("DIC")(X, as.numeric(y),
                              coefs, sigma, pb)

  ## Return model fit
  return(
    do.call(cbind.data.frame,
            r)
  )

}

# Helper functions for special cat() ----

# Special cat function for burn-in statistic
cat.burnin <- function(x) {

  # Dims
  j <- dim(x)[2]
  n <- dim(x)[1]

  # Fill setting for cat()
  f <- (nchar(crayon::green("0.00")) + 1) * j + (j-1) + nchar("chain 1 ")

  # Round
  x <- round(abs(x), digits=2)

  # Determine colors
  is_zero <- x == 0
  less_than_05 <- abs(x) < 0.05

  # Turn x into a character vector
  x_v <- as.character(round(unname(unlist(x)), digits=2))

  # For those elements of length 1, add decimals
  x_v_l1 <- nchar(x_v) == 1
  x_v_l3 <- nchar(x_v) == 3

  # Set
  x_v[x_v_l1] <- paste0(x_v[x_v_l1], ".00")
  x_v[x_v_l3] <- paste0(x_v[x_v_l3], "0")

  # Paste
  cat(paste0("        ", paste0("b", 0:(j-1), collapse="   "), "\n"))
  cat(ifelse(is_zero, crayon::green(x_v),
             ifelse(less_than_05, crayon::yellow(x_v),
                    crayon::red(x_v))), fill=f, labels=paste0("Chain ", 1:2), sep = " ")
}

# Special cat function for Gelman-Rubin statistic
cat.GR <- function(x) {

  # Round
  x <- round(x, digits=2)

  # To vector
  x <- unname(x[1,])

  # Length
  j <- length(x)

  # Fill setting for cat()
  f <- (nchar(crayon::green("0.00")) + 1) * j + (j-1) + nchar("chain 1 ")

  # Determine colors
  is_one <- x == 1

  # Turn x into a character vector
  x_v <- as.character(x)

  # For those elements of length 1, add decimals
  x_v_l1 <- nchar(x_v) == 1
  x_v_l3 <- nchar(x_v) == 3

  # Set
  x_v[x_v_l1] <- paste0(x_v[x_v_l1], ".00")
  x_v[x_v_l3] <- paste0(x_v[x_v_l3], "0")

  # Paste
  cat(paste0("        ", paste0("b", 0:(j-1), collapse="   "), "\n"))
  cat(ifelse(is_one, crayon::green(x_v), crayon::red(x_v)), fill=f, labels="Value  ", sep = " ")

}

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
                             heteroskedastic = FALSE, correlated_errors = FALSE, ...) {

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
  # Populate matrix
  if(!is.null(seed)) {
    set.seed(seed)
  }
  binary_j <- order(runif(j))[0:binary]

  for(i in 1:j) {
    if(!i %in% binary_j) {
      # Mean and sd
      m <- runif(1, 5, 10) * sign(runif(1, -1, 1))
      sd <- runif(1, 1, 3)
      X[,i] <- rnorm(n, m, sd)
    } else {
      X[,i] <- rbinom(n, 1, 0.5)
    }
  }

  # Coefficients
  if(!is.null(seed)) {
    set.seed(seed)
  }
  # Intercept + j coefs
  coef <- c(rnorm(1, 0, 20), rnorm(j, 0, 5))
  sigmaSq <- runif(1, 1, 10)

  #browser()

  # Make the example heteroskedastic
  if(heteroskedastic) {
    sigmaSq <- rep(1e-10, n)
    for(i in setdiff(1:j, binary_j)) {
      X[,i] <- X[,i] + seq(1,mean(X[,i]) + degree*sd(X[,i]), length.out = n)
      # Set sigma
      pw <- runif(1, 1, 1.6)
      # Scale power for degree
      if(degree <= 1) {
        degreepw <- 1
      } else {
        degreepw <- degree
      }
      pw <- (pw + (1-(1/degreepw)))/pw
      sigmaSq <- sigmaSq + (abs(X[,i])^pw)
    }
    sigmaSq

    # Add intercept
    X <- cbind(rep(1, n), X)

    # Generate y
    y <- rnorm(n,
               mean = X %*% matrix(coef),
               sd = sqrt(sigmaSq))

  } else if(correlated_errors) { # Make the example correlated

    # Create a time-series object
    if(!is.null(seed)) {
      set.seed(seed)
    }

    # Correlated values
    ord <- arima.sim(list(order = c(1,0,0), ar = degree), n = n, sd=sqrt(sigmaSq))

    # Add intercept
    X <- cbind(rep(1, n), X)

    # Generate y
    y <- rnorm(n,
               mean = X %*% matrix(coef) + ord,
               sd = sqrt(sigmaSq))

  } else {
    # Sigma squared
    if(!is.null(seed)) {
      set.seed(seed)
    }

    # Add intercept
    X <- cbind(rep(1, n), X)

    # Generate y
    y <- rnorm(n,
               mean = X %*% matrix(coef),
               sd = sqrt(sigmaSq))
  }

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

## Sample size helper
sample_size <- function(n, pk) {
  n / (1 + 2 * sum(pk))
}
