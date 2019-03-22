# Helper functions

# Ensure that the user passed a valid formula
#
# @param column names of the design matrix X
# @param formula formula as passed to linear_model()
#
# @return list containing the dependent variable (DV) and independent variable(s) (IV)
check_formula <- function(varnames, formula) {

  ## Coerce to character
  formula <- as.character(formula)

  ## Retrieve DV
  DV <- formula[2]

  ## Make sure DV is in variables
  if( !(DV %in% varnames) ) {

    stop("Dependent variable not found in dataset")

  }

  ## Retrieve IVs
  IV <- formula[3]

  ## If IV == ., move on
  if( IV == '.' ) {

    ## Store results in list
    res <- list(
      "DV" = DV,
      "IV" = varnames[which(varnames != DV)]
    )

  } else if( grepl("\\+", IV) ) {

    ## If '+' in IVs, then multiple IVs

    # Split the IVs
    IVs <- strsplit(IV,"(\\s)?\\+(\\s)?")[[1]]

    # Check for each if in dataset
    if( !all(IVs[which(!grepl("\\*", IVs))] %in% varnames) ) {

      stop("Not all independent variables found in dataset")

    }

    # Store results in a list
    res <- list(
      "DV" = DV,
      "IV" = IVs
    )

  } else {

    ## We have the situation that there is only one IV

    if( !(IV[which(!grepl("\\*", IVs))] %in% varnames) ) {

      stop("Independent variable not found in dataset.")

    }

    ## Store results in list
    res <- list(
      "DV" = DV,
      "IV" = IV
    )

  }

  ## Return the DV / IV list
  return(res)

}

# Perform checks for the model
#
# @param formula formula as passed to linear_model()
# @param data data as passed by the user
#
# @return a list containing variable names, the formula, the design matrix X (as a model.matrix()), the outcome variable y, the number of observations n and the number of predictors m
perform_checks <- function(formula, data, center) {

  ## Check if formula is formula or if it can be coerced to one
  if( !(is(formula)[1] == "formula") ) formula <- as.formula(formula)

  ## Check if data is data frame or matrix
  if( !(is.data.frame(data)) ) {

    stop("'data' must be a data frame")

  }

  ## Retrieve variable names from data
  varnames <- colnames(data)

  ## Check if formula correct and all variable names present
  vars <- check_formula(varnames, formula)

  ## Retrieve X & y matrix / vectors
  y <- data[,vars$DV]

  ## Assert that the outcome variable is a numeric vector
  if( !is.numeric(y) ) {

    stop("outcome vector y must be numeric")

  }

  # Center
  if(center) {
    # Index for IV
    iv_index <- setdiff(names(data), vars$DV)
    # Retrieve numeric
    numeric_vars <- iv_index[sapply(data[, iv_index], is.numeric)]
    # Center variables
    data[,numeric_vars] <- apply(data[,numeric_vars], 2, function(x) x - mean(x))
  }

  ## Create design matrix
  ## (-1 ignores the intercept)
  X <- model.matrix(as.formula(paste0(vars$DV, " ~ ",
                                      paste0(vars$IV, collapse = "+"))),
                    data)

  ## Subset y by rows (may be deleted)
  y <- y[as.numeric(row.names(X))]

  ## Number of observations
  n <- nrow(X)
  ## Number of columns
  m <- ncol(X)

  ## Put in list & return
  res <- list(
    "inputs" = list(
      "formula" = formula,
      "variables" = vars,
      "y" = y,
      "X" = X,
      "n" = n,
      "m" = m,
      "center" = center
    )
  )

  return(res)

}

# Check density parameters
# See function prior() in blm.R
check_density_params <- function(density, params, req_params, range) {

  # Parameter names
  nparam <- names(params)

  # Check if required params met
  if (!all(nparam %in% req_params)) {

    # Start building error message
    msg <- paste0("Density '", density, "' requires parameters '",
                  paste(req_params, collapse=", "),
                  "' but user passed '",
                  paste(nparam, collapse=", "), "'")

    # Error
    stop(msg)

  }

  # Check if required params have legal values
  for (param in nparam) {

    # Supplied value
    param_supplied <- params[[param]]

    # Allowed value
    param_allowed <- range[[param]]

    # In range?
    if (param_supplied < param_allowed[1] | param_supplied > param_allowed[2]) {

      msg <- paste0("Illegal value supplied for parameter '", param, "'. Parameter must be in range (",
                    paste(param_allowed, collapse=", "), ")")

      # Raise error
      stop(msg)

    }

  }

}

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

# Helper function that calls the Julia Gibbs sampler
gibbs_sampler <- function(X, y, initial_values, iterations, priors, burn) {

  # Unroll initial values
  w <- initial_values$w
  sigma <- initial_values$sigma

  # TODO: ensure that user passes valid iterations / priors (integers)
  r <- .blm$julia$eval("gibbs_sampler")(X, y, w, sigma, as.integer(iterations),
                                        unname(priors))

  # Burn
  return(r[-1:-burn,])

}

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
  # Populate matrix
  if(!is.null(seed)) {
    set.seed(seed)
  }
  binary_j <- order(runif(j))[1:binary]

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
