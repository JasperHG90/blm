## Checks for params / allowed values passed by user go here

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

  ## Store results in list
  res <- list(
    "DV" = DV,
    "IV" = strsplit(IV, "\\+")[[1]] %>% trimws(),
    "IV_in_data" = varnames[which(DV != varnames)]
  )

  ## Return the DV / IV list
  return(res)

}

# Perform checks for the model
#
# @param formula formula as passed to linear_model()
# @param data data as passed by the user
#
# @return a list containing variable names, the formula, the design matrix X (as a model.matrix()), the outcome variable y, the number of observations n and the number of predictors m
check_init_blm <- function(formula, data) {

  ## Check if formula is formula or if it can be coerced to one
  if( !(is(formula)[1] == "formula") ) formula <- as.formula(formula)

  ## Check if data is data frame or matrix
  if( !(is.data.frame(data)) ) {

    stop("'data' must be a data frame")

  }

  ## Strip row names
  row.names(data) <- 1:nrow(data)

  ## Retrieve variable names from data
  varnames <- colnames(data)

  ## Retrieve DV / IV
  vars <- check_formula(varnames, formula)

  ## Retrieve X & y matrix / vectors
  y <- data[,vars$DV]

  ## Assert that the outcome variable is a numeric vector
  if( !is.numeric(y) ) {

    stop("outcome vector y must be numeric")

  }

  ## Create design matrix
  ## (-1 ignores the intercept)
  ## model.matrix automatically throws error if formula misspecified
  X <- model.matrix(as.formula(paste0(vars$DV, " ~ ",
                                      paste0(vars$IV, collapse = "+"))),
                    data)

  ## Subset y by rows (may be deleted)
  ## Apply listwise deletion
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
      "m" = m
    )
  )

  return(res)

}

# Check inputs for sampling posterior
check_sampling_inputs <- function(iterations, chains, thinning, burn) {

  # All inputs should be integers > 0
  if(!is.integer(iterations) | !(iterations > 0)) {
    stop("Argument 'iterations' must be an integer > 0")
  } else if(!is.integer(chains) | !(chains > 0)) {
    stop("Argument 'chains' must be an integer > 0")
  } else if(!is.integer(thinning) | !(thinning > 0)) {
    stop("Argument 'thinning' must be an integer > 0")
  } else if(!is.integer(burn) | !(burn > 0)) {
    stop("Argument 'burn' must be an integer > 0")
  }

  # Iterations cannot be less than burn-in
  if (iterations < burn) {
    stop("blm cannot sample fewer iterations than burn-in period")
  }

  # Raise warnings if iterations to burn ratio is small
  if (ceiling(iterations / burn) < 2) {
    warning("The number of iterations is very low. This will likely yield unstable results.")
  }

  # Raise warning if burn less than 1.000
  if (burn < 1000) {
    warning("You have specified the burn-in samples to be less than 1.000. This will likely yield unstable results.")
  }

  # Return invisible
  return(invisible(""))

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
