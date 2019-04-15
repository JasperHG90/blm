# Core blm functions

# Export ----

#' Create a blm object
#'
#' The blm() function sets up the bayesian linear model
#'
#' @param formula linear regression formula
#' @param data data frame containing the variables in the linear regression formula
#'
#' @return blm object containing roadmap for bayesian linear modeling
#' @examples
#' data("attitudes") # Load simulated data
#' bfit <- blm("attitude ~ .", data = attitudes)
#' print(bfit)
#'
#' @export
blm <- function(formula, data) {

  # Check formula against data
  inputs <- check_init_blm(formula, data)

  # Retrieve X
  X <- inputs$inputs$X
  # Number of params + 1 (for intercept)
  j <- ncol(X)
  # Variable names
  vn <- colnames(X)

  # Define priors
  priors_init <- priors(j, vn)

  # Set samplers to Gibbs
  # j coefficients + sigma
  samplers <- lapply(1:(j+1), function(x) list("Gibbs"))
  # Set names
  names(samplers) <- c(paste0("b", 0:(j-1)), "sigma")

  # Set up sampler settings
  # Note these are defaults
  samplers_init <- sampler(chains = 1, iterations = 10000, burn = 1000,
                           thinning = 1, vn, priors_init, samplers)

  # Collect data, add class & return
  final <- list(
    "input" = inputs$inputs,
    "sampling_settings" = samplers_init,
    "priors" = priors_init
  )

  # Add structure
  class(final) <- "blm"

  # Return
  return(final)

}

# Internal use only -----

#' Priors class
#'
#' @param j number of parameters
#' @param vn variable names of the j parameters
priors <- function(j, vn) {

  priors <- list()
  # Add  + 1 (for sigma)
  for (param in 1:(j + 1)) {

    if (param == max(j+1)) {

      priors[["sigma"]] <- prior("gamma", alpha=0.01, beta=0.01)

    } else {

      # Coef name
      coef_name <- paste0("b", param - 1)

      # Set prior
      if ((param-1) > 0) {

        priors[[coef_name]] <- prior("normal", mu=0, sd=1000, varname = vn[param])

      } else {

        priors[[coef_name]] <- prior("normal", mu=0, sd=1000, varname = "Intercept")

      }

    }

  }

  # Add structure
  class(priors) <- "priors"

  # Return
  return(priors)

}

#' Sets up a prior class
#'
#' Every blm object contains an S3 class called 'priors'. This object is used to store information about priors and has its own methods (mainly used for internal purposes). The user can specify priors by means of the function 'set_priors()'
#'
#' @param density either Normal or Gamma density
#'
#' @return prior object
prior <- function(density, ...) {

  # Get options
  dparams <- getOption("blm_dparams")
  dparam_values <- getOption("blm_allowed_dparam_values")

  # Check if density in list of allowed
  if (!density %in% names(dparams)) {
    stop(paste0("Density '", density, "' not known"))
  }

  # Coerce optional arguments to list
  opts <- list(...)

  # Remove varname if present
  if ("varname" %in% names(opts)) {

    varname <- opts$varname
    opts$varname <- NULL

  }
  # Is the prior informative?
  if ("informative" %in% names(opts)) {

    informative <- TRUE
    opts$informative <- NULL

  } else {
    informative <- FALSE
  }

  # For each density, check if required arguments supplied and allowed
  check_density_params(density, opts, dparams[[density]], dparam_values[[density]])

  # Add varname back
  if ("varname" %in% ls()) {
    opts$varname <- varname
  }

  # Add informative
  opts$informative <- informative

  # Add class
  class(opts) <- c(density, "prior")

  # Return
  return(opts)

}

#' Sets up a sampler class
#'
#' @return sampler object
sampler <- function(chains, iterations, burn, thinning, vars, priors, samplers) {

  # Initialize settings
  chains_init <- vector("list", chains)

  # Fill chain
  for(i in seq_along(chains_init)) {
    chains_init[[i]] <- chain(priors, iterations, burn, thinning, vars, samplers)
  }

  # Names
  names(chains_init) <- paste0("chain_", 1:chains)

  # Structure
  class(chains_init) <- "sampler"

  # Return
  return(chains_init)

}

#' Set up a single chain
chain <- function(priors, iterations, burn, thinning, vars, samplers, ...) {

  pts <- list(...)

  # Check if inits in function
  if("inits" %in% names(pts)) {
    inits_user_defined <- pts$inits
  } else {
    inits_user_defined <- FALSE
  }

  # Set up the object
  opts <- list(
    "iterations" = as.integer(iterations),
    "burn" = as.integer(burn),
    "thinning" = as.integer(thinning),
    "params" = paste0("b", 0:length(vars)),
    "varnames" = vars,
    "samplers" = samplers,
    "initial_values" = initialize_chain_values(priors),
    "inits_user_defined" = inits_user_defined
  )

  # Return
  class(opts) <- c("chain")
  return(opts)

}

#' Sets up posterior class
#'
#' @param samples posterior samples drawn using MCMC algorithm
#' @param burn burn-in samples to be removed
#' @param accepts vector of length (# of coefficients + 2) indicating the number of accepted draws for each coefficient (this includes sigma and the intercept). For Gibbs, this will always be equal to the number of iterations.
#'
#' @return list containing matrices of chain values (dims: iterations x coef + 1) of object posterior
posterior <- function(samples, burn, accepts) {

  # New object
  post <- list(
    "burn" = burn,
    "samples" = samples,
    "accepted" = accepts
  )

  # Add structure
  class(post) <- "posterior"

  # Return
  return(post)

}
