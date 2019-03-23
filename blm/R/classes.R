# Core blm functions

#' Create a blm object
#'
#' The blm() function sets up the bayesian linear model
#'
#' @param formula linear regression formula
#' @param data data frame containing the variables in the linear regression formula
#' @param center TRUE if you want to grand-mean center the continuous variables in the data
#'
#' @return blm object containing roadmap for bayesian linear modeling
#' @export
blm <- function(formula, data, center = FALSE) {

  # Check formula against data
  inputs <- perform_checks(formula, data, center)

  # Retrieve X
  X <- inputs$inputs$X
  # Number of params + 1 (for intercept)
  j <- ncol(X)
  # Variable names
  vn <- colnames(X)

  # Retrieve summary statistics

  # Define priors
  priors <- list()
  # Add  + 1 (for sigma)
  for (param in 1:(j + 1)) {

    if (param == max(j+1)) {

      priors[["sigma"]] <- prior("gamma", alpha=0.01, beta=0.01)
      priors[["sigma"]]$informative <- FALSE

    } else {

      # Coef name
      coef_name <- paste0("b", param - 1)

      # Set prior
      if ((param-1) > 0) {

        priors[[coef_name]] <- prior("normal", mu=0, sd=1000, varname = vn[param])

      } else {

        priors[[coef_name]] <- prior("normal", mu=0, sd=1000)

      }

      priors[[coef_name]]$informative <- FALSE

    }

  }

  # Collect data, add class & return
  final <- list(
    "input" = inputs$inputs,
    "sampling_settings" = list(
      "chains" = 1,
      "iterations" = 10000,
      "burn" = 1000,
      "thinning" = 1
    ),
    "priors" = priors
  )

  # Add class
  class(final) <- "blm"

  # Return
  return(final)

}

# Sets up a prior object
#
# @param density either Normal or Gamma density
#
# @return prior object
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

  # For each density, check if required arguments supplied and allowed
  check_density_params(density, opts, dparams[[density]], dparam_values[[density]])

  # Add varname back
  if ("varname" %in% ls()) {
    opts$varname <- varname
  }

  # Add class
  class(opts) <- c(density, "prior")

  # Return
  return(opts)

}

