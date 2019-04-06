## Methods for blm class

# GENERIC METHODS ------

# Print method
#' @export
print.blm <- function(x) {

  # Print information
  msg <- paste0(
    crayon::bold("Bayesian Linear Model (BLM) object:\n\n"),
    "Data:\n",
    "\tPredictors: ", x$input$m - 1, "\n",
    "\tOutcome: ", x$input$variables$DV, "\n",
    "\tCentered: ", x$input$center, "\n\n"
  )

  # Cat to console
  cat(msg)
  print(get_value(x, "sampling_settings"))
  print(get_value(x, "priors"))

}

# Determine the sampling options
#' @export
#' @importFrom magrittr '%>%'
set_sampling_options.blm <- function(x, chains = 1, iterations = 10000,
                                     thinning = 1, burn = 1000) {

  # Unroll data
  opts <- get_value(x, "sampling_settings")

  # Update sampling options
  if("posterior" %in% names(x)) {
    # TODO: refuse to update sampler values if posterior already sampled!
    stop("Posterior already sampled. Cannot update sampling settings. Use update_posterior() to sample more values from the posterior distribution or use delete_posterior() to remove it from your results. After deletion, you can update the sampler settings.")
  }

  # If missing, then take from current values
  if(missing(chains)) {
    # Update # chains
    chains <- as.integer(length(opts))
  } else {
    chains <- as.integer(chains)
  }
  if (missing(iterations)) {
    iterations <- opts$chain_1$iterations
  } else {
    iterations <- as.integer(iterations)
  }
  if (missing(burn)) {
    burn <- opts$chain_1$burn
  } else {
    burn <- as.integer(iterations)
  }
  if (missing(thinning)) {
    thinning <- opts$chain_1$thinning
  } else {
    thinning <- as.integer(thinning)
  }

  # Checks
  check_sampling_inputs(iterations, chains, thinning, burn)

  # Retrieve varnames
  vn <- opts$chain_1$varnames

  # Call
  x <- get_value(x, "priors") %>%
    # Update options
    set_options(opts, chains, iterations, burn, thinning, .) %>%
    # Set results as new sampling objects
    set_value(x, "sampling_settings", .)

  # Return
  return(x)

}

# Set sampler !! Gibbs or MH
set_sampler.blm <- function(x, par, type = c("Gibbs", "MH")) {

  # Match arg
  type <- match.arg(type)

  # Check if par in data
  if(!(par %in% names(get_value(x, "priors")))) {

    stop("Passed parameter names (", par,") not found in data.")

  }

  # Update sampler
  x <- get_value(x, "sampling_settings") %>%
    set_sampler(., par, type = type) %>%
    set_value(x, "sampling_settings", .)

  # Return
  return(x)

}

# Set initial values!!
set_initial_values.blm <- function(x, ...) {

  # Opts
  opts <- list(...)

  # Which chain?
  chain_names <- get_value(x, "sampling_settings") %>%
    names()

  # If names(opts) != names(chain_names)
  if(!all(names(opts) %in% chain_names)) {
    stop("Supplied chain names that are not yet created. Current chain names are '", paste0(chain_names, collapse = ", "), "'.")
  }
  # For every element in the passed user list of initial values, check if the names are proper
  for(val in opts) {
    # If not b or tau in values ...
    if(!("b" %in% names(val))) {
      stop("'val' must be a list with names 'b' for coefficients and 'sigma' for sigma")
    }
    if(!("sigma") %in% names(val)) {
      stop("'val' must be a list with names 'b' for coefficients and 'sigma' for sigma")
    }
    if(length(val[["sigma"]]) != 1) {
      stop("Length of supplied initial values for sigma should equal 1")
    }
    if(val[["sigma"]] <= 0) {
      stop("'sigma' must be larger than 0.")
    }
    # If wrong length
    if(length(val[["b"]]) != nrow(ex$w)) {
      stop("Length of supplied initial values does not correspond to the number of parameters")
    }

  }

  # Retrieve a single example from the blm object
  ex <- get_value(x, "sampling_settings") %>%
    get_value(., "chain_1") %>%
    get_value(., "initial_values")

  # For each chain, set initial values
  for(i in names(opts)) {
    # Set coefficients
    x[["sampling_settings"]][[i]][["initial_values"]][["w"]][1:3,1] <- opts[[i]]$b
    # Set sigma
    x[["sampling_settings"]][[i]][["initial_values"]][["sigma"]] <- opts[[i]]$sigma
    # Set user defined starting postions
    x[["sampling_settings"]][[i]][["inits_user_defined"]] <- TRUE
  }

  # Return
  return(x)

}

# Class method to change priors for available parameters
#' @export
#' @importFrom magrittr '%>%'
set_prior.blm <- function(x, par, ...) {

  # Get opts
  new_prior_values <- list(...)

  # Names in data?
  if(!(par %in% names(get_value(x, "priors")))) {

    stop("Passed parameter names (", par,") not found in data.")

  }

  # Call update priors
  x <- get_value(x, "priors") %>%
    # Update priors with new values
    set_prior(., par, new_prior_values) %>%
    # Assign new values to x
    set_value(x, "priors", .)

  # Update initial values in each sampling chain
  x <- get_value(x, "sampling_settings") %>%
    # Draw new initial values
    set_initial_values(., get_value(x, "priors")) %>%
    # Assign updated sampling settings to x
    set_value(x, "sampling_settings", .)

  # Return
  return(x)

}

# SAMPLING -----

# Execute a blm plan
#' @export
#' @importFrom magrittr '%>%'
sample_posterior.blm <- function(x) {

  # unroll data
  X <- x$input$X
  y <- x$input$y

  # sample the posterior
  x <- get_value(x, "sampling_settings") %>%
        # Sample method for class 'sampler'
        sample(., X, y, get_value(x, "priors")) %>%
        # Add posterior samples to blm object
        set_value(x, "posterior", .)

  # Return blm results
  return(x)

}

# Update an already sampled blm plan
#' @export
#' @importFrom magrittr '%>%'
update_posterior.blm <- function(x, iterations = 1000) {

  # unroll data
  X <- x$input$X
  y <- x$input$y

  # Posterior
  posterior_samples <- get_value(x, "posterior")

  # Sampling settings and set each value for iterations (i.e. each chain) to the value of additional iterations
  sampling_settings <- lapply(get_value(x, "sampling_settings"),
                              function(y) y$iterations <- iterations)

  # All sampling settings are now fixed except iterations. Draw new samples.
  posterior_updated_samples <- get_value(x, "sampling_settings") %>%
    # Sample method for class 'sampler'
    sample(., X, y, get_value(x, "priors"))

  # Combine the posterior samples
  for(i in seq_along(1:length(posterior_samples))) {

    posterior_samples[[i]] <- rbind(posterior_samples[[i]], posterior_updated_samples[[i]])

  }

  # Add to blm object
  x <- set_value(x, "posterior", posterior_samples)

  # Return
  return(x)

}

# Remove the posterior
#' @export
delete_posterior.blm <- function(x) {

  # Set posterior to NULL
  return(set_value(x, "posterior", NULL))

}



# CONVERGENCE -----

# EVALUATION ----
