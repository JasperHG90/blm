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
  print(get_value(bfit, "sampling_settings"))
  print(get_value(bfit, "priors"))

}

# Do gibbs sampling on blm
# This should be a class method
# Initial weight correction ==> divide coefficient weights by sqrt(sample size) to avoid them getting too large
# when using uninformative priors
#' @export
#' @importFrom magrittr '%>%'
set_sampling_options.blm <- function(blm, chains = 1, iterations = 10000,
                                     thinning = 1, burn = 1000) {

  # Checks
  check_sampling_inputs(iterations, chains, thinning, burn)

  # Unroll data
  opts <- get_value(blm, "sampling_options")

  # Update sampling options
  if("posterior" %in% names(blm)) {
    # TODO: refuse to update sampler values if posterior already sampled!
    stop("Posterior already sampled. Cannot update sampling settings. Use update_posterior() to sample more values from the posterior distribution or use delete_posterior() to remove it from your results. After deletion, you can update the sampler settings.")
  }

  # If missing, then take from current values
  if(missing(chains)) {
    # Update # chains
    chains <- length(opts)
  }
  if (missing(iterations)) {
    iterations <- opts$chain_1$iterations
  }
  if (missing(burn)) {
    burn <- opts$chain_1$burn
  }
  if (missing(thinning)) {
    thinning <- opts$chain_1$thinning
  }

  # Retrieve varnames
  vn <- opts$chain_1$varnames

  # Call
  x <- get_value(blm, "priors") %>%
    # Update options
    set_options(opts, chains, iterations, burn, thinning, .) %>%
    # Set results as new sampling objects
    set_value(x, "sampling_options", .)

  # Return
  return(blm)

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
