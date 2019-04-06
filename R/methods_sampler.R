## Methods for sampler

#' Print method for sampler class
#'
print.sampler <- function(x) {

  # Paste message

  # Note: currently all iterations / burn / thinning settings are the SAME for all chains
  #  This makes life easier for the time being.

  msg <- paste0(
    "Sampler:\n",
    "\tChains: ", length(x) , "\n",
    "\tIterations: ", get_value(x, "chain_1")$iterations , "\n",
    "\tThinning: ", get_value(x, "chain_1")$thinning, "\n",
    "\tBurn: ", get_value(x, "chain_1")$burn , "\n\n"
  )

  # Cat
  cat(msg)

}

#' Sample the posterior distribution
sample.sampler <- function(x, X, y, priors) {

  # For each chain, sample posterior sequentially
  posterior_samples <- lapply(x, function(z) {

    # Unroll data
    iterations <- z$iterations
    thinning <- z$thinning
    initial_values_current <- z$initial_values
    samplers <- z$samplers

    # Call mc sampler
    r <- mc_sampler(X, y, initial_values_current, iterations, thinning, priors, samplers)

    # Add names
    colnames(r) <- z$varnames

    # Return
    return(z)

  })

  # Add structure
  posterior_samples <- posterior(posterior_samples)

  # Return
  return(posterior_samples)

}

#' Updating sampling options
set_options.sampler <- function(x, chains = 1, iterations = 10000,
                                burn = 1000, thinning = 1, priors) {

  # Create chains
  x <- sampler(chains = chains, iterations = iterations, burn = burn,
               thinning = thinning, vn, priors)

  # Return
  return(x)

}

#' Updating the type of sampler (e.g. gibbs or mh)
#' @importFrom magrittr '%>%'
set_sampler.sampler <- function(x, label, type = c("Gibbs", "MH")) {

  # Match arg
  type <- match.arg(type)

  # Retrieve labels
  labs <- get_value(x, "chain_1") %>%
    get_value(., "params")

  # Check if label in labs
  if(!(label %in% labs)) {
    stop(paste0("Label '", label, "' not found. (available labels: ", paste0(labs, collapse=", "), ")"))
  }

  # For each chain, change the value of the sampler belonging to label
  for( i in seq_along(x) ) {

    # Set new value
    x[[i]][["samplers"]][which(labs == label)] <- type

  }

  # Return
  return(x)

}

# Re-draw initial values for each chain (this is necessary when priors are set)
set_initial_values.sampler <- function(x, vpriors) {

  # For each chain, do ...
  for( i in seq_along(x) ) {

    if(!x[[i]][["inits_user_defined"]]) {

      x[[i]] <- set_value(x[[i]], "initial_values", initialize_chain_values(vpriors))

    }

  }

  # Return
  return(x)

}

# Let user set specific initial values
#user_specific_inits.sampler <- function(x, chain, values) {

  # Update initial values
  iv <- get_value(x[[chain]], "initial_values")

  # Set initial values


}
