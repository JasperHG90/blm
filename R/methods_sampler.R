## Methods for sampler

# Print method for sampler class
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

# Sample the posterior distribution
#' @importFrom maggritr '%>%'
postsamp.sampler <- function(x, X, y, priors) {

  # Set burn
  b <- get_value(x, "chain_1") %>%
    get_value(., "burn")

  # For each chain, sample posterior sequentially
  posterior_samples <- lapply(x, function(z) {

    # Unroll data
    iterations <- get_value(z, "iterations")
    thinning <- get_value(z, "thinning")
    initial_values_current <- get_value(z, "initial_values")
    samplers <- get_value(z, "samplers")
    zeta <- get_value(z, "zeta")

    # Call mc sampler
    r <- mc_sampler(X, y, initial_values_current, iterations, thinning, priors, samplers, zeta)

    # Add names
    colnames(r) <- c(z$varnames, "sigma")

    # Return
    return(r)

  })

  # Add structure
  posterior_samples <- posterior(posterior_samples, b)

  # Return
  return(posterior_samples)

}

# Updating sampling options
set_options.sampler <- function(x, chains = 1, iterations = 10000,
                                burn = 1000, thinning = 1, varnames, priors) {

  # Create chains
  x <- sampler(chains = chains, iterations = iterations, burn = burn,
               thinning = thinning, varnames, priors)

  # Return
  return(x)

}

# Updating the type of sampler (e.g. gibbs or mh)
#' @importFrom magrittr '%>%'
set_sampler.sampler <- function(x, label, type = c("Gibbs", "MH"), zeta=0.25) {

  # Match arg
  type <- match.arg(type)

  # Check zeta
  if(zeta < 0) {
    stop("'zeta' cannot be negative.")
  }

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

    # Add zeta
    x[[i]][["samplers"]]

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

