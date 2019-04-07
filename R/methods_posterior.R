## Methods for the posterior class

#' Burn samples from a posterior sample dataset
#'
burn.posterior <- function(x) {

  # Return posterior without burn-in samples
  x <- get_value(x, "samples") %>%
    lapply(., function(y) {
      y[-1:-get_value(x, "burn"),]
    })

  # Return
  return(x)

}

#' Collapse all chains into one big sample
bind.posterior <- function(x) {

  # Burn
  burn(x) %>%
    # Bind
    do.call(rbind.data.frame, .)

}

#' Append samples to a posterior distribution
append_samples.posterior <- function(x, updates) {

  samples <- get_value(x, "samples")
  new_samples <- get_value(updates, "samples")

  # Combine the posterior samples
  for(i in seq_along(x)) {

    # Make an append method for posterior class
    x[["samples"]][[i]] <- rbind(samples[[i]], new_samples[[i]])

  }

  # Return
  return(x)

}
