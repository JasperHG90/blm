## Methods for prior class

#' @export
print.priors <- function(x) {

  # Set up priors matrix
  pr_coef <- matrix(0L, ncol=length(x)-1, nrow = 2,
                    dimnames=list(
                      c("mu", "tau"),
                      names(x)[-length(x)]
                    ))
  pr_sigma <- matrix(0L, ncol=1, nrow=2,
                     dimnames = list(
                       c("rate", "scale"),
                       c("sigma")
                     ))

  # Add values
  for(i in seq_along(x)) {
    if(names(x)[i] == "sigma") {
      pr_sigma[1,1] <- get_value(x, "sigma") %>%
                          get_value(., "alpha")
      pr_sigma[2,1] <- get_value(x, "sigma") %>%
                          get_value(., "beta")
    } else {
      pr_coef[1,i] <- get_value(x, names(x)[i]) %>%
                        get_value(., "mu")
      pr_coef[2,i] <- get_value(x, names(x)[i]) %>%
                        get_value(., "sd")
    }

  }

  # Cat to console
  print.listof(list("Priors (Coefficients)" = pr_coef))
  print.listof(list("Priors (Residuals)" = pr_sigma))

}

# Set new prior values on a 'priors' object containing multiple 'prior' object
set_prior.priors <- function(x, par, ...) {

  # Accept a named list (as already been checked by the top-level function)
  x <- get_value(x, par) %>%
    set_prior(., list(...)[[1]]) %>%
    set_value(x, par, .)

  # Return
  return(x)

}

# Set a single prior object to new value
set_prior.prior <- function(x, ...) {

  opts <- list(...)[[1]]

  # Call method
  if(is(x) == "normal") {
    if(!all(c("mu", "sd") %in% names(opts))) {
      stop("Either or both 'mu' and 'sd' not passed to function.")
    }
    # Set prior
    return(set_prior_normal(x, mu = opts$mu, sd = opts$sd, varname = get_value(x, "varname")))
  } else if(is(x) == "gamma") {
    if(!all(c("alpha", "beta") %in% names(opts))) {
      stop("Either or both 'alpha' and 'beta' not passed to function.")
    }
    # Set prior
    return(set_prior_gamma(x, alpha = opts$alpha, beta = opts$beta))
  }

}

# Set a single normal prior object to new value
set_prior_normal <- function(x, mu, sd, varname) {

  # User sets an informative prior
  informative <- TRUE
  # Which is a normal
  density <- "normal"

  # Create a new prior object
  x <- prior(density, mu = mu, sd = sd, varname = varname, informative = TRUE)

  # Return
  return(x)

}

# Set a single gamma prior object to a new value
set_prior_gamma <- function(x, alpha, beta) {

  # Return
  return(prior("gamma", alpha = alha, beta=beta, informative = TRUE))

}

