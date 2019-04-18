## Methods for blm class

# Generic functions (print, summary etc.)

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

# Summary method
#' @export
#' @importFrom crayon bold
#' @importFrom crayon green
#' @importFrom crayon red
summary.blm <- function(x) {

  var_names <- c(colnames(x$input$X), "sigma")

  ### Formula
  form <- as.character(x$input$formula)
  formchar <- paste0(form[2], " ", form[1], " ", form[3])

  ## Construct individual parts
  general <- paste0(
    "Formula: '", crayon::italic(formchar), "'"
  )

  ## User has already sampled or not
  has_sampled <- paste0(
    "Sampled: ", ifelse("posterior" %in% names(x),
                        crayon::green('TRUE'),
                        crayon::red('FALSE'))
  )

  ## num. observations + predictors
  obs <- list(
    "Sampling settings" = matrix(c(x$input$n,x$input$m - 1,length(x$sampling_settings),
                                   x$sampling_settings$chain_1$iterations,
                                   x$sampling_settings$chain_1$thinning,
                                   x$sampling_settings$chain_1$burn),
                                 nrow = 1,
                                 dimnames = list(
                                   c(""),
                                   c("Obs.",
                                     "Predictors",
                                     "Chains",
                                     "Iterations",
                                     "Thinning",
                                     "Burn")
                                 )))

  ### If not sampled yet ...
  if(!"posterior" %in% names(x)) {
    cat(crayon::bold("Model results for blm object:"))
    cat("\n\n")
    cat(general)
    cat("\n\n")
    cat(has_sampled)
    cat("\n\n")
    print.listof(obs)
    return(cat(""))
  }

  ### Statistics

  # Calculate MAP for each chain & combine
  MAPV <- t(do.call(rbind, MAP(get_value(x, "posterior"))))
  # Add MC error
  MAPV <- cbind(MAPV, MAPV[,2] / sqrt(x$sampling_settings$chain_1$iterations))
  # Add TS error using effective sample size
  MAPV <- cbind(MAPV, MAPV[,2] / sqrt(get_value(x, "posterior") %>% effss()))
  # Round
  MAPV <- round(MAPV, digits=3)

  # Amend names
  colnames(MAPV) <- c("Est. (mean)", "SD", "NAIVE MCERR.", "TS MCERR.")

  # Calculate CI
  CIV <- t(round(CCI(get_value(x, "posterior")), digits = 3))

  # Print MAP & SE
  cat(crayon::bold("Model results for blm object:"))
  cat("\n\n")
  cat(general)
  cat("\n\n")
  cat(has_sampled)
  cat("\n\n")
  print.listof(obs)
  print.listof(list("Maximum a posteriori (MAP) estimates" = MAPV))
  print.listof(list("95% credible interval" = CIV))

}

# Plot method
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_density
#' @importFrom tidyr gather
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_replace_all
plot.blm <- function(x, type=c("history",
                               "autocorrelation",
                               "density"),
                     ...) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Retrieve posterior data
  pd <- get_value(x, "posterior")
  # Burn
  pd <- set_value(pd, "samples", burn(pd))

  ## If autocorrelation plot
  if(type == "autocorrelation") {

    # Get opts
    opts <- list(...)

    # If not chain specified and multiple chains, choose chain 1 but raise error
    if(!"chain" %in% names(opts) & length(get_value(x, "sampling_settings")) > 1) {

      warning("Choosing chain 1 for autocorrelation plot. You can choose a different chain by passing 'chain = <number>' to the plot() function.")

      chain <- 1

    } else {

      chain <- opts$chain

    }

    # Get data from posterior
    qd <- get_value(pd, "samples")[[paste0("chain_", chain)]]

    # Lag data
    qd <- qd %>%
      as.data.frame() %>%
      lapply(., function(x) autocor(x, n=40)) %>%
      do.call(rbind.data.frame, .)

    # Add variable name
    qd$id <- stringr::str_replace_all(row.names(qd), "\\.[0-9]{1,2}", "")

    # Plot
    ggplot2::ggplot(qd, ggplot2::aes(x=lag, y=correlation, group=id)) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::scale_x_continuous(breaks= scales::pretty_breaks()) +
      ggplot2::scale_y_continuous(limits = c(-1,1)) +
      ggplot2::facet_wrap(id ~ .)

  } else {

    samples <- get_value(pd, "samples")
    # Bind data
    for(i in seq_along(samples)) {
      df <- data.frame(samples[[i]])
      df$chain <- i
      df$iteration <- (x$sampling_settings$chain_1$burn + 1):x$sampling_settings$chain_1$iterations
      samples[[i]] <- df
    }

    # To long format
    samples <- samples %>%
      do.call(rbind.data.frame, .) %>%
      tidyr::gather(key = parameter, value = value,
                    -chain, -iteration)

    #### Plot

    # Make into factor
    pd$chain <- as.factor(pd$chain)

    ## History plot

    if(type == "history") {

      ggplot2::ggplot(samples, ggplot2::aes(x = iteration,
                                            y=value,
                                            color=as.factor(chain),
                                            group = parameter)) +
        ggplot2::geom_line(alpha=0.4) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "None") +
        #geom_smooth(method="lm") +
        ggplot2::facet_wrap("parameter ~ .", scales = "free_y",
                            ncol=2)

      ## Density plot

    } else if(type == "density") {

      ggplot2::ggplot(samples, ggplot2::aes(x=value,
                                            fill = chain)) +
        ggplot2::geom_density(alpha=0.4) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "None") +
        ggplot2::facet_wrap("parameter ~ .", scales = "free")

    } else { ## Unknown! Raise error

      ## Raise error
      stop(paste0("Type '", type, "' not a allowed."))

    }

  }

}

# Get coefficients from blm object
#' @export
coef.blm <- function(x, type = c("mean", "mode", "median")) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Match arg if user does not specify type
  type <- match.arg(type)

  # Bind posterior data
  pb <- get_posterior_samples(x)
  # Remove sigma
  pb <- pb[,-ncol(pb)]

  # Compute
  r <- switch(type,
              "mode" = apply(pb, 2, calc_mode),
              "mean" = apply(pb, 2, mean),
              "median" = apply(pb, 2, median))

  # Return
  return(r)

}

# Call coef method using coefficients
#' @export
coefficients.blm <- function(x, type = c("mean", "mode", "median")) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Match arg
  type <- match.arg(type)

  # Call coef
  return(coef(x, type))

}

# Predict method
#' @export
predict.blm <- function(x, type = c("mean", "mode", "median")) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Match arg
  type <- match.arg(type)

  # Get coefficients
  w <- matrix(coef(x, type = type), ncol=1)

  # Predict
  return(
    x$input$X %*% w
  )

}

# Residuals
#' @export
residuals.blm <- function(x, type = c("mean", "mode", "median")) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Match arg
  type <- match.arg(type)

  # Predict
  pred <- predict(x, type=type)

  # Subtract
  return(x$input$y - pred)

}

# Residuals
#' @export
resid.blm <- function(x, type = c("mean", "mode", "median")) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Match arg
  type <- match.arg(type)

  return(residuals(x, type = type))

}

# Functions that build up the blm object with settings etc. ------

# Get variable names and mappings to parameters
#' @export
get_parameter_names.blm <- function(x) {

  cat(crayon::bold("Parameter / variable names for blm object:"))
  cat("\n\n")

  # Retrieve parameter names
  p <- unlist(lapply(x$priors, function(x) x$varname))[-length(x$priors)-1]

  # To df
  par_names <- data.frame("variable" = names(p))
  row.names(par_names) <- unname(p)
  colnames(par_names) <- ""

  # Print
  print.listof(list("Mapping" = par_names))

  }

#' @param chains number of chains to initialize
#' @param iterations number of iterations to run the MC sampler
#' @param thinning set thinning parameter
#' @param burn number of samples that will be discarded
#'
#' @examples
#' bfit <- set_sampling_options(bfit, chains=2, iterations=20000, thinning=3, burn = 2000)
#'
#' @importFrom magrittr '%>%'
#' @rdname set_sampling_options
#' @export
set_sampling_options.blm <- function(x, chains = 1, iterations = 10000,
                                     thinning = 1, burn = 1000) {

  # Unroll data
  opts <- get_value(x, "sampling_settings")

  # Update sampling options
  if("posterior" %in% names(x)) {
    if(!missing(chains) | !missing(iterations) | !missing(thinning)) {
      # TODO: refuse to update sampler values if posterior already sampled!
      stop("Posterior already sampled. Cannot update sampling settings. Use update_posterior() to sample more values from the posterior distribution or use delete_posterior() to remove it from your results. After deletion, you can update the sampler settings.")
    } else if(!missing(burn)) {
      # Update the burn settings in the posterior
      x$posterior$burn <- burn #set_value(x$posterior, "burn", burn)
    }
  }

  # If missing, then take from current values
  if(missing(chains)) {
    # Update # chains
    chains <- as.integer(length(opts))
  } else {
    chains <- as.integer(chains)
  }
  if (missing(iterations)) {
    iterations <- as.integer(opts$chain_1$iterations)
  } else {
    iterations <- as.integer(iterations)
  }
  if (missing(burn)) {
    burn <- as.integer(opts$chain_1$burn)
  } else {
    burn <- as.integer(burn)
  }
  if (missing(thinning)) {
    thinning <- as.integer(opts$chain_1$thinning)
  } else {
    thinning <- as.integer(thinning)
  }

  # Checks
  check_sampling_inputs(iterations, chains, thinning, burn)

  # Retrieve varnames
  vn <- opts$chain_1$varnames

  # Retrieve samplers
  samplers <- opts$chain_1$samplers

  # Call
  x <- get_value(x, "priors") %>%
    # Update options
    set_options(opts, chains, iterations, burn, thinning, vn, ., samplers) %>%
    # Set results as new sampling objects
    set_value(x, "sampling_settings", .)

  # Return
  return(x)

}

#' @param par name of the coefficient, given as b0 (intercept) or b1, ... for coefficients.
#' @param type name of the sampler to use. Either 'Gibbs' or 'MH'. Defaults to 'Gibbs'.
#' @param lambda tuning parameter for the variance of the proposal distribution in metropolis-hastings
#'
#' @examples
#' bfit <- set_sampler(bfit, "b0", type="MH")
#'
#' @importFrom magrittr '%>%'
#' @rdname set_sampler
#' @export
set_sampler.blm <- function(x, par, type = c("Gibbs", "MH"), lambda=0.25) {

  # Match arg
  type <- match.arg(type)

  # Check if par in data
  if(!(par %in% names(get_value(x, "priors")))) {

    stop("Passed parameter names (", par,") not found in data.")

  } else if(par == "sigma" & type == "MH") {

    stop("MH algorithm not implemented for sigma.")

  }

  # Update sampler
  x <- get_value(x, "sampling_settings") %>%
    set_sampler(., par, type = type, zeta=lambda) %>%
    set_value(x, "sampling_settings", .)

  # Return
  return(x)

}

#' @param ... named argument containing (1) a chain name (e.g. 'chain_1') and (2) a list of initial values. See examples.
#'
#' @importFrom magrittr '%>%'
#' @rdname set_initial_values
#' @examples
#' bfit <- set_initial_values(bfit, chain_1 = list("b" = c(2,5,3,4,5), sigma=0.1))
#'
#' @export
set_initial_values.blm <- function(x, ...) {

  # Opts
  opts <- list(...)

  # m (number of coefs)
  m <- x$input$m

  # Which chain?
  chain_names <- get_value(x, "sampling_settings") %>%
    names()

  # Retrieve a single example from the blm object
  ex <- get_value(x, "sampling_settings") %>%
    get_value(., "chain_1") %>%
    get_value(., "initial_values")

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

  # For each chain, set initial values
  for(i in names(opts)) {
    # Set coefficients
    x[["sampling_settings"]][[i]][["initial_values"]][["w"]][1:m,1] <- opts[[i]]$b
    # Set sigma
    x[["sampling_settings"]][[i]][["initial_values"]][["sigma"]] <- opts[[i]]$sigma
    # Set user defined starting postions
    x[["sampling_settings"]][[i]][["inits_user_defined"]] <- TRUE
  }

  # Return
  return(x)

}

#' Class method to change priors for available parameters
#'
#' @param par name of the coefficient, given as b0 (intercept) or b1, ... for coefficients.
#' @param ... other arguments passed to function. In the case of intercept and coefficients, these are 'mu' and 'sd' (prior means and variances). In the case of the residual variance, the rate and shape/scale parameter must be passed as 'alpha' and 'beta'
#'
#' @examples
#' bfit <- set_prior(bfit,"b0", mu=40, sd=7)
#' bfit <- set_prior(bfit, "sigma", alpha=2, beta=1)
#'
#' @export
#' @rdname set_prior
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
#' @rdname sample_posterior
#' @importFrom magrittr '%>%'
sample_posterior.blm <- function(x) {

  # unroll data
  X <- x$input$X
  y <- x$input$y

  # sample the posterior
  x <- get_value(x, "sampling_settings") %>%
        # Sample method for class 'sampler'
        postsamp(., X, y, get_value(x, "priors")) %>%
        # Add posterior samples to blm object
        set_value(x, "posterior", .)

  # Return blm results
  return(x)

}

# Update an already sampled blm plan
#' @param iterations number of additional samples to draw from the posterior.
#'
#' @examples
#' bfit <- update_posterior(blm, iterations=5000)
#'
#' @return blm object with updated count of posterior samples.
#' @export
#' @rdname update_posterior
#' @importFrom magrittr '%>%'
update_posterior.blm <- function(x, iterations = 1000) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # unroll data
  X <- x$input$X
  y <- x$input$y

  # Posterior
  posterior_samples <- get_value(x, "posterior")
  # Dims
  n <- nrow(posterior_samples[["samples"]][[1]])
  m <- ncol(posterior_samples[["samples"]][[1]])

  # Sampling settings and set each value for iterations (i.e. each chain) to the value of additional iterations
  # Set starting values to the last obtained values in the chain
  sampling_settings <- get_value(x, "sampling_settings")
  for(i in seq_along(sampling_settings)) {
    sampling_settings[[i]]["iterations"] <- as.integer(iterations)
    sampling_settings[[i]][["initial_values"]][["w"]][,1] <- unname(posterior_samples[["samples"]][[i]][n,1:(m-1)])
    sampling_settings[[i]][["initial_values"]][["sigma"]] <- unname(posterior_samples[["samples"]][[i]][n,m])
  }

  # All sampling settings are now fixed except iterations. Draw new samples.
  x <- sampling_settings %>%
    # Sample method for class 'sampler'
    postsamp(., X, y, get_value(x, "priors")) %>%
    # Append samples to original
    append_samples(posterior_samples, .) %>%
    # Add to blm object
    set_value(x, "posterior", .)

  # Update number of iterations on the sample
  x <- get_value(x, "sampling_settings") %>%
    lapply(., function(z) set_value(z, "iterations",
                                    z$iterations + iterations)) %>%
    set_value(x, "sampling_settings", .)

  # Class is removed by lapply()
  class(x$sampling_settings) <- "sampler"

  # Return
  return(x)

}

# Remove the posterior
#' @export
delete_posterior.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Set posterior to NULL
  return(set_value(x, "posterior", NULL))

}

# Retrieve posterior samples
#' @export
get_posterior_samples.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Bind method on posterior samples
  get_value(x, "posterior") %>%
    bind()

}

# CONVERGENCE -----

# Convergence diagnostics
#' @export
convergence_diagnostics.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  var_names <- c(colnames(x$input$X), "sigma")

  ### Formula
  form <- as.character(x$input$formula)
  formchar <- paste0(form[2], " ", form[1], " ", form[3])

  ## Construct individual parts
  general <- paste0(
    "Formula: '", crayon::italic(formchar), "'"
  )

  ## User has already sampled or not
  has_sampled <- paste0(
    "Sampled: ", ifelse("posterior" %in% names(x),
                        crayon::green('TRUE'),
                        crayon::red('FALSE'))
  )

  ### If not sampled yet ...
  if(!"posterior" %in% names(x)) {
    cat(crayon::bold("Convergence diagnostics for blm object:"))
    cat("\n\n")
    cat(general)
    cat("\n\n")
    cat(has_sampled)
    return(cat(""))
  } else {

    # Burn-in diagnostics
    pd <- get_value(x, "posterior")
    # Burn
    pd <- set_value(pd, "samples", burn(pd))
    # Burn-in
    burning_diag <- pd %>% burnin_diagnostic()

    # If multiple chains
    if(length(x$sampling_settings) > 1) {
      # Calculate Gelman-Rubin
      GRS <- pd %>%
        GR(., x$sampling_settings$chain_1$iterations)
    }

    # Cat
    cat(crayon::bold("Convergence diagnostics for blm object:"))
    cat("\n\n")
    cat(general)
    cat("\n\n")
    cat(has_sampled)
    cat("\n\n")
    cat("Burn-in diagnostic:\n")
    cat.burnin(burning_diag)
    cat("\n")
    cat("Gelman-Rubin statistic:\n")
    cat.GR(GRS)
  }

}

# EVALUATION ----

# posterior predictive checks
# This returns a SEPARATE object ==> all the simulations are heavy on the memory.
#' @param iterations number of iterations to run for posterior predictive checks. This number will be tagged on to the burn parameter specified under \link[blm]{set_sampling_options}.
#' @export
#' @rdname evaluate_ppc
evaluate_ppc.blm <- function(x, iterations = 2000) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Get burn value
  burn <- x$sampling_settings$chain_1$burn
  # Update iterations
  iterations <- burn + iterations

  # Set up posterior predictive check by sampling from the posterior
  inputs <- list(
    settings = list(
      iterations = iterations,
      burn = burn
    )
  )

  # Get priors etc.
  priors <- x$priors
  thinning <- x$sampling_settings$chain_1$thinning
  samplers <- x$sampling_settings$chain_1$samplers
  iv <- x$sampling_settings$chain_1$initial_values
  chains <- 1
  X <- x$input$X
  y <- x$input$y

  # Check values
  check_sampling_inputs(as.integer(iterations), as.integer(chains), as.integer(thinning),
                        as.integer(burn))

  # Call the gibbs sampler, simulate y values and compute residuals
  r <- ppc_julia(X, y, iv, iterations, priors, thinning, burn, samplers)

  # Add results
  inputs$data$initial_values <- iv
  inputs$data$X <- X
  inputs$data$sim_y <- r$sim_y
  inputs$data$residuals <- r$residuals

  # Add class to input
  class(inputs) <- "ppc"

  browser()

  # Add the results to the data
  inputs <- normality_check(inputs, r$skewness)
  inputs <- homoskedast_check(inputs, r$heteroskedasticity)
  inputs <- independence_check(inputs, r$independence)

  # Return
  return(inputs)

}

# Model fit
#' @export
evaluate_model_fit.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Model fit
  mfit <- DIC(x$input$X, x$input$y,
              get_value(x, "posterior"))

  # Bind models
  final <- round(mfit, digits=3)

  # Row names
  row.names(final) <- c("(Model)")

  # Print
  cat(crayon::bold("Model fit for blm object:\n\n"))
  print.listof(
    list("Model DIC"=final)
  )

}

# Effective sample size
#' @export
#' @rdname evaluate_effective_sample_size
evaluate_effective_sample_size.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Effective sample size
  eff <- get_value(x, "posterior") %>%
    effss(order=15) %>%
    round()

  # To df
  df <- data.frame(n = eff)
  row.names(df) <- c(x$sampling_settings$chain_1$varnames, "sigma")

  # Names
  cat(crayon::bold("Effective sample size for blm object:\n\n"))
  print.listof(
    list("Effective sample size"=df)
  )

}

# Accepted draws
#' @export
#' @importFrom magrittr '%>%'
#' @rdname evaluate_accepted_draws
evaluate_accepted_draws.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Get value from posterior
  acc <- get_value(x, "posterior") %>%
    get_value(., "accepted")

  # Bind row-wise & transpose
  acc_b <- t(do.call(rbind.data.frame, acc))

  # Row & column names
  row.names(acc_b) <- c(x$sampling_settings$chain_1$varnames, "sigma")
  colnames(acc_b) <- names(acc)

  # To data frame
  acc_b <- as.data.frame(acc_b)

  # Print
  cat(crayon::bold("Accepted draws for blm object:\n\n"))
  print.listof(
    list("Accepted draws"=acc_b)
  )

}

# R-squared
#' @export
#' @param iterations number of samples to draw from the posterior distribution
#' @rdname evaluate_R2
evaluate_R2.blm <- function(x, iterations = 4000) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Get burn value
  burn <- x$sampling_settings$chain_1$burn
  # Update iterations
  iterations <- burn + iterations

  # Set up R-squared by sampling from the posterior
  inputs <- list(
    settings = list(
      iterations = iterations,
      burn = burn
    )
  )

  # Get priors etc.
  priors <- x$priors
  thinning <- x$sampling_settings$chain_1$thinning
  samplers <- x$sampling_settings$chain_1$samplers
  iv <- x$sampling_settings$chain_1$initial_values
  chains <- 1
  X <- x$input$X
  y <- x$input$y

  # Check values
  check_sampling_inputs(as.integer(iterations), as.integer(chains), as.integer(thinning),
                        as.integer(burn))

  # Call the gibbs sampler, simulate y values and compute R2
  r <- bayes_R2(X, y, iv, iterations, priors, thinning, burn, samplers)

  # Add
  inputs$rsquared <- r$rsquared
  inputs$posterior_draws <- r$posterior_draws

  # Add class to input
  class(inputs) <- "R2"

  # Return
  return(inputs)

}
