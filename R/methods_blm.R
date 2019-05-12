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
    "\tOutcome: ", x$input$variables$DV, "\n\n"
  )

  # Cat to console
  cat(msg)
  print(get_value(x, "sampling_settings"))
  print(get_value(x, "priors"))

  # If hypotheses, cat
  if(contains(x, "hypotheses")) {
    print(get_value(x, "hypotheses"))
  }

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
  CIV <- t(round(CCI(get_value(x, "posterior")), digits = 4))

  # Print MAP & SE
  cat(crayon::bold("Model results for blm object:"))
  cat("\n\n")
  cat(general)
  cat("\n\n")
  cat(has_sampled)
  cat("\n\n")
  print.listof(obs)
  cat("---------------------\n\n")
  print.listof(list("Maximum a posteriori (MAP) estimates" = MAPV))
  print.listof(list("95% credible interval" = CIV))
  cat("---------------------\n\n")

  # If R-squared, cat value
  if("rsq" %in% names(x)) {
    summary(get_value(x, "rsq"))
  }

  # If bayes factor for the model, cat value
  if(contains(x, "model_BF")) {
    # Cat model DIC
    summary(get_value(x, "DIC"),
            x %>%
              get_value("null_model") %>%
              get_value("DIC"))
    summary(get_value(x, "model_BF"))
    cat("(BF > 1 is evidence for the user-defined model against the intercept-only model)\n\n")
  } else {
    # Cat model DIC
    summary(get_value(x, "DIC"))
  }

  # if hypotheses
  if(contains(x, "hypotheses")) {
    # Cat
    summary(get_value(x, "hypotheses"))
  }

}

#' Plot a blm object
#'
#' @param x blm object
#' @param type one of 'history' (trace plot), 'autocorrelation', 'density' or 'nullmodel'. The latter option plots all plots on one grid for the intercept-only model if it is computed. See also \code{compute_null_model()}
#' @param ... other options passed for specific plots. Accepted arguments are:
#' \describe{
#'     \item{chain}{Integer >= 1 specifying which chain to use for autocorrelation plot. If not chosen, chain 1 is automatically used.}
#' }
#'
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
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 labs
#' @importFrom tidyr gather
#' @importFrom gridExtra grid.arrange
#' @importFrom scales pretty_breaks
#' @importFrom stringr str_replace_all
plot.blm <- function(x, type=c("history",
                               "autocorrelation",
                               "density",
                               "nullmodel"),
                     ...) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # If null model
  if(type == "nullmodel") {

    # Retrieve null model
    nullmodel <- get_value(x, "null_model")

    # Plot trace plot
    p1 <- plot(nullmodel, "history") + theme(legend.position = "none")

    # Plot density plot
    p2 <- plot(nullmodel, "density") + theme(legend.position = "none")

    # Plot autocorrelation
    p3 <- plot(nullmodel, "autocorrelation", chain=1)

    # Plot on grid
    return(grid.arrange(p1, p2, p3, nrow=3))

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
    if(!"chain" %in% names(opts)) {

      # Emit warning if chains > 1
      if(length(get_value(x, "sampling_settings")) > 1) {
        warning("Choosing chain 1 for autocorrelation plot. You can choose a different chain by passing 'chain = <number>' to the plot() function.")
      }

      # Choose default
      chain <- 1

    } else {

      chain <- opts$chain

    }

    # Call autocor plot
    autocorrelation_plot(pd, chain)

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

      history_plot(samples)

    } else if(type == "density") {

      ## Density plot

      density_plot(samples)

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

#' @export
set_sampling_options.blm <- function(x, chains = 1, iterations = 10000,
                                     thinning = 1, burn = 1000) {

  # Unroll data
  opts <- get_value(x, "sampling_settings")

  # Update sampling options
  if("posterior" %in% names(x)) {
    if(!missing(chains) | !missing(iterations) | !missing(thinning)) {
      # refuse to update sampler values if posterior already sampled!
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

# Add a hypothesis to a blm object
#' @export
set_hypothesis.blm <- function(x, name, hypothesis_user) {

  # Hypothesis in blm object?
  if(!contains(x, "hypotheses")) {
    x <- set_value.blm(x, "hypotheses", list())
    # Add structure
    class(x$hypotheses) <- "hypotheses"
  }

  # If name not a string, reject
  if(!is.character(name)) {
    stop("'name' of the hypothesis must be a character")
  }

  # Get parameters
  pars <- names(x$priors)
  pars <- pars[pars != "sigma"]

  # Does the hypothesis already exist?
  current_hyps <- lapply(x$hypotheses, function(x) x$hypothesis) %>%
    unlist()
  # Simply return
  if(hypothesis_user %in% current_hyps) {
    return(x)
  }

  # Add hypothesis
  x[["hypotheses"]][[name]] <- hypothesis(hypothesis_user, pars)

  # Set class
  class(x$hypotheses) <- "hypotheses"

  # Return
  return(x)

}

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

#' @export
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

#' @export
compute_null_model.blm <- function(x, iterations=10000, chains=1, thinning=1, burn = 1000) {

  # Checks
  iterations <- as.integer(iterations)
  chains <- as.integer(chains)
  burn <- as.integer(burn)
  thinning <- as.integer(thinning)

  # Checks
  check_sampling_inputs(iterations, chains, thinning, burn)

  # Set up intercept-only model
  form <- as.character(x$input$formula)
  # Reconstruct
  form_null <- paste0(form[2], " ~ 1")

  # Get data
  dat <- as.data.frame(x$input$y)
  colnames(dat) <- x$input$variables$DV

  # Construct null model
  int_only <- blm(form_null, data=dat) %>%
    set_sampling_options(chains=chains, iterations=iterations,
                         thinning=thinning, burn=burn)

  # add to object and return
  x$null_model <- int_only
  # Return
  return(x)

}

# SAMPLING -----

# Execute a blm plan
#' @export
sample_posterior.blm <- function(x) {

  # unroll data
  X <- x$input$X
  y <- x$input$y

  # sample the posterior
  x <- get_value(x, "sampling_settings") %>%
    # Sample method for class 'sampler'
    postsamp(., X, y, get_value(x, "priors")) %>%
    # Add posterior samples to blm object
    set_value(x, "posterior", .) %>%
    # Calculate R-squared
    evaluate_R2() %>%
    # Evaluate the null model
    # (this only happens if the null model is specified using 'compute_null_model()')
    # (which is good because otherwise we'd get stuck in an infinite loop :-/)
    # Null model calls null model calls null model calls null model ......
    sample_null_model() %>%
    # Calculate model DIC
    evaluate_model_fit() %>%
    # Calculate model bayes factor (if null model is sampled)
    evaluate_model_BF() %>%
    # Evaluate hypotheses (if exist)
    evaluate_hypotheses() %>%
    # Compute outliers
    evaluate_outliers()

  # Return blm results
  return(x)

}

# Update an already sampled blm plan
#' @export
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
    set_value(x, "posterior", .) %>%
    # Calculate model DIC
    evaluate_model_fit() %>%
    # Calculate R-squared
    evaluate_R2() %>%
    # Calculate model bayes factor (if null computed)
    evaluate_model_BF() %>%
    # Evaluate hypotheses (if exist)
    evaluate_hypotheses() %>%
    # Compute outliers
    evaluate_outliers()

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

# CONVERGENCE & EVALUATION -----

# Convergence diagnostics
#' @export
evaluate_convergence_diagnostics.blm <- function(x) {

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

# R-squared
#' @export
evaluate_R2.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Set up R-squared by sampling from the posterior
  inputs <- list()

  # Get inputs
  X <- x$input$X
  y <- x$input$y

  # Get posterior samples
  postsamps <- x %>%
    get_value("posterior") %>%
    bind() %>%
    as.matrix()

  # Call compute R2
  r <- bayes_R2(X, y, postsamps)

  # Add
  inputs$rsquared <- r

  # Add class to input
  class(inputs) <- "R2"

  # Return
  x$rsq <- inputs
  return(x)

}

# posterior predictive checks
#' @export
evaluate_ppc.blm <- function(x, return_all=TRUE, p=1) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Check if p between 0 and 1
  if(!(p >= 0 & p <= 1)) {
    stop("'p' must be between 0 and 1")
  }

  # Get y from inputs
  X <- x$input$X
  y <- x$input$y

  # Get posterior samples
  postsamps <- x %>%
    get_value("posterior") %>%
    bind() %>%
    as.matrix()

  # Subset
  selected <- runif(nrow(postsamps))
  postsamps <- postsamps[selected <= p, ]

  # Results
  inputs <- list()

  # Call the gibbs sampler, simulate y values and compute residuals
  r <- ppc_julia(X, y, postsamps)

  # Add results
  inputs$data$sim_y <- r$sim_y
  inputs$data$residuals <- r$residuals

  # Add class to input
  class(inputs) <- "ppc"

  # Add the results to the data
  inputs <- normality_check(inputs, r$skewness)
  inputs <- homoskedast_check(inputs, r$heteroskedasticity)
  inputs <- independence_check(inputs, r$independence)

  # Check if keep samples
  if(return_all) {

    x$ppc <- inputs
    return(x)

  } else {

    inputs$data$sim_y <- NULL
    inputs$data$residuals <- NULL
    x$ppc <- inputs
    return(x)

  }

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

  # Structure
  res <- list(
    "DIC" = final
  )
  class(res) <- "DIC"

  # Add to object
  x$DIC <- res
  # Return
  return(x)

}

# Effective sample size
#' @export
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

# Evaluate hypotheses
#' @export
evaluate_hypotheses.blm <- function(x) {

  # Silently exit if no hypotheses
  if(!contains(x, "hypotheses")) {

    return(x)

  }

  # Error if no posterior
  if(!contains(x, "posterior")) {
    stop("Posterior not yet sampled.")
  }

  # Otherwise, compute fit / complexity for each hypothesis
  hypres <- x %>%
    get_value("hypotheses")

  # Get posterior samples
  psamp <- x %>%
    get_posterior_samples()
  # Drop sigma
  psamp <- psamp[,-ncol(psamp)]

  # Get priors
  priors <- x %>%
    get_value("priors")

  # Set colnames of posterior to parameter names
  vn <- lapply(priors, function(x) x$varname) %>%
    unlist() # Always in the same order as posterior
  colnames(psamp) <- names(vn)

  # SD of y
  y_sd <- sd(x$input$y)

  # For each hypothesis, evaluate
  for(hypothesis_i in seq_along(hypres)) {

    # Retrieve hypothesis
    chyp <- hypres[[hypothesis_i]]

    # Compute complexity
    c_i <- compute_hypothesis_complexity(chyp, priors, y_sd)

    # Compute fit
    f_i <- compute_hypothesis_fit(chyp, psamp, y_sd)

    # Make results
    hres <- list(
      "hypothesis" = list(
        "complexity" = c_i,
        "fit" = f_i
      ),
      "complement" = list(
        "complexity" = 1-c_i,
        "fit" = 1-f_i
      ),
      # The Bayes Factor against complement (c) or the unconstrained hypothesis (u) ==> hoijtink p37
      "BF_c" = ifelse((1-c_i) == 0, 0, (f_i / c_i) / ((1-f_i) / (1-c_i))),
      "BF_u" = ifelse(c_i == 0, 0, f_i / c_i)
    )

    # Add to hypothesis
    hypres[[hypothesis_i]]$result <- hres

  }

  # Set value for hypotheses
  x <- x %>%
    set_value("hypotheses", hypres)

  # Return
  return(x)

}

# Evaluate model Bayes' Factor
# Based on Wagemakers (2007)
# Compare the model against the Null, meaning
#  H0: intercept-only fits better
#  HA: model containing predictors fits better
#' @export
evaluate_model_BF.blm <- function(x) {

  # If the model does not contain a null model ...
  if(!contains(x, "null_model")) {

    return(x)

  } else {

    nm <- get_value(x, "null_model")

    # To list
    mod_BF <- list(
      "inputs" = list(
        "model0" = list(
          "n" = nm$input$n,
          "p" = nm$DIC$DIC$`Eff. P`,
          "Neg. LL" = nm$DIC$DIC$LL,
          "BIC" = BIC_blm(nm$input$n, round(nm$DIC$DIC$`Eff. P`, digits=2), nm$DIC$DIC$LL)
          ),
        "model1" = list(
          "n" = x$input$n,
          "p" = x$DIC$DIC$`Eff. P`,
          "Neg. LL" = x$DIC$DIC$LL,
          "BIC" = BIC_blm(x$input$n, round(x$DIC$DIC$`Eff. P`, digits=2), x$DIC$DIC$LL)
        )
      )
    )

    # Add bf
    mod_BF$BF <- BF(mod_BF$inputs$model0$BIC,
                    mod_BF$inputs$model1$BIC)

    # Add structure
    class(mod_BF) <- "BF"

    # Add to model
    x$model_BF <- mod_BF

    # Return
    return(x)

  }

}

# Evaluate outliers
#' @export
evaluate_outliers.blm <- function(x) {

  # Check if posterior in blm object
  if(!"posterior" %in% names(x)) {
    stop("Posterior not yet constructed.")
  }

  # Results list
  inputs <- list()

  # Get inputs
  X <- x$input$X
  y <- x$input$y

  # Get posterior samples
  postsamps <- x %>%
    get_value("posterior") %>%
    bind() %>%
    as.matrix()

  # Compute outliers
  r <- julia_outliers(X, y, postsamps)

  # Add
  inputs$results <- r

  # Add class to input
  class(inputs) <- "outliers"

  # Return
  x$outliers <- inputs
  return(x)

}

# Not exported ----

# Sample the null model if it exists in the data!
# (this method is not exported)
sample_null_model.blm <- function(x) {

  # Check if null model in object
  if(!contains(x, "null_model")) {

    # If not, simply return the object
    return(x)

  } else {

    # Sample the null
    nm <- x$null_model %>%
      sample_posterior()

    # Add to object
    x$null_model <- nm

    # Return
    return(x)
  }

}
