# S3 Methods

# Class 'blm' methods -----

# Class method to change priors for available parameters
#' @export
set_priors.blm <- function(blm, ...) {

  # Get opts
  priors_new <- list(...)

  # Get current priors
  priors <- blm$priors

  # Empty vector to keep track of warnings
  warning_vars <- c()

  # Set new priors (based on names)
  for(i in names(priors_new)) {

    if(i %in% names(priors)) {

      # If varname in priors, then add to priors_new
      if("varname" %in% names(priors[[i]])) {
        priors_new[[i]]$varname <- priors[[i]]$varname
      }

      # Update priors
      priors[[i]] <- priors_new[[i]]

      # Add marker
      priors[[i]]$informative <- TRUE

    } else {

      # Keep track of warnings
      warning_vars <- c(warning_vars, i)

    }

  }

  # Emit warnings
  #  ...

  # Set new priors and return
  blm$priors <- priors

  # Return
  return(blm)

}

# Do gibbs sampling on blm
# This should be a class method
# TODO: add thinning
# Initial weight correction ==> divide coefficient weights by sqrt(sample size) to avoid them getting too large
# when using uninformative priors
#' @export
sampling_options.blm <- function(blm, chains = 1, iterations = 10000, thinning = 3,
                                 burn = 1000) {

  ##### Checks ######

  if (chains < 1) {
    stop("'chains' cannot be less than 1")
  }

  if (iterations < burn) {
    stop("blm cannot sample fewer iterations than burn-in period")
  }

  if (ceiling(iterations / burn) < 2) {
    warning("The number of iterations is very low. This will likely yield unstable results.")
  }

  if (burn < 1000) {
    warning("You have specified the burn-in samples to be less than 1.000. This will likely yield unstable results.")
  }

  ##### End checks #####

  # Update settings
  if(!missing(chains)) {
    # Update # chains
    blm$sampling_settings$chains <- chains
  }
  if (!missing(iterations)) {
    blm$sampling_settings$iterations <- iterations
  }
  if (!missing(burn)) {
    blm$sampling_settings$burn <- burn
  }
  if (!missing(thinning)) {
    blm$sampling_settings$thinning <- thinning
  }

  # Return
  return(blm)

}

# Execute a blm plan
#' @export
sample_posterior.blm <- function(blm) {

  # number iterations & chains
  chains <- blm$sampling_settings$chains
  iterations <- blm$sampling_settings$iterations
  priors <- blm$priors
  burn <- blm$sampling_settings$burn

  # unroll data
  X <- blm$input$X
  y <- blm$input$y

  # Draw initial values for each chain
  initial_values <- lapply(1:chains, function(x) {

    # Draw initial values
    iv <- initialize_chain_values(priors)

    # Return
    return(iv)

  })

  # Add to sampling options
  blm$sampling_settings$initial_values <- initial_values
  # (names)
  names(blm$sampling_settings$initial_values) <- paste0("chain", 1:blm$sampling_settings$chains)

  # Run foreach (possibly parallel for each chain)
  posterior <- vector("list", chains)

  # For each list
  for(i in seq_along(1:chains)) {

    # Initial values for the current chain
    initial_values_current <- initial_values[[i]]

    # Call the gibbs sampler (helpers.R)
    posterior[[i]] <- gibbs_sampler(X, y, initial_values_current, iterations, priors, burn)

  }

  # Name output list
  names(posterior) <- paste0("chain_", 1:chains)

  # Name output data columns
  posterior <- lapply(posterior, function(x) {

    # Retrieve varname if listed in prior
    varnames <- unname(unlist(lapply(priors, function(x) x[["varname"]])))

    # Construct names
    cnames <- c("Intercept", varnames, "sigma")

    # Set colnames
    colnames(x) <- cnames

    # Return
    return(x)

  })

  # Append posterior
  blm$posterior <- posterior

  # Return results
  return(blm)

}

# Print method
#' @export
print.blm <- function(blm) {

  # Print information
  msg <- paste0(
    "Bayesian Linear Model (BLM) object:\n\n",
    "Data:\n",
      "\tPredictors: ", blm$input$m - 1, "\n",
      "\tOutcome: ", blm$input$variables$DV, "\n",
      "\tCentered: ", blm$input$center, "\n\n",
    "Sampler:\n",
      "\tChains: ", blm$sampling_settings$chains , "\n",
      "\tIterations: ", blm$sampling_settings$iterations , "\n",
      "\tBurn: ", blm$sampling_settings$burn , "\n\n")

  # Set up priors matrix
  pr_coef <- matrix(0L, ncol=length(blm$priors)-1, nrow = 2,
               dimnames=list(
                 c("mu", "tau"),
                 names(blm$priors)[-length(blm$priors)]
               ))
  pr_sigma <- matrix(0L, ncol=1, nrow=2,
                     dimnames = list(
                       c("rate", "scale"),
                       c("sigma")
                     ))

  # Add values
  for(i in seq_along(blm$priors)) {
    if(names(blm$priors)[i] == "sigma") {
      pr_sigma[1,1] <- blm$priors[[i]]$alpha
      pr_sigma[2,1] <- blm$priors[[i]]$beta
    } else {
      pr_coef[1,i] <- blm$priors[[i]]$mu
      pr_coef[2,i] <- blm$priors[[i]]$sd
    }

  }

  cat(msg)
  print.listof(list("Priors (Coefficients)" = pr_coef))
  print.listof(list("Priors (Residuals)" = pr_sigma))

}

# Summary method
#' @export
summary.blm <- function(blm) {

  var_names <- c(colnames(blm$input$X), "sigma")

  ### Formula
  form <- as.character(blm$input$formula)
  formchar <- paste0(form[2], " ", form[1], " ", form[3])

  ## Construct individual parts
  general <- paste0(
    "Formula: '", formchar, "'"
  )

  ## User has already sampled or not
  has_sampled <- paste0(
    "Sampled: ", "posterior" %in% names(blm)
  )

  ## num. observations + predictors
  obs <- list(
    "Sampling settings" = matrix(c(blm$input$n,blm$input$m - 1,blm$sampling_settings$chains,
                  blm$sampling_settings$iterations,blm$sampling_settings$burn),
                nrow = 1,
                dimnames = list(
                  c(""),
                  c("Obs.",
                    "Predictors",
                    "Chains",
                    "Iterations",
                    "Burn")
                )))

  ### If not sampled yet ...
  if(!"posterior" %in% names(blm)) {
    cat("Bayesian Linear Model (BLM) results:")
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
  MAPV <- do.call(rbind.data.frame, MAP(blm$posterior))
  # Add MC error
  MAPV[3,] <- MAPV[2,] / sqrt(blm$sampling_settings$iterations)
  # Round
  MAPV <- round(MAPV, digits = 3)

  # Amend names
  row.names(MAPV) <- c("Est.", "SE", "MCERR.")
  colnames(MAPV) <-var_names

  # Calculate CI
  CIV <- list(
    "95% credible interval" = round(CI(blm$posterior), digits = 3)
  )

  # Burn-in diagnostics
  burning_diag <- round(burnin_diagnostic(blm$posterior), digits=3)

  # Print MAP & SE
  cat("Bayesian Linear Model (BLM) results:")
  cat("\n\n")
  cat(general)
  cat("\n\n")
  cat(has_sampled)
  cat("\n\n")
  print.listof(obs)
  print.listof(list("Maximum a posteriori (MAP) estimates" = MAPV))
  print.listof(CIV)
  print.listof(list("Burn-in diagnostics" = burning_diag))

  # If multiple chains
  if(blm$sampling_settings$chains > 1) {
    # Calculate Gelman-Rubin
    GRS <- list(
      "Gelman-Rubin Statistic" = round(GR(blm$posterior,
                                          blm$sampling_settings$iterations), digits = 3)
    )
    print.listof(GRS)
  }

}

# Plot method
#' @export
plot.blm <- function(blm, type=c("history",
                                 "autocorrelation",
                                 "density"),
                     ...) {



  # Retrieve posterior data
  pd <- blm$posterior

  ## If autocorrelation plot
  if(type == "autocorrelation") {

    # Get opts
    opts <- list(...)

    # If not chain specified and multiple chains, choose chain 1 but raise error
    if(!"chain" %in% names(opts) & blm$sampling_settings$chains > 1) {

      warning("Choosing chain 1 for autocorrelation plot. You can choose a different chain by passing 'chain = <number>' to the plot() function.")

      chain <- 1

    } else {

      chain <- opts$chain

    }

    # Get data from posterior
    qd <- pd[[paste0("chain_", chain)]]

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

    # Bind data
    for(i in seq_along(pd)) {
      df <- data.frame(pd[[i]])
      df$chain <- i
      df$iteration <- (blm$sampling_settings$burn + 1):blm$sampling_settings$iterations
      pd[[i]] <- df
    }

    # To long format
    pd <- pd %>%
      do.call(rbind.data.frame, .) %>%
      tidyr::gather(key = parameter, value = value,
                    -chain, -iteration)

    #### Plot

    # Make into factor
    pd$chain <- as.factor(pd$chain)

    ## History plot

    if(type == "history") {

      ggplot2::ggplot(pd, ggplot2::aes(x = iteration,
                     y=value,
                     color=as.factor(chain),
                     group = parameter)) +
        ggplot2::geom_line(alpha=0.4) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "None") +
        #geom_smooth(method="lm") +
        ggplot2::facet_wrap("parameter ~ .", scales = "free_y",
                   ncol=1)

    ## Density plot

    } else if(type == "density") {

      ggplot2::ggplot(pd, ggplot2::aes(x=value,
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
coef.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg if user does not specify type
  type <- match.arg(type)

  # Bind posterior data
  pb <- do.call(rbind.data.frame, blm$posterior)
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
coefficients.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg
  type <- match.arg(type)

  # Call coef
  return(coef(blm, type))

}

# posterior predictive checks
# This returns a SEPARATE object ==> all them simulations are heavy on the memory.
#' @export
posterior_predictive_checks.blm <- function(blm,
                                            iterations = 2000, burn = 1000) {

  # Set up posterior predictive check by sampling from the posterior
  inputs <- list(
    settings = list(
      iterations = iterations,
      burn = burn
    )
  )

  # Get priors etc.
  priors <- blm$priors
  X <- blm$input$X
  y <- blm$input$y

  # Draw values
  iv <- initialize_chain_values(priors)

  # Call the gibbs sampler, simulate y values and compute residuals
  r <- ppc_julia(X, y, iv, iterations, priors, burn)

  # Add results
  inputs$data$initial_values <- iv
  inputs$data$X <- X
  inputs$data$sim_y <- r$sim_y
  inputs$data$residuals <- r$residuals

  # Add class to input
  class(inputs) <- "ppc"

  # Calculate statistics
  inputs <- normality_check(inputs)
  inputs <- homoskedast_check(inputs)

  # Return
  return(inputs)

}

# Model fit
#' @export
model_fit.blm <- function(blm) {

  # Model fit
  mfit <- DIC(blm$input$X, blm$input$y, blm$posterior)

  # Intercept-only model
  mfitIO <- DIC(matrix(blm$input$X[,1],ncol=1),
                blm$input$y, lapply(blm$posterior, function(x){
                  x[,c(1,ncol(x))]
  }))

  # Bind models
  final <- round(do.call(rbind.data.frame,
                         list(mfit, mfitIO)), digits=3)

  # Row names
  row.names(final) <- c("(Model)", "(Null model)")

  # Print
  print.listof(
    list("Model DIC"=final)
  )

}

# Evaluate method
#  (plots, statistics etc.)
#  use functions from other packages to do this

# Class 'Prior' methods -----

# Drawing from a normal distribution
#' @export
draw_value.normal <- function(prior) {

  rnorm(1, prior$mu, prior$sd)

}

# Drawing from a gamma distribution
# Add a small value to prevent sigma being 0 (bad stuff happens)
#' @export
draw_value.gamma <- function(prior) {

  rgamma(1, prior$alpha, prior$beta) + runif(1, 1e-10, 1e-08)

}

# Class 'ppc' methods -----

# Check normality of errors assumption
#' @export
normality_check.ppc <- function(ppc) {

  # Unroll data
  iterations <- ppc$settings$iterations
  burn <- ppc$settings$burn
  resids <- ppc$data$residuals

  # Effective iterations
  eff <- (iterations - burn)

  # 3. Calculate the skewness for each sample
  skewed <- matrix(0L, nrow = eff, ncol=2)
  colnames(skewed) <- c("simulated", "observed")

  # For each desired sample, calculate
  for(i in 1:eff) {

    # Compute skewness
    skewed[i,] <- c(e1071::skewness(resids[i,,1]),
                    e1071::skewness(resids[i,,2]))

  }

  # Compute where sim > obs (bayesian p-value)
  bpv <- mean(skewed[,1] > skewed[,2])

  # To list
  ppc$data$skewness <- skewed

  # Add results
  if(!"results" %in% names(ppc)) {
    ppc$results <- list(
      "normality" = bpv
    )
  } else {
    ppc$results$normality <- bpv
  }

  # Return
  return(ppc)

}

# Check homoskedasticity assumption
#' @export
homoskedast_check.ppc <- function(ppc) {

  # Unroll data
  iterations <- ppc$settings$iterations
  burn <- ppc$settings$burn
  resids <- ppc$data$residuals
  X <- ppc$data$X

  # Effective iterations
  eff <- (iterations - burn)

  # 3. Calculate the skewness for each sample
  homosked <- matrix(0L, nrow = eff, ncol=2)
  colnames(homosked) <- c("simulated", "observed")

  # For each desired sample, calculate
  for(i in 1:eff) {

    # Square residuals (no negative)
    ys <- resids[i,,1]^2
    yo <- resids[i,,2]^2

    # Calculate homoskedasticity
    homosked[i,1] <- test_for_homoskedasticity(ys, X)
    homosked[i,2] <- test_for_homoskedasticity(yo, X)

  }

  # Compute where sim > obs (bayesian p-value)
  bpv <- mean(homosked[,1] > homosked[,2])

  # To list
  ppc$data$homosked <- homosked

  # Add results
  if(!"results" %in% names(ppc)) {
    ppc$results <- list(
      "homosked" = bpv
    )
  } else {
    ppc$results$homosked <- bpv
  }

  # Return
  return(ppc)

}

# Summary
#' @export
summary.ppc <- function(ppc) {

  msg <- paste0(
    "Posterior Predictive Checks (PPC) for blm object:",
    "\n\n"
  )

  # Bind results
  bpr <- round(do.call(cbind.data.frame, post_pc$results), digits=3)
  colnames(bpr) <- c("Normality", "Heteroskedasticity")
  row.names(bpr) <- c("p")

  res <- list(
    "Bayesian p-value" = bpr
  )

  # Cat
  cat(msg)
  print.listof(res)

}

# Plot method
#' @export
plot.ppc <- function(ppc, type=c("normality", "heteroskedasticity", "rmse")) {

  # Get data
  data <- switch(
    type,
    "normality" = ppc$data$skewness,
    "heteroskedasticity" = ppc$data$homosked,
    "rmse" = ppc$data$rmse
  )

  # Plot data
  data %>%
    as.data.frame() %>%
    tidyr::gather(data, value) %>%
    ggplot2::ggplot(., ggplot2::aes(x=value, fill=data)) +
      ggplot2::geom_density(alpha=0.6) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste0("Observed and simulated results for test '", type, "'"))

}
