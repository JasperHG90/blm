## Methods for the posterior class

# Burn samples from a posterior sample dataset
burn.posterior <- function(x) {

  # Return posterior without burn-in samples
  x <- get_value(x, "samples") %>%
    lapply(., function(y) {
      y[-1:-get_value(x, "burn"),]
    })

  # Return
  return(x)

}

# Collapse all chains into one big sample
bind.posterior <- function(x) {

  # Burn
  burn(x) %>%
    # Bind
    do.call(rbind.data.frame, .)

}

# Append samples to a posterior distribution
append_samples.posterior <- function(x, updates) {

  # Get samples
  samples <- get_value(x, "samples")
  new_samples <- get_value(updates, "samples")

  # Add the number of accepted values
  accepts <- get_value(x, "accepted")
  new_accepts <- get_value(updates, "accepted")

  # Combine the posterior samples
  for(i in seq_along(x[["samples"]])) {

    # Append data
    x[["samples"]][[i]] <- rbind(samples[[i]], new_samples[[i]])
    x[["accepted"]][[i]] <- accepts[[i]] + new_accepts[[i]]

  }

  # Return
  return(x)

}

# MAP estimates for posterior distribution
MAP.posterior <- function(x) {

  # Join the data row-wise
  D <- bind(x)

  # Get posterior MAP values
  MAP_R <- apply(D, 2, mean)
  SE_R <- apply(D, 2, sd)

  # Return
  return(
    list(
      "MAP" = MAP_R,
      "SE" = SE_R
    )
  )

}

# Calculate the 95% CCI
CCI.posterior <- function(x) {

  # Bind data row-wise
  D <- bind(x)

  # Calculate 95% CI
  post_credint <- apply(D, 2,
                        function(y) quantile(y, c(0.025, 0.25, 0.5, 0.75, 0.975)))

  # Return
  post_credint

}

# Calculate the effective sample size for a MCMC sample
#' @importFrom stats acf
effss.posterior <- function(x, order=30) {

  # Bind data
  x <- bind(x)

  ## Compute estimates of autocorrelation
  autocorr <- stats::acf(x, lag.max = order, type =c("correlation"),
                         plot = FALSE, na.action = na.fail, demean = TRUE)

  ## Get n
  n <- nrow(x)

  ## Results matrix
  ss <- rep(0, ncol(x))

  ## Populate
  for(i in 1:ncol(x)) {
    ## SUBSET: 2 - order (input), ith variable, ith variable.
    ss[i] <- round(min(sample_size(n=n, pk = autocorr$acf[2:order, i, i]), n))
  }

  ## Return sample size
  return(ss)

}

# Gelman-rubin statistic
#https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
# Compare the within / between variances
GR.posterior <- function(x, iterations) {

  ## Burn
  x <- set_value(x, "samples", burn(x))

  ## For each chain, and each statistic, calculate the between and within
  ## variance
  M <- get_value(x, "samples") %>%
    length()
  N <- iterations
  m <- get_value(x, "samples")[["chain_1"]] %>%
    ncol()

  ## Open up results matrix
  RM <- matrix(0L, ncol = m,
               nrow = 1,
               dimnames = list(
                 c(""),
                 colnames(get_value(x, "samples")[["chain"]] %>%
                            colnames())
               ))

  ## For each parameter
  for(j in seq_along(1:m)) {

    # Get chain means
    cm <- matrix(unname(
      sapply(x[["samples"]], function(y) mean(y[,j]))
    ), ncol=1)

    # Get chain vars
    cv <- matrix(unname(
      sapply(x[["samples"]], function(y) var(y[,j]))
    ), ncol=1)

    # Subtract parameter mean over all chains
    cmm <- cm - mean(cm[,1])

    # Between variance
    B <- (N / (M - 1)) * t(cmm) %*% cmm

    # Within variance
    W <- (1/M) * sum(cv)

    # Pooled variance
    V <- (((N-1) / N) * W) + (((M+1) / (M*N)) * B)

    # GR stat
    RM[1, j] <- V / W

  }

  # Return
  return(RM)

}

# Burn-in period diagnostics
# Run a linear regression on squared coefficient draws from posterior
burnin_diagnostic.posterior <- function(x) {

  # Burn
  x <- set_value(x, "samples", burn(x))

  out_post <- lapply(seq_along(x[["samples"]]), function(chain_it) {

    # Get current chain
    chain <- x[["samples"]][[chain_it]]

    # For each column, compute linear coef
    out <- lapply(seq_along(1:ncol(chain)), function(y) {

      df <- data.frame(
        "y" = chain[,y]^2,
        "index"= (1:nrow(chain))
      )

      # Linear reg
      linr <- lm("y ~ index", data=df)

      # Get coef
      linr$coefficients["index"]

    })

    # Name out
    names(out) <- colnames(chain)

    # Bind
    out <- do.call(cbind.data.frame, out)
    row.names(out) <- paste0("chain ", chain_it)

    # Return
    return(out)

  })

  # Bind
  diagnostics <- do.call(rbind.data.frame, out_post)

  # Return
  return(diagnostics)

}
