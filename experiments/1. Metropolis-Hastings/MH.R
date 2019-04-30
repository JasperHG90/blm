## MH script

# Correlation ----

x <- scale(rnorm(1377, mean=0, sd=1), center=TRUE, scale=TRUE)[,1]
y <- scale(rnorm(1377, mean=0, sd=1), center=TRUE, scale=TRUE)[,1]
#y <- scale(x + runif(1377, -2, 2), center=TRUE, scale=TRUE)[,1]

# Logged conditional posterior distribution
lnpost <- function(r) {
  -690*log(1-r^2) -.5*(sum(x^2)-2*r*sum(x*y)+sum(y^2))/(1-r^2)
}

# Results matrix
k <- 10000
corr <- matrix(0, k)
# Accepted draws
accepted <- 0

# Loop over the k draws
# (note that we start at 2 [==>] this means that the starting value is implicitly set to 0)
for(i in 2:k) {

  # Proposal distribution
  # Update the previous value with a uniform proposal density
  corr[i] <- corr[i-1] + runif(1, min = -.07, max=.07)

  # If the absolute value of the current draw is larger than one, reject the sample
  # (this is a check to see whether the current draw is bounded by -1 and 1)
  if( abs(corr[i]) > 1 ) {

    corr[i] <- corr[i-1]

  }

  # Evaluate the log-likelihood of the function
  ll_corr_current <- lnpost(corr[i])
  ll_corr_previous <- lnpost(corr[i-1])

  # Draw a uniform and log (since the likelihood is logged also)
  u <- log(runif(1, min=0, max=1))

  # Accept or reject based on the difference between corr_current and corr_previous
  # (why difference? Because we logged the values ==> log(a/b) = log(a) - log(b))
  ll_diff <- ll_corr_current - ll_corr_previous

  # Reject the sample if ll_diff < u
  if( ll_diff < u ) {

    # Set current value to previous value
    corr[i] <- corr[i-1]

    # Continue
    next

  }

  # Set accepted + 1
  accepted <- accepted + 1

  # Every 100 draws, print progress to console
  if ( (i %% 100) == 0 ) {
    print(c(i,corr[i],accepted/i))
  }

}

hist(corr)
summary(corr)

# Linear regression using Gibbs only ----

## Plan: model b[0] + b[1]*x + e

# Draw x
x <- rnorm(100, mean=0, sd=3)
# Coef (intercept)
b0 <- 4.2
b1 <- 1.87
# Generate y
res_variance <- 2.1^2
y <- rnorm(100, mean = (b0 + b1*x), sd = sqrt(res_variance))

# Create data frame for data
dat <- data.frame(x=x, y=y)

## Create a model matrix
X <- model.matrix(y ~ x, data=dat)
N <- nrow(X)

# Number of replications
k <- 20000
# Results
posterior_sampled <- matrix(0, nrow=k, ncol=3)

# First row of the posterior_sampled matrix are the initial values
posterior_sampled[1,] <- c(2, -2, 0.01)

# Posterior mean helper function
gibbs_posterior_mu <- function(X, xj, y, w, sigma, mu0, tau0) {

  # Numerator
  num1 <- sum(xj * (y - X %*% w)) / sigma
  num2 <- mu0 / tau0

  # Denominator
  den1 <- t(xj) %*% xj / sigma
  den2 <- 1 / tau0

  # Return
  return( ( num1 + num2 ) / ( den1 + den2 ) )

}

# Posterior tau helper function
gibbs_posterior_tau <- function(xj, sigma, tau0) {

  return( 1 / ( (t(xj) %*% xj/sigma) + (1/tau0) ) )

}

## Posterior rate for IG gibbs
gibbs_posterior_rate <- function(N, prior_rate) {

  return( (N / 2) + prior_rate )

}

## Posterior scale for IG gibbs
gibbs_posterior_scale <- function(X, y, w, prior_scale) {

  return(
    (sum((y - X %*% w)^2) / 2) + prior_scale
  )

}

## Priors
pr_mu_b0 <- 0
pr_tau_b0 <- 1000
pr_mu_b1 <- 0
pr_tau_b1 <- 1000
pr_rate_sig <- 0.001
pr_scale_sig <- 0.001

## Iterate
for( i in 2:k ) {

  ## Retrieve w (coefficients)
  # Ignore the last value ==> this is for sigma
  w <- posterior_sampled[i-1, -3]
  sigma <- posterior_sampled[i-1, 3]

  ## For each parameter, retrieve posterior samples using gibbs only

  #### b[0]

  w[1] <- rnorm(1,
                mean = gibbs_posterior_mu(matrix(X[,-1], ncol=1), matrix(X[,1], ncol=1),
                                          y, w[-1], sigma, pr_mu_b0, pr_tau_b0),
                sd = sqrt(gibbs_posterior_tau(X[,1], sigma, pr_tau_b0)))

  #### b[1]

  w[2] <- rnorm(1,
                mean = gibbs_posterior_mu(matrix(X[,-2], ncol=1), matrix(X[,2], ncol=1),
                                          y, w[-2], sigma, pr_mu_b0, pr_tau_b0),
                sd = sqrt(gibbs_posterior_tau(X[,2], sigma, pr_tau_b0)))

  #### sigma

  sigma <- 1 / rgamma(1,
                      gibbs_posterior_rate(N, pr_rate_sig),
                      gibbs_posterior_scale(X, y, w, pr_scale_sig))

  # Update results matrix
  posterior_sampled[i,] <- c(w, sigma)

}

# square-root the posterior sigma
posterior_sampled[,3] <- sqrt(posterior_sampled[,3])

# MAP values
apply(posterior_sampled, 2, mean)

# Plots
library(ggplot2)
library(dplyr)
library(tidyr)
density_plot <- function(x) {

  ## To data frame and names
  x <- as_data_frame(x)
  colnames(x) <- c("b0", "b1", "residual_variance")

  ## Burn 1.000
  x <- x[-2:-2000,]

  ## To long format
  x %>%
    gather(parameter, value) %>%
    # Plot
    ggplot(., aes(x=value, group=parameter)) +
    geom_density() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    facet_wrap("parameter ~ .", scales="free_y")
}

# Trace plot
trace_plot <- function(x) {

  ## To data frame and names
  x <- as_data_frame(x)
  colnames(x) <- c("b0", "b1", "residual_variance")

  ## To long
  x %>%
    mutate(index = 1:k) %>%
    gather(parameter, value, -index) %>%
    ## Plot trace
    ggplot(., aes(x=index, y=value, group=parameter)) +
      geom_line() +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      facet_grid("parameter ~ .", scales = "free_y")

}

# Autocorrelation
# Helper function for autocorrelation
autocor <- function(x, n=10) {

  # Results
  res <- rep(0, n)

  # Lag for each n and calculate correlation
  for(i in 1:n) {
    res[i] <- cor(x, c(rep(NA, i), x[1:(length(x)-i)]),
                  use="complete.obs")
  }

  # Return
  return(
    data.frame(
      "lag" = 1:n,
      "correlation" = res
    )
  )

}

# Plot
autocor_plot <- function(x) {

  ## To data frame and names
  x <- as_data_frame(x)
  colnames(x) <- c("b0", "b1", "residual_variance")

  ## Compute autocorrelation
  qd <- x %>%
    as.data.frame() %>%
    lapply(., function(y) autocor(y, n=40)) %>%
    do.call(rbind.data.frame, .)

  # Add variable name
  qd$id <- stringr::str_replace_all(row.names(qd), "\\.[0-9]{1,2}", "")

  # Plot
  ggplot2::ggplot(qd, ggplot2::aes(x=lag, y=correlation, group=id)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_x_continuous(breaks= scales::pretty_breaks()) +
    ggplot2::scale_y_continuous(limits = c(-1,1)) +
    ggplot2::facet_wrap(id ~ .)

}

## Plot
density_plot(posterior_sampled)
trace_plot(posterior_sampled)
autocor_plot(posterior_sampled)

# Linear regression using Gibbs/MH mix ----

# Plan: sample b[0] using MH, sample b[1] and sigma using Gibbs

# Posterior density
# i.e. the posterior density proportional to a normal up to a constant
# It will be logged to simplify working with the exponents
MH_posterior_coef <- function(bj, xj, X, y, w, sigma, prior_mu, prior_tau) {

  ## Left side
  left <- -(bj^2) * ( ((t(xj) %*% xj)/(2*sigma)) + (1 / (2*prior_tau)))

  ## Right side
  right <- bj * ( ( (sum(xj * (y - X %*% w))) / (sigma) ) + (prior_mu / prior_tau) )

  ## Add left and right
  return(left + right)

}

# One step of MH algorithm
MH_one_step <- function(b_previous, xj, X, y, w, sigma, prior_mu, prior_tau, zeta) {

  ## Proposal density
  b_proposed <- rnorm(1, mean=b_previous, sd=zeta)

  ## Draw from the posterior
  g_proposed <- MH_posterior_coef(b_proposed, xj, X, y, w, sigma, prior_mu, prior_tau)
  g_previous <- MH_posterior_coef(b_previous, xj, X, y, w, sigma, prior_mu, prior_tau)

  ## Difference (because logged)
  g_diff <- g_proposed - g_previous

  ## Draw random variate from uniform
  u <- log(runif(1, min=0, max=1))

  ## If gdiff < u ==> set b_current to b_previous
  if( g_diff < u ) {

    b_current <- b_previous

  } else { # Accepted the proposed value

    b_current <- b_proposed

  }

  # Return b
  return(b_current)

}

# We are not resetting results matrix etc.
posterior_sampled_mh <- matrix(0, nrow=k, ncol=3)
posterior_sampled_mh[1,] <- posterior_sampled[1,]

# Keep track of accepted draws for MH
accepted <- 0

# Tuning parameter for proposal density
zeta <- 0.25

## Iterate
for( i in 2:k ) {

  ## Retrieve w (coefficients)
  # Ignore the last value ==> this is for sigma
  w <- posterior_sampled_mh[i-1, -3]
  sigma <- posterior_sampled_mh[i-1, 3]

  ## For each parameter, retrieve posterior samples using a mix of gibbs / MH only

  #### b[0] ==> Metropolis-Hastings

  b0_prev <- w[1]
  w[1] <- MH_one_step(w[1], matrix(X[,1], ncol=1), matrix(X[,-1], ncol=1),
                      y, w[-1], sigma, pr_mu_b0, pr_tau_b0, zeta)

  # If accepted, increment
  if(w[1] != b0_prev) {
    accepted <- accepted+1
  }

  #### b[1] ==> Gibbs

  w[2] <- rnorm(1,
                mean = gibbs_posterior_mu(matrix(X[,-2], ncol=1), matrix(X[,2], ncol=1),
                                          y, w[-2], sigma, pr_mu_b0, pr_tau_b0),
                sd = sqrt(gibbs_posterior_tau(X[,2], sigma, pr_tau_b0)))

  #### sigma ==> Gibbs

  sigma <- 1 / rgamma(1,
                      gibbs_posterior_rate(N, pr_rate_sig),
                      gibbs_posterior_scale(X, y, w, pr_scale_sig))

  # Update results matrix
  posterior_sampled_mh[i,] <- c(w, sigma)

}

# Accepted values
accepted / k

# square-root the posterior sigma
posterior_sampled_mh[,3] <- sqrt(posterior_sampled_mh[,3])

# MAP values
apply(posterior_sampled_mh, 2, mean)
apply(posterior_sampled, 2, mean)

# Plots
density_plot(posterior_sampled_mh)
trace_plot(posterior_sampled_mh)
autocor_plot(posterior_sampled_mh)

# Compare densities
data_frame(
  "Gibbs" = posterior_sampled[,1],
  "MH" = posterior_sampled_mh[,1]
) %>%
  gather(sampler, value) %>%
  ggplot(., aes(x=value, fill=sampler)) +
    geom_density(alpha=0.3)
