## Debugging julia code

# Functions -----

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


# Posterior tau helper function
gibbs_posterior_tau <- function(xj, sigma, tau0) {

  return( 1 / ( (t(xj) %*% xj/sigma) + (1/tau0) ) )

}

# Posterior rate for IG gibbs
gibbs_posterior_rate <- function(N, prior_rate) {

  return( (N / 2) + prior_rate )

}

# Posterior scale for IG gibbs
gibbs_posterior_scale <- function(X, y, w, prior_scale) {

  return(
    (sum((y - X %*% w)^2) / 2) + prior_scale
  )

}

# Debug -----

library(JuliaCall)

# Set up blm
.blm <- new.env(parent = emptyenv())
# Set up Julia
.blm$julia <- JuliaCall::julia_setup()
# Install Distributions/statistics package if needed
.blm$julia$install_package_if_needed("Distributions")
# Load Distributions/statistics package
.blm$julia$library("Distributions")
# Source gibbs sampler julia functions
.blm$julia$source("inst/julia/blm.jl")

# Set up blm object
#tjulia <- blm("Alcohol ~ .", data=wine2)
##w <- tjulia$sampling_settings$chain_1$initial_values$w
##w[,1] <- c(-1, 1.2, 3)
#sigma <- tjulia$sampling_settings$chain_1$initial_values$sigma

# Tuning parameter for proposal
zeta <- 0.25

# Draw x
x <- rnorm(500, mean=0, sd=3)
# Coef (intercept)
b0 <- 4.2
b1 <- 1.87
# Generate y
res_variance <- 2.1^2
y <- rnorm(500, mean = (b0 + b1*x), sd = sqrt(res_variance))

# Create data frame for data
dat <- data.frame(x=x, y=y)

library(blm)
blm_setup()
fit <- blm("y ~ x", data=dat) %>%
  set_sampler("b1", "MH", lambda=zeta) %>%
  set_sampling_options(burn = 2000, iterations = 50000, chains=2) %>%
  #set_initial_values(chain_1 = list('b'=c(-10,20), "sigma"=1)) %>%
  sample_posterior()

plot(fit, "history")
plot(fit, "autocorrelation")

fit <- fit %>%
  update_posterior(iterations = 30000) %>%
  set_sampling_options(burn = 30000)

fit %>%
  evaluate_effective_sample_size()

fit %>%
  summary()

# Priors
pr_mu_b0 <- fit$priors$b0$mu
pr_tau_b0 <- fit$priors$b0$sd
pr_mu_b1 <- fit$priors$b1$mu
pr_tau_b1 <- fit$priors$b1$sd
pr_rate_sig <- fit$priors$sigma$alpha
pr_scale_sig <- fit$priors$sigma$beta

# Create a model matrix
X <- fit$input$X
N <- nrow(X)

# Number of replications
k <- 20000
# Results
posterior_sampled <- matrix(0, nrow=k, ncol=3)

# First row of the posterior_sampled matrix are the initial values
posterior_sampled[1,] <- c(fit$sampling_settings$chain_1$initial_values$w,
                           fit$sampling_settings$chain_1$initial_values$sigma)

# We are not resetting results matrix etc. but making a new results matrix
posterior_sampled_mh <- matrix(0, nrow=k, ncol=3)
posterior_sampled_mh[1,] <- posterior_sampled[1,]

# Keep track of accepted draws for MH
accepted <- 0

## Iterate
for( i in 2:k ) {

  ## Retrieve w (coefficients)
  # Ignore the last value ==> this is for sigma
  w <- posterior_sampled_mh[i-1, -3]
  sigma <- posterior_sampled_mh[i-1, 3]

  ## For each parameter, retrieve posterior samples using a mix of gibbs / MH only

  #### b[0] ==> Metropolis-Hastings

  b0_prev <- w[1]
   # w[1] <- MH_one_step(w[1], matrix(X[,1], ncol=1), matrix(X[,-1], ncol=1),
   #                     y, w[-1], sigma, pr_mu_b0, pr_tau_b0, zeta)

  r <- .blm$julia$eval("MH_one_step_coef")(w[1], # b_previous
                                              matrix(X[,1], ncol=1), # xj
                                              matrix(X[,-1], ncol=1), # X
                                              y, # y
                                              fit$sampling_settings$chain_1$initial_values$w[-1], # w
                                              fit$sampling_settings$chain_1$initial_values$sigma, # sigma
                                              fit$priors$b0$mu, # prior mu
                                              fit$priors$b0$sd, # prior tau
                                              zeta) # zeta

  w[1] <- r$samp

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

# square-root the posterior sigma
posterior_sampled_mh[,3] <- sqrt(posterior_sampled_mh[,3])

# As coda mcmc
library(coda)
pd <- as.mcmc(posterior_sampled_mh)

# Plot
plot(pd)

# Summary
summary(pd)

# Acceptance rate
accepted / k

convergence_diagnostics(fit)

plot(fit, "history")

os <- .blm$julia$eval("MCMC_one_step")(X, y,
                                       fit$sampling_settings$chain_1$initial_values$w,
                                       fit$sampling_settings$chain_1$initial_values$sigma,
                                       unname(fit$priors),
                                       unname(fit$sampling_settings$chain_1$samplers))
