## Jags file
library(rjags)

linm <- "model{
  #### Priors
  # Sigma^2
  tau ~ dgamma(0.001, 0.001)

  # Coefficients
  for(coef in 1:2) {
    b[coef] ~ dnorm(0, 1/100)
  }

  #### Likelihood ==> for each example in the data SEPARATELY
  for(i in 1:length(Compensation)) {
    # Calculate the mean value
    mu[i] <- b[1] + b[2] * Age[i]
    # Use dnorm to calculate the density at the mean value
    Compensation[i] ~ dnorm(mu[i], tau)
  }

  #### Variables of interest
  sigma2 <- 1/tau
  sigma <- sqrt(sigma2)
}"
# Specify initial values for the chains
# (this is optional)
init <- list(
  init1 <- list(b = c(7, -5), tau=1),
  init2 <- list(b = c(-5, 7), tau=2)
)
# Text connection
con <- textConnection(linm)
# Jags model
library(rjags)
d <- list(
  "Compensation" = directors_scaled$Compensation,
  "Age" = directors_scaled$Age
)
jmod <- jags.model(file=con, data=d, n.chains=2, inits = init)
close(con)
# Run the model
update(jmod, n.iter = 5000)
# Specify parameter names
params <- c("b", "sigma")
# Run the chain
res <- coda.samples(jmod, variable.names = params, n.iter=20000)
