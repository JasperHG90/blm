rm(list=ls())

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

## posterior coefficient
posterior_mh_coef <- function(bj, xj, X, y, w, sigma, prior_mu, prior_tau) {

  # Left
  left <- -(bj)^2 * ((sum(t(xj) %*% xj)/(2*sigma))+(1/(2*prior_tau)))

  # Right
  right <- bj* ((sum(xj * (y - X %*% w)) / sigma) + (prior_mu / prior_tau) )

  return(left + right)

}

# Perform one step of metropolis-hastings algorithm
mh_one_step <- function(b_previous, zeta, xj, X, y, w, sigma, prior_mu, prior_tau) {

  # Draw proposal
  b_proposal <- rnorm(1, mean = b_previous, sd = zeta)

  # Sample random value from uniform ==> log it
  u <- log(runif(1, min=0, max=1))

  # Proposal
  p_previous <- posterior_mh_coef(b_previous, xj, X, y, w, sigma, prior_mu, prior_tau)
  p_current <- posterior_mh_coef(b_proposal, xj, X, y, w, sigma, prior_mu, prior_tau)

  # Acceptance ratio r
  r <- p_current - p_previous

  # If r less than u, set b_current to b_previous. Else, set b_current to b_proposal and accept the draw
  if( r < u ) {
    b_current <- b_previous
  } else {
    b_current <- b_proposal
  }

  # Return
  return(b_current)

}

# Priors
pr_mu_b0 <- 0
pr_tau_b0 <- 1000
pr_mu_b1 <- 0
pr_tau_b1 <- 1000
pr_rate_sig <- 0.001
pr_scale_sig <- 0.001

# Draw x
x1 <- rnorm(500, mean=25, sd=5) + runif(500, -4, 4)
x2 <- rnorm(500, mean=-50, sd=15) + runif(500, -6, 6)
# Coef (intercept)
b0 <- 4.2
b1 <- 1.87
b2 <- 8.2
# Generate y
res_variance <- 500^2
y <- rnorm(500, mean = (b0 + b1*x1 +b2*x2), sd = sqrt(res_variance))

# Create data frame for data
dat <- data.frame(x1=x1,x2=x2 , y=y)
dvar <- apply(dat, 2, var)
dmean <- apply(dat, 2, mean)
dat <- as.data.frame(scale(dat,center=TRUE, scale=TRUE))
y <- dat$y

# Load mtcars
# data("mtcars")
#
# round(apply(mtcars, 2, var), digits=3)
# dat <- data.frame(y = mtcars$mpg, x1=mtcars$hp, x2=mtcars$qsec)
#
# # Create a model matrix
# dat <- as.data.frame(scale(dat, center=TRUE, scale=TRUE))
# y <- dat$y
#
# # Load wine data
# # Load data
 d <- rattle.data::wine %>% as.data.frame()
#
# # Run linear regression on each numeric variable to see which variable can be predicted from the others
 cn_num <- colnames(d)[-1]
#
# # Center data
 d2 <- d
 d2[,sapply(d2, is.numeric)] <- scale(d2[,sapply(d2, is.numeric)], center = TRUE, scale=TRUE)
#
# # Subset
 dat <- data.frame(y = d2$Proline, x1=d2$Color, x2=d2$Hue)
 y <- dat$Proline

# Create model matrix
X <- model.matrix(y ~ ., data=dat)
N <- nrow(X)

# Diamonds data
data("trees")
wine <- rattle.data::wine

hist(wine$Alcohol)
hist(wine$Ash)
hist(wine$Alcalinity)
hist(wine$Phenols)
hist(wine$Nonflavanoids)

dat <- data.frame(
  "y" = wine$Alcohol,
  "x1" = wine$Ash,
  "x2" = wine$Alcalinity
)

dat$x1 <- dat$x1 - mean(dat$x1)
dat$x2 <- dat$x2 - mean(dat$x2)

# Create model matrix
X <- model.matrix(y ~ ., data=dat)
N <- nrow(X)
y <- dat$y

# Number of replications
k <- 20000
# Results
posterior_sampled <- matrix(0, nrow=k, ncol=4)

# First row of the posterior_sampled matrix are the initial values
posterior_sampled[1,] <- c(2, -2, 5, 0.01)

accepted <- 0

# Value for tau / zeta
zeta <- 5

# Iterate
for( i in 2:k ) {

  # Retrieve w (coefficients)
  # Ignore the last value ==> this is for sigma
  w <- posterior_sampled[i-1, -4]
  sigma <- posterior_sampled[i-1, 4]

  # For each parameter, retrieve posterior samples using gibbs only

  #### b[0] ==> Metropolis-hastings algorithm

  # # Store previous value of bj in a vector
  #  b_prev <- w[1]
  # #
  # # # One step of the MH algorithm to update bj
  #  w[1] <- mh_one_step(w[1], zeta, matrix(X[,1], ncol=1),
  #                      as.matrix(X[,-1], ncol=1), y, w[-1],
  #                      sigma, pr_mu_b0, pr_tau_b0)
  #
  # # # Update acceptance tally if draw accepted
  #  if(b_prev != w[1]) {
  #    accepted <- accepted + 1
  #  }

  w[1] <- rnorm(1,
                mean = gibbs_posterior_mu(X[,-1], matrix(X[,1], ncol=1),
                                          y, w[-1], sigma, pr_mu_b0, pr_tau_b0),
                sd = sqrt(gibbs_posterior_tau(X[,1], sigma, pr_tau_b0)))

  #### b[1]

  # w[2] <- rnorm(1,
  #               mean = gibbs_posterior_mu(X[,-2], matrix(X[,2], ncol=1),
  #                                         y, w[-2], sigma, pr_mu_b0, pr_tau_b0),
  #               sd = sqrt(gibbs_posterior_tau(X[,2], sigma, pr_tau_b0)))

  # Store previous value of bj in a vector
  b_prev <- w[2]
  #
  # # One step of the MH algorithm to update bj
  w[2] <- mh_one_step(w[2], zeta, matrix(X[,2], ncol=1),
                      X[,-2], y, w[-2],
                      sigma, pr_mu_b1, pr_tau_b1)

  # # Update acceptance tally if draw accepted
  if(b_prev != w[1]) {
    accepted <- accepted + 1
  }

  #### b[2]

  w[3] <- rnorm(1,
                mean = gibbs_posterior_mu(X[,-3], matrix(X[,3], ncol=1),
                                          y, w[-3], sigma, pr_mu_b0, pr_tau_b0),
                sd = sqrt(gibbs_posterior_tau(X[,3], sigma, pr_tau_b0)))

  #### sigma

  sigma <- 1 / rgamma(1,
                      gibbs_posterior_rate(N, pr_rate_sig),
                      gibbs_posterior_scale(X, y, w, pr_scale_sig))

  # Update results matrix
  posterior_sampled[i,] <- c(w, sigma)

}

library(ggplot2)
ggplot(dat, aes(x=x1, y=y)) +
  geom_point()

p <- plot_ly(dat, x=~x1, y=~x2, z=~y, size=3) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Color'),
                      yaxis = list(title = 'Hue'),
                      zaxis = list(title = 'Proline')))
p
library(plotly)

accepted / k

posterior_sampled[,4] <- sqrt(posterior_sampled[,4])

library(coda)
posterior <- as.mcmc(posterior_sampled)
# MAP
round(apply(posterior, 2, mean), digits=5)

round(c(coef(lm("y~.", data=dat)), sd(resid(lm("y~.", data=dat)))), digits=5)

plot(as.mcmc(posterior[-1:-1000,]), density=FALSE)

effectiveSize(posterior)

autocorr.plot(posterior)
