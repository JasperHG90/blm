# READ DATA ----

rm(list=ls())
library(readr)
d <- read_csv("Appendices/boston-housing/train.csv")
y <- d$medv
X <- cbind(rep(1, nrow(d)), d$rm, d$crim, d$chas, d$lstat, d$dis)

# Generate data
source("code/helpers.R")
source("code/posterior_densities.R")
#d <- generate_dataset(center=TRUE, heteroskedastic = FALSE, degree = 20, seed=100)
#X <- d$X
#y <- d$y
#real <- d$parameters

# PREPROCESS ----

# Reference linear regression
library(magrittr)
lrd <- cbind(y,X[,-1]) %>% as.data.frame()
linreg <- lm("y ~ .", data=lrd)
summary(linreg)
plot(predict(linreg), resid(linreg))

# RUN WITH GIBBS SAMPLER -----

# Starting values
sigma <- rgamma(1, 0.001, 0.001) + runif(1, 1e-10, 1e-09) 
b0 <- rnorm(1,0,1000) / sqrt(nrow(X)) #0
b1 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b2 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b3 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b4 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b5 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
# Priors
priors <- list(
  "b0" = list("mu"=0, "tau"=1000),
  "b1" = list("mu"=0, "tau"=1000),
  "b2" = list("mu"=0, "tau"=1000),
  "b3" = list("mu"=0, "tau"=1000),
  "b4" = list("mu"=0, "tau"=1000),
  "b5" = list("mu"=0, "tau"=1000),
  "sigma" = list("alpha"=0.001, "beta"=0.001)
)
# To column vector
w <- matrix(c(b0,b1,b2,b3,b4,b5), ncol=1)
# Iterations
k <- 10000
# Matrix to store data
results <- matrix(0L, nrow=k, ncol = length(priors))
# For each
for(i in 1:(k+1)) {
  
  # Save values from previous iteration
  if(i > 1) {
    results[i-1,1:6] <- w[,1]
    results[i-1,7] <- 1/sqrt(sigma)
  }
  
  # Update the params
  for(j in 1:6) {
    
    # Get posterior
    posterior_values <- posterior_coef(X, y, w, j, priors[[paste0("b", (j-1))]], sigma)
    
    # Draw from posterior
    w[j,1] <- rnorm(1, posterior_values$mu_posterior, posterior_values$tau_posterior)
    
  }
  
  # Get posterior for sigma
  posterior_sigma_values <- posterior_sigma(X, y, w, priors[["sigma"]])
  
  # Draw
  sigma <- 1 / rgamma(1, posterior_sigma_values$alpha_posterior, 
                      posterior_sigma_values$beta_posterior)
  
}

burned <- 1500

# Remove the first 1.000 observations (burn-in period)
results <- results[-1:-burned,]

#

# EVALUATION ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Caterpilar plot
results %>%
  as_tibble() %>%
  mutate(iteration = (burned+1):(n() + burned)) %>%
  gather(variable, value, -iteration) %>%
  mutate(variable = factor(variable, labels = c(paste0("b", 1:6), "sigma"))) %>%
  ggplot(aes(x=iteration, y=value, group=variable, color=variable)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  facet_wrap("variable ~ .", scales = "free_y")

# Distribution of the posterior draws
results %>%
  as_tibble() %>%
  mutate(iteration = 1:n()) %>%
  gather(variable, value, -iteration) %>%
  mutate(variable = factor(variable, labels = c(paste0("b", 1:6), "sigma"))) %>%
  ggplot(aes(x=value, group=variable, fill=variable)) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "None") +
  facet_wrap("variable ~ .", scales = "free")

# Results
coef <- apply(results, 2, median)
se <- apply(results, 2, sd)

# For each ...
for(i in 1:7) {
  
  msg <- paste0(
    ifelse(i == 7, "sigma: ", paste0("coef: b",i-1)), 
    "\n",
    "BLM Est. ", round(ifelse(i == 7, coef[i], coef[i]), digits=4), "\n",
    "LR Est. ", round(ifelse(i == 7, 1/sqrt(var(residuals(linreg))), 
                             unname(coefficients(linreg)[i])), digits=4), "\n",
    #"Actual: ", round(real[[i]], digits=4), "\n",
    "SE: ", round(se[i], digits=5)
  )
  
  # Print
  cat(msg)
  cat("\n\n")
  
}

# Get point estimates -----
# See https://stattrek.com/matrix-algebra/sums-of-squares.aspx

# Posterior mean/median (EAP)
Exp_a_post_mean <- apply(results,2,mean)
Exp_a_post_med <- apply(results,2,median)

# Posterior standard deviation
post_sd <- apply(results,2,sd)

# Correlation matrix
post_cm <- cor(results)

# Var-covar matrix
h <- matrix(rep(1, nrow(results)), ncol=1)
v <- (results - ((h %*% t(h) %*% results) * (1/nrow(results))))
post_vc_matrix <- (t(v) %*% v) * (1/nrow(results))

# Calculate 95% credible interval
post_credint <- apply(results,2, function(x) quantile(x, c(0.025, 0.975)))

# MC error
post_mcerror <- post_sd / sqrt(k)

# Print
cat(Exp_a_post_mean)
cat(Exp_a_post_med)
print(post_credint)
cat(post_sd)
cat(post_mcerror)
post_cm

# Make a helper function 
autocor <- function(x, n=10) {
  
  # Results
  res <- rep(0, n)
  
  # Lag for each n and calculate correlation
  for(i in 1:n) {
    res[i] <- cor(x, c(rep(NA, i), x[1:(length(x)-i)]), use="complete.obs")
  }
  
  # Return
  return(
    data.frame(
      "lag" = 1:n,
      "correlation" = res
    )
  )
  
}

# Plot autocorrelation per variable 
library(purrr)
library(scales)
results %>%
  as.data.frame() %>%
  map(., function(x) autocor(x, n=40)) %>%
  dplyr::bind_rows(.id="id") %>%
  ggplot(., aes(x=lag, y=correlation, group=id)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks= pretty_breaks()) +
  scale_y_continuous(limits = c(-1,1)) + 
  facet_wrap(id ~ .)

## Bayesian p-value -----

# Testing the normality assumption by calculating skewness in the data

# 1. Sample the posterior +- 1000 times using Gibbs procedure

# Starting values 

# Jeffreys prior + small random value to rule out 0
sigma <- rgamma(1, 0.001, 0.001) + runif(1, 1e-10, 1e-09) 
b0 <- rnorm(1,0,1000) / sqrt(nrow(X)) #0
b1 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b2 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b3 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b4 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
b5 <- rnorm(1,0,1000) / sqrt(nrow(X))#0
# Iterations
ks <- 2000
# Matrix to store data
results_sim <- matrix(0L, nrow=ks, ncol = length(priors))
# For each
for(i in 1:(ks+1)) {
  
  # Save values from previous iteration
  if(i > 1) {
    results_sim[i-1,1:6] <- w[,1]
    results_sim[i-1,7] <- 1/sqrt(sigma)
  }
  
  # Update the params
  for(j in 0:5) {
    
    # Get posterior
    posterior_values <- posterior_coef(X, y, w, j, priors[[paste0("b", j)]], sigma)
    
    # Draw from posterior
    w[j+1,1] <- rnorm(1, posterior_values$mu_posterior, posterior_values$tau_posterior)
    
  }
  
  # Get posterior for sigma
  posterior_sigma_values <- posterior_sigma(X, y, w, priors[["sigma"]])
  
  # Draw
  sigma <- 1 / rgamma(1, posterior_sigma_values$alpha_posterior, 
                      posterior_sigma_values$beta_posterior)
  
}

# Burn first 1000
results_sim <- results_sim[-1:-1000,]

# Left-multiply X with sample T
# Precompute these values ==> Xw = yhat for each simulation
lincom <- X %*% t(results_sim[,-7])

# 2. Calculate the residual values for each of the samples
resids <- array(0L, dim = c(nrow(X), 1000, 2))

# 3. Calculate the skewness for each sample
skewed <- matrix(0L, nrow = 1000, ncol=2)
colnames(skewed) <- c("simulated", "observed")

# For each desired sample, calculate
for(i in 1:1000) {
  
  # Sample yhat
  yhat <- rnorm(nrow(lincom), mean=lincom[,i], sd=sqrt(results_sim[i,7]))
  
  # Calculate simulated resid
  resids[, i, 1] <- yhat - lincom[,i]
  
  # Calculate observed residual
  resids[, i, 2] <- y - lincom[,i]
  
  # Compute skewness
  skewed[i,] <- c(e1071::skewness(resids[,i,1]), 
                  e1071::skewness(resids[,i,2]))
  
}

# Compute where sim > obs
mean(skewed[,1] > skewed[,2]) 

# Plot residuals and skew

# Unroll array
data_frame(
  "sim" = resids[,,1] %>% as_tibble() %>% 
    gather(var, sim) %>% select(sim) %>% pull(),
  "obs" = resids[,,2] %>% as_tibble() %>% 
    gather(var, obs) %>% select(obs) %>% pull()
) %>%
  gather(data, value) %>%
  ggplot(., aes(x=value, fill = data)) +
    geom_density()

# Unroll skewness stats
skewed %>% as_tibble() %>%
  gather(data, value) %>%
  ggplot(., aes(x=value, fill=data)) +
    geom_density()

# Boston housing: data is a terrible fit. But this makes sense if you plot a histogram or qqnorm ==> there is non-normality in the data. 

# Simulated data: Near .5 so yay us ==> data is excellent fit Not strange given that we actually used normal samples

## Check the equal variance assumption

# Plan: cut the residuals at the median
# Randomly compare above or below median

# Compute heteroskedasticity for observed and simulated
heterosked <- matrix(0L, ncol=2, nrow=1000)
colnames(heterosked) <- c("simulated", "observed")

# Heteroskedasticity test
test_for_homoskedasticity <- function(y, X) {
  
  # Regress residuals on X
  fit <- lm("y ~ .", as.data.frame(cbind(y, X[,-1])))
  
  # Return
  return(summary(fit)$r.squared)
  
}

# Test for homoskedasticity ==> this does suppose that the data are normally distributed
for(i in 1:1000) {
  
  # Get first and last quantile for each residuals (simulated & observed)
  ys <- resids[,i,1]^2
  yo <- resids[,i,2]^2

  # Calculate R^2 for each in a linear regression where we test the explanatory power of independent variables on residuals
  heterosked[i,1] <- test_for_homoskedasticity(ys, X)
  heterosked[i,2] <- test_for_homoskedasticity(yo, X)
  
}

mean(heterosked[,1] > heterosked[,2])

# Plot
heterosked %>%
  as_tibble() %>%
  gather(grp, ratio) %>%
  ggplot(., aes(x=ratio, fill=grp)) +
    geom_density(alpha=0.4)
