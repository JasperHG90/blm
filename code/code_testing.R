## Code testing

# Clean wd
rm(list=ls())

# Load functions
source("code/blm.R")
source("code/generics.R")
source("code/gibbs.R")
source("code/helpers.R")
source("code/methods.R")
source("code/posterior_densities.R")
source("code/zzz.R")
source("code/eval.R")
source("code/ppc.R")

# Libraries
library(doFuture)
#library(doParallel)
library(magrittr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(tibble)

# Setting up a dataset
d <- generate_dataset(center = FALSE, seed = 687, n=1000, j=3)
X <- d$X[, -1]
y <- d$y
df <- as.data.frame(cbind(y, X))
df$V3 <- as.factor(df$V3)

#library(readr)
d <- read.csv("Appendices/boston-housing/train.csv")
y <- d$medv
X <- cbind(d$rm, d$crim, d$chas, d$lstat, d$dis)
df <- as.data.frame(cbind(y, X))
df$V4 <- as.factor(df$V4)
#df$V2 <- df$V2 - mean(df$V2)
#df$V3 <- df$V3 - mean(df$V3)
##df$V5 <- df$V5 - mean(df$V5)
#df$V6 <- df$V6 - mean(df$V6)

# Linear model for comparison
ffit <- lm("y ~ .", data=df)
summary(ffit)

# Set up parallel backend for doFuture package
registerDoFuture()

## Sequential (4 chains, 20.000 iterations) ==> +- 81 seconds
## Multicore (4 chains, 20.000 iterations) ==> +- 36 seconds

# Set up a BLM object
bfit <- blm("y ~ .", data=df, center = TRUE) %>%
  # Update sampling settings
  sampling_options(., chains = 2, iterations = 10000, burn = 2000,
                   initial_weight_correction = TRUE) %>%
  # Set prior for coefficient b2
  set_priors("b2" = prior("normal", mu=3, sd=4)) %>%
  # Sample posterior
  sample_posterior()

# Call summary
summary(bfit)

# Print object information
print(bfit)

# Plot
plot(bfit, "history")
plot(bfit, "density")
plot(bfit, "autocorrelation", chain=1)

# Posterior predictive checks
post_pc <- bfit %>%
  posterior_predictive_checks(iterations=5000, burn=4000)

# Call summary
summary(post_pc)

# Plot
plot(post_pc, "normality")
plot(post_pc, "heteroskedasticity")
plot(post_pc, "rmse")

# TODO: Add DIC (model fit)
# TODO: sigma in output
# TODO: burn-in period for PPC ==> always 1.000?
# TODO: for DIC ==> sample posterior for a subset
# TODO: add diagnostic check for proper burn-in period ==> that is, run a linear regression that checks whether there is a pattern in the data over time.