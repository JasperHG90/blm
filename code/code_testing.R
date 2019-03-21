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

# Libraries
library(doFuture)
#library(doParallel)
library(magrittr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)

# Setting up a dataset
d <- generate_dataset(center = FALSE, seed = 687, n=2000, j=3)
X <- d$X[, -1]
y <- d$y
df <- as.data.frame(cbind(y, X))
df$V3 <- as.factor(df$V3)

library(readr)
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

# Set up a BLM object
bfit <- blm("y ~ .", data=df, center = TRUE) %>%
  # Update sampling settings
  sampling_options(., chains = 2, iterations = 20000, burn = 12000,
                   initial_weight_correction = TRUE) #%>%
  # Set prior for coefficient b2
  #set_priors("b2" = prior("normal", mu=3, sd=4))

## Sequential (4 chains, 20.000 iterations) ==> +- 81 seconds
## Multicore (4 chains, 20.000 iterations) ==> +- 36 seconds

# Call summary
summary(bfit)

# Print object information
print(bfit)

# Run model
bfit <- bfit %>%
  # Sample posterior
  sample_posterior()

# Plot
plot(bfit, "history")
plot(bfit, "density")
plot(bfit, "autocorrelation", chain=1)

# Print summary
summary(bfit)

# Posterior predictive checks
post_pd <- bfit %>%
  posterior_predictive_checks()

# TODO: Add DIC (model fit)

