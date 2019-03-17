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

# Set up a parallel backend
#registerDoParallel(cores=2)

# Setting up a dataset
d <- generate_dataset(center = FALSE, seed = 687, n=2000, j=3)
X <- d$X[, -1]
y <- d$y
df <- as.data.frame(cbind(y, X))
df$V3 <- as.factor(df$V3)

# Linear model for comparison
ffit <- lm("y ~ .", data=df)
summary(ffit)

# Set up parallel backend for doFuture package
registerDoFuture()

# Set up a BLM object
bfit <- blm("y ~ .", data=df, center = FALSE) %>%
  # Update sampling settings
  sampling_options(., chains = 2, iterations = 10000, burn = 1000) %>%
  # Set prior for coefficient b2
  set_priors("b2" = prior("normal", mu=3, sd=4))

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

# Print summary
summary(bfit)
