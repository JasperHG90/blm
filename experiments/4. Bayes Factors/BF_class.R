## Hypothesis framework

## Allowed:
# 1. one, specific hypothesis
#    a > b
#    a < b
#    a = b
#    (a-b) < 1
#    2*a > 0
#    ...
# 2. Multiple hypotheses chained with '&' (AND)
#    a = 0 & b = 0
#    a < b & b < c
#    a < .1 & a > .1 ==> same as |a| > .1
#    (a-b) < .1 & (a-b) > .1 ==> same as |a-b| > .1
#    ...

library(stringr)
library(magrittr)

# Single hypothesis -----

rm(list=ls())
# Add hypothesis to blm object
library(blm)
library(dplyr)
# Load data
data("directors")
# Adjust
directors <- directors %>%
  mutate(Age = Age - mean(Age),
         Compensation = log(Compensation))
# Bridge to julia
blm_setup()
# Iterations
k <- 50000
# Fit the model
fit <- blm("Compensation ~ Age + Male", data=directors) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=5) %>%
  # Set reasonable priors
  set_prior("b1", mu=0, sd=.01) %>%
  set_prior("b2", mu=.05, sd=.03) %>%
  # Set hypotheses
  # H1: Males earn about the same as females
  set_hypothesis("b0 + b2 < b0") %>%
  # H2: Directors only earn more as they get older
  set_hypothesis("b1 > 0") %>%
  # H3: Sectors are ordered as: mu_services > mu_basic materials > mu_financial
  #set_hypothesis("b0 + b4 > b0 & b0 > b0 + b3") %>%
  set_hypothesis("|b0| > 0") %>%
  # Sample posterior
  sample_posterior()

# Evaluate hypotheses
fit <- fit %>%
  evaluate_hypotheses()

fit$hypotheses$hypothesis_1$result # We expect this to be true (mean males < mean females)
fit$hypotheses$hypothesis_2$result # Age > 0 is also in line of expectation
fit$hypotheses$hypothesis_3$result # Mu services > mu basic materials > mu financial (just for testing)
fit$hypotheses$hypothesis_4$result # This is also in line because of course abs(b0) > 0

library(bain)
## Compared to Bain
regr <- lm(Compensation ~ Age + Male, directors)
regr$call$formula
# Set seed
set.seed(453)
summary(regr)
z<-bain(regr,"Male1 < Male0 ; Age > 0", standardize = TRUE)
z # Why is BF 4 ==> it is: (fit_hypothesis / complexity_hypothesis) / (fit_complement_of_hypothesis / complexity_complement_of_hypothesis)

# In case of |u1 - u2| < .2sd ==> which sd are we referring to? p.35
