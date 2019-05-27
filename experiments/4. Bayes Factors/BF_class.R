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
  na.omit() %>%
  mutate(Age = Age - mean(Age),
         Compensation = log(Compensation))
         #Male = as.numeric(Male),
         #Male = Male - mean(Male))
# Bridge to julia
blm_setup()
# Iterations
k <- 50000
# Fit the model
fit <- blm("Compensation ~ Age + Male", data=directors) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=5) %>%
  # Set reasonable priors
  #set_prior("b1", mu=.01, sd=.2) %>% # 1 % increase / every year ==> 20% spread
  set_prior("b2", mu=0, sd=1) %>%
  # Set hypotheses
  # H1: Males earn about the same as females
  set_hypothesis("H1", "b2 > 0 & b2 < .01") %>%
  # H2: Directors only earn more as they get older
  set_hypothesis("H2", "b1 > 0") %>%
  #H3: Sectors are ordered as: mu_services > mu_basic materials > mu_financial
  #set_hypothesis("b0 + b4 > b0 & b0 > b0 + b3") %>%
  #set_hypothesis("|b0| > 0") %>%
  # Sample posterior
  sample_posterior()

# Print blm
print(fit)

# Add another hypothesis
fit <- fit %>%
  evaluate_hypotheses()

# Print blm
print(fit)

# Summary
summary(fit)

fit %>%
  get_value("hypotheses") %>%
  summary()

library(bain)
## Compared to Bain
regr <- lm(Compensation ~ Age + Male, directors)
regr$call$formula
# Set seed
set.seed(453)
summary(regr)
z<-bain(regr,"Male > 0 & Male < .1; Age > 0", standardize = FALSE)

z # Why is BF 4 ==> it is: (fit_hypothesis / complexity_hypothesis) / (fit_complement_of_hypothesis / complexity_complement_of_hypothesis)

# In case of |u1 - u2| < .2sd ==> which sd!! are we referring to? p.35
# Exact equality constraints ==> how do we do this? (p.36)
# How to calculate Pr_b ???

# If complexity is 0 ==> how to deal with this numerically if f / c == BF

h <- fit$hypotheses
