rm(list=ls())
# Subset
library(dplyr)
d <- read.csv2("docs/data/FinalData.csv") %>%
  # Select only independent directors
  filter(isED == 0) %>%
  # Select only these variables
  select(Recent_Comp_Thous, Sector, isM, Age) %>%
  # Remove NA
  filter(!is.na(Recent_Comp_Thous),
         !is.na(Age)) %>%
  # Log compensation
  mutate(compensation = log(Recent_Comp_Thous)) %>%
  # Remove these variables
  select(-Recent_Comp_Thous) %>%
  # Center data
  mutate(Age = Age - mean(Age),
         isM = isM - mean(isM)) %>%
  # Filter small sectors
  filter(!Sector %in% c("Money Center Banks", "Property & Casualty Insurance")) %>%
  # Filter sectors and add re-factor
  filter(Sector %in% c("Financial", "Services", "Basic Materials")) %>%
  mutate(Sector = factor(as.character(Sector),
                         levels = c("Financial", "Services", "Basic Materials")))

library(ggplot2)
ggplot(d, aes(x=Age, y=compensation, group=Sector, color=as.factor(isM))) +
 geom_point() +
 geom_smooth(method="lm") +
 facet_wrap(~ Sector)

ggplot(d %>% mutate(isM = as.factor(isM)), aes(x=compensation, group=isM, fill=isM)) +
 geom_density(alpha=0.3) +
 facet_wrap(~ Sector)

# Load blm
library(blm)

# Set up blm
blm_setup()

# Model 1: Age & Gender
dirmod1 <- blm("compensation ~ Age + isM",
                    data=d) %>%
  set_sampling_options(chains = 2, iterations = 20000,
                       burn = 2000, thinning=3) %>%
  # Sample
  sample_posterior()

# Convergence diagnostics
dirmod1 %>%
  convergence_diagnostics()

# Plots
plot(dirmod1, "history")
plot(dirmod1, "autocorrelation")

# Summary
summary(dirmod1)

# R-squared
m1r2 <- dirmod1 %>%
  evaluate_R2()

# Plot R2
plot(m1r2)

# PPC
m1_ppc <- dirmod1 %>%
  evaluate_ppc()

m1_ppc
# Violation assumption of normality and independence

# Model fit
dirmod1 %>%
  evaluate_model_fit()

## Model II -- Adding sectors & priors

# Make the object
dirmod2 <- blm("compensation ~ .",
                data=d) %>%
  # Update samplers
  #set_sampler("b5", type="MH", lambda = 0.005) %>%
  #set_sampler("b6", type="MH", lambda = 0.05) %>%
  # Add prior
  #set_prior("b1", mu=2, sd=4) %>%
  # Update sampling settings
  set_sampling_options(chains = 2, iterations = 40000,
                       burn = 4000, thinning=4) %>%
  #set_prior("b6", mu = 1.5, sd = 0.5) %>%
  #set_prior("b4", mu = 3, sd=1) %>%
  # Set samplers
  set_sampler("b4", type="MH", lambda=0.01) %>%
  # Set initial values
  set_initial_values(chain_1 = list("b" = c(4, 2, 3, 8, 10), "sigma"= 1),
                     chain_2 = list("b" = c(2, 4, 7, 2, 3), "sigma"=2)) %>%
  # Sample
  sample_posterior()

# Convergence
dirmod2 %>%
  set_sampling_options(burn = 15000) %>%
  convergence_diagnostics()

# Look at accepted samples
dirmod2 %>%
  evaluate_accepted_draws()

# Look at effective sample size
dirmod2 %>%
  evaluate_effective_sample_size()

# Plots
plot(dirmod2, "history")
plot(dirmod2, "autocorrelation")

# Summary
summary(dirmod2)

# R-squared
m2r2 <- dirmod2 %>%
  evaluate_R2()

# Plot R2
plot(m2r2)

# PPC
m2ppc <- dirmod2 %>%
  evaluate_ppc()

m2ppc
# Normality is now OK but heteroskedasticity and independence are violated
plot(predict(dirmod2), resid(dirmod2))
plot(m2ppc, "heteroskedasticity")
plot(m2ppc, "independence")

# Model fit
dirmod2 %>%
  evaluate_model_fit()

# Posterior predictive distribution
ind <- sample.int(nrow(dirmod2$input$X), 1)
x1 <- matrix(dirmod2$input$X[ind, ], nrow=1)
pp <- x1 %*% t(get_posterior_samples(dirmod2)[,-6])
x1
d[ind,]
hist(pp)
exp(d$compensation[ind])

# Log transformation makes is so that we need to interpret the model using percentages
# https://data.library.virginia.edu/interpreting-log-transformations-in-a-linear-model/
(exp(coef(dirmod2)[-1]) -1) * 100
summary(dirmod2)
