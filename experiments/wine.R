rm(list=ls())

# Load wine data
wine <- rattle.data::wine

hist(wine$Alcohol)
hist(wine$Ash)
hist(wine$Alcalinity)
hist(wine$Phenols)
hist(wine$Nonflavanoids)

# Create new dataset with desired predictors
wine2 <- data.frame(
  "Alcalinity" = wine$Alcalinity,
  "Alcohol" = wine$Alcohol,
  "Ash" = wine$Ash,
  "Phenols" = wine$Phenols,
  "Nonflavanoids" = wine$Nonflavanoids,
  "Type" = wine$Type
)

# Load blm
library(blm)

# Set up blm
blm_setup()

# Center data
wine2[,2:5] <- scale(wine2[,2:5], center=TRUE, scale=FALSE)

# SImulate data
# d <- blm:::generate_dataset(n=100, j=3, binary=1, seed=908, heteroskedastic = FALSE, degree=10)
# X <- d$X[,-1]
# y <- d$y
#
# # To data
# wine2 <- data.frame(
#   "Alcohol" = y,
#   "Ash" = X[,1],
#   "Alcalinity" = X[,2],
#   "RedWine" = X[,3]
# )
#
# # Center
# wine2$Ash <- wine2$Ash - mean(wine2$Ash)
# wine2$Alcalinity <- wine2$Alcalinity - mean(wine2$Alcalinity)

# Make the object
wine_fit <- blm("Alcalinity ~ .",
                data=wine2) %>%
  # Update samplers
  set_sampler("b1", type="MH", lambda = 0.25) %>%
  # Update sampling settings
  set_sampling_options(chains = 2, iterations = 50000,
                       burn = 10000, thinning=3) %>%
  # Sample the posterior
  sample_posterior()

# Plot history
plot(wine_fit, "history")

# Update posterior by drawing more samples
#wine_fit <- wine_fit %>%
#  # Set burn-in to 25000
#  set_sampling_options(burn = 25000) %>%
#  update_posterior(iterations = 50000)

# Accepted draws
wine_fit %>%
  evaluate_accepted_draws()

# Convergence diagnostics
wine_fit %>%
  convergence_diagnostics()

# Effective sample size
wine_fit %>%
  evaluate_effective_sample_size()

# Plot
plot(wine_fit, "history")
plot(wine_fit, "density")
plot(wine_fit, "autocorrelation", chain=1)

# Summary
summary(wine_fit)

# Posterior predictive checks
wine_ppc <- wine_fit %>%
  evaluate_ppc(10000)

# Summary
summary(wine_ppc)

# Plot
plot(wine_ppc, "independence")

# Model fit
wine_fit %>%
  evaluate_model_fit()

# Linear model
winelm <- lm("Alcohol ~ .", data = wine2)
plot(predict(winelm), resid(winelm))
summary(winelm)
