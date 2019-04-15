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
  "Alcohol" = wine$Alcohol,
  "Ash" = wine$Ash,
  "Alcalinity" = wine$Alcalinity,
  "Phenols" = wine$Phenols,
  "Nonflavanoids" = wine$Nonflavanoids,
  "Type" = wine$Type
)

# Set up blm
blm_setup()

# Center data
wine2[,2:5] <- scale(wine2[,2:5], center=TRUE, scale=FALSE)

# SImulate data
d <- blm:::generate_dataset(n=100, j=3, binary=1, seed=908, heteroskedastic = FALSE, degree=10)
X <- d$X[,-1]
y <- d$y

# To data
wine2 <- data.frame(
  "Alcohol" = y,
  "Ash" = X[,1],
  "Alcalinity" = X[,2],
  "RedWine" = X[,3]
)

# Center
wine2$Ash <- wine2$Ash - mean(wine2$Ash)
wine2$Alcalinity <- wine2$Alcalinity - mean(wine2$Alcalinity)

# Make the object
wine_fit <- blm("Alcohol ~ .",
                data=wine2) %>%
  # Update samplers
  set_sampler("b1", type="MH", lambda = 1) %>%
  # Update sampling settings
  set_sampling_options(chains = 2, iterations = 10000,
                       burn = 2000, thinning=1) %>%
  # Sample the posterior
  sample_posterior()

# Convergence diagnostics
wine_fit %>%
  convergence_diagnostics()

# Update posterior
wine_fit <- wine_fit %>%
  update_posterior(20000)

# Effective sample size
wine_fit %>%
  evaluate_effective_sample_size()

# Accepted draws
wine_fit %>%
  evaluate_accepted_draws()

# Plot
plot(wine_fit, "history")
plot(wine_fit, "density")
plot(wine_fit, "autocorrelation", chain=1)

# Summary
summary(wine_fit)

# Posterior predictive checks
wine_ppc <- wine_fit %>%
  evaluate_ppc(3000)

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
