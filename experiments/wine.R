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

# Make the object
wine_fit <- blm("Alcohol ~ .",
                data=wine2) %>%
  # Update samplers
  set_sampler("b1", type="MH") %>%
  # Update sampling settings
  set_sampling_options(., chains = 2, iterations = 10000, burn = 2000,
                       thinning=1) %>%
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
plot(wine_ppc, "heteroskedasticity")

# Model fit
wine_fit %>%
  evaluate_model_fit()

library(tidyr)
indep_corr <- data.frame(
  "index" = 1:nrow(corrs),
  "sampled" = corrs[,1],
  "observed" = corrs[,2]
) %>%
  gather(dataset, value, -index)


library(ggplot2)
p <- ggplot(indep_corr, aes(x=index, y=value, color=dataset, shape=dataset)) +
  geom_point(alpha=0.5) +
  theme(legend.position = "top")
# Add maginals
ggMarginal(p, type = "histogram", margins="y", groupColour = FALSE, alpha=0.3,
           groupFill = TRUE)

winelm <- lm("Alcohol ~ .", data = wine2)
plot(predict(winelm), resid(winelm))

# Draw more samples
wine_fit <- wine_fit %>%
  update_posterior(10000)
