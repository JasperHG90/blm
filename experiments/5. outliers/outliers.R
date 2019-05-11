## Outlier plots etc

library(ggplot2)
library(blm)

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
  # Compute intercept-only model
  compute_null_model() %>%
  # Set reasonable priors
  #set_prior("b1", mu=.01, sd=.2) %>% # 1 % increase / every year ==> 20% spread
  #set_prior("b2", mu=.05, sd=.03) %>%
  # Set hypotheses
  # H1: Males earn about the same as females
  set_hypothesis("b0 + b2 < b0") %>%
  # H2: Directors only earn more as they get older
  set_hypothesis("b1 > 0") %>%
  #H3: Sectors are ordered as: mu_services > mu_basic materials > mu_financial
  #set_hypothesis("b0 + b4 > b0 & b0 > b0 + b3") %>%
  #set_hypothesis("|b0| > 0") %>%
  # Sample posterior
  sample_posterior()

# Retrieveo utliers
out <- fit %>% get_value("outliers")
# Retrieve actual
y <- fit$input$y

# To df
out_df <- data.frame(
  "p" = out$results[,1],
  "y" = y
)

# 1. Density plot with quantiles
quants <- quantile(out_df$p, c(0.025, 0.25, 0.5, 0.75, 0.975))

# Plot
ggplot(out_df, aes(x=p)) +
  geom_density(fill = "lightgreen", alpha=.2) +
  theme_blm() +
  scale_x_continuous(name = latex2exp::TeX("$p(y_{sim} > y_{obs})$")) +
  geom_vline(xintercept = quants, linetype = "dashed") +
  labs(title=latex2exp::TeX("Proportion where $y_{sim} > y_{obs}$"))

# Plot observed against proportion
ggplot(out_df, aes(x=y, y=p)) +
  geom_point() +
  scale_y_continuous(limits=c(-.001,1.001)) +
  theme_blm() +
  geom_hline(yintercept = quants, linetype = "dashed") +
  labs(title="Observed data against proportion outliers")

# Get values of outliers for a threshold
# x is outliers object
# p is threshold for outliers (symmetric ==> p and 1-p)
which_outliers <- function(x, p=.75) {

  # Get indices
  ind <- which(x$results[,1] >= p | x$results[,1] <= (1-p))

  # Return
  return(ind)

}

# Cluster indices to search for pattern
ol <- which_outliers(out, p=.8)
# Euclidean distance
d <- dist(ol)
# Hierarchical clusters
cl <- hclust(d)
plot(cl)
