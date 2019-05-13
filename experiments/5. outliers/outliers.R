## Outlier plots etc

rm(list=ls())
library(blm)
blm_setup()
# ggplot
library(ggplot2)

# Simulated data (no violation) ----

# Simulate some data (no violation)
d <- blmsim(n=800, j=3, seed = 2001)
df <- cbind(d$X[,-1], d$y) %>% as.data.frame()
# Names
names(df) <- c("male", "attitude", "extraversion", "likeability")

# Prep data
df <- df %>%
  mutate(attitude = attitude - mean(attitude),
         extraversion = extraversion - mean(extraversion))

# Iterations
k <- 30000
# Fit the model
fit <- blm("likeability ~ male + attitude + extraversion", data=df) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=4) %>%
  # Compute intercept-only model
  compute_null_model() %>%
  # Sample posterior
  sample_posterior()

# Evaluate ppc
fit <- fit %>%
  evaluate_ppc()

# View
fit %>%
  get_value("ppc")

# Retrieveo ppd
out <- fit %>% get_value("ppd")
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
  geom_histogram(fill = "lightgreen", alpha=.2, color="black") +
  theme_blm() +
  scale_x_continuous(name = latex2exp::TeX("$p(y_{sim} > y_{obs})$")) +
  geom_vline(xintercept = quants, linetype = "dashed") +
  labs(title=latex2exp::TeX("Proportion where $y_{sim} > y_{obs}$"))

# Plot observed against proportion
yout <- ggplot(out_df, aes(x=(y - mean(y)) / sd(y), y=p)) +
  geom_point() +
  scale_y_continuous(limits=c(-.001,1.001)) +
  theme_blm() +
  labs(title=latex2exp::TeX("Observed outcome versus $p=(y_{sim} > y_{obs})$"))

yout

# Violate equal variance ----

# Simulate some data
d <- blmsim(n=800, j=3, seed = 2250, heteroskedastic = TRUE, degree=3)
df <- cbind(d$X[,-1], d$y) %>% as.data.frame()
# Names
names(df) <- c("male", "attitude", "extraversion", "likeability")

# Prep data
df <- df %>%
  mutate(attitude = attitude - mean(attitude),
         extraversion = extraversion - mean(extraversion))

# Iterations
k <- 30000
# Fit the model
fit <- blm("likeability ~ male + attitude + extraversion", data=df) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=4) %>%
  # Compute intercept-only model
  compute_null_model() %>%
  # Sample posterior
  sample_posterior()

# Retrieveo ppd
out <- fit %>% get_value("ppd")
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
  geom_histogram(fill = "lightgreen", alpha=.2, color="black") +
  theme_blm() +
  scale_x_continuous(name = latex2exp::TeX("$p(y_{sim} > y_{obs})$")) +
  geom_vline(xintercept = quants, linetype = "dashed") +
  labs(title=latex2exp::TeX("Proportion where $y_{sim} > y_{obs}$"))

# Plot observed against proportion
yout <- ggplot(out_df, aes(x=(y - mean(y)) / sd(y), y=p)) +
  geom_point() +
  scale_y_continuous(limits=c(-.001,1.001)) +
  theme_blm() +
  labs(title=latex2exp::TeX("Observed outcome versus $p=(y_{sim} > y_{obs})$"))

yout

# Correlated errors -----

# Simulate some data
d <- blmsim(n=800, j=3, seed = 3000, correlated_errors = TRUE, degree=.9)
df <- cbind(d$X[,-1], d$y) %>% as.data.frame()
# Names
names(df) <- c("male", "attitude", "extraversion", "likeability")

# Prep data
df <- df %>%
  mutate(attitude = attitude - mean(attitude),
         extraversion = extraversion - mean(extraversion))

# Iterations
k <- 30000
# Fit the model
fit <- blm("likeability ~ male + attitude + extraversion", data=df) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=4) %>%
  # Compute intercept-only model
  compute_null_model() %>%
  # Sample posterior
  sample_posterior()

# Retrieveo ppd
out <- fit %>% get_value("ppd")
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
  geom_histogram(fill = "lightgreen", alpha=.2, color="black") +
  theme_blm() +
  scale_x_continuous(name = latex2exp::TeX("$p(y_{sim} > y_{obs})$")) +
  geom_vline(xintercept = quants, linetype = "dashed") +
  labs(title=latex2exp::TeX("Proportion where $y_{sim} > y_{obs}$"))

# Plot observed against proportion
yout <- ggplot(out_df, aes(x=(y - mean(y)) / sd(y), y=p)) +
  geom_point() +
  scale_y_continuous(limits=c(-.001,1.001)) +
  theme_blm() +
  labs(title=latex2exp::TeX("Observed outcome versus $p=(y_{sim} > y_{obs})$"))

yout

# Directors data ----

# Add hypothesis to blm object
library(dplyr)
# Load data
data("directors")
# Adjust
directors <- directors %>%
  mutate(Age = Age - mean(Age),
         Compensation = log(Compensation),
         Male = as.numeric(Male),
         Male = Male - mean(Male))
# Iterations
k <- 30000
# Fit the model
fit <- blm("Compensation ~ Male + Age + Sector", data=directors) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=4) %>%
  # Compute intercept-only model
  compute_null_model() %>%
  # Set reasonable priors
  set_prior("b1", mu=.01, sd=.2) %>% # 1 % increase / every year ==> 20% spread
  set_prior("b2", mu=.05, sd=.03) %>%
  # Set hypotheses
  # H1: Males earn about the same as females
  set_hypothesis("H1", "b0 + b2 < b0") %>%
  # H2: Directors only earn more as they get older
  set_hypothesis("H2", "b1 > 0") %>%
  #H3: Sectors are ordered as: mu_services > mu_basic materials > mu_financial
  set_hypothesis("H3", "b0 + b4 > b0 & b0 > b0 + b3") %>%
  #set_hypothesis("|b0| > 0") %>%
  # Sample posterior
  sample_posterior()

# Retrieveo ppd
out <- fit %>% get_value("ppd")
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
  geom_histogram(fill = "lightgreen", alpha=.2, color="black") +
  theme_blm() +
  scale_x_continuous(name = latex2exp::TeX("$p(y_{sim} > y_{obs})$")) +
  geom_vline(xintercept = quants, linetype = "dashed") +
  labs(title=latex2exp::TeX("Proportion where $y_{sim} > y_{obs}$"))

# Plot observed against proportion
yout <- ggplot(out_df, aes(x=(y - mean(y)) / sd(y), y=p)) +
  geom_point() +
  scale_y_continuous(limits=c(-.001,1.001)) +
  theme_blm() +
  labs(title=latex2exp::TeX("Observed outcome versus $p=(y_{sim} > y_{obs})$"))

yout

# Get values of outliers for a threshold
# x is outliers object
# p is threshold for outliers (symmetric ==> p and 1-p)
which_outliers <- function(x, p=.025) {

  # Get indices
  ind <- which(x$results[,1] <= p | x$results[,1] >= 1-p)

  # Return
  return(ind)

}

# Cluster indices to search for pattern
ol <- which_outliers(out, p=.025)
grp <- rep(0, nrow(df))
grp[ol] <- 1
# Plot
ggplot(data.frame(x=df[,1], y=x$results[,1],
                  c = as.factor(grp)), aes(x=x,fill=c)) +
  geom_density(alpha=.4) +
  geom_vline(xintercept = mean(df[ol,1])) +
  geom_vline(xintercept = mean(df[-ol, 1]), linetype="dashed")
