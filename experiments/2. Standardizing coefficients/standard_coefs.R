## Part I: making standardized coefs

# Generate some data
n <- 500
x1 <- rnorm(n, 0, 3)
x2 <- rnorm(n, 0, 5)

# True coefficients
b0 <- 2.3
b1 <- 5
b2 <- 2

# Generate y
mu <- b0 + b1*x1 + b2*x2
y <- rnorm(n, mean=mu, sd=2)

# To df
df <- data.frame(
  "y" = y,
  "x1" = x1,
  "x2" = x2)

library(magrittr)
library(dplyr)

# Try using standardizing data
modl <- lm("y ~ .", data=df %>%
             mutate(x1 = (x1 - mean(x1)) / sd(x1),
                    x2 = (x2 - mean(x2)) / sd(x2)))

summary(modl)

# Standardize coefs
modl2 <- lm("y ~ .", data=df)
coefs <- coef(modl2)
coefs[2] * sd(x1) / sd(y)
coefs[3] * sd(x2) / sd(y)
coef(modl)

# Compare
library(lm.beta)
coef(lm.beta(modl2))

## Try fit and complexity on directors data -----

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
k <- 20000
# Fit the model
fit <- blm("Compensation ~ Age + Male", data=directors) %>%
  set_sampling_options(iterations=k, chains=2, thinning=5) %>%
  # Set reasonable priors
  set_prior("b1", mu=0, sd=.01) %>%
  set_prior("b2", mu=.05, sd=.03) %>%
  # Sample posterior distribution
  sample_posterior()

# Convergence
plot(fit, "autocorrelation")
plot(fit, "history")
plot(fit, "density")
fit %>% evaluate_convergence_diagnostics()

# Scale coefficients
scale_coef <- function(x, y) {
  return(x * sd(x) / sd(y))
}

# Priors
prior_b0 <- rnorm(40000, fit$priors$b0$mu, fit$priors$b0$sd) %>%
  scale_coef(fit$input$y)
prior_b1 <- rnorm(40000, fit$priors$b1$mu, fit$priors$b1$sd) %>%
  scale_coef(fit$input$y)
prior_b2 <- rnorm(40000, fit$priors$b1$mu, fit$priors$b1$sd) %>%
  scale_coef(fit$input$y)

# Posteriors
postsamp <- fit %>%
  get_posterior_samples()
post_b0 <- scale_coef(postsamp$`(Intercept)`, fit$input$y)
post_b1 <- scale_coef(postsamp$Age, fit$input$y)
post_b2 <- scale_coef(postsamp$Male1, fit$input$y)

## HYP: Male1 > Male0 ==> b2+b0 > b0 (Male + intercept larger than simply intercept, which equals female + intercept)

c <- mean((prior_b0 + prior_b2) > prior_b0)
f <- mean((post_b0 + post_b2) > post_b0)

BF1u <- f / c
BF1u

## HYP: Male1 > Male0 > Age
c <- mean(((prior_b0 + prior_b2) > prior_b0) & (prior_b0 > prior_b1))
f <- mean(((post_b0 + post_b2) > post_b0) & (post_b0 > post_b1))

## HYP: Age > 0
c <- mean((prior_b1 > 0))
f <- mean((post_b1 > 0))

library(bain)
## Compared to Bain
regr <- lm(Compensation ~ Age + Male, directors)
regr$call$formula
# Set seed
set.seed(453)
summary(regr)
z<-bain(regr,"Male1 > Male0 > Age; Age >0", standardize = TRUE)
z # Why is BF 4????

# Plotting priors and posteriors
p1 <- ggplot(dd, aes(x=b1_prior, y=b2_prior)) +
  geom_density_2d() +
  geom_abline(intercept = 0, slope=1, linetype="dashed", size=2)
p2 <- ggplot(dd, aes(x=b1_post, y=b2_post)) +
  geom_density_2d() +
  geom_abline(intercept = 0, slope=1, linetype="dashed", size=2)
library(gridExtra)
grid.arrange(p1, p2, ncol=2)


