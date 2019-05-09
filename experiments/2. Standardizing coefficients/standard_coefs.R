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
                    x2 = (x2 - mean(x2)) / sd(x2),
                    x3 = (x3 - mean(x3)) / sd(x3),
                    x4 = (x4 - mean(x4)) / sd(x4)))

summary(modl)

# Standardize coefs
modl2 <- lm("y ~ .", data=df)
coefs <- coef(modl2)
coefs[2] * sd(x1) / sd(y)
coefs[3] * sd(x2) / sd(y)
coefs[4] * sd(x3) / sd(y)
coefs[5] * sd(x4) / sd(y)
coef(modl)

# Compare
library(lm.beta)
coef(lm.beta(modl2))

# Generate a distribution for the betas
k <- 20000

# New coefs
b0 <- 2
b1 <- 8.1 # Close to prior
b2 <- 8 # Not close to prior

# Regenerate y
mu <- b0 + b1*x1 + b2*x2
y <- rnorm(n, mean=mu, sd=4)

# Set data
data <- data.frame(y=y, x1=x1 - mean(x1), x2=x2 - mean(x2))

library(blm)
library(dplyr)
data("sesamesim")
sesamesim <- sesamesim %>%
  mutate(funumb = funumb - mean(funumb),
         peabody =peabody - mean(peabody),
         prenumb = prenumb - mean(prenumb))
blm_setup()
k <- 20000
fit <- blm("postnumb ~ prenumb + funumb + peabody", data=sesamesim) %>%
  set_sampling_options(iterations=k, chains=2) %>%
  sample_posterior()

plot(fit, "autocorrelation")
plot(fit, "history")
plot(fit, "density")

# Scale coefficients
scale_coef <- function(x, y) {
  return(x * sd(x) / sd(y))
}

prior_b1 <- rnorm(40000, fit$priors$b1$mu, fit$priors$b1$sd) %>%
  scale_coef(fit$input$y)
prior_b2 <- rnorm(40000, fit$priors$b2$mu, fit$priors$b2$sd) %>%
  scale_coef(fit$input$y)
prior_b3 <- rnorm(40000, fit$priors$b3$mu, fit$priors$b3$sd) %>%
  scale_coef(fit$input$y)

postsamp <- fit %>%
  get_posterior_samples()
post_b1 <- scale_coef(postsamp$prenumb, fit$input$y)
post_b2 <- scale_coef(postsamp$funumb, fit$input$y)
post_b3 <- scale_coef(postsamp$peabody, fit$input$y)

dprior <- data.frame(
  "b1" = c(prior_b1),
  "b2" = prior_b2
)

dpost <- data.frame(
  "b1" = post_b1,
  "b2" = post_b2
)

library(ggplot2)
ggplot() +
  geom_density_2d(data=dprior, aes(x=b1, y=b2))

ggplot() +
  geom_density_2d(data=dpost, aes(x=b1, y=b2), color="red")

# HYP: b1 > b2 > b3
# Larger than 0
c <- mean((prior_b2 > prior_b3) & (prior_b1 > prior_b2))

# Posterior
f <- mean((post_b2 > post_b3) & (post_b1 > post_b2))

BF1u <- f / c

## HYP: pre > fu + pea
c <- mean((prior_b1) > (prior_b2 + prior_b3))
f <- mean((post_b1 > (post_b2 + post_b3)))
BF2u <- f / c

BF12 <- BF1u / BF2u

##### DIRECTORS

library(blm)
library(dplyr)
data("directors")
directors <- directors %>%
  mutate(Age = Age - mean(Age),
         Compensation = log(Compensation))
blm_setup()
k <- 20000
fit <- blm("Compensation ~ Age + Male", data=directors) %>%
  set_sampling_options(iterations=k, chains=2, thinning=5) %>%
  sample_posterior()

plot(fit, "autocorrelation")
plot(fit, "history")
plot(fit, "density")

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

## Compared to Bain
df <- data
regr <- lm(Compensation ~ Age + Male, directors)
regr$call$formula
# Set seed
set.seed(453)
summary(regr)
z<-bain(regr,"Male1 > Male0", standardize = TRUE)
z # Why is BF 4????

# Larger than 0
cases <- which(abs(prior_b1 - prior_b2) < (0.001 * sd(prior_b1)))
diff <- prior_b1 - prior_b2
diff[cases[1]]
c <-  abs(prior_b1 - prior_b2) < (0.001 * sd(prior_b1))
mean(c)



sd(post_b1)
sd(post_b2)

f <- mean(abs(post_b1 - post_b2) < .0028)

BFH1 <- f / c
BFHu

BFHu <- 1-f

# BF of H1 against Hu
f / BFHu

library(ggplot2)
ggplot(data.frame(pb1 = post_b1, pb2 = post_b2), aes(x=pb1, y=pb2)) +
  geom_density_2d()

# Generate posterior b
b1post <- rnorm(k, b1, .5)
b2post <- rnorm(k, b2, .8)

library(ggplot2)
prs <- data.frame("b1" = b1pr, "b2"=b2pr)
ggplot(prs, aes(x=b1, y=b2)) +
  geom_density_2d() # Unstandardized

# Standardize posterior and prior
dd <- data.frame(
  "b1_prior" = b1pr * sd(b1pr) / sd(y),
  "b2_prior" = b2pr * sd(b2pr) / sd(y),
  "b1_post" = b1post * sd(b1post) / sd(y),
  "b2_post" = b2post * sd(b2post) / sd(y)
)

p1 <- ggplot(dd, aes(x=b1_prior, y=b2_prior)) +
  geom_density_2d() +
  geom_abline(intercept = 0, slope=1, linetype="dashed", size=2)
p2 <- ggplot(dd, aes(x=b1_post, y=b2_post)) +
  geom_density_2d() +
  geom_abline(intercept = 0, slope=1, linetype="dashed", size=2)
library(gridExtra)
grid.arrange(p1, p2, ncol=2)

# Set hypothesis
# Not supported by the prior
# But supported by the data
h1 <- "b1 < b2"

# Get outer boundaries
c <- with(dd, mean(b1_prior > b2_prior))
f <- with(dd, mean(b1_post > b2_post))

# This is the bayes factor
BF <- f / c
BF

lm.beta(lm("y ~ .", data.frame(x1=x1, x2=x2, y=y)))
coef(lm("y ~ .", data.frame(x1=x1, x2=x2, y=y)))

library(bain)
?bain::bain()

