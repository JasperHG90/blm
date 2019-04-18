# Get posterior samples
w <- do.call(rbind.data.frame, fit_blm$posterior$samples)

# Split w and sigma
sigma <- w[,4]
w <- w[,-4]

# Create linear combinations
lincom <- bfit$input$X %*% t(w)

# Function to calculate R2 for each
bayes_R2 <- function(lincom_current, sigma) {

  # Variance in lincom current
  vlc <- var(lincom_current)

  # R2
  return(vlc / (vlc + sigma^2))

}

# Calculate for each linear combination
r2v <- rep(0, ncol(lincom))

for(i in 1:length(r2v)) {
  r2v[i] <- bayes_R2(lincom[,i], sigma[i])
}

library(blm)
library(dplyr)
data("mtcars")
df <- mtcars %>%
  select(mpg, wt, qsec, am) %>%
  mutate(wt = wt - mean(wt),
         qsec = qsec - mean(qsec),
         am = am - mean(am))
row.names(df) <- 1:nrow(df)
blm_setup()
fit_blm <- blm("mpg ~ wt + qsec + am", data = mtcars) %>%
  set_sampling_options(chains=1, iter=10000, burn=2000, thinning=20) %>%
  sample_posterior()

plot(fit_blm, "history")

bayes_R2 <- function(fit) {
  y_pred <- rstanarm::posterior_linpred(fit)
  var_fit <- apply(y_pred, 1, var)
  var_res <- as.matrix(fit, pars = c("sigma"))^2
  var_fit / (var_fit + var_res)
}

library(rstanarm)
(fit <- stan_lm(attitude ~ ., data=df, prior=NULL,
                # the next line is only to make the example go fast enough
                chains = 1, iter = 4000, seed = 12345))

rsq_bayes <- bayes_R2(fit)

hist(rsq_bayes)
