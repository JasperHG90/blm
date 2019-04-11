## Run simulations for posterior predictive checks

dat <- function() {
  # SImulate data
  d <- blm:::generate_dataset(n=500, j=3, binary=1, heteroskedastic = FALSE, degree=1)
  X <- d$X[,-1]
  y <- d$y

  # Bind
  df <- as.data.frame(cbind(y, X))
  colnames(df) <- c("attitude", "extraversion", "agreeableness", "male")

  return(df)
}

r1 <- rep(0, 500)
r2 <- r1
r3 <- r1
for(i in 1:100) {
  # Center
  datc <- as.data.frame(scale(dat(), center=TRUE, scale=FALSE))
  fit <- blm("attitude ~ .", data=datc) %>%
    # Update sampling options
    set_sampling_options(chains=2, iterations=10000) %>%
    sample_posterior() %>%
    evaluate_ppc()
  r1[i] <- fit$results$homosked
  r2[i] <- fit$results$normality
  r3[i] <- fit$results$independence
}
