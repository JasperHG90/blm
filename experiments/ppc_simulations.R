## Simulation study for PPC

library(blm)
blm_setup()

# iterations
k <- 1000

# Results
sim_res <- matrix(0, ncol=3, nrow=k)

# Loop
t1 <- Sys.time()
for(i in 1:k) {

  # Simulate data
  d <-blm:::generate_dataset(n = 200, j=2, binary = 1, heteroskedastic = FALSE,
                              correlated_errors = TRUE, degree=0.4)

  # Unroll
  df <- as.data.frame(cbind(d$y, d$X[,-1]))
  # Names
  colnames(df) <- c("y", paste0("x", 1:(ncol(df) - 1)))

  # Blm object
  bfit <- blm("y ~ .", data=df) %>%
    set_sampling_options(iterations = 15000, burn = 2000, chains = 2) %>%
    sample_posterior()

  # Get ppc results
  ppcr <- bfit %>%
    evaluate_ppc(iterations = 4000)

  # Store results
  sim_res[i, ] <- c(ppcr$results$normality, ppcr$results$homosked, ppcr$results$independence)

}
t2 <- Sys.time() - t1

# Store data
saveRDS(sim_res, "experiments/ppc_autocor_data.rds")
