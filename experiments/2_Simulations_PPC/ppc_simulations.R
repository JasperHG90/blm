## Simulation study for PPC
## Plan: Run the plans below ==> without any violations and violating to different degrees

library(blm)
blm_setup()

# iterations
k <- 1000

# Results
sim_res <- array(0L, c(k, 3, 3))

# Grid of values to test
# See: generate_dataset in helpers.R (line 222)
grid <- list(
  "noviolation" = list(
    "results" = array(0L, c(k, 3, 1)),
    "degrees" = 1 # Degrees will be ignored if no violation but it eases coding in the for-loop below
  ),
  "heterosked" = list(
    "degrees" = c("mild" = 1, "medium"=3, "severe" = 5),
    "results" = sim_res
  ),
  "indep" = list(
    "degrees" = c("mild"=0.1, "medium"=0.5, "severe"=0.9), # Determines 'ar' param for an arima.sim
    "results" = sim_res
  )
)

# For each test
# Progress bar
library(utils)
pb <- txtProgressBar(min = 0, max = (k * (3 * 2)) + k, style = 3) # (k * (length(grid) * length(degrees))) + k [==>] (last k is for no violation)
# keep track of total i for progress bar
i_pb <- 0
# For each assumption
for(j in seq_along(grid)) {

  # For each degree
  for(d in seq_along(grid[[j]][["degrees"]])) {

    # For each simulation
    for(i in 1:k) {

      # Increment i_pb
      i_pb <- i_pb + 1

      # Set pb
      setTxtProgressBar(pb, i_pb)

      # Unroll data
      noviolation <- ifelse(j == 1, TRUE, FALSE)
      heteroskedastic <- ifelse(j == 2, TRUE, FALSE)
      independence <- !heteroskedastic & !noviolation

      # Simulate data
      dat <-blm::blmsim(n = 100, j=2, binary = 0, heteroskedastic = heteroskedastic,
                        correlated_errors = independence, degree=grid[[j]][["degrees"]][[d]])

      # Unroll
      df <- as.data.frame(cbind(dat$y, dat$X[,-1]))
      # Names
      colnames(df) <- c("y", paste0("x", 1:(ncol(df) - 1)))

      # Blm object
      bfit <- blm("y ~ .", data=df) %>%
        set_sampling_options(iterations = 15000, burn = 2000, chains = 1) %>%
        sample_posterior() %>%
        # PPC
        evaluate_ppc(p=1)

      # Store ppc
      ppcr <- get_value(bfit, "ppc")
      # Store results
      grid[[j]][["results"]][i, , d] <- c(ppcr$results$normality, ppcr$results$homosked, ppcr$results$independence)

      # Save every 10 iterations
      if(i %% 10 == 0) {
        # Store data
        saveRDS(grid, "experiments/2_Simulations_PPC/ppc_simulated/tmp.rds")
      }

    }

  }

}

# Store data
saveRDS(grid, "experiments/2_Simulations_PPC/ppc_simulated/final.rds")
