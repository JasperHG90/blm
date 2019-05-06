## Library
library(blm)
library(dplyr)

# Set up blm
blm_setup()

# Directors data
#data("directors")

# Generate some data
d <- blm:::generate_dataset(n = 200, j=4, binary = 1)
data <- as.data.frame(cbind(d$X[,-1], d$y))

# Modify
directors <- directors %>%
  # Log compensation
  mutate(Compensation = log(Compensation)) %>%
  # Subtract mean from Age variable
  mutate(Age = Age - mean(Age),
         Male = as.numeric(Male) - mean(as.numeric(Male)))

# Build the model
fit1 <- blm("V5 ~ .", data=data) %>%
  # Set sampling options
  set_sampling_options(chains=2, iterations=40000, thinning=3, burn=2000) %>%
  # Compute the null model
  compute_null_model(iterations=25000) %>%
  # Sample the data
  sample_posterior()

# Assess convergence
plot(fit1, "nullmodel")

summary(fit1)

