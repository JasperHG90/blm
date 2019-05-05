## Generate a multilevel dataset and try out JAGS program
rm(list=ls())

# Set seed
set.seed(501)

# Sample size (employees)
n <- 400
# Groups (companies)
k <- 10
# Sample companies
companies <- sample(1:k, n, replace=TRUE)

# Specify level 2 mean and variance (company level)
# That is, the mean of level 2: (in '000 dollars)
u_expected <- 120
u_sd <- 10
# Generate the group means using these numbers
u_companies <- rnorm(k, mean = u_expected, sd = u_sd)

# Specify level 1 intercept of each employee, depending on the company they work for
intcp <- matrix(0, nrow = n, ncol = k)
# Populate with 1 if employee belongs to company k
for (i in 1:n) {
  intcp[i, companies[i]] <- 1
}
# Create true intercept for employees
b0_true <- intcp %*% u_companies
# Create true slope for age
b1_true_age <- 2.2
# Generate data for age (sorted)
age <- sort(rnorm(n, mean=40, sd=5))
# True residual variance
sd_true <- 10
# Create the outcome variable y (compensation)
compensation <- rnorm(n, mean = b0_true + b1_true_age * age, sd = 10)

# To data frame
directors <- data.frame(compensation = compensation, age=age, company=companies)

# Convert to JAGS format ----

# Jags data
jagsdata_s3 <- with(directors, list(compensation = compensation, company=company,
                                  n=n, k=k))

# Jags model
jmodl <- "model{
}"



# Initial values
init_values <- function(){
  list(tau_a = runif(1), tau = runif(1))
}

library(rjags)
# Specify model in JAGS
mod_io <- jags.model(textConnection(jmodl),
                     data = jagsdata_s3,
                     inits = init_values,
                     n.chains=2)

# Burn
update(mod_io, n.iter=20000)

# Draw samples
params <- c("E_u0","sigma", "sigma_u0")
# Run the chain
res <- coda.samples(mod_io, variable.names = params, n.iter=80000, thin = 10)

summary(res)

# Run ppc
dic.samples(mod_io, 4000, "pD")
