## Testing

source("blm.R")

# Load data & Run
data("iris")
iris <- iris[,-5]
library(magrittr)
res <- blm("Petal.Width ~ .", iris)  %>% 
  set_priors("b0" = prior("normal", mu=4, sd=8),
             "b1" = prior("normal", mu=12, sd=12)) %>%
  sample_posterior(chains = 2, iterations = 20000, burn=2000, julia=TRUE) %>%
  execute()

# (ADD THINNING)

library(doFuture)


# Load Julia
library(JuliaCall)
#julia <- julia_setup()

# Source gibbs_one_step
JuliaCall::julia_source("gos.jl")
res <- julia_eval("gibbs_sampling")(8, 1000000, 1)
