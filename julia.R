# Julia
julia <- JuliaCall::julia_setup()
library(JuliaCall)
#julia_install_package_if_needed("Distributions")
julia_library("Distributions")

init <- bfit$sampling_settings$initial_values
w <- init$chain1$w[,1]
sigma <- init$chain1$sigma
# Source gibbs_one_step
JuliaCall::julia_source("julia/posterior.jl")

# Sample from gibbs sampler
t1 <- Sys.time()
t <- julia_eval("gibbs_sampler")(bfit$input$X, bfit$input$y, w, 
                                 sigma, as.integer(30000), 
                                 as.list(unname(bfit$priors)),
                                 as.integer(1000))
Sys.time() - t1