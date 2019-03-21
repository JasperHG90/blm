# S3 generics

# Set priors for a blm object
#
# @param blm blm object
# @param ... optional arguments needed to specify a prior
# Update priors (generic)
set_priors <- function(x, ...) {
  UseMethod("set_priors")
}

# Sample the posterior distribution
#
# @param blm blm object
# @param chains number of chains to run
# @param iterations number of iterations per chain
# @param burn number of burn-in iterations
# @param julia use julia for computations?
# @return updated blm object
sampling_options <- function(x, ...) {
  UseMethod("sampling_options", x)
}

# Execute a blm plan
#
sample_posterior <- function(x) {
  UseMethod("sample_posterior", x)
}

# Set up posterior predictive checks
posterior_predictive_checks <- function(x, ...) {
  UseMethod("posterior_predictive_checks")
}

# Model fit
model_fit <- function(x) {
  UseMethod("model_fit")
}

# ---------- Functions for priors ---------------

# Draw value from a prior
draw_value <- function(x) {
  UseMethod("draw_value")
}

# ---------- Functions for posterior predictive checks ----------

normality_check <- function(x) {
  UseMethod("normality_check")
}

# Generic for homoskedasticity
homoskedast_check <- function(x) {
  UseMethod("homoskedast_check")
}

# Generic for model rmse
model_rmse <- function(x) {
  UseMethod("model_rmse")
}
