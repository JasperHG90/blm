# S3 generics

## Exported ----

#' Set sampling options on a blm object
#'
#' @param x blm object
#' @param ... other options (chains, iterations etc.)
#'
#' @return updated blm object
#' @export
set_sampling_options <- function(x, ...) {
  UseMethod("set_sampling_options", x)
}

#' Update the sampler used for a specific variable
#'
#' The use can set the sampler to either (1) 'Gibbs' (default for all coefficients) or (2) 'MH' (metropolis-hastings algorithm, only for coefficients).
#'
#' @param x a blm object
#' @param ... other arguments passed to the function
#'
#' @return updated blm object
#' @export
set_sampler <- function(x, ...) {
  UseMethod("set_sampler")
}

#' Set initial values for a chain
#'
#' @param x blm object
#' @param ... other arguments passed to the function
#'
#' @return updated blm object
#' @export
set_initial_values <- function(x, ...) {
  UseMethod("set_initial_values")
}

#' Set priors for a blm object
#'
#' @param x blm object
#' @param ... optional arguments needed to specify a prior
#'
#' @return updated blm object with new priors
#' @export
set_prior <- function(x, ...) {
  UseMethod("set_prior")
}

#' Sample the posterior distribution
#'
#' @param x blm object
#'
#' @return blm object containing sampled posterior
#' @export
sample_posterior <- function(x) {
  UseMethod("sample_posterior", x)
}

#' Draw more samples from a sampled posterior distribution.
#'
#' @param x blm object
#' @param ... other options passed to function
#'
#' @export
update_posterior <- function(x, ...) {
  UseMethod("update_posterior", x)
}

#' Remove the posterior samples from a blm object.
#'
#' Generally, there is no good reason to delete the posterior samples. The only situation in which this may occur if the user has already sampled the posterior and wants to change the sampling plan, which is locked after the posterior is sampled.
#'
#' @param x blm object
#'
#' @return blm object minus posterior samples
#' @export
delete_posterior <- function(x) {
  UseMethod("delete_posterior", x)
}

#' Retrieve the posterior samples as a data frame from a blm object.
#'
#' @param x blm object
#'
#' @return a matrix of dimensions iterations x (variables + 1)
#' @export
get_posterior_samples <- function(x) {
  UseMethod("get_posterior_samples", x)
}

#' Set up posterior predictive checks
#'
#' @param x blm object
#'
#' @return prints summary of the ppc to the R console
#' @export
evaluate_ppc <- function(x, ...) {
  UseMethod("evaluate_ppc")
}

#' Assess blm model fit
#'
#' @param x blm object
#'
#' @return prints summary of model fit (DIC) to the R console
#' @export
evaluate_model_fit <- function(x) {
  UseMethod("evaluate_model_fit")
}

#' Print diagnostics
#'
#' @param x blm object
#'
#' @return prints summary of convergence diagnostics to R console
#' @export
convergence_diagnostics <- function(x) {
  UseMethod("convergence_diagnostics")
}

## Internal use only, not exported -----

#' Get a value from an S3 object
#'
#' @param x object with a get_value() method
#' @param var string. name of the value to retrieve from the object.
get_value <- function(x, var) {
  UseMethod("get_value")
}

#' Set a value for an S3 object
#'
#' @param x object with a set_value() method
#' @param var string. name of the value you want to change
#' @param val new value for var
set_value <- function(x, var, val) {
  UseMethod("set_value")
}

# Set priors on an S3 object (blm/priors)
set_priors <- function(x, ...) {
  UseMethod("set_priors")
}

# Set options (sampler object)
set_options <- function(x, ...) {
  UseMethod("set_options")
}

# Given a sampling plan (sampler object), sample the posterior distribution (resulting in a posterior object)
postsamp <- function(x,...) {
  UseMethod("postsamp")
}

# Combine existing posterior samples with newly sampled posterior samples (posterior object)
append_samples <- function(x, ...) {
  UseMethod("append_samples")
}

# Burn n samples from a posterior distribution (posterior object)
burn <- function(x, ...) {
  UseMethod("burn")
}

# Bind all samples from all chains together (posterior object)
bind <- function(x, ...) {
  UseMethod("bind")
}

# Retrieve maximum a posteriori (MAP) estimates (posterior object)
MAP <- function(x) {
  UseMethod("MAP", x)
}

# Retrieve 95% credible interval (posterior object)
CCI <- function(x) {
  UseMethod("CCI", x)
}

# Calculate the gelman-rubin statistic (posterior object)
GR <- function(x, ...) {
  UseMethod("GR", x)
}

# Calculate burn-in diagnostic (posterior object)
burnin_diagnostic <- function(x, ...) {
  UseMethod("burnin_diagnostic", x)
}

# ---------- Functions for priors ---------------

#' Draw a value from a prior distribution
#'
#' @param x prior object
#'
#' @return random variate drawn from a Normal or Gamma distribution
draw_value <- function(x) {
  UseMethod("draw_value")
}

# ---------- Functions for posterior predictive checks ----------

#' Check a posterior sample for normality assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for normality assumption
normality_check <- function(x, ...) {
  UseMethod("normality_check")
}

#' Check a posterior sample for homoskedasticity assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for homoskedasticity assumption
homoskedast_check <- function(x, ...) {
  UseMethod("homoskedast_check")
}

#' Check a posterior sample for independence of errors assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for independence of errors assumption
independence_check <- function(x, ...) {
  UseMethod("independence_check")
}

#' Check a posterior sample for the assumption that errors are normally distributed with mean mu and sd sigma
#'
#' @param x ppc object
#'
#' @return posterior predictive check on normal distribution of errors assumption
norm_of_errors_check <- function(x, ...) {
  UseMethod("norm_of_errors_check")
}
