# S3 generics

#' Get a value from an S3 object
#'
#' @param x object with a get_value() method
#' @param var string. name of the value to retrieve from the object.
#'
#' @export
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

#' Set priors for a blm object
#'
#' @param x blm object
#' @param ... optional arguments needed to specify a prior
#'
#' @return blm object with new priors
#' @export
set_priors <- function(x, ...) {
  UseMethod("set_priors")
}

set_prior <- function(x, ...) {
  UseMethod("set_prior")
}

set_options <- function(x, ...) {
  UseMethod("set_options")
}

set_sampler <- function(x, ...) {
  UseMethod("set_sampler")
}

set_initial_values <- function(x, ...) {
  UseMethod("set_initial_values")
}

#' Sample the posterior distribution
#'
#' @param blm blm object
#' @param chains number of chains to run
#' @param iterations number of iterations per chain
#' @param burn number of burn-in iterations
#' @param julia use julia for computations?
#'
#' @return updated blm object
#' @export
set_sampling_options <- function(x, ...) {
  UseMethod("set_sampling_options", x)
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

update_posterior <- function(x) {
  UseMethod("update_posterior", x)
}

delete_posterior <- function(x) {
  UseMethod("delete_posterior", x)
}

#' Set up posterior predictive checks
#'
#' @param x blm object
#'
#' @return prints summary of the ppc to the R console
#' @export
posterior_predictive_checks <- function(x, ...) {
  UseMethod("posterior_predictive_checks")
}

#' Assess blm model fit
#'
#' @param x blm object
#'
#' @return prints summary of model fit (DIC) to the R console
#' @export
model_fit <- function(x) {
  UseMethod("model_fit")
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

# ---------- Functions for priors ---------------

#' Draw a value from a prior distribution
#'
#' @param x prior object
#'
#' @return random variate drawn from a Normal or Gamma distribution
#' @export
draw_value <- function(x) {
  UseMethod("draw_value")
}

# ---------- Functions for posterior predictive checks ----------

#' Check a posterior sample for normality assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for normality assumption
#' @export
normality_check <- function(x, ...) {
  UseMethod("normality_check")
}

#' Check a posterior sample for homoskedasticity assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for homoskedasticity assumption
#' @export
homoskedast_check <- function(x, ...) {
  UseMethod("homoskedast_check")
}

#' Check a posterior sample for independence of errors assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for independence of errors assumption
#' @export
independence_check <- function(x, ...) {
  UseMethod("independence_check")
}

#' Check a posterior sample for the assumption that errors are normally distributed with mean mu and sd sigma
#'
#' @param x ppc object
#'
#' @return posterior predictive check on normal distribution of errors assumption
#' @export
norm_of_errors_check <- function(x, ...) {
  UseMethod("norm_of_errors_check")
}
