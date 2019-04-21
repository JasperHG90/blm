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
  UseMethod("set_sampler", x)
}

#' Set initial values for a chain
#'
#' @param x blm object
#' @param ... other arguments passed to the function
#'
#' @return updated blm object
#' @export
set_initial_values <- function(x, ...) {
  UseMethod("set_initial_values", x)
}

#' Set priors for a blm object
#'
#' @param x blm object
#' @param ... optional arguments needed to specify a prior
#'
#' @return updated blm object with new priors
#' @export
set_prior <- function(x, ...) {
  UseMethod("set_prior", x)
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
  UseMethod("evaluate_ppc", x)
}

#' Assess blm model fit
#'
#' @param x blm object
#'
#' @return prints summary of model fit (DIC) to the R console
#' @export
evaluate_model_fit <- function(x) {
  UseMethod("evaluate_model_fit", x)
}

#' Print diagnostics
#'
#' @param x blm object
#'
#' @return prints summary of convergence diagnostics to R console
#' @export
evaluate_convergence_diagnostics <- function(x) {
  UseMethod("evaluate_convergence_diagnostics", x)
}

#' Calculate the effective sample size
#'
#' @param x a blm object
#'
#' @return prints effective sample size information to screen
#'
#' @export
evaluate_effective_sample_size <- function(x) {
  UseMethod("evaluate_effective_sample_size", x)
}

#' View the number of accepted posterior draws for each parameter and chain
#'
#' If a user chooses to use the Gibbs sampler only, then the number of accepted draws is equal to the number of iterations in each chain. For the Metropolis-Hastings algorithm, this depends on the value of lambda set by the user in the \link[blm]{set_sampler} function.
#'
#' @param x a blm object
#'
#' @return prints the number of accepted draws for each parameter in each chain to the console
#' @export
evaluate_accepted_draws <- function(x) {
  UseMethod("evaluate_accepted_draws", x)
}

#' Calculate the Bayesian R-squared value
#'
#' @param x a blm object
#'
#' @return object of class "R2" containing input parameters, posterior draws and r-squared values computed on the posterior draws
#'
#' @seealso Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2018). R-squared for Bayesian regression models. The American Statistician, (just-accepted), 1-6.
#' @export
evaluate_R2 <- function(x, ...) {
  UseMethod("evaluate_R2", x)
}

#' Retrieve mapping from variable to parameter names
#'
#' blm objects will mostly emit parameter names (b0, b1, etc) instead of variable names. This convenience function returns the mapping from parameter names to variable names.
#'
#' @param x a blm object
#'
#' @return prints mapping from parameter to variable names to the console
#'
#' @export
get_parameter_names <- function(x) {
  UseMethod("get_parameter_names", x)
}

#' Get a value from an S3 object
#'
#' @param x object with a get_value() method
#' @param var string. name of the value to retrieve from the object.
#' @export
get_value <- function(x, var) {
  UseMethod("get_value", x)
}

## Internal use only, not exported -----

#' Set a value for an S3 object
#'
#' @param x object with a set_value() method
#' @param var string. name of the value you want to change
#' @param val new value for var
set_value <- function(x, var, val) {
  UseMethod("set_value", x)
}

# Set priors on an S3 object (blm/priors)
set_priors <- function(x, ...) {
  UseMethod("set_priors", x)
}

# Set options (sampler object)
set_options <- function(x, ...) {
  UseMethod("set_options", x)
}

# Given a sampling plan (sampler object), sample the posterior distribution (resulting in a posterior object)
postsamp <- function(x,...) {
  UseMethod("postsamp", x)
}

# Combine existing posterior samples with newly sampled posterior samples (posterior object)
append_samples <- function(x, ...) {
  UseMethod("append_samples", x)
}

# Burn n samples from a posterior distribution (posterior object)
burn <- function(x, ...) {
  UseMethod("burn", x)
}

# Bind all samples from all chains together (posterior object)
bind <- function(x, ...) {
  UseMethod("bind", x)
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

# Calculate effective sample size (posterior object)
effss <- function(x, ...) {
  UseMethod("effss", x)
}

# ---------- Functions for priors ---------------

#' Draw a value from a prior distribution
#'
#' @param x prior object
#'
#' @return random variate drawn from a Normal or Gamma distribution
draw_value <- function(x) {
  UseMethod("draw_value", x)
}

# ---------- Functions for posterior predictive checks ----------

#' Check a posterior sample for normality assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for normality assumption
normality_check <- function(x, ...) {
  UseMethod("normality_check", x)
}

#' Check a posterior sample for homoskedasticity assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for homoskedasticity assumption
homoskedast_check <- function(x, ...) {
  UseMethod("homoskedast_check", x)
}

#' Check a posterior sample for independence of errors assumption
#'
#' @param x ppc object
#'
#' @return posterior predictive check for independence of errors assumption
independence_check <- function(x, ...) {
  UseMethod("independence_check", x)
}

#' Check a posterior sample for the assumption that errors are normally distributed with mean mu and sd sigma
#'
#' @param x ppc object
#'
#' @return posterior predictive check on normal distribution of errors assumption
norm_of_errors_check <- function(x, ...) {
  UseMethod("norm_of_errors_check", x)
}
