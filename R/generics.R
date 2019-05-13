# S3 generics

## Exported ----

#' Set sampling options on a blm object
#'
#' @param x blm object
#' @param chains number of chains to initialize
#' @param iterations number of iterations to run the MC sampler
#' @param thinning set thinning parameter
#' @param burn number of samples that will be discarded
#'
#' @examples
#' bfit <- set_sampling_options(bfit, chains=2, iterations=20000, thinning=3, burn = 2000)
#'
#' @importFrom magrittr '%>%'
#' @return updated blm object
#' @export
set_sampling_options <- function(x, chains, iterations, thinning, burn) {
  UseMethod("set_sampling_options", x)
}

#' Update the sampler used for a specific variable
#'
#' The use can set the sampler to either (1) 'Gibbs' (default for all coefficients) or (2) 'MH' (metropolis-hastings algorithm, only for coefficients).
#'
#' @param x a blm object
#' @param par name of the coefficient, given as b0 (intercept) or b1, ... for coefficients.
#' @param type name of the sampler to use. Either 'Gibbs' or 'MH'. Defaults to 'Gibbs'.
#' @param lambda tuning parameter for the variance of the proposal distribution in metropolis-hastings
#'
#' @examples
#' bfit <- set_sampler(bfit, "b0", type="MH")
#'
#' @importFrom magrittr '%>%'
#'
#' @return updated blm object
#' @export
set_sampler <- function(x, ...) {
  UseMethod("set_sampler", x)
}

#' Set initial values for a chain
#'
#' @param x blm object
#' @param ... named argument containing (1) a chain name (e.g. 'chain_1') and (2) a list of initial values. See examples.
#'
#' @importFrom magrittr '%>%'
#'
#' @examples
#' bfit <- set_initial_values(bfit, chain_1 = list("b" = c(2,5,3,4,5), sigma=0.1))
#'
#' @return updated blm object
#' @export
set_initial_values <- function(x, ...) {
  UseMethod("set_initial_values", x)
}

#' Change priors for available parameters
#'
#' @param par name of the coefficient, given as b0 (intercept) or b1, ... for coefficients.
#' @param ... other arguments passed to function. In the case of intercept and coefficients, these are 'mu' and 'sd' (prior means and variances). In the case of the residual variance, the rate and shape/scale parameter must be passed as 'alpha' and 'beta'
#'
#' @examples
#' bfit <- set_prior(bfit,"b0", mu=40, sd=7)
#' bfit <- set_prior(bfit, "sigma", alpha=2, beta=1)
#' @importFrom magrittr '%>%'
#'
#' @return updated blm object with new priors
#' @export
set_prior <- function(x, par, ...) {
  UseMethod("set_prior", x)
}

#' Add a hypothesis to a blm model
#'
#' To evaluate informative hypotheses, users may add these using this function. See details for more information.
#'
#' @details
#' This function accepts two kind of hypotheses. The first are 'simple' or single hypotheses, such as:
#' \itemize{
#'    \item{a > b}
#'    \item{a < b}
#'    \item{a = b}
#'    \item{a-b < 1}
#'    \item{2*a > 0}
#'    \item{...}
#' }
#' The second kind are 'complex' or 'multiple' hypotheses which are chained with '&' (AND). Examples are:
#' \itemize{
#'    \item{a = 0 & b = 0}
#'    \item{a < b & b < c}
#'    \item{a < .1 & a > .1 OR |a| > .1}
#'    \item{a-b < .1 & a-b > .1 OR |a-b| > .1}
#'    \item{...}
#' }
#'
#' @param x a blm object
#' @param name name for the hypothesis (e.g. 'H1')
#' @param hypothesis_user a hypothesis.
#'
#' @return a blm object with a new or updated 'hypotheses' object
#'
#' @seealso cite Hoitink
#' @export
set_hypothesis <- function(x, name, hypothesis_user) {
  UseMethod("set_hypothesis", x)
}

#' Sample the posterior distribution
#'
#' Calling the \code{sample_posterior()} function executes a blm sampling plan. The user should set priors, initial values etc. before executing this function, because by executing it the user will 'lock' the sampling plan. That is, no changes may be made except for changing the 'burn' parameter.
#'
#' @details
#' For details on the implementation of the Gibbs sampler, see https://tinyurl.com/y5vk9x35. For details on the implementation of the Metropolis-Hastings sampler, see https://tinyurl.com/y3zp3l3w.
#'
#' @param x blm object
#'
#' @return blm object containing sampled posterior
#'
#' @importFrom magrittr '%>%'
#'
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    sample_posterior()
#'
#' @export
sample_posterior <- function(x) {
  UseMethod("sample_posterior", x)
}

#' Draw more samples from a sampled posterior distribution.
#'
#' @param x blm object
#' @param iterations number of additional samples to draw from the posterior.
#'
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    sample_posterior()
#' # Update posterior samples
#' fit <- fit %>% update_posterior(iterations=2000)
#'
#' @return blm object with updated count of posterior samples.
#' @importFrom magrittr '%>%'
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
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    sample_posterior()
#' # Delete posterior samples
#' fit <- fit %>% delete_posterior
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
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    sample_posterior()
#' # Retrieve samples (burned automatically)
#' psamples <- get_posterior_samples(fit)
#'
#' @return a matrix of dimensions iterations x (variables + 1)
#' @export
get_posterior_samples <- function(x) {
  UseMethod("get_posterior_samples", x)
}

#' Compute the intercept-only model for comparison purposes
#'
#' @param x blm object
#' @param iterations number of iterations used to sample the model
#' @param chains number of separate chains used to sample the model
#' @param thinning thinning parameter
#' @param burn number of examples to burn
#'
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    compute_null_model(iterations=10000, burn=1000) %>%
#'    sample_posterior()
#' # Plot the null model
#' plot(fit, "nullmodel")
#'
#' @return blm object containing instructions on how to sample the null model
#'
#' @export
compute_null_model <- function(x, iterations=10000, chains=1, thinning=1,  burn=1000) {
  UseMethod("compute_null_model", x)
}

#' Set up posterior predictive checks
#'
#' @param x blm object
#' @param return_all logical. If TRUE, then the function will return the discrepancy statistics for the simulated and observed data.
#' @param p proportion of samples used to compute posterior predictive checks. Defaults to p=1.
#'
#' @return prints summary of the ppc to the R console
#' @export
evaluate_ppc <- function(x, return_all=TRUE, p=1) {
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

#' Evaluate hypotheses specified by the user
#'
#' @param x blm object
#'
#' @return returns the blm object containing Bayes' Factors
#'
#' @seealso
#' \itemize{
#'   \item{Hoijtink, H., Mulder, J., Van Lissa, C. J., & Gu, X. (2019). A tutorial on testing hypotheses using the Bayes factor.}
#'   \item{Hoijtink, H., Gu, X., & Mulder, J. (2018). Bayesian evaluation of informative hypotheses for multiple populations. British Journal of Mathematical and Statistical Psychology.}
#'   \item{Gu, X., Mulder, J., & Hoijtink, H. (2018). Approximated adjusted fractional Bayes factors: A general method for testing informative hypotheses. British Journal of Mathematical and Statistical Psychology, 71(2), 229-261.}
#'   \item{Gu, X., Hoijtink, H., & Mulder, J. (2016). Error probabilities in default Bayesian hypothesis testing. Journal of Mathematical Psychology, 72, 130-143.}
#' }
#'
#' @export
evaluate_hypotheses <- function(x) {
  UseMethod("evaluate_hypotheses", x)
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
#' @importFrom magrittr '%>%'
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
evaluate_R2 <- function(x) {
  UseMethod("evaluate_R2", x)
}

#' Evaluate the posterior predictive distributions for each observation
#'
#' This function simulates outcome variables for each of the posterior parameters. It then returns the proportion where:
#'
#' \deqn{p = y_{\text{simulated}} > y_{\text{observed}}}
#'
#' This allows the user to look for patterns against e.g. the predictors or the outcome. See examples and the reference to Lynch below.
#'
#' @param x a blm object
#'
#' @return object of class "outliers" is added to the blm object
#'
#' @seealso Lynch, S. M. (2007). Introduction to applied Bayesian statistics and estimation for social scientists. Springer Science & Business Media. pp.178-182
#' @export
evaluate_ppd <- function(x) {
  UseMethod("evaluate_ppd", x)
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
#'
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    sample_posterior()
#' # Get DIC and print
#' fit %>% get_value("DIC") %>% summary()
#'
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

# Sample the null model
sample_null_model <- function(x) {
  UseMethod("sample_null_model", x)
}

#' Evaluate the model Bayes Factor
#'
#' The model Bayes Factor computes the evidence of the model versus the evidence found in the intercept-only model. To compute the model Bayes Factor, the user needs to explicitly state that they want to compute the null model (see examples). This Bayes Factor value is approximated using Wagemakers' approximation (see below for citation).
#'
#' @seealso Wagenmakers, E. J. (2007). A practical solution to the pervasive problems of p values. Psychonomic bulletin & review, 14(5), pp. 796-799.
#'
#' @examples
#' # Load the directors data
#' data("directors")
#' # Fit a blm model
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    # Compute the null model next to the desired model
#'    compute_null_model() %>%
#'    # Sample posterior data
#'    sample_posterior()
#' # Check the convergence of the null model
#' plot(fit, "nullmodel")
#' # The Bayes Factor is in the summary
#' summary(fit)
#'
#' @return blm object with additional slot called 'model_BF' containing Bayes Factor information
#' @export
evaluate_model_BF <- function(x) {
  UseMethod("evaluate_model_BF", x)
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
