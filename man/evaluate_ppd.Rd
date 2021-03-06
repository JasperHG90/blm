% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{evaluate_ppd}
\alias{evaluate_ppd}
\title{Evaluate the posterior predictive distributions for each observation}
\usage{
evaluate_ppd(x)
}
\arguments{
\item{x}{a blm object}
}
\value{
object of class "outliers" is added to the blm object
}
\description{
This function simulates outcome variables for each of the posterior parameters. It then returns the proportion where:
}
\details{
\deqn{p = y_{\text{simulated}} > y_{\text{observed}}}

This allows the user to look for patterns against e.g. the predictors or the outcome. See examples and the reference to Lynch below.
}
\seealso{
Lynch, S. M. (2007). Introduction to applied Bayesian statistics and estimation for social scientists. Springer Science & Business Media. pp.178-182
}
