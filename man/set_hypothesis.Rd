% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{set_hypothesis}
\alias{set_hypothesis}
\title{Add a hypothesis to a blm model}
\usage{
set_hypothesis(x, name, hypothesis_user)
}
\arguments{
\item{x}{a blm object}

\item{name}{name for the hypothesis (e.g. 'H1')}

\item{hypothesis_user}{a hypothesis.}
}
\value{
a blm object with a new or updated 'hypotheses' object
}
\description{
To evaluate informative hypotheses, users may add these using this function. See details for more information.
}
\details{
This function accepts two kind of hypotheses. The first are 'simple' or single hypotheses, such as:
\itemize{
   \item{a > b}
   \item{a < b}
   \item{a = b}
   \item{a-b < 1}
   \item{2*a > 0}
   \item{...}
}
The second kind are 'complex' or 'multiple' hypotheses which are chained with '&' (AND). Examples are:
\itemize{
   \item{a = 0 & b = 0}
   \item{a < b & b < c}
   \item{a < .1 & a > .1 OR |a| > .1}
   \item{a-b < .1 & a-b > .1 OR |a-b| > .1}
   \item{...}
}
}
\seealso{
cite Hoitink
}
