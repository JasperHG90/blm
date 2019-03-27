# Bayesian Linear Model (`blm`)

The blm library contains tools to run Bayesian linear regression models in R. Some of its features are:

- Lazy evaluation of linear model objects
- Automatic generation of starting values
- Automatically center continuous variables
- Core functions are fast thanks to Julia code
- Gibbs & Metropolis-Hastings samplers

It further contains the following features to assess the model:

- Model fit (DIC)
- Posterior predictive checks for linear regression assumptions (normality, independence and homoskedasticity)
- Sampling history, autocorrelation, and density plots
- Gelman-rubin statistic
- Bayes factor

## Installation

* Julia (>= v1.0.0)
* R
* blm (install gh)

