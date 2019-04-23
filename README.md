# Bayesian Linear Model (`blm`)

I created this R library to implement some core Bayesian ideas taught in the course "Introduction to Bayesian Statitics" at Utrecht University. Secondly, it also provided me with an opportunity to further practice Julia [LINK]. Julia code is fast, but needs to compile on the first run. Hence, the highly repetitive nature of the code for Bayesian statistics (small functions that are repeatedly used) makes Julia an ideal programming language to run Markov Chain Monte Carlo (MCMC) samplers.

Thanks to the R library JuliaCall [LINK], it is possible to create a near-seamless bridge between R code and Julia code (much like the Reticulate [LINK] library does for Python). This bridge is used as follows:

1. All model-building aspects are constructed in R. Users specify a model much like they would with the `lm()` function. blm imports the magrittr library [LINK] to facilitate model-building, and allows the user to build a model using a 'tidy' workflow. The functions associated with this part of the model begin with `set_` (e.g. `set_priors()`). 
2. After the model is specified by the user, the sampling plan is executed (much like the `collect()` command from the dbplyr library [LINK]). This part calls the Julia code to draw posterior samples. The functions associated with this part of the model end with `_posterior` (e.g. `sample_posterior()`, `update_posterior()`, `delete_posterior()`).
3. Once the posterior distribution is sampled, all plotting functionality and summary statistics of the posterior is executed in R. If additional samples or computations need to be carried out (e.g. posterior predictive checks or DIC), these are executed in Julia. The functions associated with this part start with `evaluate_` (e.g. `evaluate_model_fit()`).

All core functions return S3 objects. Data embedded in these objects can be retrieved using the `get_value()` function (exported to facilitate the tidy workflow).

The following documents contain specific information or implementation notes:

- Course summary
- Implementation notes Gibbs sampler
- Implementation notes Metropolis-Hastings sampler

## Installation

Installing the R library from GitHub is straightforward

```r
devtools::install_github("JasperHG90/blm")
```

R dependencies will be installed automatically. You will also need an installation of Julia (>= 1.0.0). [MORE INFO].

To test whether R can find the Julia installation, execute the following in R:

```r
library(blm)
# Create the bridge between R and Julia
blm_setup()
```

If this succeeds, you should see the following message:

```shell
Julia version 1.0.3 at location /home/jasper/julia-1.0.3/bin will be used.
Loading setup script for JuliaCall...
Finish loading setup script for JuliaCall.
```

This means you are good to go.

## Short example

To build a Bayesian Linear Model (blm) object, start by executing the following:

```r
# Load data
data("directors")
# Log the compensation variable ('000 GBR)
directors$Compensation <- log(directors$Compensation)
# Build the model
dirfit <- blm("Compensation ~ Age", data=directors)
```

This creates a blm object for the 'directors' data. At this point, we have a model with uninformative priors (mu=0, sd=1000) and initial values that are drawn from these priors:

```r
print(dirfit)
```

```shell
Bayesian Linear Model (BLM) object:

Data:
	Predictors: 1
	Outcome: Compensation

Sampler:
	Chains: 1
	Iterations: 10000
	Thinning: 1
	Burn: 1000

Priors (Coefficients) :
      b0   b1
mu     0    0
tau 1000 1000

Priors (Residuals) :
      sigma
rate   0.01
scale  0.01
```

Furthermore, the sampling settings are set to just $1$ chain and $10.000$ iterations. We can update the model as follows:

```r
# Set the following values for the dirfit model
dirfit <- dirfit %>%
  # Sampling settings
  set_sampling_options(chains=2, iterations=20000, burn=3000, thinning=2) %>%
  # Set priors for Age (we're on a log scale)
  set_prior("b1", mu=0.2, sd = 0.001) %>%
  # Set initial values
  set_initial_values(chain_1 = list("b" = c(-1, 1), "sigma" = 1),
                     chain_2 = list("b" = c(1, -1), "sigma" = 2))
# Print object
print(dirfit)
```

```shell
Bayesian Linear Model (BLM) object:

Data:
	Predictors: 1
	Outcome: Compensation

Sampler:
	Chains: 2
	Iterations: 20000
	Thinning: 2
	Burn: 3000

Priors (Coefficients) :
      b0    b1
mu     0 0.200
tau 1000 0.001

Priors (Residuals) :
      sigma
rate   0.01
scale  0.01
```

View the full example here [LINK].
