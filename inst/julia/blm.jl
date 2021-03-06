#=

BLM (Bayesian Linear Model) helper functions in Julia

Written by: Jasper Ginn
Course: Introduction to Bayesian Statistics
          @Utrecht University
          @instructor1: Prof. Herbert Hoijtink
          @instructor2: Fayette Klaassen

=#

#=
PART I: MCMC sampler in in Julia
=#

function MCMC_sampler(X::Array{Float64}, y::Array{Float64}, w,
                       sigma::Float64, iterations::Int, thinning::Int, priors,
                       samplers)

    #=
    Run the MCMC sampler to obtain the conditional posterior coefficients / sigma values

    :param X: Design matrix containing j+1 coefficients (first coefficient must equal intercept variable). The design matrix is created in R using model.matrix()
    :param y: Numeric array with outcome variables. Length(y) == nrow(X)
    :param w: Coefficient matrix containing initial values. nrow(w) == ncol(X)
    :param sigma: Single initial value for the variance.
    :param iterations: Number of sampling iteratios in each chain
    :param thinning: thinning parameter. The value for thinning determines the observation that will be used as a draw from the posterior.
                        That is, if thinning = 3, then the first two draws will be discarded and only the third draw will be returned as a draw
                        from the posterior distribution.
    :param priors: Dict containing information on the priors of each coefficient and sigma
    :param samplers: Dict containing information on which sampling strategy to follow for each coefficient (gibbs or MH)

    :return: Tensor of shape N * J (iterations * parameters).
    :seealso:
        - Lynch, S. M. (2007). Introduction to applied Bayesian statistics and estimation for social scientists. Springer Science & Business Media.
        - Gelman, A., Stern, H. S., Carlin, J. B., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian data analysis. Chapman and Hall/CRC.
    =#

    # Check if w is an Array as required by this function.
    # This is an issue of conversion between julia & R
    if isa(w,AbstractArray) == false

        w = [w]

    end;

    # Open up results for accepted draws
    accepted = zeros(size(X)[2] + 1)

    # Open up a results matrix
    res = Array{Float64}(undef, iterations, (size(X)[2] + 1))

    # For each iteration
    for i in 1:iterations

        for j in 1:thinning

          # One step of the MCMC sampler
          r = MCMC_one_step(X, y, w, sigma, priors, samplers)

          # Set new values for weights and sigma
          w = r["w"]
          sigma = r["sigma"]
          # If j equal to thinning, add accepts
          if j == thinning
              accepted += r["accepted"]
            end;

          end;

        # Add to results
        res[i,1:end-1] = w
        res[i,end] = sigma

        end;

    # Square-root the sigma values
    res[:,end] = sqrt.(res[:,end])

    # Return results
    return(
        Dict(
            "posterior" => res,
            "accepted" => accepted
        )
    )

    end;

function MCMC_one_step(X::Array{Float64}, y::Array{Float64}, w::Array{Float64},
                        sigma::Float64, priors, samplers)

    #=
    Run a single iteration of the MCMC sampler

    :param X: Design matrix containing j+1 coefficients (first coefficient must equal intercept variable). The design matrix is created in R using model.matrix()
    :param y: Numeric array with outcome variables. Length(y) == nrow(X)
    :param w: Coefficient matrix containing initial values. nrow(w) == ncol(X)
    :param sigma: Single initial value for the variance.
    :param priors: Dict containing information on the priors of each coefficient and sigma
    :param samplers: Dict containing information on which sampling strategy to follow for each coefficient (gibbs or MH)

    :return: Dict containing draws for each conditional posterior distribution.
    =#

    # Local value for accepted draws
    accepted = zeros(size(X)[2] + 1)

    # For each prior / parameter, do
    for j in 1:length(priors)

        # If j == nparam, then sigma
        if j == length(priors)

            sigma = gibbs_one_step_resvar(X, w, y, priors[j][:alpha], priors[j][:beta])
            # Increment accepted
            accepted[j] += 1

        # Else, coefficient
        # See the Appendix in 'docs' folder (section 5) for an elaboration on the calculation of posterior values
        else

            ## Check sampler
            ## Use Gibbs sampling? ==> sample conditional posteriors
            if samplers[j][1] == "Gibbs"

                w[j] = gibbs_one_step_coef(X, y, w, sigma, priors[j][:mu], priors[j][:sd], j)
                # Increment accepted
                accepted[j] += 1

            ## Use MH sampling? ==> sample conditional posteriors that are proportional up to a constant
        elseif samplers[j][1] == "MH"

                local b_current::Dict

                # Call Metropolis-Hastings algorithm
                b_current = MH_one_step_coef(w[j], X[:, j], X[:, 1:end .!= j], y,
                                            w[1:end .!= j], sigma, priors[j][:mu],
                                            priors[j][:sd], samplers[j][2])

                # Set new value for w
                w[j] = b_current["samp"]

                # Set accepted
                accepted[j] = b_current["accepted"]

            else

                # Unknown sampler. ERROR OUT

                end;

            end;

        end;

    # Return the coefficient matrix
    return(
        Dict(
            "w" => w,
            "sigma" => sigma,
            "accepted" => accepted
        )
    )

    end;

function MH_one_step_coef(b_previous::Float64, xj::Array{Float64}, X::Array{Float64}, y::Array{Float64},
                         w::Array{Float64}, sigma::Float64, prior_mu::Float64, prior_tau::Float64,
                         zeta::Float64)

    #=
    Perform one draw (iteration) of the Metropolis-Hastings sampler

    :param b_previous: Scalar. Previous value of the jth coefficient.
    :param xj: Array. values of the jth column of the design matrix X as a column vector
    :param X: Array. design matrix X excluding the jth column
    :param y: Array. outcome values y
    :param w: Array. parameter matrix w excluding the jth row (i.e. the jth coefficient)
    :param sigma: Scalar. variance parameter
    :param prior_mu: Scalar. prior mean
    :param prior_tau: Scalar. prior variance
    :para, zeta: Scalar. tuning parameter. Controls the variance of the proposal distribution (a normal).

    :return: A dict containing:
               - samp: either a new value or b_proposed
               - accepted: indicator (1/0) if proposal was accepted

    :seealso:
        - 'MH.pdf' in '/docs' folder for implementation notes
    =#

    # Draw proposal
    b_proposed = rand(Normal(b_previous, zeta))

    # Draw from the posterior (current coef and previous coef) & ompare difference to logged random uniform variate
    if (posterior_coef(b_proposed, xj, X, y, w, sigma, prior_mu, prior_tau) - posterior_coef(b_previous, xj, X, y, w, sigma, prior_mu, prior_tau)) < log(rand(Uniform(0, 1)))

        return(
            Dict(
                "samp" => b_previous,
                "accepted" => 0
            )
        )

    else

        return(
            Dict(
                "samp" => b_proposed,
                "accepted" => 1
            )
        )

        end;

    end;

function gibbs_one_step_coef(X::Array{Float64}, y::Array{Float64}, w::Array{Float64},
                             sigma::Float64, prior_mu::Float64, prior_sd::Float64, j::Int)

     #=
     Perform one draw (iteration) of the Gibbs sampler

     :param X: Array. design matrix X
     :param y: Array. outcome values y
     :param w: Array. parameter matrix w
     :param sigma: Scalar. variance parameter
     :param prior_mu: Scalar. prior mean
     :param prior_tau: Scalar. prior variance
     :param j: Scalar. indicator value for the jth coefficient. (e.g. if this is the second coefficient then j=2)

     :return: single draw of the posterior distribution
     =#

    return(
        rand(Normal(posterior_mu(X[:,1:end .!= j], X[:,j], y, w[1:end .!= j], sigma,
                                        prior_mu, prior_sd),
                    sqrt(posterior_tau(X[:,j], sigma, prior_sd))))
    )

    end;

function gibbs_one_step_resvar(X::Array{Float64}, w::Array{Float64}, y::Array{Float64},
                              prior_alpha::Float64, prior_beta::Float64)

      #=
      Perform one draw (iteration) of the Gibbs sampler (residual variance)

      :param X: Array. design matrix X
      :param y: Array. outcome values y
      :param w: Array. parameter matrix w
      :param prior_alpha: Scalar. prior scale
      :param prior_beta: Scalar. prior rate

      :return: single draw of the posterior distribution
      =#


    return(
        rand(InverseGamma(posterior_rate(size(X)[1], prior_alpha),
                          posterior_scale(X, w, y, prior_beta)))
    )

    end;

## Posterior densities --------

## Posterior densities
function posterior_mu(X::Array{Float64}, xj::Array{Float64}, y::Array{Float64},
                      w::Array{Float64}, sigma::Float64, mu0::Float64, tau0::Float64)

    #=
    Calculate the posterior mean for an intercept / slope coefficient

    :param X: The design matrix containing the data, where the jth column is removed
    :param xj: Column vector containing the data corresponding to the jth coefficient
    :param y: Data for the outcome variable
    :param w: Column vector containing the coefficients, where the jth row is removed
    :param sigma: Variance parameter. Updated iteratively.
    :param mu0: Prior mean
    :param tau0: Prior variance

    :return: posterior mean for jth coefficient
    :seealso: The file 'conditionalposterior.pdf' in the 'docs' folder for an elaboration on the derivation of the
                conditional posterior values for each parameter.
    =#

    # Numerator
    return(((sum(xj .* (y - X * w)) / sigma) + (mu0/tau0)) / ((transpose(xj) * xj / sigma) + (1 / tau0)))

    end;

function posterior_tau(xj::Array{Float64}, sigma::Float64, tau0::Float64)

    #=
    Calculate the posterior variance for an intercept / slope coefficient

    :param xj: Column vector containing the data corresponding to the jth coefficient
    :param sigma: Variance parameter. Updated iteratively.
    :param tau0: Prior variance

    :return: posterior variance for jth coefficient
    :seealso: The file 'conditionalposterior.pdf' in the 'docs' folder for an elaboration on the derivation of the
                conditional posterior values for each parameter.
    =#

    return( 1 / ( (transpose(xj) * (xj ./ sigma)) + (1/tau0)) )

    end;

function posterior_rate(n::Int, a0::Float64)

    #=
    Calculate the posterior rate for the sigma value

    :param n: number of examples in the data
    :param a0: prior rate value

    :return: posterior rate
    :seealso: The file 'conditionalposterior.pdf' in the 'docs' folder for an elaboration on the derivation of the
                conditional posterior values for each parameter.
    =#

    return( (n/2) + a0 )

    end;

function posterior_scale(X::Array{Float64}, w::Array{Float64}, y::Array{Float64}, b0::Float64)

    #=
    Calculate the posterior scale for the sigma value

    :param X: The design matrix containing the data, where the jth column is removed
    :param y: Data for the outcome variable
    :param w: Column vector containing the coefficients, where the jth row is removed
    :param b0: prior scale value

    :return: posterior scale
    :seealso: The file 'conditionalposterior.pdf' in the 'docs' folder for an elaboration on the derivation of the
                conditional posterior values for each parameter.
    =#

    return( ( sum( ( y - (X * w) ).^2 ) / 2 ) + b0 )

    end;

function posterior_coef(bj::Float64, xj::Array{Float64}, X::Array{Float64}, y::Array{Float64},
                        w::Array{Float64}, sigma::Float64, prior_mu::Float64, prior_tau::Float64)

    #=
    Sample the logged non-normalized posterior for a coefficient.

    :param b_previous: Scalar. Previous value of the jth coefficient.
    :param xj: Array. values of the jth column of the design matrix X as a column vector
    :param X: Array. design matrix X excluding the jth column
    :param y: Array. outcome values y
    :param w: Array. parameter matrix w excluding the jth row (i.e. the jth coefficient)
    :param sigma: Scalar. variance parameter
    :param prior_mu: Scalar. prior mean
    :param prior_tau: Scalar. prior variance
    :para, zeta: Scalar. tuning parameter. Controls the variance of the proposal distribution (a normal).

    :return: Scalar. single draw from the logged and unnormalized posterior distribution

    :seealso:
        - 'MH.pdf' in '/docs' folder for implementation notes
    =#

    left = -(bj)^2 * (( sum( transpose(xj) * xj ) / (2*sigma) ) + ( 1 / ( 2*prior_tau ) ))
    right = bj * ((sum( xj .* (y .- X * w) ) / sigma) + (prior_mu / prior_tau) )

    return(left + right)

    end;

#=
PART II: Posterior Predictive Checks
=#

function ppc_draws(X::Array{Float64}, y::Array{Float64}, sampled::Array{Float64})

    #=
    Perform posterior predictive checking.

    This function is used to draw data from the posterior using the gibbs sampler and to simulate values of the
    outcome variable y to compute posterior predictive checks. The function only simulates / computes the data. The
    statistic is calculated in R. The checks that are computed are:

        1. Normality of errors (using skewness).
        2. Homogeneity of variances (using adj. R^2 of a linear regression on squared residuals).
        3. Independence of errors (using correlation of lagged residuals).

    :param X: Design matrix containing j+1 coefficients (first coefficient must equal intercept variable). The design matrix is created in R using model.matrix()
    :param y: Numeric array with outcome variables. Length(y) == nrow(X)
    :param sampled: posterior samples for the parameters

    :return: Dictionary containing:
        1. sim_y              [==>] 2D Tensor ([iterations - burn] * 1) of simulated y values
        2. residuals          [==>] 3D Tensor ([iterations - burn] * n. examples * 2) containing residuals for simulated and observed data
        3. heteroskedasticity [==>] 2D Tensor ([iterations - burn] * 2) containing test values (simulatd and observed data) on the test for heteroskedasticity
        4. skewness           [==>] 2D Tensor ([iterations - burn] * 2) containing test values on the test for normality of errors
        5. independence       [==>] 2D Tensor ([iterations - burn] * 2) containing test values on the test for independence of errors
    =#

    # 1. Precompute linear combinations as a function of X and the posterior set of params theta = [b0,b1,...,bn]
    lincom = X * transpose(sampled[: , 1:end-1]) # lincom dims: rows(X) * (iterations-burn)

    # 2. Set up a results matrix for simulated yhat values
    res = Array{Float64}(undef, size(sampled)[1], size(X)[1])

    # 3. Set up a results matrix for residuals, skewness stats, heteroskedasticity stats and correlation stats
    resids = Array{Float64}(undef, size(sampled)[1], size(X)[1], 2)
    skewed = Array{Float64}(undef, size(sampled)[1], 2)
    heterosked = Array{Float64}(undef, size(sampled)[1], 2)
    correlate = Array{Float64}(undef, size(sampled)[1], 2)

    # 4. Draw from a normal distribution given each specific y value and populate results matrix
    for j in 1:(size(sampled)[1])

        # Simulate y
        res[j,:] = simulate_y(lincom[:,j], sampled[j,end])

        # Compute residuals (simulated)
        resids[j, :, 1] = res[j, :] .- lincom[:, j]

        # Compute residuals (observed)
        resids[j, :, 2] = y .- lincom[:, j]

        # Adjusted R^2 values for simulated & observed data
        heterosked[j, 1] = adjR2(X, resids[j, :, 1].^2)
        heterosked[j, 2] = adjR2(X, resids[j, :, 2].^2)

        # Skewness for simulated & observed data
        skewed[j, 1] = skewness(resids[j, :, 1])
        skewed[j, 2] = skewness(resids[j, :, 2])

        # Correlation for simulated & observed data
        correlate[j, 1] = independence(resids[j, :, 1], 1)
        correlate[j, 2] = independence(resids[j, :, 2], 1)

        end;

    # 5. Return
    return(Dict(
      "heteroskedasticity" => heterosked,
      "skewness" => skewed,
      "independence" => correlate
    ))

    end;

# Simulate values from a random
function simulate_y(computed_y::Array{Float64}, sd::Float64)

    #=
    Simulate a vector of outcome variable values for y, drawn from a normal distribution

    :param computed_y: single-dimensional array of values (these will be used as means in the random draw)
    :param sd: scalar. Single standard deviation value.

    :return: single-dimensional array with simulated values for y
    =#

    # Set up single-dimensional array
    res = Array{Float64}(undef, size(computed_y)[1], 1)

    # For each value in computed_y, generate a normal
    for i in 1:size(computed_y)[1]

        res[i,1] = rand(Normal(computed_y[i,1], sd))

        end;

    # Return
    return(res)

    end;

function predict_y(X::Array{Float64}, y::Array{Float64})

  #=
  Predict outcome variable y

  :param X: Array. Design matrix
  :param y: Array. Vector containing dependent variable values.

  :return: predicted value for y computed using OLS estimates
  =#

  # Linear model, prediction and return
  return(X * ( inv( transpose(X) * X ) * transpose(X) * y))

  end;

function adjR2(X::Array{Float64}, y::Array{Float64})

  #=
  Calculate sums of squares & R2

  :param X: Array. Design matrix
  :param y: Array. Vector containing dependent variable values.

  :return: Scalar. Adjusted R^2 value
  =#

  local TSS::Float64
  local RSS::Float64
  local R2::Float64

  # TSS
  TSS = sum((y .- mean(y)).^2)

  # RSS
  RSS = sum( (predict_y(X,y) .- y).^2 )

  # R2
  R2 = 1 - (RSS / TSS)

  # Adjusted R2
  return( 1 - ((1 - R2) * (size(X)[1] - 1)) / (size(X)[1] - (size(X)[2]-1) - 1 - 1) ) # Last 1 == intercept

  end;

function skewness(x::Array{Float64})

  #=
  Calculate skewness in standardized vector x

  :param x: Array. Column vector of length n

  :return: Scalar. Skewness statistic
  =#

  return((sum(x.^3) / size(x)[1]) / ((sum(x.^2) / size(x)[1])^1.5))

  end;

function center(x::Array{Float64})

  #=
  Convenience function used to center an Array by its mean

  :param x: Array. Column vector of length n

  :return: Array with same dimensions of X, such that each column has been grand-mean centered
  =#

  return(x .- (sum(x) / size(x)[1]))

  end;

function cor(x::Array{Float64}, y::Array{Float64})

  #=
  Calculate the pearson correlation coefficient

  :param x: Array. Column vector of length n
  :param y: Array. Column vector of length n

  :return: estimate of the pearson correlation coefficient
  =#

  return( (sum(center(x) .* center(y)) ) / ( sqrt(sum(center(x).^2)) * sqrt(sum(center(y).^2)) ) )

  end;

function independence(x::Array{Float64}, k::Int)

  #=
  Check correlation of x against x lagged by one

  :param x: Array. Column vector of length n
  :param k: Scalar. Indicates the lag value. (i.e. lag of 1 lags x by one)

  :return: correlation of x with x lagged by k
  =#

  return(cor(x[k+1:end,:], x[1:end-k,:]))

  end;

#=
PART III: Model Fit
=#

# Convenience function for density from pdf (dnorm)
function dnorm(x::Float64, mu::Float64, sigma::Float64)

  #=
  Compute the probability density function of x for a normal with mean mu and sd sigma

  :param x: Scalar. value for which you want the pdf
  :param mu: mean of the normal distribution.
  :param sigma: standard deviation of the normal distribution.

  :return: pdf of x
  =#

  pdf(Normal(mu, sigma), x)

  end;

# Convenience function for LL
function logLikelihood(y::Array{Float64}, pred_y::Array{Float64}, sd::Float64)

    #=
    Compute the log likelihood for an outcome given the predicted value

    :param y: Array. Column vector with actual outcome values
    :param pred_y : Array. Column vector of predicted outcome values.
    :sd: standard deviation of the posterior distribution at sample i

    :return: Scalar. Sum of log likelihood for each value y.
    =#

  # Results
  R = Array{Float64}(undef, size(y)[1])

  # For each i in length y/pred_y
  for i in 1:size(y)[1]

    R[i] = dnorm(y[i], pred_y[i], sd)

    end;

  # Sum of log and return
  return(sum(log.(R)))

  end;

# Compute DIC
function DIC(X::Array{Float64}, y::Array{Float64}, w, sigma::Float64, posterior::Array{Float64})

  #=
  Compute the Deviance Information Criterion (DIC)

  :param X: Array. Design matrix X
  :param y: Array. Column vector with outcome variables
  :param sigma: Scalar. MAP value of the posterior standard deviation (residual standard error)
  :param posterior: posterior samples drawn using MCMC.

  :return: Dict containing:
    - DIC: the DIC value
    - Eff. P: estimate of the effective number of parameters
    - LL : summed log likelihood
  =#

  # Check if w is an Array as required by this function.
  # This is an issue of conversion between julia & R
  if isa(w,AbstractArray) == false

      w = [w]

  end;

  # Linear combinations
  predy = X * w

  # Log likelihood of MAP values
  LL = logLikelihood(y, predy, sigma)

  # Linear combinations for each sample
  lincom = X * transpose(posterior[:,1:end-1])

  # Results matrix for P(..)
  P = Array{Float64}(undef, size(lincom)[1])

  # For each sample
  for i in 1:size(lincom)[1]

    P[i] = logLikelihood(y, lincom[:,i], posterior[i,end])

    end;

  # Return
  return(
    Dict(
      "DIC" => -2*LL + 2*(LL-mean(P)),
      "Eff. P" => 2*(LL-mean(P)),
      "LL" => LL
    )
  )

  end;

#=
PART IV: R-squared
=#

function bayes_R2(X::Array{Float64}, y::Array{Float64}, sampled::Array{Float64})

    #=
    Calculate Bayesian R-squared value

    This function calculates the Bayesian R-squared value. For more information, see the article below.

    :param X: Design matrix containing j+1 coefficients (first coefficient must equal intercept variable). The design matrix is created in R using model.matrix()
    :param y: Numeric array with outcome variables. Length(y) == nrow(X)
    :param sampled: posterior samples for the parameters

    :return: Array. Column vector containing r-squared value for each of the samples.

    :seealso:
        - Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2018). R-squared for Bayesian regression models. The American Statistician, (just-accepted), 1-6.
    =#

    # 1. Precompute linear combinations as a function of X and the posterior set of params theta = [b0,b1,...,bn]
    lincom = X * transpose(sampled[: , 1:end-1]) # lincom dims: rows(X) * (iterations-burn)

    # 2. Set up a results matrix for simulated yhat values
    res = Array{Float64}(undef, size(sampled)[1], 1)

    # 3. Compute R-squared
    for j in 1:(size(sampled)[1])

        # R2 = var(pred_y) / (var(pred_y) + sigma^2)
        res[j,1] = var(lincom[:,j]) / (var(lincom[:,j]) + sampled[j, end]^2)

        end;

    # 4. Return
    return(res)

    end;

#=
PART V: Compute outliers
=#

function compute_ppd(X::Array{Float64}, y::Array{Float64}, sampled::Array{Float64})

    #=
    Calculate the proportion of cases where simulated y values exceed the observed y values

    :param X: Design matrix containing j+1 coefficients (first coefficient must equal intercept variable). The design matrix is created in R using model.matrix()
    :param y: Numeric array with outcome variables. Length(y) == nrow(X)
    :param sampled: posterior samples for the parameters

    :return: Array. Column vector proportions where y_simulated > y_observed

    :seealso:
        - Lynch, S. M. (2007). Introduction to applied Bayesian statistics and estimation for social scientists. Springer Science & Business Media. pp.178-182
    =#

    # 1. Precompute linear combinations as a function of X and the posterior set of params theta = [b0,b1,...,bn]
    lincom = X * transpose(sampled[: , 1:end-1]) # lincom dims: rows(X) * (iterations-burn)

    # 2. Set up a results matrix for simulated yhat values
    res = Array{Float64}(undef, size(sampled)[1], size(X)[1])
    # Prop
    psa = Array{Float64}(undef, size(y)[1], 1)

    # 3. Draw from a normal distribution given each specific y value and populate results matrix
    for j in 1:(size(sampled)[1])

        # Simulate y
        res[j,:] = simulate_y(lincom[:,j], sampled[j,end])

        end;

    # 4. Calculate proportion sim > actual for each individual y-outcome
    for i in 1:(size(psa)[1])

        psa[i,1] = (sum(res[:,i] .> y[i,1]) / size(res)[1])

        end;

    # 4. Return proportions
    return(psa)

    end;
