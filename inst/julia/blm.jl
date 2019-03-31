#=

BLM (Bayesian Linear Model) helper functions in Julia

Written by: Jasper Ginn
Course: Introduction to Bayesian Statistics
          @Utrecht University
          @instructor1: Prof. Herbert Hoijtink
          @instructor2: Fayette Klaassen

=#

#=
PART I: Gibbs sampler in Julia
=#

function gibbs_sampler(X::Array{Float64}, y::Array{Float64}, w::Array{Float64},
                       sigma::Float64, iterations::Int, thinning::Int, priors)

    #=
    Run the Gibbs sampler to obtain the conditional posterior coefficients / sigma values

    :param X:
    :param y:
    :param w:
    :param sigma:
    :param iterations:
    :param priors:

    :return:
    :seealso:
    =#

    # Open up a results matrix
    res = Array{Float64}(undef, iterations, (size(X)[2] + 1))

    # For each iteration
    for i in 1:iterations

        for j in 1:thinning

          # One step of the gibbs sampler
          r = gibbs_one_step(X, y, w, sigma, priors)

          # Set new values for weights and sigma
          w = r["w"]
          sigma = r["sigma"]

          end;

        # Add to results
        res[i,1:end-1] = w
        res[i,end] = sigma

        end;

    # Square-root the sigma values
    res[:,end] = sqrt.(res[:,end])

    # Return results
    return(res)

    end;

function gibbs_one_step(X::Array{Float64}, y::Array{Float64}, w::Array{Float64},
                        sigma::Float64, priors)

    #=
    Run a single iteration of the gibbs sampler

    :param X:
    :param y:
    :param w:
    :param sigma:
    :param priors:

    :return:
    :seealso:
    =#

    # For each, do
    for j in 1:length(priors)

        # If j == nparam, then sigma
        if j == length(priors)

            sigma = rand(InverseGamma(posterior_rate(size(X)[1], priors[j][:alpha]),
                                      posterior_scale(X, w, y, priors[j][:beta])))

        # Else, coefficient
        else

            w[j] = rand(Normal(posterior_mu(X[:,1:end .!= j], X[:,j], y, w[1:end .!= j], sigma,
                                            priors[j][:mu], priors[j][:sd]),
                               sqrt(posterior_tau(X[:,j], sigma, priors[j][:sd]))))

            end;

        end;

    # Return the coefficient matrix
    return(
        Dict(
            "w" => w,
            "sigma" => sigma
        )
    )

    end;

## Posterior densities --------

## Posterior densities
function posterior_mu(X::Array{Float64}, xj::Array{Float64}, y::Array{Float64},
                      w::Array{Float64}, sigma::Float64, mu0::Float64, tau0::Float64)

    # Numerator
    return(((sum(xj .* (y - X * w)) / sigma) + (mu0/tau0)) / ((transpose(xj) * xj / sigma) + (1 / tau0)))

    end;

function posterior_tau(xj::Array{Float64}, sigma::Float64, tau0::Float64)

    return( 1 / ( (transpose(xj) * (xj ./ sigma)) + (1/tau0)) )

    end;

function posterior_rate(n::Int, a0::Float64)

    return( (n/2) + a0 )

    end;

function posterior_scale(X::Array{Float64}, w::Array{Float64}, y::Array{Float64}, b0::Float64)

    return( ( sum( ( y - (X * w) ).^2 ) / 2 ) + b0 )

    end;

#=
PART II: Posterior Predictive Checks
=#

function ppc_draws(X::Array{Float64}, y::Array{Float64}, w::Array{Float64},
                       sigma::Float64, iterations::Int, thinning::Int, burn::Int, priors)

    #=
    This function is used to draw data from the posterior using the gibbs sampler and to simulate values of the
    outcome variable y to compute posterior predictive checks.

    :param effective_iterations: number of iterations minus burn-in period
    =#

    # 1. Call the gibbs sampler
    sampled = gibbs_sampler(X, y, w, sigma, iterations, thinning, priors)

    # 2. Burn
    sampled = sampled[burn+1:end, :]

    # 3. Precompute linear combinations as a function of X and the posterior set of params theta = [b0,b1,...,bn]
    lincom = X * transpose(sampled[: , 1:end-1]) # lincom dims: rows(X) * (iterations-burn)

    # 4. Set up a results matrix for simulated yhat values
    res = Array{Float64}(undef, iterations - burn, size(X)[1])

    # 5. Set up a results matrix for residuals, skewness stats, heteroskedasticity stats and correlation stats
    resids = Array{Float64}(undef, iterations - burn, size(X)[1], 2)
    skewed = Array{Float64}(undef, iterations - burn, 2)
    heterosked = Array{Float64}(undef, iterations - burn, 2)
    correlate = Array{Float64}(undef, iterations - burn, 2)

    # 5. Draw from a normal distribution given each specific y value and populate results matrix
    for j in 1:(iterations - burn)

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

    # 6. Return
    return(Dict(
      "sim_y" => res,
      "residuals" => resids,
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
  =#

  # Linear model, prediction and return
  return(X * ( inv( transpose(X) * X ) * transpose(X) * y))

  end;

function adjR2(X::Array{Float64}, y::Array{Float64})

  #=
  Calculate sums of squares & R2
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
  =#

  return((sum(x.^3) / size(x)[1]) / ((sum(x.^2) / size(x)[1])^1.5))

  end;

function center(x::Array{Float64})

  #=
  Convenience function used to center a vector by its mean
  =#

  return(x .- (sum(x) / size(x)[1]))

  end;

function cor(x::Array{Float64}, y::Array{Float64})

  #=
  Pearson correlation coefficient
  =#

  return( (sum(center(x) .* center(y)) ) / ( sqrt(sum(center(x).^2)) * sqrt(sum(center(y).^2)) ) )

  end;

function independence(x::Array{Float64}, k::Int)

  #=
  Check correlation of x against x lagged by one
  =#

  return(cor(x[k+1:end,:], x[1:end-k,:]))

  end;

#=
PART III: Model Fit
=#

# Convenience function for density from pdf (dnorm)
function dnorm(x, mu, sigma)

  pdf(Normal(mu, sigma), x)

  end;

# Convenience function for LL
function logLikelihood(y, pred_y, sd)

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
function DIC(X, y, w, sigma, posterior)

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
