## Gibbs sampler -------

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

    # For each iterations
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
