## Gibbs sampler -------

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
    
    # Number of params
    nparam = length(priors)
    
    # Number of examples
    n = size(X)[1]
    
    # For each, do
    for j in 1:nparam
        
        # If j == nparam, then sigma
        if j == nparam
            
            sigma = 1 / rand(Gamma(posterior_rate(n, priors[j][:alpha]),
                                   posterior_scale(X, w, y, priors[j][:beta])))
        
        # Else, coefficient    
        else
            
            w[j] = rand(Normal(posterior_mu(X[:,1:end .!= j], X[:,j], y, w[1:end .!= j], sigma, 
                                            priors[j][:mu], priors[j][:sd]),
                               posterior_tau(X[:,j], sigma, priors[j][:sd])))
            
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
    
function gibbs_sampler(X::Array{Float64}, y::Array{Float64}, w::Array{Float64}, 
                       sigma::Float64, iterations::Int, priors, burn::Int)
    
    # Open up a results matrix
    res = Array{Float64}(undef, iterations, (size(X)[2] + 1))
    
    # For each iterations
    for i in 1:iterations
        
        # One stepp of the gibbs sampler
        r = gibbs_one_step(X, y, w, sigma, priors)
        
        # Set new values for weights and sigma
        w = r["w"]
        sigma = r["sigma"]
        
        # Add to results
        res[i,1:end-1] = w
        res[i,size(res)[2]] = sigma
        
        end;
    
    # Return results
    return(res)
    
    end;

## Posterior densities --------

## Posterior densities
function posterior_mu(X::Array{Float64}, xj::Array{Float64}, y::Array{Float64}, 
                      w::Array{Float64}, sigma::Float64, mu0, tau0)
    
    # Numerator
    return(((sum(xj .* (y - X * w)) / sigma) + (mu0/tau0)) / ((transpose(xj) * xj / sigma) + (1 / tau0)))
    
    end;

function posterior_tau(xj::Array{Float64}, sigma::Float64, tau0)
    
    return( 1 / ( (transpose(xj) * (xj ./ sigma)) + (1/tau0)) )
    
    end;

function posterior_rate(n::Int, a0::Float64)
    
    return( (n/2) + a0 )
    
    end;

function posterior_scale(X::Array{Float64}, w::Array{Float64}, y::Array{Float64}, b0::Float64)
    
    return( ( sum( ( y - (X * w) ).^2 ) / 2 ) + b0 )
    
    end;