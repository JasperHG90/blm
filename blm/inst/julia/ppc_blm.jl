# Julia functions for posterior predictive checks ==> simulation of y samples

function ppc_draws(X::Array{Float64}, y::Array{Float64}, w::Array{Float64},
                       sigma::Float64, iterations::Int, burn::Int, priors)

    #=
    This function is used to draw data from the posterior using the gibbs sampler and to simulate values of the
    outcome variable y to compute posterior predictive checks.

    :param effective_iterations: number of iterations minus burn-in period
    =#

    # 1. Call the gibbs sampler
    sampled = gibbs_sampler(X, y, w, sigma, iterations, 1, priors)

    # 2. Burn
    sampled = sampled[burn+1:end, :]

    # 3. Precompute linear combinations as a function of X and the posterior set of params theta = [b0,b1,...,bn]
    lincom = X * transpose(sampled[: , 1:end-1]) # lincom dims: rows(X) * (iterations-burn)

    # 4. Set up a results matrix for simulated yhat values
    res = Array{Float64}(undef, iterations - burn, size(X)[1])

    # 5. Set up a results matrix for residuals
    resids = Array{Float64}(undef, iterations - burn, size(X)[1], 2)

    # 5. Draw from a normal distribution given each specific y value and populate results matrix
    for j in 1:(iterations - burn)

        # Simulate y
        res[j,:] = simulate_y(lincom[:,j], sampled[j,end])

        # Compute residuals (simulated)
        resids[j, :, 1] = res[j, :] .- lincom[:, j]

        # Compute residuals (observed)
        resids[j, :, 2] = y .- lincom[:, j]

        end;

    # 6. Return
    return(Dict(
      "sim_y" => res,
      "residuals" => resids
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
