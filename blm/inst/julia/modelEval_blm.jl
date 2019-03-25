# Model evaluation

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
