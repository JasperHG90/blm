# Baby julia function

function gibbs_one_step(b0=0)

  return(b0)

  end;
  
function gibbs_sampling(b0=1, iterations=10, chains=1)

  res = []

  for i = 1:chains
  
    for j = 1:iterations
    
      push!(res, gibbs_one_step(b0))
    
      end;
  
    end;

  return(res)

  end;