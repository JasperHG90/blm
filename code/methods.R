# S3 Methods

# Class 'blm' methods -----

# Class method to change priors for available parameters 
set_priors.blm <- function(blm, ...) {
  
  # Get opts
  priors_new <- list(...)
  
  # Get current priors
  priors <- blm$priors
  
  # Empty vector to keep track of warnings
  warning_vars <- c()
  
  # Set new priors (based on names)
  for(i in names(priors_new)) {
    
    if(i %in% names(priors)) {
      
      # If varname in priors, then add to priors_new
      if("varname" %in% names(priors[[i]])) {
        priors_new[[i]]$varname <- priors[[i]]$varname
      }
      
      # Update priors
      priors[[i]] <- priors_new[[i]]
      
    } else {
      
      # Keep track of warnings
      warning_vars <- c(warning_vars, i)
      
    }
    
  }
  
  # Emit warnings
  #  ...
  
  # Set new priors and return
  blm$priors <- priors
  
  # Return
  return(blm)
  
}

# Do gibbs sampling on blm
# This should be a class method
# TODO: add thinning
sampling_options.blm <- function(blm, chains = 1, iterations = 10000, 
                                 burn = 1000) {
  
  ##### Checks ######
  
  if (chains < 1) {
    stop("'chains' cannot be less than 1")
  }
  
  if (iterations < burn) {
    stop("blm cannot sample fewer iterations than burn-in period")
  }
  
  if (ceiling(iterations / burn) < 5) {
    warning("The number of iterations is very low. This will likely yield unstable results.")
  }
  
  if (burn < 1000) {
    warning("You have specified the burn-in samples to be less than 1.000. This will likely yield unstable results.")
  }
  
  ##### End checks #####
  
  # Update settings
  if(!missing(chains)) {
    # Update # chains
    blm$sampling_settings$chains <- chains
    # Update parallel processing
    blm$sampling_settings$parallel = ifelse(chains > 1, "multisession", "sequential")
  } 
  if (!missing(iterations)) {
    blm$sampling_settings$iterations <- iterations
  } 
  if (!missing(burn)) {
    blm$sampling_settings$burn <- burn
  }
  
  # Return
  return(blm)
  
}

# Execute a blm plan
sample_posterior.blm <- function(blm) {

  # Set processing strategy (sequential or multiprocessing)
  strategy <- paste0("future::", blm$sampling_settings$parallel)
  plan(strategy)
  
  # number iterations & chains
  chains <- blm$sampling_settings$chains
  iterations <- blm$sampling_settings$iterations
  priors <- blm$priors
  burn <- blm$sampling_settings$burn
  
  # unroll data
  X <- blm$input$X
  y <- blm$input$y
  
  # Draw initial values for each chain
  initial_values <- lapply(1:chains, function(x) initialize_chain_values(priors))
  
  # Run foreach (possibly parallel for each chain)
  posterior <- foreach(k = 1:chains, .export = c("initial_values", "iterations", "priors", "burn", "X", "y",
                                                 "gibbs_sampler", "gibbs_one_iteration", "posterior_sigma",
                                                 "posterior_coef", "posterior_mu", "posterior_tau")) %dopar% {
  #posterior <- list()
  #for(k in seq_along(1:chains)) {
    
    # Initial values for the current chain
    initial_values_current <- initial_values[[k]]
    
    # Call the gibbs sampler
    gibbs_sampler(X, y, initial_values_current, iterations, priors, burn) #posterior[[k]] <- 
    
  }
  
  # Name output list
  names(posterior) <- paste0("chain_", 1:chains)
  
  # Name output data columns
  posterior <- lapply(posterior, function(x) {
    
    # Retrieve varname if listed in prior
    varnames <- unname(unlist(lapply(priors, function(x) x[["varname"]])))
    
    # Construct names 
    cnames <- c("Intercept", varnames, "sigma")
    
    # Set colnames
    colnames(x) <- cnames

    # Return
    return(x)
        
  })
  
  # Append posterior
  blm$posterior <- posterior
  
  # Return results
  return(blm)
  
}

# Print method
print.blm <- function(blm) {
  
  # Print information
  msg <- paste0(
    "Bayesian Linear Model (BLM) object:\n\n",
    "Data:\n",
      "\tPredictors: ", blm$input$m - 1, "\n",
      "\tOutcome: ", blm$input$variables$DV, "\n",
      "\tCentered: ", blm$input$center, "\n\n",
    "Sampler:\n",
      "\tChains: ", blm$sampling_settings$chains , "\n",
      "\tIterations: ", blm$sampling_settings$iterations , "\n",
      "\tBurn: ", blm$sampling_settings$burn , "\n\n",
    "Priors:\n")
  
  # Set up priors matrix
  pr_coef <- matrix(0L, ncol=length(blm$priors)-1, nrow = 2,
               dimnames=list(
                 c("mu", "tau"),
                 names(blm$priors)[-length(blm$priors)]
               ))
  pr_sigma <- matrix(0L, ncol=1, nrow=2,
                     dimnames = list(
                       c("rate", "scale"),
                       c("sigma")
                     ))
  
  # Add values
  for(i in seq_along(blm$priors)) {
    if(names(blm$priors)[i] == "sigma") {
      pr_sigma[1,1] <- blm$priors[[i]]$alpha
      pr_sigma[2,1] <- blm$priors[[i]]$beta
    } else {
      pr_coef[1,i] <- blm$priors[[i]]$mu
      pr_coef[2,i] <- blm$priors[[i]]$sd
    }
     
  }
  
  cat(msg)
  cat("\t")
  print.listof(list("Coefficients" = pr_coef))
  cat("\t")
  print.listof(list("Variance" = pr_sigma))

}

# Summary method
summary.blm <- function(blm) {
  
  var_names <- c(colnames(blm$input$X), "sigma")
  
  ### Formula
  form <- as.character(blm$input$formula)
  formchar <- paste0(form[2], " ", form[1], " ", form[3])
  
  ## Construct individual parts
  general <- paste0(
    "Formula: '", formchar, "'"
  )
  
  ## User has already sampled or not
  has_sampled <- paste0(
    "Sampled: ", "posterior" %in% names(blm)
  )
  
  ## num. observations + predictors
  obs <- list(
    "Sampling settings" = matrix(c(blm$input$n,blm$input$m - 1,blm$sampling_settings$chains,
                  blm$sampling_settings$iterations,blm$sampling_settings$burn),
                nrow = 1,
                dimnames = list(
                  c(""),
                  c("Obs.",
                    "Predictors",
                    "Chains",
                    "Iterations",
                    "Burn")
                )))
  
  ### If not sampled yet ...
  if(!"posterior" %in% names(blm)) {
    cat("Bayesian Linear Model (BLM) results:")
    cat("\n\n")
    cat(general)
    cat("\n\n")
    cat(has_sampled)
    cat("\n\n")
    print.listof(obs)
    return(cat(""))
  }
  
  ### Statistics
  
  # Calculate MAP for each chain & combine
  MAPV <- round(do.call(rbind.data.frame, MAP(blm$posterior)), digits = 3)
  # Add MC error
  MAPV[3,] <- MAPV[2,] / sqrt(blm$sampling_settings$iterations)
  
  # Amend names
  row.names(MAPV) <- c("Est.", "SE", "MCERR.")
  colnames(MAPV) <-var_names
  
  # Calculate CI
  CIV <- list(
    "95% credible interval" = round(CI(blm$posterior), digits = 3)
  )
  
  # Print MAP & SE
  cat("Bayesian Linear Model (BLM) results:")
  cat("\n\n")
  cat(general)
  cat("\n\n")
  cat(has_sampled)
  cat("\n\n")
  print.listof(obs)
  print.listof(list("Maximum a posteriori (MAP) estimates" = MAPV))
  print.listof(CIV)
  
}

# Evaluate method
#  (plots, statistics etc.)
#  use functions from other packages to do this

# Class 'Prior' methods -----

# Drawing from a normal distribution
draw_value.normal <- function(prior) {
  
  rnorm(1, prior$mu, prior$sd)
  
}

# Drawing from a gamma distribution
draw_value.gamma <- function(prior) {
  
  rgamma(1, prior$alpha, prior$beta) + 1e-10
  
}
