# Helper functions

# Exported functions ----

#' Check if a list / object contains a value
#'
#' @param x object or list
#' @param val character value of the name of the object/list element to check
#'
#' @examples
#' x <- list("a"=1, "b"=2)
#' contains(x, "a")
#'
#' @return TRUE if the value is present, FALSE if not.
#' @export
contains <- function(x, val) {

  return(val %in% names(x))

}

# Not exported -----

# Functions used to parse/check hypotheses -----

# Scale a vector of coefficients to unit measurement
# This creates scaled beta coefficients
scale_coef <- function(b, y_sd) {
  b * sd(b) / y_sd
}

# Parse hypothesis
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom stringr str_extract
parse_hypothesis <- function(hypothesis_this) {

  # Get from the global environment
  operators_allowed <- getOption("blm_hypothesis_operators")$operators

  # Check if operations allowed
  detect_operator <- stringr::str_detect(hypothesis_this, operators_allowed)

  # If all false, then error
  if(all(detect_operator == FALSE)) {
    stop(paste0("No allowed operators found in hypothesis. Allowed operators are '",
                paste0(operators_allowed, collapse=", "), "'"))
  }

  # Which operators found?
  # Split into: primary operators (>, <, =)
  # Group operators (|, &)
  pops <- operators_allowed[1:3][detect_operator[1:3]]
  gops <- operators_allowed[4:5][detect_operator[4:5]]

  # If multiple primary operators found but not & then ==> error
  if((length(gops) == 0) & (length(pops) > 1)) {
    stop("Multiple primary operators ('=', '<', '>') passed without a group operator ('&')")
  }

  # Check for '&' (multiple hypotheses)
  if(any(gops == "\\&")) {
    hyps <- stringr::str_split(hypothesis_this, "\\&") %>%
      # Trim whitespace
      lapply(trimws) %>%
      # To list
      unlist() %>%
      as.list()
  } else {
    hyps <- list(hypothesis_this)
  }

  # For each operator, split at left and right side
  hyps_split <- hyps %>%
    lapply(function(x) {

      # Detect operator
      pop_current <- stringr::str_extract(x, pops) %>%
        na.omit(.) %>%
        as.character()

      # Split string
      sp <- stringr::str_split(x, pop_current)[[1]] %>%
        trimws()

      # Check length
      if(length(sp) > 2) {
        stop(paste0("Too many elements in hypothesis: '", paste0(pop_current), "'"))
      }

      # Any numeric?
      nums <- c()
      for(i in sp) nums <- c(nums, suppressWarnings(!is.na(as.numeric(i))))
      # Left and right parts
      left <- ifelse(nums[1], as.numeric(sp[1]), sp[1])
      right <- ifelse(nums[2], as.numeric(sp[2]), sp[2])

      # Result
      res <- list(
        "operator" = pop_current,
        "left" = list(
          "expression" = left,
          "abs_values" = stringr::str_detect(left, "\\|"),
          "algebra" = stringr::str_detect(left, "\\*|\\-|\\+|\\/"),
          "is_numeric" = nums[1]
        ) ,
        "right" = list(
          "expression" = right,
          "abs_values" = stringr::str_detect(right, "\\|"),
          "algebra" = stringr::str_detect(right, "\\*|\\-|\\+|\\/"),
          "is_numeric" = nums[2]
        ))

      # Parse parameter names
      res$left$params <- parse_parameters(res$left$expression, res$left$abs_values, res$left$algebra,
                                          res$left$is_numeric)
      res$right$params <- parse_parameters(res$right$expression, res$right$abs_values, res$right$algebra,
                                           res$right$is_numeric)

      # If absolute values, replace | | for R command
      if(res$left$abs_values) {
        res$left$expression <- stringr::str_replace_all(res$left$expression, "\\|", "") %>%
          paste0("abs(", ., ")")
      }

      if(res$right$abs_values) {
        res$right$expression <- stringr::str_replace_all(res$right$expression, "\\|", "") %>%
          paste0("abs(", ., ")")
      }

      # Add algebra operators
      if(res$left$algebra) {
        res$left$algebra_operator <- stringr::str_extract(res$left$expression, "\\*|\\-|\\+|\\/")
      }
      if(res$right$algebra) {
        res$right$algebra_operator <- stringr::str_extract(res$right$expression, "\\*|\\-|\\+|\\/")
      }

      # Return
      return(res)

    })

  # Names of hypotheses
  names(hyps_split) <- paste0("group_", 1:length(hyps_split))

  # Add order of evaluation by reversing groups
  hyps_split <- hyps_split[length(hyps_split):1]

  # Return
  return(hyps_split)

}

# Check if all parsed parameters in parameters in data
check_hypothesis <- function(hypothesis_parsed, parameters) {

  # Retrieve parameters
  pars_in_hyp <- lapply(hypothesis_parsed, function(x) list(x$left$params, x$right$params)) %>%
    unlist(recursive = TRUE) %>%
    unique()

  # If pars not in parameters of data, throw error
  if(any(!pars_in_hyp %in% parameters)) {
    stop(paste0("User passed parameters ('",
                paste0(pars_in_hyp[!pars_in_hyp %in% parameters], collapse=", "),
                "') in hypotheses that are not present in data"))
  }

}

# Parse parameters from an expression
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
parse_parameters <- function(exp, abs_values, algebra, is_numeric) {
  # Assign r
  r <- exp
  # Check abs values
  if(abs_values) {
    # Strip
    r <- stringr::str_replace_all(r, "\\|", "")
  }
  if(algebra) {
    r <- stringr::str_split(r, "\\*|\\-|\\+|\\/")[[1]] %>%
      trimws()
  }
  if(is_numeric) {
    r <- c()
  }

  # Is any a number?
  # If empty
  nnums <- c()
  for(i in seq_along(r)){
    if(suppressWarnings(is.na(as.numeric(r[i])))) {
      nnums <- c(nnums, i)
    }
  }
  r <- r[nnums]

  return(r)
}

# Compute the complexity of a hypothesis
# @param hypothesis_i the current hypothesis under evaluation
# @param priors object of class priors containing prior information
# @param y_sd standard deviation of the outcome variable y
# CITE HOITINK
compute_hypothesis_complexity <- function(hypothesis_i, priors, y_sd) {

  # Get parsed version of the hypothesis
  hpar <- hypothesis_i$parsed

  # The order has already been reversed if multiple hypotheses!

  # Get an overview of the parameters mentioned
  params_in_hyp <- lapply(hpar, function(x) c(x$left$params, x$right$params)) %>%
    unlist(recursive = TRUE) %>%
    unique()

  # Draw samples from the priors
  prior_samples <- lapply(params_in_hyp, function(x) {

    # Prior information
    prinf <- priors[[x]]

    # Draw
    rnorm(50000, prinf$mu, prinf$sd) %>%
      scale_coef(y_sd)

  })

  # Names
  names(prior_samples) <- params_in_hyp

  # Assign to environment
  for(i in seq_along(prior_samples)) {

    # Assign value to local scope. i.e. this has the effect of b0 <- value
    assign(names(prior_samples)[i], prior_samples[[i]])

  }

  # Evaluate the parameters
  evaluated_hyps <- rep(0, 50000)

  for(i in seq_along(hpar)) {

    # Current
    cur <- hpar[[i]]

    # Evaluate the expression
    evaluated_hyps <- evaluated_hyps + eval(parse(text=paste(cur$left$expression, cur$operator, cur$right$expression)))

  }

  # Only accept evaluated_hyps if it is TRUE in all cases. i.e. value for pos i == # of evaluted hypotheses
  evaluated_hyps <- ifelse(evaluated_hyps == length(hpar), 1, 0)

  # Return proportion (complexity)
  return(mean(evaluated_hyps))

}

# Compute the fit of a hypothesis
# @param hypothesis_i the current hypothesis under evaluation
# @param posterior posterior samples
# @param y_sd standard deviation of the outcome variable y
# CITE HOITINK
compute_hypothesis_fit <- function(hypothesis_i, posterior, y_sd) {

  # Get parsed version of the hypothesis
  hpar <- hypothesis_i$parsed

  # The order has already been reversed if multiple hypotheses!

  # Get an overview of the parameters mentioned
  params_in_hyp <- lapply(hpar, function(x) c(x$left$params, x$right$params)) %>%
    unlist(recursive = TRUE) %>%
    unique()

  # Extract posterior from the posterior samples
  post_samples <- vector("list", length(params_in_hyp))

  # Names
  names(post_samples) <- params_in_hyp

  # Assign to environment
  for(i in seq_along(post_samples)) {

    # Assign value to local scope. i.e. this has the effect of b0 <- value
    assign(params_in_hyp[i], posterior[,params_in_hyp[i]] %>% scale_coef(y_sd))

  }

  # Evaluate the parameters
  evaluated_hyps <- rep(0, nrow(posterior))

  for(i in seq_along(hpar)) {

    # Current
    cur <- hpar[[i]]

    # Evaluate the expression
    evaluated_hyps <- evaluated_hyps + eval(parse(text=paste(cur$left$expression,
                                                             cur$operator,
                                                             cur$right$expression)))

  }

  # Only accept evaluated_hyps if it is TRUE in all cases. i.e. value for pos i == # of evaluted hypotheses
  evaluated_hyps <- ifelse(evaluated_hyps == length(hpar), 1, 0)

  # Return proportion (fit)
  return(mean(evaluated_hyps))

}


# Helper functions for sampling ----

# Initiate initial values for a Gibbs chain
initialize_chain_values <- function(priors) {

  # For each prior, draw initial values and construct a weight matrix
  w <- matrix(0L, nrow = length(priors)-1, ncol=1)

  # To grid
  for(i in seq_along(priors)) {
    if(names(priors)[i] == "sigma") {
      sigma <- draw_value(priors[[i]])
    } else {
      w[i,1] <- draw_value(priors[[i]])
    }
  }

  # Return
  return(
    list(
      "w" = w,
      "sigma" = sigma
    )
  )

}

# Helper function that calls the Julia MC sampler
mc_sampler <- function(X, y, initial_values, iterations, thinning, priors, samplers) {

  # Unroll initial values
  w <- initial_values$w
  # If intercept-only model, then w is a scalar. Coerce to array
  # JuliaCall screws this up in conversion!!!
  if(is.vector(w)) {
    w <- matrix(w, ncol=1, nrow=1)

  }

  sigma <- initial_values$sigma

  # Call MCMC sampler in Julia
  r <- .blm$julia$eval("MCMC_sampler")(X, as.numeric(y), w, sigma, as.integer(iterations),
                                       as.integer(thinning), unname(priors), unname(samplers))

  # Burn
  return(r)

}

# Helper functions for evaluation -----

# Autocorrelation plot
# @param pd posterior samples
# @param chain chain to show autocorrelation plot for
autocorrelation_plot <- function(pd, chain) {

  # Get data from posterior
  qd <- get_value(pd, "samples")[[paste0("chain_", chain)]]

  # Lag data
  qd <- qd %>%
    as.data.frame() %>%
    lapply(., function(x) autocor(x, n=40)) %>%
    do.call(rbind.data.frame, .)

  # Add variable name
  qd$id <- stringr::str_replace_all(row.names(qd), "\\.[0-9]{1,2}", "")

  # Plot
  ggplot2::ggplot(qd, ggplot2::aes(x=lag, y=correlation, group=id)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::scale_x_continuous(breaks= scales::pretty_breaks()) +
    ggplot2::scale_y_continuous(limits = c(-1,1)) +
    theme_blm() +
    ggplot2::facet_wrap(id ~ .) +
    ggplot2::labs(title="Autocorrelation plot")


}

# History plot
history_plot <- function(samples) {

  ggplot2::ggplot(samples, ggplot2::aes(x = iteration,
                                        y=value,
                                        color=as.factor(chain),
                                        group = parameter)) +
    ggplot2::geom_line(alpha=0.4) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "chain",
                                labels = paste0("chain ", 1:length(unique(samples$chain)))) +
    theme_blm() +
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::facet_wrap("parameter ~ .", scales = "free_y",
                        ncol=2) +
    ggplot2::labs(title = "Trace plot", subtitle = "each chain is indicated by its own color")

}

# Density plot
density_plot <- function(samples) {

  ggplot2::ggplot(samples, ggplot2::aes(x=value,
                                        fill = as.factor(chain))) +
    ggplot2::geom_density(alpha=0.4) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "chain",
                                labels = paste0("chain ", 1:length(unique(samples$chain)))) +
    theme_blm() +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   legend.position="none") +
    ggplot2::facet_wrap("parameter ~ .", scales = "free") +
    ggplot2::labs(title = "Posterior densities", subtitle = "each chain is indicated by its own color")

}

# Helper function for autocorrelation
autocor <- function(x, n=10) {

  # Results
  res <- rep(0, n)

  # Lag for each n and calculate correlation
  for(i in 1:n) {
    res[i] <- cor(x, c(rep(NA, i), x[1:(length(x)-i)]),
                  use="complete.obs")
  }

  # Return
  return(
    data.frame(
      "lag" = 1:n,
      "correlation" = res
    )
  )

}

# Calculate mode
calc_mode <- function(x) {
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

# Posterior predictive checks in julia
ppc_julia <- function(X, y, posterior_samples) {

  # Call Julia function for posterior predictive checks
  return(
    .blm$julia$eval("ppc_draws")(X, as.numeric(y), posterior_samples)
  )

}

# R-squared calculation in julia
bayes_R2 <- function(X, y, posterior_samples) {

  # Call Julia function for posterior predictive checks
  return(
    .blm$julia$eval("bayes_R2")(X, as.numeric(y), posterior_samples)
  )

}

# Model DIC
# See: http://kylehardman.com/BlogPosts/View/6
DIC <- function(X, y, posterior_samples) {

  ## Map values
  map <- MAP(posterior_samples)$MAP

  ## Bind data
  pb <- bind(posterior_samples)
  ## To matrix (expected by julia)
  pb <- as.matrix(pb)

  ## Coef separate from sigma
  sigma <- unname(map["sigma"])
  coefs <- matrix(map[-length(map)], ncol=1)

  ## Two parts to DIC ==> (1) sum of log of likelihood P(y|theta_MAP)

  # Call Julia implementation
  # NOTE: sigma here is the standard deviation. We use posterior values. These have been squared-rooted
  #  after sampling.
  r <- .blm$julia$eval("DIC")(X, as.numeric(y),
                              coefs, sigma, pb)

  ## Return model fit
  return(
    do.call(cbind.data.frame,
            r)
  )

}

# Calculate the model BIC
BIC_blm <- function(n, p, LL) {
  log(n)*p - 2*LL
}

# Calculate the Bayes factor of the model against the null model
# See Wagemakers 2007
BF <- function(BIC0, BIC1) {
  exp((BIC0 - BIC1)/2)
}

## Sample size helper
sample_size <- function(n, pk) {
  n / (1 + 2 * sum(pk))
}

# Helper functions for special cat() ----

# Special cat function for burn-in statistic
cat.burnin <- function(x) {

  # Dims
  j <- dim(x)[2]
  n <- dim(x)[1]

  # Fill setting for cat()
  f <- (nchar(crayon::green("0.00")) + 1) * j + (j-1) + nchar("chain 1 ")

  # Round
  x <- round(abs(x), digits=2)

  # Determine colors
  is_zero <- x == 0
  less_than_05 <- abs(x) < 0.05

  # Turn x into a character vector
  x_v <- as.character(round(unname(unlist(x)), digits=2))

  # For those elements of length 1, add decimals
  x_v_l1 <- nchar(x_v) == 1
  x_v_l3 <- nchar(x_v) == 3

  # Set
  x_v[x_v_l1] <- paste0(x_v[x_v_l1], ".00")
  x_v[x_v_l3] <- paste0(x_v[x_v_l3], "0")

  # Paste
  cat(paste0("        ", paste0("b", 0:(j-1), collapse="   "), "\n"))
  cat(ifelse(is_zero, crayon::green(x_v),
             ifelse(less_than_05, crayon::yellow(x_v),
                    crayon::red(x_v))), fill=f, labels=paste0("Chain ", 1:2), sep = " ")
}

# Special cat function for Gelman-Rubin statistic
cat.GR <- function(x) {

  # Round
  x <- round(x, digits=2)

  # To vector
  x <- unname(x[1,])

  # Length
  j <- length(x)

  # Fill setting for cat()
  f <- (nchar(crayon::green("0.00")) + 1) * j + (j-1) + nchar("chain 1 ")

  # Determine colors
  is_one <- x == 1

  # Turn x into a character vector
  x_v <- as.character(x)

  # For those elements of length 1, add decimals
  x_v_l1 <- nchar(x_v) == 1
  x_v_l3 <- nchar(x_v) == 3

  # Set
  x_v[x_v_l1] <- paste0(x_v[x_v_l1], ".00")
  x_v[x_v_l3] <- paste0(x_v[x_v_l3], "0")

  # Paste
  cat(paste0("        ", paste0("b", 0:(j-1), collapse="   "), "\n"))
  cat(ifelse(is_one, crayon::green(x_v), crayon::red(x_v)), fill=f, labels="Value  ", sep = " ")

}

# Misc ----

# See documentation below.
generate_dataset <- function(n = 2000, j = 5, binary = 1, seed=NULL,
                             heteroskedastic = FALSE, correlated_errors = FALSE, ...) {

  opts <- list(...)
  if("degree" %in% names(opts)) {
    degree <- opts$degree
    if(!is.numeric(degree)) {
      stop("'degree' must be numeric")
    } else if(degree <= 0) {
      stop("'degree' cannot be less than or equal to 0")
    } else {
      # Continue
      degree <- degree
    }
  } else{
    degree <- 1
  }

  # Empty matrix
  X <- matrix(0L, ncol=j, nrow=n)
  # Determine which variables will be binary
  # Populate matrix
  if(!is.null(seed)) {
    set.seed(seed)
  }
  binary_j <- order(runif(j))[0:binary]

  for(i in 1:j) {
    if(!i %in% binary_j) {
      # Mean and sd
      m <- runif(1, 5, 10) * sign(runif(1, -1, 1))
      sd <- runif(1, 1, 3)
      X[,i] <- rnorm(n, m, sd)
    } else {
      X[,i] <- rbinom(n, 1, 0.5)
    }
  }

  # Intercept + j coefs
  coef <- c(rnorm(1, 0, 20), rnorm(j, 0, 5))
  sigmaSq <- runif(1, 1, 10)

  #browser()

  # Make the example heteroskedastic
  if(heteroskedastic) {
    sigmaSq <- rep(1e-10, n)
    for(i in setdiff(1:j, binary_j)) {
      X[,i] <- X[,i] + seq(1,mean(X[,i]) + degree*sd(X[,i]), length.out = n)
      # Set sigma
      pw <- runif(1, 1, 1.6)
      # Scale power for degree
      if(degree <= 1) {
        degreepw <- 1
      } else {
        degreepw <- degree
      }
      pw <- (pw + (1-(1/degreepw)))/pw
      sigmaSq <- sigmaSq + (abs(X[,i])^pw)
    }
    sigmaSq

    # Add intercept
    X <- cbind(rep(1, n), X)

    # Generate y
    y <- rnorm(n,
               mean = X %*% matrix(coef),
               sd = sqrt(sigmaSq))

  } else if(correlated_errors) { # Make the example correlated

    # Correlated values
    ord <- arima.sim(list(order = c(1,0,0), ar = degree), n = n, sd=sqrt(sigmaSq))

    # Add intercept
    X <- cbind(rep(1, n), X)

    # Generate y
    y <- rnorm(n,
               mean = X %*% matrix(coef) + ord,
               sd = sqrt(sigmaSq))

  } else {

    # Add intercept
    X <- cbind(rep(1, n), X)

    # Generate y
    y <- rnorm(n,
               mean = X %*% matrix(coef),
               sd = sqrt(sigmaSq))
  }

  # List of real values (for later)
  real <- list(
    coef = coef,
    sigma = sigmaSq
  )

  # Return
  return(
    list(
      "X" = X,
      "y" = y,
      "parameters" = real,
      "seed" = ifelse(is.null(seed), NA, seed)
    )
  )

}

#' Generate a dataset with normal variables and outcome for testing purposes
#'
#' @param n number of examples
#' @param j number of variables
#' @param binary number of columns with 0/1 outcomes
#' @param seed seed for pseudo-random number generator
#'
#' @return list containing true values and data frame with n rows and j columns with simulated data.
#' @export
blmsim <- generate_dataset
