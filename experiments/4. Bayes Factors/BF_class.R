## Hypothesis framework

## Allowed:
# 1. one, specific hypothesis
#    a > b
#    a < b
#    a = b
#    (a-b) < 1
#    2*a > 0
#    ...
# 2. Multiple hypotheses chained with '&' (AND)
#    a = 0 & b = 0
#    a < b & b < c
#    a < .1 & a > .1 ==> same as |a| > .1
#    (a-b) < .1 & (a-b) > .1 ==> same as |a-b| > .1
#    ...

library(stringr)
library(magrittr)

# Single hypothesis -----

# Operators that are allowed
operators_allowed <- c("=", ">", "<", "\\&", "\\|")

# Pass a hypothesis and parse it. THen return object of class 'hypothesis'
hypothesis <- function(hypothesis_user, parameters) {

  # Parse hypothesis
  hypothesis_parsed <- parse_hypothesis(hypothesis_user)

  # Check if variable names in data
  check_hypothesis(hypothesis_parsed, parameters)

  # Results
  hyp <- list(
    "hypothesis" = hypothesis_user,
    "parsed" = hypothesis_parsed
  )

  # Add structure
  class(hyp) <- "hypothesis"

  # Return parsed hypothesis
  return(hyp)

}

# Parse hypothesis
parse_hypothesis <- function(hypothesis_this) {

  # Check if operations allowed
  detect_operator <- str_detect(hypothesis_this, operators_allowed)

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
    hyps <- str_split(hypothesis_this, "\\&") %>%
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
      pop_current <- str_extract(x, pops) %>%
        na.omit(.) %>%
        as.character()

      # Split string
      sp <- str_split(x, pop_current)[[1]] %>%
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
          "abs_values" = str_detect(left, "\\|"),
          "algebra" = str_detect(left, "\\*|\\-|\\+|\\/"),
          "is_numeric" = nums[1]
        ) ,
        "right" = list(
          "expression" = right,
          "abs_values" = str_detect(right, "\\|"),
          "algebra" = str_detect(right, "\\*|\\-|\\+|\\/"),
          "is_numeric" = nums[2]
      ))

      # Parse parameter names
      res$left$params <- parse_parameters(res$left$expression, res$left$abs_values, res$left$algebra,
                                          res$left$is_numeric)
      res$right$params <- parse_parameters(res$right$expression, res$right$abs_values, res$right$algebra,
                                           res$right$is_numeric)

      # Add algebra operators
      if(res$left$algebra) {
        res$left$algebra_operator <- str_extract(res$left$expression, "\\*|\\-|\\+|\\/")
      }
      if(res$right$algebra) {
        res$right$algebra_operator <- str_extract(res$right$expression, "\\*|\\-|\\+|\\/")
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
parse_parameters <- function(exp, abs_values, algebra, is_numeric) {
  # Assign r
  r <- exp
  # Check abs values
  if(abs_values) {
    # Strip
    r <- str_replace_all(r, "\\|", "")
  }
  if(algebra) {
    r <- str_split(r, "\\*|\\-|\\+|\\/")[[1]] %>%
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

# Methods

# Print method for class hypothesis
print.hypothesis <- function(x) {

  # Hypothesis
  usr_hyp <- x$hypothesis

  # Number of elements
  elements <- length(x$parsed)

  # cat
  msg <- paste0(
    crayon::bold("Bayesian Linear Model (BLM) hypothesis:\n\n"),
    "Hypothesis: ", usr_hyp, "\n",
    "Elements: ", elements
  )
  cat(msg)

}


# Hypotheses class ----

#' Add a hypothesis to a blm model
#'
#' To evaluate informative hypotheses, users may add these using this function. The following kind of hypotheses are allowed:
#' 1. one, specific hypothesis, e.g
#'    a > b
#'    a < b
#'    a = b
#'    a-b < 1
#'    2*a > 0
#'    ...
#' 2. Multiple hypotheses chained with '&' (AND), e.g,
#'    a = 0 & b = 0
#'    a < b & b < c
#'    a < .1 & a > .1 ==> same as |a| > .1
#'    a-b < .1 & a-b > .1 ==> same as |a-b| > .1
#'    ...
#'
#' @param x a blm object
#' @param hypothesis_user a hypothesis.
#'
#' @return a blm object with a new or updated 'hypotheses' object
#'
#' @seealso cite Hoitink
#' @export
set_hypothesis <- function(x, hypothesis_user) {
  UseMethod("set_hypothesis", x)
}

# Add a hypothesis to a blm object
set_hypothesis.blm <- function(x, hypothesis_user) {

  # Hypothesis in blm object?
  if(!blm::contains(x, "hypotheses")) {
    x <- blm:::set_value.blm(x, "hypotheses", list())
    # Add structure
    class(x$hypotheses) <- "hypotheses"
  }

  # Get parameters
  pars <- names(x$priors)
  pars <- pars[pars != "sigma"]

  # Does the hypothesis already exist?
  current_hyps <- lapply(x$hypotheses, function(x) x$hypothesis) %>%
    unlist()
  # Simply return
  if(hypothesis_user %in% current_hyps) {
    return(x)
  }

  # Add hypothesis
  x$hypotheses <- append(x$hypotheses, list(hypothesis(hypothesis_user, pars)))

  # Set names
  names(x$hypotheses) <- paste0("hypothesis_", 1:length(x$hypotheses))

  # Return
  return(x)

}

# Add hypothesis to blm object
library(blm)
library(dplyr)
# Load data
data("directors")
# Adjust
directors <- directors %>%
  mutate(Age = Age - mean(Age),
         Compensation = log(Compensation))
# Bridge to julia
blm_setup()
# Iterations
k <- 20000
# Fit the model
fit <- blm("Compensation ~ Age + Male", data=directors) %>%
  # Specify MCMC options
  set_sampling_options(iterations=k, chains=2, thinning=5) %>%
  # Set reasonable priors
  set_prior("b1", mu=0, sd=.01) %>%
  set_prior("b2", mu=.05, sd=.03) %>%
  # Set hypotheses
  set_hypothesis("b0 > 0") %>%
  set_hypothesis("b0 + b2 > b0") %>%
  set_hypothesis("|b1| > .05 & b2 > .1")
  # Sample posterior distribution
  sample_posterior()

