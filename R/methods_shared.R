## Shared methods

#' @export
get_value.blm <- function(x, var) {

  ## Get this value
  return(x[[var]])

}

# Get a value from a blm or related object (internal objects)
get_value.sampler <-
  get_value.priors <-
  get_value.prior <-
  get_value.posterior <-
  get_value.chain <-
  get_value.R2 <-
  get_value.ppc <-
  get_value.DIC <- function(x, var) {

    ## Get this value
    return(x[[var]])

  }

# Set a value to a new value
set_value.sampler <-
  set_value.priors <-
  set_value.prior <-
  set_value.posterior <-
  set_value.blm <-
  set_value.chain <- function(x, var, val) {

    ## Set this value
    x[[var]] <- val

    ## Return
    return(x)

  }
