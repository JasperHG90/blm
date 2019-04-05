## Shared methods

#' @export
get_value.sampler <-
  get_value.priors <-
  get_value.prior <-
  get_value.posterior <-
  get_value.blm <-
  get_value.chain <- function(x, var) {

    ## Get this value
    return(x[[var]])

  }

#' @export
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
