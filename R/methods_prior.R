## Methods for prior class

#' @export
print.priors <- function(x) {

  # Set up priors matrix
  pr_coef <- matrix(0L, ncol=length(x)-1, nrow = 2,
                    dimnames=list(
                      c("mu", "tau"),
                      names(x)[-length(x)]
                    ))
  pr_sigma <- matrix(0L, ncol=1, nrow=2,
                     dimnames = list(
                       c("rate", "scale"),
                       c("sigma")
                     ))

  # Add values
  for(i in seq_along(x)) {
    if(names(x)[i] == "sigma") {
      pr_sigma[1,1] <- get_value(x, "sigma") %>%
                          get_value(., "alpha")
      pr_sigma[2,1] <- get_value(x, "sigma") %>%
                          get_value(., "beta")
    } else {
      pr_coef[1,i] <- get_value(x, names(x)[i]) %>%
                        get_value(., "mu")
      pr_coef[2,i] <- get_value(x, names(x)[i]) %>%
                        get_value(., "sd")
    }

  }

  # Cat to console
  print.listof(list("Priors (Coefficients)" = pr_coef))
  print.listof(list("Priors (Residuals)" = pr_sigma))

}

# IF UPDATE PRIORS, DRAW NEW INITIAL VALUES !!!!
update_priors.priors <- function(x) {



}

