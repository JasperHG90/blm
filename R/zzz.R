# These settings are set up when the package loads

.onLoad <- function(libname = find.package("blm"), pkgname="blm") {

  # Density params
  options(
    blm_dparams = list(
      "normal" = c("mu", "sd"),
      "beta" = c("alpha", "beta"),
      "gamma" = c("alpha", "beta")
    )
  )

  # Legal param values
  options(
    blm_allowed_dparam_values = list(
      "normal" = list(
        "mu" = c(-Inf, Inf),
        "sd" = c(0, Inf)
      ),
      "beta" = list(
        "alpha" = c(0, Inf),
        "beta" = c(0, Inf)
      ),
      "gamma" = list(
        "alpha" = c(0, Inf),
        "beta" = c(0, Inf)
      )
  ))

  # Operators allowed when setting hypotheses
  options(
    blm_hypothesis_operators = list(
      "operators" = c("=", ">", "<", "\\&", "\\|")
    )
  )

}

