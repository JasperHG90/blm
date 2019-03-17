# Settings

# Density params
dparams <- list(
  "normal" = c("mu", "sd"),
  "beta" = c("alpha", "beta"),
  "gamma" = c("alpha", "beta")
)

# Legal param values
dparam_values <- list(
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
)