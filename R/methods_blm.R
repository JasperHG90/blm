## Methods for blm class

# GENERIC METHODS ------

# Print method
#' @export
print.blm <- function(x) {

  # Print information
  msg <- paste0(
    crayon::bold("Bayesian Linear Model (BLM) object:\n\n"),
    "Data:\n",
    "\tPredictors: ", x$input$m - 1, "\n",
    "\tOutcome: ", x$input$variables$DV, "\n",
    "\tCentered: ", x$input$center, "\n\n"
  )

  # Cat to console
  cat(msg)
  print(get_value(bfit, "sampling_settings"))
  print(get_value(bfit, "priors"))

}

# CONVERGENCE -----

# EVALUATION ----
