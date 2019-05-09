## Methods for hypothesis / hypotheses classes

# Print method for class hypothesis
#' @export
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
