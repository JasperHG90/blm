## Methods for sampler



#' Print method for sampler class
#'
print.sampler <- function(x) {

  # Paste message

  # Note: currently all iterations / burn / thinning settings are the SAME for all chains
  #  This makes life easier for the time being.

  msg <- paste0(
    "Sampler:\n",
    "\tChains: ", length(x) , "\n",
    "\tIterations: ", get_value(x, "chain_1")$iterations , "\n",
    "\tThinning: ", get_value(x, "chain_1")$thinning, "\n",
    "\tBurn: ", get_value(x, "chain_1")$burn , "\n\n"
  )

  # Cat
  cat(msg)

}
