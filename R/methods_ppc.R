## Methods for posterior predictive checks

# Check normality of errors assumption
normality_check.ppc <- function(ppc, skewed) {

  # Compute where sim > obs (bayesian p-value)
  bpv <- mean(skewed[,1] > skewed[,2])

  # To list
  ppc$data$skewness <- skewed

  # Add results
  if(!"results" %in% names(ppc)) {
    ppc$results <- list(
      "normality" = bpv
    )
  } else {
    ppc$results$normality <- bpv
  }

  # Return
  return(ppc)

}

# Check homoskedasticity assumption
homoskedast_check.ppc <- function(ppc, heterosked) {

  # Legacy
  homosked <- heterosked

  # Compute where sim > obs (bayesian p-value)
  bpv <- mean(homosked[,1] > homosked[,2])

  # To list
  ppc$data$homosked <- homosked

  # Add results
  if(!"results" %in% names(ppc)) {
    ppc$results <- list(
      "homosked" = bpv
    )
  } else {
    ppc$results$homosked <- bpv
  }

  # Return
  return(ppc)

}

# Check independence of errors assumption
independence_check.ppc <- function(ppc, independence) {

  # Compute where sim > obs (bayesian p-value)
  bpv <- mean(independence[,1] > independence[,2])

  # To list
  ppc$data$independence <- independence

  # Add results
  if(!"results" %in% names(ppc)) {
    ppc$results <- list(
      "independence" = bpv
    )
  } else {
    ppc$results$independence <- bpv
  }

  # Return
  return(ppc)
}

# Summary
#' @export
summary.ppc <- function(ppc) {

  msg <- paste0(
    crayon::bold("Posterior Predictive Checks (PPC) for blm object:"),
    "\n\n"
  )

  # Bind results
  bpr <- round(do.call(cbind.data.frame, ppc$results), digits=3)
  colnames(bpr) <- c("Normality", "Heteroskedasticity", "Independence")
  row.names(bpr) <- c("p")

  res <- list(
    "Bayesian p-value" = bpr
  )

  # Cat
  cat(msg)
  print.listof(res)

}

# Print method == summary method for ppc
#' @export
print.ppc <- function(ppc) {
  summary(ppc)
}

# Plot method
#' @export
plot.ppc <- function(ppc, type=c("normality", "heteroskedasticity", "independence")) {

  type <- match.arg(type)

  # Get data
  data <- switch(
    type,
    "normality" = ppc$data$skewness,
    "heteroskedasticity" = ppc$data$homosked,
    "independence" = ppc$data$independence
  )

  # Plot data
  data <- data %>%
    as.data.frame() %>%
    tidyr::gather(dataset, value)

  # Update labels
  data$dataset <- ifelse(data$dataset == "V1", "Simulated", "Observed")

  # Plot
  data %>%
    ggplot2::ggplot(., ggplot2::aes(x=value, fill=dataset)) +
    ggplot2::geom_density(alpha=0.6) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Observed and simulated results for test '", type, "'"))

}
