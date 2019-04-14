## Methods for posterior predictive checks

# Exported -----

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
#' @importFrom magrittr '%>%'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggExtra ggMarginal
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
    as.data.frame()

  # Add index
  data$index <- 1:nrow(data)

  # Gather data
  data <- tidyr::gather(data, dataset, value, -index)

  # Update labels
  data$dataset <- ifelse(data$dataset == "V1", "Simulated", "Observed")

  # Plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x=index, y=value, color=dataset)) +
    ggplot2::geom_point(alpha=0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::ggtitle(paste0("Observed and simulated results for test '", type, "'"))

  # Add marginal histogram
  ggExtra::ggMarginal(p, type = "histogram", margins="y", groupColour = FALSE, alpha=0.3,
                      groupFill = TRUE)

}

# Not exported -----

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
