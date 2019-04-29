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
  bpr <- round(do.call(cbind.data.frame, ppc$results), digits=4)
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

#' Plot posterior predictive samples
#'
#' Plotting is only possible if you choose to return 'return_samples=TRUE' in \code{evaluate_ppc()}
#'
#' @param ppc a posterior predictive check object (ppc)
#' @param type one of 'normality', 'heteroskedasticity' or 'independence'
#'
#' @importFrom magrittr '%>%'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggExtra ggMarginal
#' @importFrom latex2exp TeX
#'
#' @examples
#' data("directors")
#' fit <- blm("Compensation ~ Age", data=directors) %>%
#'    sample_posterior()
#' # Calculate PPC
#' fit <- fit %>%
#'    evaluate_ppc(return_samples=TRUE)
#' # Plot
#' plot(fit %>% get_value('ppc'), 'normality')
#'
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

  # Make y label for plot
  if(type == "normality") {
    ylab <- "Skewness"
  } else if(type == "independence") {
    ylab <- "Pearson correlation coefficient"
  } else {
    ylab <- latex2exp::TeX("Adj. R^2 (residuals)")
  }

  # Plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x=index, y=value, color=dataset)) +
    ggplot2::geom_point(alpha=0.5) +
    ggplot2::labs(title = "Observed (red) and simulated (blue) results",
                  subtitle = paste0("test: ", type)) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "Type") +
    ggplot2::scale_y_continuous(name = ylab) +
    theme_blm() +
    ggplot2::theme(legend.position = "none")

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
