## Methods for miscellaneous classes

#' Summary method for DIC
#'
#' @param x DIC object from blm library
#' @param ... DIC of the intercept-only model (if available)
#' @export
summary.DIC <- function(x, ...) {

  names(x) <- "Model fit"

  # Get opts
  opts <- list(...)
  # If element in opts
  if(length(opts) > 0) {
    x[[1]] <- rbind(opts[[1]][[1]], x[[1]])
    row.names(x[[1]]) <- c("Null model", "User model")
  }

  print.listof(x)

}

# Summary method for rsquared
#' @export
summary.R2 <- function(x) {

  # List for summary
  tab <- matrix(0L, nrow = 1, ncol = 5)
  rownames(tab) <- "(Model)"
  colnames(tab) <- c("2.5%", "25%", "50%", "75%", "97.5%")
  tab[1,] <- round(quantile(x$rsquared, c(0.025, 0.25, 0.5, 0.75, 0.975)), digits=4)

  # Print
  print.listof(list("Model R-squared" = tab))

}

# Summary method for Model Bayes Factor
#' @export
summary.BF <- function(x) {

  printBF <- list(
    "Wagemakers' Approx. Bayes Factor" = data.frame(
      "BF" = round(x$BF, digits=4)
    )
  )

  # Print
  print.listof(printBF)

}

# R2 methods
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_histogram
#' @importFrom latex2exp TeX
#' @export
plot.R2 <- function(x) {

  # Turn draws to data frame
  df <- data.frame(y=x$rsquared)

  # Tex
  title <- latex2exp::TeX("\\textbf{Bayesian $R^2$ posterior}")
  subt <- "vertical line is the median of the distribution"

  # Plot
  ggplot2::ggplot(df, ggplot2::aes(x=y)) +
    ggplot2::geom_histogram(color="black", fill="lightgrey") +
    ggplot2::labs(title = title, subtitle = subt) +
    ggplot2::geom_vline(xintercept = median(x$rsquared), linetype="dashed") +
    theme_blm() +
    ggplot2::theme(axis.title = element_blank())

}
