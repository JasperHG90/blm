## Methods for miscellaneous classes

# R2 methods
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_histogram
#' @export
plot.R2 <- function(x) {

  # Turn draws to data frame
  df <- data.frame(y=x$rsquared)

  ggplot(df, aes(x=y)) +
    geom_histogram(color="black")

}
