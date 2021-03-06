## GGplot theme

#' Custom theme for blm plots
#'
#' @param text_size text size for blm plots
#' @param line_size size of the axis lines
#' @param rect_size sizs of the background lines
#'
#' @return adds blm theme to a plot generated by ggplot2
#'
#' @seealso https://www.r-bloggers.com/custom-themes-in-ggplot2/
#'
#' @export
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 '%+replace%'
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 rel
#' @importFrom grDevices rgb
theme_blm <- function(text_size = 12,
                      line_size = text_size / 170,
                      rect_size = text_size / 170){
  # Take theme minimal ...
  theme_minimal(base_size = text_size) %+replace%
  # ... and replace these values
    theme(
      plot.title = element_text(
        color = rgb(25, 43, 65, maxColorValue = 255),
        face = "bold",
        hjust = 0,
        size = rel(1)),
      plot.subtitle = element_text(
        color = rgb(25, 43, 65, maxColorValue = 255),
        face = "italic",
        hjust = 0,
        size = rel(0.9)),
      legend.title = element_text(
        color = rgb(105, 105, 105, maxColorValue = 255),
        size = rel(0.8)),
      legend.text = element_text(
        color = rgb(105, 105, 105, maxColorValue = 255),
        size = rel(0.6)),
      axis.title = element_text(
        color = rgb(105, 105, 105, maxColorValue = 255),
        size = rel(1.2)),
      axis.text = element_text(
        color = rgb(105, 105, 105, maxColorValue = 255),
        size = rel(1)),
      axis.line = element_line(color = rgb(105, 105, 105, maxColorValue = 255), size = 0.5),
      complete = TRUE
    )
}
