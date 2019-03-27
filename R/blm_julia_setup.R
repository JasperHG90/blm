## Julia setup
.blm <- new.env(parent = emptyenv())

#' Initial setup for blm package.
#'
#' \code{blm_setup} runs the julia setup for the blm package
#'
#' @param ... arguments passed to \code{JuliaCall::julia_setup}.
#'
#' @examples
#' \dontrun{
#' # blm_setup()
#' }
#'
#' @export
blm_setup <- function(...) {
  # Set up Julia
  .blm$julia <- JuliaCall::julia_setup(...)
  # Install Distributions/statistics package if needed
  .blm$julia$install_package_if_needed("Distributions")
  # Load Distributions/statistics package
  .blm$julia$library("Distributions")
  # Source gibbs sampler julia functions
  .blm$julia$source(system.file("julia/blm.jl", package = "blm"))
}
