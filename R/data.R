#' Compensation of 336 independent directors in S&P500 firms
#'
#' This dataset contains a cross-section of the compensation in GBR of 336 independent directors at 52 S&P500 firms in 3 sectors. The data were collected from Boardex (http://corp.boardex.com/) in 2014.
#'
#' @format A data frame containing 336 observations and 4 variables:
#' \describe{
#'    \item{Sector (factor)}{Sector to which the industry belongs.}
#'    \item{Industry (factor)}{Industry to which the company belongs.}
#'    \item{Company (factor)}{Company the independent directors wors at.}
#'    \item{Male (factor)}{Binary dummy indicating if the director is male.}
#'    \item{Age (numeric)}{Age of the director.}
#'    \item{Compensation (numeric)}{Compensation, in thousands GBR.}
#' }
#'
#' @source Boardex, 2014
"directors"
