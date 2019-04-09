# S3 Methods

#update.blm

# Class 'blm' methods -----



# Summary method
#' @export
#' @importFrom crayon bold
#' @importFrom crayon green
#' @importFrom crayon red
summary.blm <- function(blm) {

  var_names <- c(colnames(blm$input$X), "sigma")

  ### Formula
  form <- as.character(blm$input$formula)
  formchar <- paste0(form[2], " ", form[1], " ", form[3])

  ## Construct individual parts
  general <- paste0(
    "Formula: '", crayon::italic(formchar), "'"
  )

  ## User has already sampled or not
  has_sampled <- paste0(
    "Sampled: ", ifelse("posterior" %in% names(blm),
                        crayon::green('TRUE'),
                        crayon::red('FALSE'))
  )

  ## num. observations + predictors
  obs <- list(
    "Sampling settings" = matrix(c(blm$input$n,blm$input$m - 1,blm$sampling_settings$chains,
                  blm$sampling_settings$iterations,blm$sampling_settings$thinning, blm$sampling_settings$burn),
                nrow = 1,
                dimnames = list(
                  c(""),
                  c("Obs.",
                    "Predictors",
                    "Chains",
                    "Iterations",
                    "Thinning",
                    "Burn")
                )))

  ### If not sampled yet ...
  if(!"posterior" %in% names(blm)) {
    cat(crayon::bold("Model results for blm object:"))
    cat("\n\n")
    cat(general)
    cat("\n\n")
    cat(has_sampled)
    cat("\n\n")
    print.listof(obs)
    return(cat(""))
  }

  ### Statistics

  # Calculate MAP for each chain & combine
  MAPV <- do.call(cbind.data.frame, MAP(blm$posterior))
  # Add MC error
  MAPV[,3] <- MAPV[,2] / sqrt(blm$sampling_settings$iterations)
  # Round
  MAPV <- round(MAPV, digits = 3)

  # Amend names
  colnames(MAPV) <- c("Est. (mean)", "SD", "MCERR.")
  rownames(MAPV) <- c(paste0("b", 0:(length(var_names)-2)), "sigma")

  # Calculate CI
  CIV <- t(round(CI(blm$posterior), digits = 3))
  # Rownames
  rownames(CIV) <- c(paste0("b", 0:(length(var_names)-2)), "sigma")

  # Print MAP & SE
  cat(crayon::bold("Model results for blm object:"))
  cat("\n\n")
  cat(general)
  cat("\n\n")
  cat(has_sampled)
  cat("\n\n")
  print.listof(obs)
  print.listof(list("Maximum a posteriori (MAP) estimates" = MAPV))
  print.listof(list("95% credible interval" = CIV))

}



# Get coefficients from blm object
#' @export
coef.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg if user does not specify type
  type <- match.arg(type)

  # Bind posterior data
  pb <- do.call(rbind.data.frame, blm$posterior)
  # Remove sigma
  pb <- pb[,-ncol(pb)]

  # Compute
  r <- switch(type,
              "mode" = apply(pb, 2, calc_mode),
              "mean" = apply(pb, 2, mean),
              "median" = apply(pb, 2, median))

  # Return
  return(r)

}

# Call coef method using coefficients
#' @export
coefficients.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg
  type <- match.arg(type)

  # Call coef
  return(coef(blm, type))

}

# Predict method
#' @export
predict.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg
  type <- match.arg(type)

  # Get coefficients
  w <- matrix(coef(blm, type = type), ncol=1)

  # Predict
  return(
    blm$input$X %*% w
  )

}

# Residuals
#' @export
residuals.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg
  type <- match.arg(type)

  # Predict
  pred <- predict(blm, type=type)

  # Subtract
  return(blm$input$y - pred)

}

# Residuals
#' @export
resid.blm <- function(blm, type = c("mean", "mode", "median")) {

  # Match arg
  type <- match.arg(type)

  return(residuals(blm, type = type))

}



# Evaluate method
#  (plots, statistics etc.)
#  use functions from other packages to do this

# Class 'Prior' methods -----



# Class 'ppc' methods -----

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
