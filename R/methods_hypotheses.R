## Methods for hypothesis / hypotheses classes

# Print method for class hypothesis
#' @export
print.hypothesis <- function(x) {

  # Hypothesis
  usr_hyp <- x$hypothesis

  # Number of elements
  elements <- length(x$parsed)

  # Name
  name <- x$name

  # cat
  msg <- paste0(
    crayon::bold(name), "\n\n",
    "Hypothesis: ", usr_hyp, "\n",
    "Elements: ", elements
  )
  cat(msg)

}

# Summary for hypotheses
#' @export
summary.hypotheses <- function(x) {

  # Grab each hypothesis
  hyps <- lapply(x, function(x) x$hypothesis) %>%
    unlist() %>%
    unname()

  # Append H number
  for(i in seq_along(hyps)) hyps[i] <- paste0("H", i, ": ", hyps[i])

  # Add H_u
  hyps <- c(hyps, "Hu:")

  # Make results
  res <- matrix(0, ncol=5, nrow=length(hyps))
  colnames(res) <- c("complexity", "fit", "BF_c", "Pr_a", "Pr_b")
  row.names(res) <- hyps

  # Get complexity
  comp <- lapply(x, function(x) x$result$hypothesis$complexity) %>% unlist() %>% unname()
  # Get fit
  fit <- lapply(x, function(x) x$result$hypothesis$fit) %>% unlist() %>% unname()
  # Get bf
  BF_u <- lapply(x, function(x) x$result$BF_u) %>% unlist() %>% unname()
  BF_c <- lapply(x, function(x) x$result$BF_c)%>% unlist() %>% unname()
  # Set pr_a NOT SURE
  Pr_a <- BF_u / sum(BF_u)
  # Set pr_b NOT SURE
  Pr_b <- BF_u / (1 + sum(BF_u))

  # For Hu
  Pr_bu <- 1 - sum(Pr_b)

  # Add to matrix
  res[,1] <- c(comp, 0)
  res[,2] <- c(fit, 0)
  res[,3] <- c(BF_c, 0)
  res[,4] <- c(Pr_a, 0)
  res[,5] <- c(Pr_b, Pr_bu)

  # Round
  res <- round(res, digits=3)

  # Set H0 values to 0
  H0_row <- nrow(res)
  res[H0_row, -5] <- NA

  # To df
  df <- as.data.frame(res)
  NA_vals <- is.na(df)
  # Do not print NA values but print .
  df <- format(df)
  df[NA_vals] <- "."

  # Cat
  print.listof(
    list("Hypothesis evaluation"=df),
    na.print="."
  )

}

# Print method
#' @export
print.hypotheses <- function(x) {

  x[[1]]$name <- names(x)[1]
  print(x[[1]])
  for(h in 2:length(x)) {
    x[[h]]$name <- names(x)[h]
    cat("\n\n")
    print(x[[h]])
  }

}
