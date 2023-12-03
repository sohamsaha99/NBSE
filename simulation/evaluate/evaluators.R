# one file for all evaluators
# these will work on replicated outputs
# not on single output.
frobenius_diff <- function(result, truth_global) {
  diff_vector <- rep(NA, length(result))
  k <- 0
  for (Phat in result) {
    k <- k + 1
    diff_vector[k] <- sqrt(mean((Phat - truth_global)^2))
  }
  mean(diff_vector)
}

d2inf_diff <- function(result, truth_global) {
  diff_vector <- rep(NA, length(result))
  k <- 0
  for (Phat in result) {
    k <- k + 1
    diff_vector[k] <- max(sqrt(rowSums((Phat- truth_global)^2)) / sqrt(nrow(Phat)))
  }
  mean(diff_vector)
}