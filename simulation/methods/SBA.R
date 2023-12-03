# Function for SBA
library(graphon)
SBA <- function(dat, ...) {
  A <- dat$adjacency
  A[A == -1] <- 0
  est.SBA(A, delta = 0.5)$P
}