# Function for EDS
library(graphon)
EDS <- function(dat, ...) {
  A <- dat$adjacency
  A[A == -1] <- 0
  est.LG(A, K = 2)$P
}