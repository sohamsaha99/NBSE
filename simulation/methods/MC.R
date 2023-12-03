# Function for matrix completion
library(graphon)
MC <- function(dat, ...) {
  A <- dat$adjacency
  A[A == -1] <- 0
  est.completion(A)
}