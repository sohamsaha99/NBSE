# Function for USVT
library(graphon)
USVT <- function(dat, ...) {
  A <- dat$adjacency
  A[A == -1] <- 0
  est.USVT(A, eta = 0.01)$P
}