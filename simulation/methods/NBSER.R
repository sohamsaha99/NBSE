# Needed to use distance estimation function from NBS
source("NBSE.R", chdir = TRUE, local = TRUE)

# Function for NBSER
NBSER <- function(dat, num_iter = 8, ...) {
  A <- dat$adjacency
  V <- dat$V
  D <- est.distance.nbser(A, V)
  # Matrix for indicator whether observed
  N <- nrow(A)
  observed <- matrix(0, nrow = N, ncol = N)
  for (k in seq_len(length(V))) {
    observed[V[[k]], V[[k]]] <- 1
  }
  Atemp <- A
  kernel_mat <- NULL
  for (iter in seq_len(num_iter)) {
    res <- est.nbdsmoothsingle(Atemp, D, kernel_mat)
    P <- res$P
    kernel_mat <- res$kernel_mat
    Atemp <- observed * A + (1 - observed) * P
  }
  # (P + t(P)) / 2
  P
}
# Extend distance matrix to two overlapping subgraphs
# Input : D1, D2. Distance matrices of two graphs V1, V2, with string dimnames
#       : Assumption - D1, D2 have dimnames, those are V1, V2 respectively.
# Output: Extend D to V_1 \cup V_2
est.distance.extendr <- function(D1, D2) {
  # Set of indices
  V1 <- rownames(D1)
  V2 <- rownames(D2)
  V <- union(V1, V2)
  n <- length(V)
  D <- matrix(NA, nrow = n, ncol = n)
  rownames(D) <- paste(V)
  colnames(D) <- rownames(D)
  # Fill D with values from D1, D2
  D[V1, V1] <- D1[V1, V1]
  D[V2, V2] <- D2[V2, V2]
  # Find intersection. On intersection, distance is average.
  inter <- paste(intersect(V1, V2))
  D[inter, inter] <- (D1[inter, inter] + D2[inter, inter]) / 2
  # For i \in V1 \ V2, j \in V_2 \ V_1
  diff12 <- paste(setdiff(V1, V2))
  diff21 <- paste(setdiff(V2, V1))
  for (i in diff12) {
    for (j in diff21) {
      # Find D[i, j]
      s <- mean(D1[i, inter] + D2[j, inter])
      D[i, j] <- s
      D[j, i] <- D[i, j]
    }
  }
  D
}

# Estimate distance matrix for NBSER
# Input : A. Adjacency matrix
#       : V. List of vectors to indicate subgraphs. V[[k]] is index set of V_k
# Output: Estimated distance matrix
est.distance.nbser <- function(A, V) {
  # Size
  n <- nrow(A)
  # Get traversal
  path <- get_traversal(find_mst(create_supergraph(V)))
  # Create distance matrix for first two subgraphs in path and combine
  v <- V[[1]]
  D1 <- aux_nbdsmooth(A[v, v], length(v))
  colnames(D1) <- rownames(D1) <- paste(v)
  v <- V[[2]]
  D2 <- aux_nbdsmooth(A[v, v], length(v))
  colnames(D2) <- rownames(D2) <- paste(v)
  D <- est.distance.extendr(D1, D2)
  # Extend D matrix along path
  for (p in path[c(-1, -2)]) {
    v <- V[[p]]
    D2 <- aux_nbdsmooth(A[v, v], length(v))
    colnames(D2) <- rownames(D2) <- paste(v)
    D <- est.distance.extendr(D, D2)
  }
  D
}
