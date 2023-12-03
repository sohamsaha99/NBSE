# file for generating data based on parameters
# source("../utils/utils.R", chdir = TRUE, local = TRUE)

# params contain parameters that denote truth and additional parameters
# returns a list with entries: data, truth_local, knowledge_local
generate_data <- function(params) {
  P <- params[["P"]]
  # Use the following when overlap fraction is given
  overlap <- params[["overlap"]]
  N <- nrow(P)
  m0 <- floor(N * overlap / 2)
  V <- list(1:(0.5*N + m0), 0.5*N - m0 + 1:(0.5*N + m0))
  # Use the following when subgraphs are given
  # V <- params[["V"]]
  A <- generate_full_adjacency(P)
  A <- mask_adjacency(A, V)
  # No oracle knowledge
  list(
    data = list(adjacency = A, V = V)
  )
}

# params contain parameters that denote truth and additional parameters
# returns a list with entries: data, truth_local, knowledge_local
generate_data_given_subgraphs <- function(params) {
	P <- params[["P"]]
	V <- params[["V"]]
	A <- generate_full_adjacency(P)
	A <- mask_adjacency(A, V)
	# generating oracle knowledge
	list(
		data = list(adjacency = A, V = V)
	)
}

# Create random adjacency matrix from P
# Input : P. Probability matrix. Symmetric
# Output: Random adjacency matrix from this model
generate_full_adjacency <- function(P) {
  N <- nrow(P)
  A <- P * 0.0
  # Create adjacency matrix
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      A[i, j] <- rbinom(1, 1, P[i, j])
      A[j, i] <- A[i, j]
      # Diagonals remain zero
    }
  }
  A
}

# Create adjacency matrix from overlapping subgraph structure
# Input : A. Full adjacency matrix
#       : V. List of vectors to indicate subgraphs. V[[k]] is index set of V_k
# Output: Matrix with same dimension as A. Unobserved entries are -1
mask_adjacency <- function(A, V) {
  N <- nrow(A)
  observed <- matrix(0, nrow = N, ncol = N)
  for (k in seq_len(length(V))) {
    observed[V[[k]], V[[k]]] <- 1
  }
  observed * A + (1 - observed) * (-1)
}
