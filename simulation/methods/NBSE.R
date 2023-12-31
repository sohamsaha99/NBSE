# Needed to use distance estimation function from NBS
source("NBS.R", chdir = TRUE, local = TRUE)

# Function for NBSE
NBSE <- function(dat, num_iter = 8, ...) {
  A <- dat$adjacency
  V <- dat$V
  D <- est.distance.nbse(A, V)
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

# Function for generating supergraph based on given subgraphs
# Input : list of K vectors. k^th vector contains indices for V_k
#       : weights. Possible values "intersection", "random", "binary"
# Output: K \times K symmetric adjacency matrix of the supergraph
#       : diagonals are 0. Other entries are according to "weights" parameter.
create_supergraph <- function(indices_list, weights = "intersection") {
  K <- length(indices_list)
  if (K == 1) return(matrix(0, nrow = 1, ncol = 1))
  adj <- matrix(0, nrow = K, ncol = K)
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      # Create edge weights based on intersection size
      adj[i, j] <- sum(indices_list[[i]] %in% indices_list[[j]])
      if (adj[i, j] > 0) {
        adj[i, j] <- 1 / adj[i, j]
      }
      adj[j, i] <- adj[i, j]
    }
  }
  if (weights == "binary") {
    adj <- 1 * (adj > 0)
  } else if (weights == "random") {
    adj <- 1 * (adj > 0)
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        adj[i, j] <- adj[i, j] * runif(1)
        adj[j, i] <- adj[i, j]
      }
    }
  }
  adj
}

# Prim's algorithm to find minimal cost spanning tree
# Input : Adjacency matrix, weighted
# Output: T \times 3 matrix, each row represents an edge:
#         edge start, edge end, edge cost
find_mst <- function(adj) {
  # number of vertices
  K <- nrow(adj)
  # Replace 0s with infinity to denote no edge
  adj[adj == 0] <- Inf
  # Initialize matrix to return
  mst <- matrix(NA, nrow = K - 1, ncol = 3)
  # Indicator for each vertex whether included in tree
  included <- rep(FALSE, K)
  # Initiate with one vertex
  included[1] <- TRUE # Possibly make it random
  # Keep running until all vertices are included
  for (num_included in seq_len(K-1)) {
    # Find adj[u, v] for u included, v not included
    adj_temp <- adj[included, !included, drop = FALSE]
    # Find minimum cost edge
    index <- arrayInd(which.min(adj_temp), dim(adj_temp))
    # Add the new node to tree
    mst[num_included, 1] <- which(included)[index[1, 1]]
    mst[num_included, 2] <- which(!included)[index[1, 2]]
    mst[num_included, 3] <- adj_temp[index]
    # Change flag for vertex and increment count
    included[mst[num_included, 2]] <- TRUE
  }
  mst
}

# Get traversal based on sorting of edge weights
# Input : tree. A matrix with 3 columns. Each column represents an edge.
# Output: Sequence of nodes for traversal
get_traversal <- function(tree) {
  # How to get traversal?
  n <- nrow(tree) + 1
  seq_len(n)
  # Change third column of mst and re-run Prim's
}

# Extend distance matrix to two overlapping subgraphs
# Input : D1, D2. Distance matrices of two graphs V1, V2, with string dimnames
#       : Assumption - D1, D2 have dimnames, those are V1, V2 respectively.
# Output: Extend D to V_1 \cup V_2
est.distance.extend <- function(D1, D2) {
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
      s1 <- min(D1[i, inter] + D2[j, inter])
      s2 <- max(abs(D1[i, inter] - D2[j, inter]))
      D[i, j] <- (s1 + s2) / 2
      D[j, i] <- D[i, j]
    }
  }
  D
}

# Estimate distance matrix for NBSE
# Input : A. Adjacency matrix
#       : V. List of vectors to indicate subgraphs. V[[k]] is index set of V_k
# Output: Estimated distance matrix
est.distance.nbse <- function(A, V) {
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
  D <- est.distance.extend(D1, D2)
  # Extend D matrix along path
  for (p in path[c(-1, -2)]) {
    v <- V[[p]]
    D2 <- aux_nbdsmooth(A[v, v], length(v))
    colnames(D2) <- rownames(D2) <- paste(v)
    D <- est.distance.extend(D, D2)
  }
  D
}
