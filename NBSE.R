# Function for generating supergraph based on given subgraphs
# Input : list of K vectors. k^th vector contains indices for V_k
#       : weights. Possible values "intersection", "random", "binary"
# Output: K \times K symmetric adjacency matrix of the supergraph
#       : diagonals are 0. Other eies are according to "weights" parameter.
create_supergraph <- function(indices_list, weights = "intersection") {
  K <- length(indices_list)
  if (K == 1) return(matrix(0, nrow = 1, ncol = 1))
  adj <- matrix(0, nrow = K, ncol = K)
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      # Create edge weights based on intersection size
      adj[i, j] <- sum(indices_list[[i]] %in% indices_list[[j]])
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
  included[1] <- TRUE
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

# Get traversal based on sorting on edge weights
# Input : tree. A matrix with 3 columns. Each column represents an edge.
# Output: Sequence of nodes for traversal
get_traversal <- function(tree) {
  # How to get traversal?
  n <- nrow(tree) + 1
  seq_len(n)
}

# Estimate distance matrix of NBS. From package 'graphon'
# Input : A. Adjacency matrix with 0, 1
#       : N. Number of nodes.
# Output: Distance matrix of same dimension as described in Zhang et al.
aux_nbdsmooth <- function(A, N) {
  # definitions
  D = array(0,c(N,N))
  A_sq = (A%*%A)/N
  
  # compute D : dissimilarity
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      # val = max(abs(A_sq[i,]-A_sq[j,])) # this is my original work
      tgtvec = abs(A_sq[i,]-A_sq[j,])
      tgtvec[i] = 0
      tgtvec[j] = 0
      val = max(tgtvec) # tgtvec2 is Li Chen's
      D[i,j] = val
      D[j,i] = val
    }
  }
  
  # return result
  return(D)
}

# NBS, Zhang et al. From package 'graphon'
# Input : A. Adjacency matrix
# Output: Estimated probability matrix of same dimesion
est.nbdsmoothsingle <- function(A){
  # 1. size
  N = nrow(A)
  h = sqrt(log(N)/N)
  
  # 2. compute dissimilarity measure
  D = aux_nbdsmooth(A, N)
  
  # 3. quantiled as logical
  kernel_mat = matrix(0,N,N)
  for (i in 1:N){
    kernel_mat[i,] = as.double(D[i,]<quantile(D[i,],h))
  }
  
  # 4. L1 normalization of each row
  kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
  
  # 5. Compute P
  P = kernel_mat %*% A;
  P = (P+t(P))/2;
  
  ## (3) outputs
  res = vector("list")
  res$h = h
  res$P = P
  return(res)
}

# Example adjacency matrix
# A = matrix(rbinom(10000, 1, 1 - 2^(-0.5)), nrow = 100)
# A = A + t(A)
# A = 1 * (A > 0)
# diag(A) = 0

# NBSE
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

# D1, D2
d1 <- matrix(runif(36), nrow = 6)
d1 <- d1 + t(d1)
rownames(d1) <- colnames(d1) <- paste(c(2, 4, 1, 5, 6, 7))
d2 <- matrix(runif(49), ncol = 7)
d2 <- d2 + t(d2)
rownames(d2) <- colnames(d2) <- paste(c(3, 8, 1, 5, 4, 9, 10))

# Estimate distance matrix for NBSE
# Input : A. Adjacency matrix
#       : V. List of vectors to indicate subgraphs. V[[k]] is index set of V_k
# Output: Estimated distance matrix
est.distance.nbse <- function(A, V) {
  # Size
  n <- nrow(A)
  # Get traversal
  path <- get_traversal(create_supergraph(V))
  # Create distance matrix for first two subgraphs in path and combine
  v <- V[[1]]
  D1 <- aux_nbdsmooth(A[v, v], length(v))
  colnames(D1) <- rownames(D1) <- v
  v <- V[[2]]
  D2 <- aux_nbdsmooth(A[v, v], length(v))
  colnames(D2) <- rownames(D2) <- v
  D <- est.distance.extend(D1, D2)
  # Extend D matrix along path
  for (p in path[c(-1, -2)]) {
    v <- V[[p]]
    D2 <- aux_nbdsmooth(A[v, v], length(v))
    colnames(D2) <- rownames(D2) <- v
    D <- est.distance.extend(D, D2)
  }
  D
}

# TODO: Create adjacency matrix from graphon

