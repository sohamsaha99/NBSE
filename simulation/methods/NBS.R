# Function for NBS
library(Rcpp)
NBS <- function(dat, ...) {
  A <- dat$adjacency
  A[A == -1] <- 0
  D <- aux_nbdsmooth(A, nrow(A))
  est.nbdsmoothsingle(A, D)$P
  # P <- est.nbdsmoothsingle(A, D)$P
  # (P + t(P)) / 2
}

# C++ function to estimate D from A^2
cppFunction('NumericMatrix D_from_Asq(NumericMatrix Asq, int N) {
           NumericMatrix D(N, N);
           NumericVector tgtvec(N);
            for (int i=0; i < N-1; i++) {
              for (int j=i+1; j < N; j++) {
                tgtvec = abs(Asq(i, _) - Asq(j, _));
                tgtvec[i] = 0;
                tgtvec[j] = 0;
                double val = max(tgtvec);
                D(i, j) = val;
                D(j, i) = val;
              }
            }
            return D;
}')

# Estimate distance matrix of NBS. From package 'graphon'
# Input : A. Adjacency matrix with 0, 1
#       : N. Number of nodes.
# Output: Distance matrix of same dimension as described in Zhang et al.
aux_nbdsmooth <- function(A, N) {
  # definitions
  D = array(0,c(N,N))
  # A_sq = (A%*%A)/N
  A_sq = (crossprod(A))/N
  
  # compute D : dissimilarity
  D = D_from_Asq(A_sq, N)
  # for (i in 1:(N-1)){
  #   for (j in (i+1):N){
  #     # val = max(abs(A_sq[i,]-A_sq[j,])) # this is my original work
  #     tgtvec = abs(A_sq[i,]-A_sq[j,])
  #     tgtvec[i] = 0
  #     tgtvec[j] = 0
  #     val = max(tgtvec) # tgtvec2 is Li Chen's
  #     D[i,j] = val
  #     D[j,i] = val
  #   }
  # }
  
  # return result
  return(D)
}

# NBS, Zhang et al. From package 'graphon'. Slightly modified
# Input : A. Adjacency matrix
# Output: Estimated probability matrix of same dimesion
est.nbdsmoothsingle <- function(A, D, kernel_mat = NULL){
  # 1. size
  N = nrow(A)
  h = sqrt(log(N)/N)
  if (is.null(kernel_mat)) {
    # 2. compute dissimilarity measure
    # D is function argument
    # D = aux_nbdsmooth(A, N)
    
    # 3. quantiled as logical
    kernel_mat = matrix(0,N,N)
    for (i in 1:N){
      kernel_mat[i,] = as.double(D[i,]<quantile(D[i,],h))
    }
    
    # 4. L1 normalization of each row
    kernel_mat = kernel_mat/(outer(rowSums(kernel_mat),rep(1,N))+1e-10)
    # kernel_mat = (kernel_mat + t(kernel_mat)) / 2
  }
  # 5. Compute P
  P = kernel_mat %*% A
  P = (P+t(P))/2
  
  ## (3) outputs
  res = vector("list")
  res$h = h
  res$P = P
  res$kernel_mat = kernel_mat
  return(res)
}