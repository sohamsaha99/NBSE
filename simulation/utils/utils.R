# Input : N. Number of nodes
#         graphon. Function with inputs x, y \in [0, 1]
# Output: N*N probability matrix from the graphon model
generate_P <- function(N, graphon) {
  # Sample \xi_1, ..., \xi_n iid from U(0, 1)
  xi <- runif(N, min = 0, max = 1)
  # xi <- sort(xi)
  # Create N*N matrix of probabilities
  P <- outer(xi, xi, FUN = graphon)
  list(P = P, xi = xi)
}

# Example graphon function
graphon <- function(x, y) {
  # 1 / (1 + exp(- x - y))
  # sin(5*pi*(x+y - 1) + 1)/2 + 0.5
  # 1 - 0.5 * (x + y + abs(x - y))
  (x^2 + y^2) * cos(1 / (x^2 + y^2)) / 3 + 0.15
}
