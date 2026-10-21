# Gauss-Hermite rule for weight exp(-z^2), computed in base R (replaces
# fastGHQuad::gaussHermiteData in the research code).
#
# Nodes: eigenvalues of the symmetric Jacobi matrix (Golub-Welsch), polished
# by Newton steps on the orthonormal Hermite recurrence. Weights: the
# Christoffel form w_l = 1 / sum_{k<n} p_k(z_l)^2 with p_k the orthonormal
# polynomials for exp(-z^2). This keeps full RELATIVE accuracy in the tiny
# outer weights, which matters because adaptive quadrature multiplies them by
# exp(z^2).
#
# Returns list(z = nodes (increasing), logw = log(w) + z^2), the form used by
# the adaptive E-step.
.gauss_hermite <- function(n) {
  stopifnot(n >= 1)
  if (n == 1) return(list(z = 0, logw = log(sqrt(pi))))
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off
  J[cbind(2:n, 1:(n - 1))] <- off
  z <- sort(eigen(J, symmetric = TRUE, only.values = TRUE)$values)
  # orthonormal Hermite polynomials p_0..p_n at z (weight exp(-z^2)):
  #   p_0 = pi^(-1/4), p_1 = sqrt(2) z p_0,
  #   p_{k+1} = sqrt(2/(k+1)) z p_k - sqrt(k/(k+1)) p_{k-1}
  ortho <- function(z) {
    P <- matrix(0, length(z), n + 1)
    P[, 1] <- pi^(-1 / 4)
    P[, 2] <- sqrt(2) * z * P[, 1]
    for (k in seq_len(n - 1)) P[, k + 2] <- sqrt(2 / (k + 1)) * z * P[, k + 1] - sqrt(k / (k + 1)) * P[, k]
    P
  }
  for (it in 1:3) {                      # Newton polish: p_n(z) = 0, p_n' = sqrt(2 n) p_{n-1}
    P <- ortho(z)
    z <- z - P[, n + 1] / (sqrt(2 * n) * P[, n])
  }
  P <- ortho(z)
  logw <- -log(rowSums(P[, 1:n, drop = FALSE]^2))
  list(z = z, logw = logw + z^2)
}
