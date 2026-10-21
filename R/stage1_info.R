# Stage-1 score and observed information (ported from score_stage1() in
# chen-item-cov-ext/R/stage1_fit.R, info_fast() in R/stage1_fast.R and
# stage1_person_scores() in R/stage1_fit.R). All take the quadrature nodes
# Theta (N x n) and posterior weights q (N x n) of the current E-step.
# Packed layout: beta (K), eta (K), a (J), b (J), g_{2..J} row-major ((J-1) K).

# Analytic score of the marginal log-likelihood (Fisher's identity: the
# posterior expectation of the complete-data score).
.score_stage1 <- function(a, b, g, x, dat, beta, eta, q, Theta) {
  N <- nrow(dat); J <- ncol(dat); K <- ncol(x); n <- ncol(Theta)
  mu_i   <- as.vector(x %*% beta)
  sig2_i <- as.vector(exp(x %*% eta))
  obs    <- !is.na(dat)
  y0     <- ifelse(obs, dat, 0)
  lin0   <- rep(1, N) %*% t(b) + x %*% t(g)
  da <- numeric(J); db <- numeric(J); dg <- matrix(0, J, K)
  s_beta <- numeric(N); s_eta <- numeric(N)
  for (l in seq_len(n)) {
    th  <- Theta[, l]
    P   <- stats::plogis(outer(th, a) + lin0)
    R   <- (y0 - P) * obs
    qR  <- q[, l] * R
    da  <- da + colSums(th * qR)
    db  <- db + colSums(qR)
    dg  <- dg + t(qR) %*% x
    r   <- th - mu_i
    s_beta <- s_beta + q[, l] * r / sig2_i
    s_eta  <- s_eta  + q[, l] * (-0.5 + 0.5 * r^2 / sig2_i)
  }
  c(as.vector(t(x) %*% s_beta), as.vector(t(x) %*% s_eta), da, db, as.vector(t(dg[-1, , drop = FALSE])))
}

# Observed information by Louis' identity, anchor removed:
#   I = sum_i E[-H_i | y_i] - sum_i Var(s_i | y_i),
# with the posterior variance of the complete-data score accumulated node by
# node as cross-products.
.info_stage1 <- function(a, b, g, x, dat, beta, eta, q, Theta) {
  N <- nrow(dat); J <- ncol(dat); K <- ncol(x); n <- ncol(Theta)
  D <- 2 * K + 2 * J + (J - 1) * K
  mu_i  <- as.vector(x %*% beta); sig2 <- as.vector(exp(x %*% eta))
  obs <- !is.na(dat); y0 <- ifelse(obs, dat, 0)
  lin0 <- rep(1, N) %*% t(b) + x %*% t(g)
  ix <- .stage1_index(J, K)
  gpos <- function(k) 2 * K + 2 * J + (seq_len(J - 1) - 1) * K + k

  score_l <- function(l) {                 # complete-data score of every person at node l
    th <- Theta[, l]; r <- th - mu_i
    R  <- (y0 - stats::plogis(outer(th, a) + lin0)) * obs
    S  <- matrix(0, N, D)
    S[, ix$beta] <- x * (r / sig2)
    S[, ix$eta]  <- x * (-0.5 + 0.5 * r^2 / sig2)
    S[, ix$a]    <- th * R
    S[, ix$b]    <- R
    for (k in seq_len(K)) S[, gpos(k)] <- x[, k] * R[, -1, drop = FALSE]
    S
  }
  Sbar <- matrix(0, N, D)
  for (l in seq_len(n)) Sbar <- Sbar + q[, l] * score_l(l)
  VarS <- matrix(0, D, D)
  for (l in seq_len(n)) VarS <- VarS + crossprod(sqrt(q[, l]) * (score_l(l) - Sbar))

  E_H  <- matrix(0, D, D)                  # expected negative complete-data Hessian
  E_r  <- rowSums(q * (Theta - mu_i))
  E_r2 <- rowSums(q * (Theta - mu_i)^2)
  E_H[ix$beta, ix$beta] <- crossprod(x / sig2, x)
  E_H[ix$eta,  ix$eta]  <- crossprod(x * (0.5 * E_r2 / sig2), x)
  E_H[ix$beta, ix$eta]  <- crossprod(x * (E_r / sig2), x)
  E_H[ix$eta,  ix$beta] <- t(E_H[ix$beta, ix$eta])
  m <- 2 + K
  pairs <- which(upper.tri(diag(m), diag = TRUE), arr.ind = TRUE)
  Hs <- matrix(0, nrow(pairs), J)
  for (l in seq_len(n)) {
    th <- Theta[, l]
    P  <- stats::plogis(outer(th, a) + lin0)
    V  <- P * (1 - P) * obs
    Z  <- cbind(th, 1, x)
    Hs <- Hs + crossprod(q[, l] * Z[, pairs[, 1], drop = FALSE] * Z[, pairs[, 2], drop = FALSE], V)
  }
  for (j in seq_len(J)) {
    Hj <- matrix(0, m, m); Hj[pairs] <- Hs[, j]; Hj[pairs[, 2:1]] <- Hs[, j]
    pos <- c(ix$a[j], ix$b[j], if (j > 1) sapply(seq_len(K), gpos)[j - 1, ])
    keep <- if (j > 1) 1:m else 1:2
    E_H[pos, pos] <- Hj[keep, keep]
  }
  E_H - VarS
}

# Per-person scores of the marginal log-likelihood at a fit's estimates
# (N x D, packed layout, item order as fitted: anchor first). Uses the fit's
# final posterior weights and nodes. Column sums equal .score_stage1().
.stage1_person_scores <- function(fit, dat, x) {
  a <- fit$alpha.vec; b <- fit$beta.vec; g <- fit$gamma.mat; q <- fit$post; Theta <- fit$nodes
  N <- nrow(dat); J <- ncol(dat); K <- ncol(x); n <- ncol(Theta)
  D <- 2 * K + 2 * J + (J - 1) * K
  mu_i <- as.vector(x %*% fit$beta); sig2 <- as.vector(exp(x %*% fit$eta))
  obs <- !is.na(dat); y0 <- ifelse(obs, dat, 0)
  lin0 <- rep(1, N) %*% t(b) + x %*% t(g)
  ix <- .stage1_index(J, K)
  gpos <- function(k) 2 * K + 2 * J + (seq_len(J - 1) - 1) * K + k
  S <- matrix(0, N, D)
  for (l in seq_len(n)) {
    th <- Theta[, l]; r <- th - mu_i; ql <- q[, l]
    R  <- (y0 - stats::plogis(outer(th, a) + lin0)) * obs
    S[, ix$beta] <- S[, ix$beta] + ql * x * (r / sig2)
    S[, ix$eta]  <- S[, ix$eta]  + ql * x * (-0.5 + 0.5 * r^2 / sig2)
    S[, ix$a]    <- S[, ix$a]    + ql * th * R
    S[, ix$b]    <- S[, ix$b]    + ql * R
    for (k in seq_len(K)) S[, gpos(k)] <- S[, gpos(k)] + ql * x[, k] * R[, -1, drop = FALSE]
  }
  S
}
