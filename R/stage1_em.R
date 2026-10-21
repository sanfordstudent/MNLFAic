# Stage 1 building blocks with adaptive Gauss-Hermite quadrature (ported from
# chen-item-cov-ext/R/stage1_aq.R and the Theta-based parts of
# R/stage1_fast.R; the fixed-quadrature engine and the original EM of Chen et
# al. are not included).
#
# Model (anchored; item 1 = the anchor in the fitted order, its DIF fixed at 0):
#   logit P(y_ij = 1 | theta) = a_j theta + b_j + x_i' g_j,
#   theta_i ~ N(x_i' beta, exp(x_i' eta)).
# Adaptive quadrature centers the nodes at each person's posterior mode m_i and
# scales them by the posterior curvature s_i:
#   int f_i(theta) dtheta ~= sqrt(2) s_i sum_l w_l exp(z_l^2) f_i(m_i + sqrt(2) s_i z_l),
# with f_i(theta) = P(y_i | theta) N(theta; mu_i, sigma_i^2). Posterior weights
# q_il are proportional to w_l exp(z_l^2) f_i(theta_il). Everything downstream
# (item M-step, score, Louis information) needs only the node matrix Theta and
# the weights q. Missing responses (NA) contribute nothing.

# Posterior mode and curvature scale of theta for every person, by Newton's
# method (the log-posterior is strictly concave; steps capped at 1).
.posterior_mode <- function(a, lin0, x, obs, y0, beta, eta, start = NULL,
                            tol = 1e-8, max_iter = 100) {
  mu <- as.vector(x %*% beta); s2 <- as.vector(exp(x %*% eta))
  m  <- if (is.null(start)) mu else start
  for (it in seq_len(max_iter)) {
    P    <- stats::plogis(outer(m, a) + lin0)
    grad <- as.vector(((y0 - P) * obs) %*% a) - (m - mu) / s2
    info <- as.vector((P * (1 - P) * obs) %*% (a^2)) + 1 / s2
    step <- pmax(pmin(grad / info, 1), -1)
    m    <- m + step
    if (max(abs(step)) < tol) break
  }
  P    <- stats::plogis(outer(m, a) + lin0)
  info <- as.vector((P * (1 - P) * obs) %*% (a^2)) + 1 / s2
  list(mode = m, scale = 1 / sqrt(info), iterations = it)
}

# Adaptive E-step. Returns list(q = N x n posterior weights, loglik, loglik_i,
# Theta = N x n nodes, mode = posterior modes for warm starts).
.estep_aq <- function(a, b, g, x, dat, beta, eta, n_nodes = 21, start = NULL) {
  N <- nrow(dat)
  rule <- .gauss_hermite(n_nodes)
  obs  <- !is.na(dat); y0 <- ifelse(obs, dat, 0)
  lin0 <- rep(1, N) %*% t(b) + x %*% t(g)
  pm   <- .posterior_mode(a, lin0, x, obs, y0, beta, eta, start)
  mu   <- as.vector(x %*% beta); sd_i <- sqrt(as.vector(exp(x %*% eta)))
  Theta <- pm$mode + sqrt(2) * pm$scale %o% rule$z
  L <- matrix(0, N, n_nodes)                                 # log f_i(theta_il)
  for (l in seq_len(n_nodes)) {
    e <- outer(Theta[, l], a) + lin0
    L[, l] <- rowSums((y0 * stats::plogis(e, log.p = TRUE) + (1 - y0) * stats::plogis(-e, log.p = TRUE)) * obs) +
              stats::dnorm(Theta[, l], mu, sd_i, log = TRUE)
  }
  Lw  <- L + rep(1, N) %*% t(rule$logw)
  mx  <- apply(Lw, 1, max)
  ew  <- exp(Lw - mx)
  tot <- rowSums(ew)
  ll_i <- log(sqrt(2) * pm$scale) + mx + log(tot)
  list(q = ew / tot, loglik = sum(ll_i), loglik_i = ll_i, Theta = Theta, mode = pm$mode)
}

# Expected complete-data log-likelihood of each item (J-vector).
.q_items <- function(a, b, g, x, dat, Theta, q) {
  N <- nrow(dat); obs <- !is.na(dat); y0 <- ifelse(obs, dat, 0)
  lin0 <- rep(1, N) %*% t(b) + x %*% t(g)
  out <- numeric(ncol(dat))
  for (l in seq_len(ncol(Theta))) {
    eta_l <- outer(Theta[, l], a) + lin0
    out <- out + colSums(q[, l] * obs * (y0 * eta_l - log1p(exp(-abs(eta_l))) - pmax(eta_l, 0)))
  }
  out
}

# Item M-step: `n_newton` Newton steps on each item's expected complete-data
# log-likelihood, all items at once, with per-item step halving so that no
# item's objective decreases. Item 1 (the anchor) has g fixed at 0.
.mstep_items <- function(a, b, g, x, dat, q, Theta, n_newton = 2) {
  N <- nrow(dat); J <- ncol(dat); K <- ncol(x); m <- 2 + K
  obs <- !is.na(dat); y0 <- ifelse(obs, dat, 0)
  pairs <- which(upper.tri(diag(m), diag = TRUE), arr.ind = TRUE)   # (u, v), u <= v
  Qold  <- .q_items(a, b, g, x, dat, Theta, q)
  for (it in seq_len(n_newton)) {
    lin0 <- rep(1, N) %*% t(b) + x %*% t(g)
    G <- matrix(0, m, J)                   # gradient, rows = (a, b, g_1..g_K)
    H <- matrix(0, nrow(pairs), J)         # negative Hessian, upper triangle
    for (l in seq_len(ncol(Theta))) {
      th <- Theta[, l]
      P  <- stats::plogis(outer(th, a) + lin0)
      R  <- (y0 - P) * obs
      V  <- P * (1 - P) * obs
      Z  <- cbind(th, 1, x)
      G  <- G + crossprod(q[, l] * Z, R)
      Zp <- q[, l] * Z[, pairs[, 1], drop = FALSE] * Z[, pairs[, 2], drop = FALSE]
      H  <- H + crossprod(Zp, V)
    }
    step <- matrix(0, m, J)
    for (j in seq_len(J)) {
      Hj <- matrix(0, m, m); Hj[pairs] <- H[, j]; Hj[pairs[, 2:1]] <- H[, j]
      free <- if (j == 1) 1:2 else 1:m
      step[free, j] <- tryCatch(solve(Hj[free, free], G[free, j]), error = function(e) 0)
    }
    frac <- rep(1, J)
    for (h in 0:10) {
      a_n <- a + frac * step[1, ]; b_n <- b + frac * step[2, ]
      g_n <- g + t(step[3:m, , drop = FALSE]) * frac
      Qnew <- .q_items(a_n, b_n, g_n, x, dat, Theta, q)
      bad <- Qnew < Qold - 1e-10
      if (!any(bad) || h == 10) break
      frac[bad] <- frac[bad] / 2
    }
    frac[bad] <- 0                          # give up on items that still decrease
    a <- a + frac * step[1, ]; b <- b + frac * step[2, ]
    g <- g + t(step[3:m, , drop = FALSE]) * frac
    g[1, ] <- 0
    Qold <- .q_items(a, b, g, x, dat, Theta, q)
  }
  list(a = a, b = b, g = g)
}

# Structural M-step on fixed nodes: maximize sum_i sum_l q_il log N(theta_il;
# x_i' beta, exp(x_i' eta)). beta given eta: weighted least squares; eta given
# beta: Newton on a concave objective; alternated.
.mstep_struct <- function(x, Theta, q, beta, eta, n_outer = 3, n_inner = 10) {
  E1 <- rowSums(q * Theta)
  for (o in seq_len(n_outer)) {
    s2   <- as.vector(exp(x %*% eta))
    beta <- as.vector(solve(crossprod(x / s2, x), crossprod(x / s2, E1)))
    mu   <- as.vector(x %*% beta)
    V    <- rowSums(q * (Theta - mu)^2)
    for (it in seq_len(n_inner)) {
      r    <- V * exp(-as.vector(x %*% eta))
      grad <- crossprod(x, 0.5 * (r - 1))
      hess <- crossprod(x * (0.5 * r), x)
      step <- as.vector(solve(hess, grad))
      eta  <- eta + step
      if (max(abs(step)) < 1e-10) break
    }
  }
  list(beta = beta, eta = eta)
}

# EM with adaptive quadrature. Stops when the largest parameter change is
# below tol, or (with stall_ratio) once EM enters its slow linear crawl: after
# min_iter iterations, the change shrank by less than stall_ratio. EM is only a
# warm start for Newton.
.em_aq <- function(a, b, g, x, dat, beta, eta, n_nodes = 21, tol = 0.005,
                   max_iter = 500, stall_ratio = NULL, min_iter = 3, verbose = FALSE) {
  t0 <- proc.time()[3]
  J <- ncol(dat); K <- ncol(x)
  g <- if (is.null(dim(g))) matrix(g, J, K) else g
  g[1, ] <- 0
  mode <- NULL; prev_change <- Inf; stop_reason <- "max_iter"
  for (it in seq_len(max_iter)) {
    old <- c(a, b, g, beta, eta)
    es  <- .estep_aq(a, b, g, x, dat, beta, eta, n_nodes, start = mode)
    mode <- es$mode
    ip  <- .mstep_items(a, b, g, x, dat, es$q, es$Theta)
    st  <- .mstep_struct(x, es$Theta, es$q, beta, eta)
    a <- ip$a; b <- ip$b; g <- ip$g; beta <- st$beta; eta <- st$eta
    change <- max(abs(c(a, b, g, beta, eta) - old))
    if (verbose) message(sprintf("  EM %3d: loglik %.4f  max change %.2e  (%.1f s)",
                                 it, es$loglik, change, proc.time()[3] - t0))
    if (change < tol) { stop_reason <- "tol"; break }
    if (!is.null(stall_ratio) && it >= min_iter && change > stall_ratio * prev_change) {
      stop_reason <- "stalled"; break
    }
    prev_change <- change
  }
  list(beta = beta, eta = eta, alpha.vec = a, beta.vec = b, gamma.mat = g,
       iterations = it, stop_reason = stop_reason, mode = mode)
}
