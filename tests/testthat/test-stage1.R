# Stage 1 (adaptive quadrature, EM + Newton). Ported from
# chen-item-cov-ext/tests/test_stage1_aq.R and test_stage1_newton.R on the
# baseline data (N = 1000, J = 20, K = 3). The comparisons with the fixed-
# quadrature engine are not ported (that engine is not in the package); the
# fit is compared with Chen et al.'s original EM fit instead.

baseline <- function() {
  fx <- fixture("fixture_chen_baseline")
  x <- fx$data$x; dat <- fx$data$dat
  N <- nrow(dat)
  list(fx = fx, x = x, dat = dat, tr = fx$truth, em = fx$em, J = ncol(dat), K = ncol(x), N = N,
       obs = !is.na(dat), y0 = ifelse(is.na(dat), 0, dat),
       lin0 = rep(1, N) %*% t(fx$truth$b) + x %*% t(fx$truth$g))
}

# fine-grid log posterior kernel for person i
grid_person <- function(b, i, n_grid = 8001) {
  tr <- b$tr; x <- b$x
  mu <- sum(x[i, ] * tr$beta); s <- sqrt(exp(sum(x[i, ] * tr$eta)))
  u <- seq(-9, 9, length.out = n_grid); th <- mu + s * u
  e <- outer(th, tr$a) + rep(1, n_grid) %*% t(b$lin0[i, ])
  lp <- as.vector(stats::plogis(e, log.p = TRUE) %*% b$dat[i, ] + stats::plogis(-e, log.p = TRUE) %*% (1 - b$dat[i, ])) +
        stats::dnorm(th, mu, s, log = TRUE)
  mx <- max(lp)
  list(mode = th[which.max(lp)], loglik = mx + log(sum(exp(lp - mx)) * s * (u[2] - u[1])))
}

test_that("posterior modes: zero gradient, agree with a fine grid", {
  b <- baseline(); tr <- b$tr
  pm <- .posterior_mode(tr$a, b$lin0, b$x, b$obs, b$y0, tr$beta, tr$eta)
  P  <- stats::plogis(outer(pm$mode, tr$a) + b$lin0)
  grad <- as.vector(((b$y0 - P) * b$obs) %*% tr$a) -
          (pm$mode - as.vector(b$x %*% tr$beta)) / exp(as.vector(b$x %*% tr$eta))
  expect_lt(max(abs(grad)), 1e-6)
  gm <- sapply(1:30, function(i) grid_person(b, i)$mode)
  expect_lt(max(abs(gm - pm$mode[1:30])), 5e-3)
})

test_that("adaptive quadrature log-likelihood matches a fine grid", {
  b <- baseline(); tr <- b$tr
  es <- .estep_aq(tr$a, tr$b, tr$g, b$x, b$dat, tr$beta, tr$eta, n_nodes = 21)
  gl <- sapply(1:100, function(i) grid_person(b, i)$loglik)
  expect_lt(max(abs(es$loglik_i[1:100] - gl)), 1e-6)
  expect_lt(max(abs(rowSums(es$q) - 1)), 1e-12)
})

test_that("structural M-step is the exact maximizer on fixed nodes", {
  b <- baseline(); tr <- b$tr; K <- b$K; x <- b$x
  es <- .estep_aq(tr$a, tr$b, tr$g, x, b$dat, tr$beta, tr$eta)
  st <- .mstep_struct(x, es$Theta, es$q, beta = rep(0, K), eta = rep(0, K), n_outer = 50)
  negQ <- function(p) { mu <- as.vector(x %*% p[1:K]); s2 <- exp(as.vector(x %*% p[(K + 1):(2 * K)]))
    -sum(es$q * (-0.5 * log(s2) - 0.5 * (es$Theta - mu)^2 / s2)) }
  opt <- stats::optim(c(st$beta, st$eta) + 0.05, negQ, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_lt(max(abs(opt$par - c(st$beta, st$eta))), 1e-4)
  expect_lte(negQ(c(st$beta, st$eta)), opt$value + 1e-8)
})

test_that("analytic score matches the numerical gradient in all blocks", {
  b <- baseline(); em <- b$em; J <- b$J; K <- b$K
  set.seed(1)
  p0 <- .pack_stage1(em$beta, em$eta, em$alpha.vec, em$beta.vec, em$gamma.mat)
  p1 <- p0 + stats::rnorm(length(p0), 0, 0.05); u1 <- .unpack_stage1(p1, J, K)
  e1 <- .estep_aq(u1$a, u1$b, u1$g, b$x, b$dat, u1$beta, u1$eta)
  s_an <- .score_stage1(u1$a, u1$b, u1$g, b$x, b$dat, u1$beta, u1$eta, e1$q, e1$Theta)
  expect_length(s_an, nrow(b$fx$info))
  ll_p <- function(p) { u <- .unpack_stage1(p, J, K); .estep_aq(u$a, u$b, u$g, b$x, b$dat, u$beta, u$eta)$loglik }
  idx <- if (identical(Sys.getenv("NOT_CRAN"), "true")) seq_along(p1) else c(1:(2 * K), sample(seq_along(p1)[-(1:(2 * K))], 15))
  s_num <- sapply(idx, function(k) { e <- replace(numeric(length(p1)), k, 1e-5); (ll_p(p1 + e) - ll_p(p1 - e)) / 2e-5 })
  blk <- list(impact = which(idx <= 2 * K), items = which(idx > 2 * K))
  for (nm in names(blk)) {
    err <- max(abs(s_an[idx[blk[[nm]]]] - s_num[blk[[nm]]])) / max(abs(s_num[blk[[nm]]]))
    expect_lt(err, 1e-5, label = paste(nm, "block relative error"))
  }
})

test_that("Newton fit converges, is start-independent, and improves on the original EM fit", {
  skip_on_cran()
  b <- baseline(); tr <- b$tr; em <- b$em; J <- b$J; K <- b$K
  st1 <- list(a = tr$a, b = tr$b, g = matrix(0, J, K), beta = tr$beta, eta = tr$eta)
  st2 <- list(a = rep(1, J), b = rep(0, J), g = matrix(0, J, K), beta = rep(0, K), eta = rep(0, K))
  fa  <- .fit_stage1(b$dat, b$x, st1, anchor = 1)
  fa2 <- .fit_stage1(b$dat, b$x, st2, anchor = 1)
  pa  <- .pack_stage1(fa$beta, fa$eta, fa$alpha.vec, fa$beta.vec, fa$gamma.mat)
  pa2 <- .pack_stage1(fa2$beta, fa2$eta, fa2$alpha.vec, fa2$beta.vec, fa2$gamma.mat)
  p0  <- .pack_stage1(em$beta, em$eta, em$alpha.vec, em$beta.vec, em$gamma.mat)
  se  <- sqrt(diag(solve(fa$info)))
  expect_true(fa$convergence$converged && fa2$convergence$converged)
  s_fa <- .score_stage1(fa$alpha.vec, fa$beta.vec, fa$gamma.mat, b$x, b$dat, fa$beta, fa$eta, fa$post, fa$nodes)
  expect_lt(drop(t(s_fa) %*% solve(fa$info, s_fa)) / 2, 1e-6)           # Newton decrement
  expect_lt(max(abs(pa2 - pa) / se), 5e-3)                              # same optimum from a generic start
  expect_true(all(diff(fa$convergence$loglik_trace) >= -1e-6))
  expect_true(all(fa$gamma.mat[1, ] == 0))
  expect_gt(min(eigen(fa$info, symmetric = TRUE, only.values = TRUE)$values), 0)
  # the original tol = 1e-4 EM fit stopped ~0.15 SE short of the optimum
  ll_em <- .estep_aq(em$alpha.vec, em$beta.vec, em$gamma.mat, b$x, b$dat, em$beta, em$eta)$loglik
  expect_gte(fa$loglik, ll_em - 1e-8)
  expect_lt(max(abs(pa - p0) / se), 0.5)
})

test_that("non-convergence returns with a warning instead of an error", {
  b <- baseline(); J <- b$J; K <- b$K
  st <- list(a = rep(1, J), b = rep(0, J), g = matrix(0, J, K), beta = rep(0, K), eta = rep(0, K))
  expect_warning(f <- .fit_stage1(b$dat, b$x, st, anchor = 1, max_newton = 1, max_fallback = 0), "did not converge")
  expect_false(f$convergence$converged)
})
