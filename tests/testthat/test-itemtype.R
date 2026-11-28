# Item types: 1PL items share one slope.

onepl_setup <- function(slope = 1.3, N = 1200) {
  set.seed(91); J <- 10
  x <- cbind(g1 = stats::rbinom(N, 1, 0.5), z = as.vector(scale(stats::rnorm(N))))
  eps <- matrix(0, J, 2); eps[2, 1] <- 0.6
  sim <- .simulate_itemcov_data(x, rep(slope, J), stats::rnorm(J, 0, 0.8), NULL, NULL, eps, c(0.3, 0.4), c(0, 0.1),
                                seed = 92)
  st <- list(a = rep(1, J), b = stats::qlogis(colMeans(sim$dat)), g = matrix(0, J, 2), beta = c(0, 0), eta = c(0, 0))
  list(J = J, N = N, x = x, dat = sim$dat, st = st)
}

test_that("restriction matrix and model covariance", {
  Tm <- .restriction_matrix(3, 1, c(TRUE, FALSE, TRUE))
  D <- 2 + 6 + 2
  expect_equal(dim(Tm), c(D, D - 1))
  expect_equal(Tm[3, ], Tm[5, ])                         # slopes of items 1 and 3 share a column
  expect_null(.restriction_matrix(3, 1, c(FALSE, FALSE, FALSE)))
  A <- crossprod(matrix(stats::rnorm(D * D), D)) + diag(D)
  expect_identical(.inv_info(A), solve(A))
  V <- .inv_info(A, Tm)
  expect_equal(V[3, ], V[5, ])
})

test_that("1PL fit: equal slopes, zero restricted score, anchor invariance, recovery", {
  skip_on_cran()
  s <- onepl_setup()
  f1 <- .fit_stage1(s$dat, s$x, s$st, onepl = rep(TRUE, s$J))
  f2 <- .fit_stage1(s$dat, s$x, s$st)
  expect_true(f1$convergence$converged)
  expect_equal(unname(f1$alpha.vec), rep(f1$alpha.vec[[1]], s$J), tolerance = 1e-12)
  expect_lt(f1$convergence$max_abs_score, 1e-4)
  expect_lte(f1$loglik, f2$loglik + 1e-6)
  expect_lt(abs(f1$alpha.vec[[1]] - 1.3), 3 * sqrt(f1$vcov[5, 5]))     # slope position 2K + 1
  fa <- .fit_stage1(s$dat, s$x, s$st, anchor = 7, onepl = rep(TRUE, s$J))
  expect_equal(fa$alpha.vec, f1$alpha.vec, tolerance = 1e-5)
  ia <- .stage1_index(s$J, 2)$a
  expect_equal(sqrt(diag(f1$vcov))[ia], rep(sqrt(f1$vcov[ia[1], ia[1]]), s$J), tolerance = 1e-10)
})

test_that("restricted information is the Jacobian of the restricted score", {
  skip_on_cran()
  s <- onepl_setup(N = 600)
  onepl <- rep(TRUE, s$J)
  f1 <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, onepl = onepl)
  Tm <- .restriction_matrix(s$J, 2, onepl)
  p0 <- .pack_stage1(f1$beta, f1$eta, f1$alpha.vec, f1$beta.vec, f1$gamma.mat)
  rscore <- function(p) {
    u <- .unpack_stage1(p, s$J, 2)
    es <- .estep_aq(u$a, u$b, u$g, s$x, s$dat, u$beta, u$eta, 21)
    as.vector(crossprod(Tm, .score_stage1(u$a, u$b, u$g, s$x, s$dat, u$beta, u$eta, es$q, es$Theta)))
  }
  I_r <- crossprod(Tm, f1$info %*% Tm)
  h <- 1e-5
  for (c in c(3, 5, 10)) {                              # eta, the common slope, an intercept
    dp <- h * Tm[, c]
    jac <- -(rscore(p0 + dp) - rscore(p0 - dp)) / (2 * h)
    expect_equal(jac, unname(I_r[, c]), tolerance = 1e-4)
  }
})

test_that("mnlfa_ic(itemtype =): all, mixed, named; inference paths; checks", {
  skip_on_cran()
  s <- onepl_setup(N = 800)
  pd <- as.data.frame(s$x)
  ctl <- mnlfa_control(M = 200, n_pilot = 100, seed = 1)
  colnames(s$dat) <- paste0("i", 1:s$J)
  f2  <- mnlfa_ic(s$dat, pd, ~ g1 + z, control = ctl)
  f2b <- mnlfa_ic(s$dat, pd, ~ g1 + z, itemtype = rep("2PL", s$J), control = ctl)
  expect_identical(f2b$stage1$gamma.mat, f2$stage1$gamma.mat)
  expect_identical(f2b$inference$residual, f2$inference$residual)
  f1 <- mnlfa_ic(s$dat, pd, ~ g1 + z, itemtype = "1PL", control = ctl)
  expect_output(print(f1), "item type: 1PL")
  sl <- confint(f1, "slopes")
  expect_equal(sl$se, rep(sl$se[1], s$J), tolerance = 1e-10)
  it <- stats::setNames(rep(c("1PL", "2PL"), each = 5), paste0("i", 1:s$J))[c(10:1)]   # named, any order
  fm <- mnlfa_ic(s$dat, pd, ~ g1 + z, itemtype = it, inference = FALSE, control = ctl)
  expect_equal(unname(fm$stage1$alpha.vec[1:5]), rep(fm$stage1$alpha.vec[[1]], 5), tolerance = 1e-12)
  expect_gt(stats::sd(fm$stage1$alpha.vec[6:10]), 0.01)
  expect_identical(fm$coding$itemtype, rep(c("1PL", "2PL"), each = 5))
  expect_output(print(fm), "5 1PL (one common slope), 5 2PL", fixed = TRUE)
  cl <- rep(1:40, length.out = s$N)
  fc <- mnlfa_ic(s$dat, pd, ~ g1 + z, itemtype = "1PL", cluster = cl, control = ctl)
  fw <- infer(fc, cluster_method = "wild")
  expect_equal(fw$inference$slopes$se, rep(fw$inference$slopes$se[1], s$J), tolerance = 1e-10)
  w <- stats::runif(s$N, 1, 2); school <- rep(1:40, each = s$N / 40)
  Rw <- ilsa_replicates(w, (school + 1) %/% 2, (school + 1) %% 2)
  fr <- mnlfa_ic(s$dat, pd, ~ g1 + z, itemtype = "1PL", weights = w, replicate_weights = Rw, control = ctl)
  expect_identical(fr$inference$settings$vcov_type, "replicate")
  expect_error(mnlfa_ic(s$dat, pd, ~ g1, itemtype = "3PL"), "1PL")
  expect_error(mnlfa_ic(s$dat, pd, ~ g1, itemtype = c("1PL", "2PL")), "one entry per item")
  s2 <- .restore_stage1_fit(.compact_stage1_fit(f1$stage1))
  expect_equal(s2$vcov, f1$stage1$vcov, tolerance = 1e-10)
})
