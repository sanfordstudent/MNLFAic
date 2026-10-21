# Stage 2 (generalized ML1). Ported from chen-item-cov-ext/tests/test_stage2.R
# (sections 1, 3-5; logit metric pinned where the checks were written for it)
# and tests/test_difficulty_metric.R (sections 1-3).

flagged <- function(eps, tol = 1e-6) which(abs(eps) > tol)
toy_pos <- function() .build_item_design(data.frame(
  pos = factor(c(rep("noun", 10), rep("adj", 6), rep("verb", 4)), levels = c("noun", "adj", "verb"))),
  categorical = "pos", omit = list(pos = "noun"))$W

test_that("no item covariates: reproduces Chen et al.'s ML1 step exactly", {
  fx <- fixture("fixture_chen_baseline")
  r0 <- .ml1_itemcov(fx$em$gamma.mat, fx$em$alpha.vec, alpha_hat = fx$em$beta, metric = "logit")
  expect_equal(unname(r0$eps), unname(fx$ml1$gamma_ml1), tolerance = 1e-12)
  expect_equal(unname(r0$delta), fx$ml1$c_vec, tolerance = 1e-12)
  expect_equal(unname(r0$alpha), fx$inference$beta_ml1, tolerance = 1e-12)
  expect_equal(r0$nu_total, r0$eps, tolerance = 1e-12)
  expect_false(any(r0$nonunique))
  rd <- .ml1_itemcov(fx$em$gamma.mat, fx$em$alpha.vec, alpha_hat = fx$em$beta, metric = "difficulty")
  expect_equal(rd$alpha, r0$alpha, tolerance = 1e-12)   # same impact in both metrics without W
})

test_that("toy A: item covariates recover delta, gamma and the residual-DIF items; no-W does not", {
  set.seed(1); J <- 20; lam <- stats::runif(J, 0.8, 1.6); W <- toy_pos()
  eps <- rep(0, J); eps[c(1, 11)] <- c(0.8, -0.8)
  nu <- as.vector(W %*% c(0.6, 0.6)) + eps + lam * 0.3
  rA  <- .ml1_itemcov(cbind(nu), lam, W = W, metric = "logit")
  rA0 <- .ml1_itemcov(cbind(nu), lam, metric = "logit")
  expect_equal(unname(rA$delta), 0.3, tolerance = 1e-8)
  expect_equal(unname(rA$gamma[[1]]), c(0.6, 0.6), tolerance = 1e-8)
  expect_identical(flagged(rA$eps[, 1]), c(1L, 11L))
  expect_gt(abs(rA0$delta - 0.3), 0.1)
  expect_true(.ml1_condition_check(lam, W, support = c(1, 11), signs = c(1, -1), metric = "logit")$holds)
})

test_that("toy B: residual DIF on 3 of 4 verbs breaks identification with item covariates", {
  set.seed(1); J <- 20; lam <- stats::runif(J, 0.8, 1.6); W <- toy_pos()
  eps <- rep(0, J); eps[17:19] <- 0.8
  nu <- eps + lam * 0.3
  rB  <- .ml1_itemcov(cbind(nu), lam, W = W, metric = "logit")
  rB0 <- .ml1_itemcov(cbind(nu), lam, metric = "logit")
  expect_equal(unname(rB0$delta), 0.3, tolerance = 1e-8)
  expect_identical(flagged(rB0$eps[, 1]), 17:19)
  expect_gt(abs(rB$delta - 0.3), 0.1)
  cB <- .ml1_condition_check(lam, W, support = 17:19, signs = 1, metric = "logit")
  expect_false(cB$holds); expect_gte(cB$rho, 1)
  expect_true(.ml1_condition_check(lam, NULL, support = 17:19, signs = 1, metric = "logit")$holds)
})

test_that("toy C: centering matters in the logit metric but not in the difficulty metric", {
  set.seed(1); J <- 20; lam <- stats::runif(J, 0.8, 1.6)
  set.seed(2); aoa <- stats::rnorm(J, 6, 1.5)
  eps <- rep(0, J); eps[c(2, 5, 9)] <- c(1, -1, 1)
  Wc <- cbind(aoa = aoa - mean(aoa)); Wu <- cbind(aoa = aoa)
  nu <- 0.4 * (aoa - mean(aoa)) + eps + lam * 0.3
  rCc <- .ml1_itemcov(cbind(nu), lam, W = Wc, metric = "logit")
  rCu <- .ml1_itemcov(cbind(nu), lam, W = Wu, metric = "logit")
  expect_equal(unname(rCc$delta), 0.3, tolerance = 1e-8)
  expect_equal(unname(rCc$gamma[[1]]), 0.4, tolerance = 1e-8)
  expect_identical(flagged(rCc$eps[, 1]), c(2L, 5L, 9L))
  expect_false(identical(flagged(rCu$eps[, 1]), c(2L, 5L, 9L)))
  nud <- lam * 0.4 * (aoa - mean(aoa)) + eps + lam * 0.3
  rDc <- .ml1_itemcov(cbind(nud), lam, W = Wc, metric = "difficulty")
  rDu <- .ml1_itemcov(cbind(nud), lam, W = Wu, metric = "difficulty")
  expect_equal(rDc$eps, rDu$eps, tolerance = 1e-8)
  expect_identical(flagged(rDc$eps[, 1]), c(2L, 5L, 9L))
  expect_equal(rDc$gamma[[1]], rDu$gamma[[1]], tolerance = 1e-8)
  expect_equal(unname(rDu$delta), unname(rDc$delta) - 0.4 * mean(aoa), tolerance = 1e-8)
})

test_that("results are invariant to the Stage-1 anchor (both metrics)", {
  fx <- fixture("fixture_chen_baseline")
  nu1 <- fx$em$gamma.mat; lam1 <- fx$em$alpha.vec; al1 <- fx$em$beta
  set.seed(3)
  items20 <- data.frame(pos = sample(c("noun", "adj", "verb"), 20, replace = TRUE, prob = c(.5, .25, .25)),
                        aoa = stats::rnorm(20))
  W20 <- .build_item_design(items20, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"))$W
  for (anchor in c(5, 13)) {
    shift <- nu1[anchor, ] / lam1[anchor]
    nu2 <- nu1 - outer(lam1, shift); al2 <- al1 + shift
    for (m in c("logit", "difficulty")) {
      a <- .ml1_itemcov(nu1, lam1, al1, W = W20, metric = m)
      b <- .ml1_itemcov(nu2, lam1, al2, W = W20, metric = m)
      expect_equal(a$eps, b$eps, tolerance = 1e-8)
      expect_equal(unlist(a$gamma), unlist(b$gamma), tolerance = 1e-8)
      expect_equal(a$alpha, b$alpha, tolerance = 1e-8)
      expect_equal(a$nu_total, b$nu_total, tolerance = 1e-8)
    }
  }
})

test_that("different item covariates per person covariate", {
  fx <- fixture("fixture_chen_baseline")
  nu1 <- fx$em$gamma.mat; lam1 <- fx$em$alpha.vec; al1 <- fx$em$beta
  set.seed(3)
  items20 <- data.frame(pos = sample(c("noun", "adj", "verb"), 20, replace = TRUE, prob = c(.5, .25, .25)),
                        aoa = stats::rnorm(20))
  W20 <- .build_item_design(items20, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"))$W
  r0 <- .ml1_itemcov(nu1, lam1, al1, metric = "logit")
  rl <- .ml1_itemcov(nu1, lam1, al1, W = list(W20, NULL, W20[, "aoa", drop = FALSE]), metric = "logit")
  expect_equal(rl$eps[, 2], r0$eps[, 2], tolerance = 1e-12)
  expect_length(rl$gamma[[2]], 0)
  expect_equal(rl$eps[, 1], .ml1_itemcov(nu1, lam1, al1, W = W20, metric = "logit")$eps[, 1], tolerance = 1e-12)
  expect_identical(names(rl$gamma[[3]]), "aoa")
  expect_error(.ml1_itemcov(nu1, lam1, al1, W = list(W20)), "list of length C")
  expect_error(.ml1_itemcov(nu1, lam1[-1], al1), "length J")
})

# ── difficulty metric (tests/test_difficulty_metric.R, sections 1-3) ──────────
difficulty_setup <- function() {
  set.seed(31)
  J <- 60
  items <- data.frame(aoa = stats::rnorm(J, 7, 1.5),
                      pos = sample(c("noun", "verb", "adjective"), J, TRUE, prob = c(.5, .25, .25)))
  des <- .build_item_design(items, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"),
                            standardize = TRUE)
  lam <- stats::runif(J, 0.7, 2.2)
  gam <- cbind(female = c(0.05, -0.12, 0.04), age = c(-0.03, 0.02, 0.06))
  eps <- matrix(0, J, 2, dimnames = list(NULL, c("female", "age")))
  eps[c(5, 12, 20, 33, 41, 50), 1] <- c(0.5, -0.4, 0.6, 0.45, -0.5, 0.4)
  eps[c(8, 17, 26, 44, 55), 2] <- c(-0.3, 0.35, 0.4, -0.3, 0.3)
  alpha <- c(female = 0.1, age = 0.6)
  nu <- -lam * (des$W %*% gam + eps)
  list(J = J, des = des, W = des$W, lam = lam, gam = gam, eps = eps, alpha = alpha, nu = nu,
       anchor_at = function(a) { sh <- nu[a, ] / lam[a]; list(nu = nu - outer(lam, sh), alpha = alpha + sh) })
}

test_that("difficulty metric: relation to the logit metric without item covariates", {
  s <- difficulty_setup(); A <- s$anchor_at(1)
  rd <- .ml1_itemcov(A$nu, s$lam, A$alpha, metric = "difficulty")
  rl <- .ml1_itemcov(A$nu, s$lam, A$alpha, metric = "logit")
  expect_equal(rd$alpha, rl$alpha, tolerance = 1e-12)
  expect_equal(rd$delta, rl$delta, tolerance = 1e-12)
  expect_equal(rd$eps, -rl$eps / s$lam, tolerance = 1e-12)
  expect_equal(rd$nu_total, -rl$nu_total / s$lam, tolerance = 1e-12)
})

test_that("difficulty metric: exact recovery of the truth from any anchor", {
  s <- difficulty_setup()
  for (a in c(1, 5, 8)) {
    A <- s$anchor_at(a)
    r <- .ml1_itemcov(A$nu, s$lam, A$alpha, s$W)
    expect_equal(unname(r$alpha), unname(s$alpha), tolerance = 1e-10)
    expect_equal(unname(unlist(r$gamma)), as.vector(s$gam), tolerance = 1e-10)
    expect_equal(unname(r$eps), unname(s$eps), tolerance = 1e-10)
    expect_equal(unname(r$nu_total), unname(s$W %*% s$gam + s$eps), tolerance = 1e-10)
  }
  for (c in 1:2) {
    S <- which(s$eps[, c] != 0)
    expect_lt(.ml1_condition_check(s$lam, s$W, S, sign(s$eps[S, c]))$rho, 1)
  }
})

test_that("difficulty metric = slope-weighted median regression of difficulty DIF on [1, W]", {
  s <- difficulty_setup(); A <- s$anchor_at(5)
  r <- .ml1_itemcov(A$nu, s$lam, A$alpha, s$W)
  for (c in 1:2) {
    b <- -A$nu[, c] / s$lam
    X <- cbind(1, s$W)
    f <- quantreg::rq.fit(s$lam * X, s$lam * b, tau = 0.5, method = "br")
    expect_equal(unname(f$coefficients[-1]), unname(r$gamma[[c]]), tolerance = 1e-9)
    expect_equal(as.vector(b - X %*% f$coefficients), unname(r$eps[, c]), tolerance = 1e-9)
  }
})

test_that("all-zero DIF input is handled without calling the solver", {
  r <- .ml1_itemcov(matrix(0, 10, 1), stats::runif(10, 0.5, 2))
  expect_equal(unname(r$delta), 0); expect_true(all(r$eps == 0))
})
