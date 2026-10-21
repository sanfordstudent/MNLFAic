test_that("Gauss-Hermite rule integrates polynomials exactly", {
  for (n in c(1, 2, 5, 21, 35)) {
    gh <- .gauss_hermite(n)
    w <- exp(gh$logw - gh$z^2)
    expect_equal(sum(w), sqrt(pi), tolerance = 1e-13)
    for (k in seq_len(n - 1)) {                   # exact up to degree 2n - 1
      expect_equal(sum(w * gh$z^(2 * k)), gamma(k + 0.5), tolerance = 1e-12)
      expect_lt(abs(sum(w * gh$z^(2 * k - 1))), 1e-10 * gamma(k + 0.5))
    }
  }
})

test_that("Gauss-Hermite rule: known two-point rule, sorted symmetric nodes", {
  gh <- .gauss_hermite(2)
  expect_equal(gh$z, c(-1, 1) / sqrt(2), tolerance = 1e-14)
  expect_equal(exp(gh$logw - gh$z^2), rep(sqrt(pi) / 2, 2), tolerance = 1e-14)
  g21 <- .gauss_hermite(21)
  expect_false(is.unsorted(g21$z))
  expect_equal(g21$z, -rev(g21$z), tolerance = 1e-13)
})

test_that(".with_seed reproduces draws and restores the caller's stream", {
  a <- .with_seed(3, stats::rnorm(4)); b <- .with_seed(3, stats::rnorm(4))
  expect_identical(a, b)
  set.seed(1); r1 <- stats::runif(2)
  set.seed(1); invisible(.with_seed(99, stats::rnorm(10))); r2 <- stats::runif(2)
  expect_identical(r1, r2)
  set.seed(2); x1 <- .with_seed(NULL, stats::runif(1)); set.seed(2)
  expect_identical(x1, stats::runif(1))
})
