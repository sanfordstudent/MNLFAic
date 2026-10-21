# Automatic Stage-1 anchor. Ported from chen-item-cov-ext/tests/test_anchor.R.
# The anchor is bookkeeping: the maximum likelihood, the anchor-invariant
# parameters and their SEs, and every Stage-2 result must agree with an
# anchor-1 fit whenever both reach the optimum.

anchor_setup <- function() {
  set.seed(51)
  J <- 25; K <- 2; N <- 2000
  items <- data.frame(aoa = stats::rnorm(J, 7, 1.5), pos = sample(c("noun", "verb", "adjective"), J, TRUE))
  des <- .build_item_design(items, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"), standardize = TRUE)
  x <- cbind(female = stats::rbinom(N, 1, 0.5), age = as.vector(scale(stats::rnorm(N))))
  lam <- stats::runif(J, 0.8, 2); nu0 <- stats::rnorm(J, 0, 1); nu0[1] <- 3.5            # item 1: easy ...
  eps <- matrix(0, J, K); eps[1, 1] <- -1.2; eps[c(6, 12), 1] <- c(0.5, -0.5); eps[c(9, 20), 2] <- c(0.4, -0.3)
  gam <- cbind(c(0.05, -0.1, 0.05), c(0, 0.05, -0.05))                                     # ... and much easier for girls
  sim <- .simulate_itemcov_data(x, lam, nu0, des$W, gam, eps, c(0.1, 0.5), c(-0.1, 0.1), seed = 52)
  st <- list(a = rep(1, J), b = stats::qlogis(pmin(pmax(colMeans(sim$dat), .01), .99)), g = matrix(0, J, K),
             beta = c(0, 0), eta = c(0, 0))
  f <- x[, "female"] == 1
  cells <- pmin(colSums(sim$dat[f, ] == 0), colSums(sim$dat[!f, ] == 0), colSums(sim$dat[f, ] == 1), colSums(sim$dat[!f, ] == 1))
  list(J = J, K = K, des = des, x = x, sim = sim, st = st, f = f, cells = cells)
}

test_that("automatic anchor avoids a nearly separated item and weakly discriminating items", {
  s <- anchor_setup()
  a <- .choose_stage1_anchor(s$sim$dat, s$x)
  r_rest <- .item_rest_correlation(s$sim$dat)
  elig <- r_rest >= stats::median(r_rest)
  expect_false(a == 1)
  expect_true(elig[a])
  expect_equal(s$cells[a], max(s$cells[elig]))
  # A balanced but uninformative item (the Kim et al. case) is screened out.
  dat2 <- s$sim$dat; dat2[, 2] <- 0
  set.seed(53)
  for (grp in list(which(s$f), which(!s$f))) dat2[sample(grp, floor(length(grp) / 2)), 2] <- 1
  cells2 <- pmin(colSums(dat2[s$f, ] == 0), colSums(dat2[!s$f, ] == 0), colSums(dat2[s$f, ] == 1), colSums(dat2[!s$f, ] == 1))
  expect_equal(which.max(cells2), 2)                # the count rule alone would pick it
  expect_false(.choose_stage1_anchor(dat2, s$x) == 2)
  dna <- s$sim$dat; set.seed(54); dna[sample(length(dna), 500)] <- NA
  expect_false(anyNA(.item_rest_correlation(dna)))
})

test_that("anchor 1 vs. automatic anchor: same optimum, same invariant quantities", {
  skip_on_cran()
  s <- anchor_setup()
  f1 <- suppressWarnings(.fit_stage1(s$sim$dat, s$x, s$st, anchor = 1))
  fa <- .fit_stage1(s$sim$dat, s$x, s$st)
  ix <- .stage1_index(s$J, s$K)
  se <- function(ff) sqrt(diag(solve(ff$info)))
  s1 <- se(f1); sa <- se(fa)
  expect_identical(f1$anchor, 1L)
  expect_true(fa$convergence$converged && fa$convergence$em_fallbacks == 0)
  expect_lt(abs(f1$loglik - fa$loglik), 1e-4)
  expect_true(all(fa$gamma.mat[fa$anchor, ] == 0))
  expect_equal(f1$alpha.vec, fa$alpha.vec, tolerance = 1e-4)
  expect_equal(f1$beta.vec, fa$beta.vec, tolerance = 1e-4)
  expect_lt(max(abs(sa[c(ix$a, ix$b)] / s1[c(ix$a, ix$b)] - 1)), 1e-3)
  expect_equal(f1$eta, fa$eta, tolerance = 1e-4)
  expect_lt(max(abs(sa[ix$eta] / s1[ix$eta] - 1)), 1e-3)
  p1 <- .ml1_itemcov(f1$gamma.mat, f1$alpha.vec, f1$beta, s$des$W)
  pa <- .ml1_itemcov(fa$gamma.mat, fa$alpha.vec, fa$beta, s$des$W)
  expect_equal(p1$alpha, pa$alpha, tolerance = 1e-4)
  expect_equal(unlist(p1$gamma), unlist(pa$gamma), tolerance = 1e-4)
  expect_equal(p1$eps, pa$eps, tolerance = 1e-4)
  i1 <- .dif_inference_itemcov(f1, s$des$W, s$des, M = 3000, seed = 9)
  ia <- .dif_inference_itemcov(fa, s$des$W, s$des, M = 3000, seed = 9)
  r <- c(ia$impact$se / i1$impact$se, ia$category_effects$se / i1$category_effects$se)
  expect_true(all(abs(r - 1) < 0.05))
})
