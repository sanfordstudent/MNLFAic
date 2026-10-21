# Cluster-robust Stage-1 covariance. Ported from chen-item-cov-ext/tests/test_cluster.R.

cluster_setup <- function(icc = 0.15) {
  set.seed(71)
  J <- 20; G <- 40; m <- 50; N <- G * m
  cl <- rep(seq_len(G), each = m)
  items <- data.frame(pos = sample(c("noun", "verb", "adjective"), J, TRUE))
  des <- .build_item_design(items, categorical = "pos", omit = list(pos = "noun"))
  lam <- stats::runif(J, 0.8, 2); nu0 <- stats::rnorm(J, 0, 1)
  treat <- rep(rep(0:1, length.out = G), each = m)                 # constant within cluster
  x <- cbind(treat = treat, age = as.vector(scale(stats::rnorm(N))))
  eps <- matrix(0, J, 2); eps[c(3, 9), 1] <- c(0.5, -0.4); eps[c(5, 14), 2] <- c(0.3, -0.3)
  gam <- cbind(c(0.05, -0.1), c(0, 0.05))
  sim0 <- .simulate_itemcov_data(x, lam, nu0, des$W, gam, eps, c(0.2, 0.5), c(0, 0), seed = 72)
  sim  <- .simulate_itemcov_data(x, lam, nu0, des$W, gam, eps, c(0.2, 0.5), c(0, 0), seed = 74, cluster = cl, icc = icc)
  st <- function(d) list(a = rep(1, J), b = stats::qlogis(pmin(pmax(colMeans(d), .01), .99)), g = matrix(0, J, 2),
                         beta = c(0, 0), eta = c(0, 0))
  list(J = J, G = G, N = N, cl = cl, des = des, x = x, sim0 = sim0, sim = sim, st = st)
}

test_that("simulator with icc keeps the marginal latent distribution", {
  s <- cluster_setup()
  res_th <- s$sim$theta - as.vector(s$x %*% c(0.2, 0.5))
  expect_lt(abs(stats::sd(res_th) - 1), 0.1)
  expect_gt(stats::var(tapply(res_th, s$cl, mean)) / stats::var(res_th), 0.05)   # clustered
})

test_that("input checks", {
  s <- cluster_setup()
  expect_error(.fit_stage1(s$sim$dat, s$x, s$st(s$sim$dat), cluster = replace(s$cl, 1, NA)), "missing values")
  expect_error(.fit_stage1(s$sim$dat, s$x, s$st(s$sim$dat), cluster = s$cl[-1]), "one value per person")
})

test_that("per-person scores; one person per cluster gives the usual sandwich", {
  skip_on_cran()
  s <- cluster_setup()
  f1 <- .fit_stage1(s$sim0$dat, s$x, s$st(s$sim0$dat), anchor = 1)
  S  <- .stage1_person_scores(f1, s$sim0$dat, s$x)
  sc <- .score_stage1(f1$alpha.vec, f1$beta.vec, f1$gamma.mat, s$x, s$sim0$dat, f1$beta, f1$eta, f1$post, f1$nodes)
  expect_equal(colSums(S), sc, tolerance = 1e-8)
  expect_lt(max(abs(colSums(S)) / sqrt(colSums(S^2))), 1e-3)
  expect_identical(f1$vcov_type, "model")
  expect_equal(f1$vcov, solve(f1$info), tolerance = 1e-10)
  fo <- .fit_stage1(s$sim0$dat, s$x, s$st(s$sim0$dat), anchor = 1, cluster = seq_len(s$N))
  r <- sqrt(diag(fo$vcov)) / sqrt(diag(solve(fo$info)))
  expect_lt(abs(stats::median(r) - 1), 0.05)
  expect_gte(mean(abs(r - 1) < 0.15), 0.9)
  expect_equal(fo$gamma.mat, f1$gamma.mat, tolerance = 1e-12)
})

test_that("clustered data: robust SEs grow for the cluster-level covariate; Stage 3 uses them", {
  skip_on_cran()
  s <- cluster_setup(); G <- s$G
  fc <- .fit_stage1(s$sim$dat, s$x, s$st(s$sim$dat), cluster = s$cl)
  fm <- .fit_stage1(s$sim$dat, s$x, s$st(s$sim$dat))
  ix <- .stage1_index(s$J, 2)
  rr <- sqrt(diag(fc$vcov)) / sqrt(diag(solve(fc$info)))
  expect_true(fc$vcov_type == "cluster" && fc$n_clusters == G && all(dim(fc$cluster_scores) == c(G, nrow(fc$info))))
  expect_gt(rr[ix$beta[1]], 1.5)
  expect_gt(rr[ix$beta[1]], rr[ix$beta[2]])
  expect_lte(qr(fc$vcov, tol = 1e-8)$rank, G - 1)
  fc1 <- .fit_stage1(s$sim$dat, s$x, s$st(s$sim$dat), anchor = 1, cluster = s$cl)
  idx <- c(ix$a, ix$b, ix$eta)
  expect_lt(max(abs(sqrt(diag(fc$vcov))[idx] / sqrt(diag(fc1$vcov))[idx] - 1)), 1e-3)
  ic  <- .dif_inference_itemcov(fc, s$des$W, s$des, M = 2000, seed = 5)
  im  <- .dif_inference_itemcov(fc, s$des$W, s$des, M = 2000, seed = 5, vcov = "model")
  im0 <- .dif_inference_itemcov(fm, s$des$W, s$des, M = 2000, seed = 5)
  icn <- .dif_inference_itemcov(fc, s$des$W, s$des, M = 2000, seed = 5, cluster_t = FALSE)
  expect_true(ic$settings$vcov_type == "cluster" && ic$settings$n_clusters == G && im$settings$vcov_type == "model")
  expect_gt(ic$impact$se[1], 1.5 * im$impact$se[1])
  expect_equal(im$impact$se, im0$impact$se, tolerance = 1e-12)
  expect_equal(im$residual$p, im0$residual$p, tolerance = 1e-12)
  expect_true(ic$settings$draw_df == G - 1 && is.na(icn$settings$draw_df) && is.na(im$settings$draw_df))
  rt <- stats::median(ic$impact$se / icn$impact$se)
  expect_lt(abs(rt / sqrt((G - 1) / (G - 3)) - 1), 0.05)
  expect_gt(mean(ic$residual$p), mean(icn$residual$p))
  s2 <- .restore_stage1_fit(.compact_stage1_fit(fc))
  expect_identical(s2$vcov_type, "cluster")
  expect_equal(s2$vcov, fc$vcov, tolerance = 1e-10)
  expect_warning(.fit_stage1(s$sim$dat, s$x, s$st(s$sim$dat), cluster = rep(1:10, length.out = s$N)), "clusters")
})
