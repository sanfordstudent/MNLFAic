# Stage 3 (bootstrap inference). Ported from chen-item-cov-ext/tests/
# test_stage3.R (sections 1-3) and test_difficulty_metric.R (sections 5-6).

baseline_fit <- function() {
  fx <- fixture("fixture_chen_baseline")
  fit <- c(fx$em, list(info = fx$info, N = fx$N))
  colnames(fit$gamma.mat) <- paste0("x", seq_len(ncol(fit$gamma.mat)))
  list(fx = fx, fit = fit, J = nrow(fit$gamma.mat), K = ncol(fit$gamma.mat))
}
baseline_design <- function(J) {
  set.seed(3)
  items20 <- data.frame(pos = sample(c("noun", "adj", "verb"), J, replace = TRUE, prob = c(.5, .25, .25)),
                        aoa = stats::rnorm(J))
  .build_item_design(items20, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"))
}

test_that("null centering without item covariates reproduces Chen et al.'s inference exactly", {
  b <- baseline_fit(); ic <- b$fx$inference_call; o <- b$fx$inference
  inf0 <- .dif_inference_itemcov(b$fit, W = NULL, center = "null", metric = "logit",
                                 M = ic$M, sig = ic$sig, seed = ic$seed)
  expect_equal(inf0$residual$est, as.vector(o$gamma_ml1), tolerance = 1e-12)
  expect_equal(inf0$impact$est, o$beta_ml1, tolerance = 1e-12)
  expect_equal(inf0$residual$p, as.vector(o$p_values), tolerance = 1e-12)
  expect_equal(inf0$residual$ci_low, as.vector(o$CI_low), tolerance = 1e-12)
  expect_equal(inf0$residual$ci_high, as.vector(o$CI_high), tolerance = 1e-12)
})

test_that("internal consistency with item covariates (three centerings)", {
  b <- baseline_fit(); J <- b$J; K <- b$K; des <- baseline_design(J)
  M <- 600
  res <- lapply(c(estimate = "estimate", threshold = "threshold", null = "null"), function(ctr)
    .dif_inference_itemcov(b$fit, W = des$W, design = des, center = ctr, metric = "logit",
                           M = M, seed = 11, keep_draws = TRUE))
  r <- res$threshold
  expect_identical(res$estimate$residual$est, res$threshold$residual$est)
  expect_identical(res$estimate$gamma$est, res$null$gamma$est)
  expect_identical(res$estimate$impact$est, res$null$impact$est)
  expect_true(nrow(r$residual) == J * K && nrow(r$total) == J * K && nrow(r$gamma) == 3 * K && nrow(r$impact) == K)
  expect_true(all(c(r$residual$se, r$gamma$se, r$impact$se, r$total$se) > 0))
  ps <- c(r$residual$p, r$gamma$p, r$impact$p)
  expect_true(all(ps >= 0 & ps <= 1))
  expect_true(all(r$residual$ci_low <= r$residual$ci_high) && all(r$gamma$ci_low <= r$gamma$ci_high) &&
              all(r$impact$ci_low <= r$impact$ci_high))
  ce1 <- r$category_effects[r$category_effects$covariate == "x1", ]
  expect_equal(ce1$est[ce1$term == "pos_noun"], -sum(ce1$est[ce1$term %in% c("pos_adj", "pos_verb")]), tolerance = 1e-12)
  expect_equal(ce1$se[ce1$term == "pos_noun"],
               stats::sd(-(r$draws$e_gamma[, "x1.pos_adj"] + r$draws$e_gamma[, "x1.pos_verb"])), tolerance = 1e-10)
  expect_true(nrow(r$joint) == K && all(r$joint$df == 2) && all(r$joint$p_boot >= 0 & r$joint$p_boot <= 1))
  expect_equal(r$total$est, as.vector(des$W %*% sapply(r$point$gamma, identity)) + r$residual$est, tolerance = 1e-10)
  expect_true(r$settings$n_residuals_zeroed > 0 && r$settings$n_residuals_zeroed < J * K)
  expect_identical(r$impact$ci_low, r$impact$ci_sym_low)
  expect_identical(r$gamma$ci_high, r$gamma$ci_sym_high)
  expect_identical(r$category_effects$ci_low, r$category_effects$ci_sym_low)
  expect_identical(r$residual$ci_low, r$residual$ci_eq_low)
  expect_identical(r$total$ci_high, r$total$ci_eq_high)
  expect_equal(r$impact$ci_sym_high - r$impact$est, r$impact$est - r$impact$ci_sym_low, tolerance = 1e-12)
  r_eq <- .dif_inference_itemcov(b$fit, W = des$W, design = des, metric = "logit", M = M, seed = 11,
                                 ci_shape_effects = "equal")
  expect_identical(r_eq$impact$ci_low, r_eq$impact$ci_eq_low)
  expect_identical(r_eq$impact$est, r$impact$est)
  expect_identical(r_eq$impact$ci_sym_low, r$impact$ci_sym_low)
  r3 <- .dif_inference_itemcov(b$fit, W = des$W, design = des, metric = "logit", M = M, seed = 11, threshold_z = 3)
  expect_lte(r$settings$n_residuals_zeroed, r3$settings$n_residuals_zeroed)   # lower threshold zeroes no more
})

test_that("threshold rules", {
  N <- 1000; J <- 20; K <- 3
  expect_equal(.resolve_threshold("sqrt2loglogN", N, J, K), sqrt(2 * log(log(N))))
  expect_equal(.resolve_threshold("sqrtlogN", N, J, K), sqrt(log(N)))
  expect_equal(.resolve_threshold("bonferroni", N, J, K), sqrt(2 * log(J * K)))
  expect_identical(.resolve_threshold(2.5, N, J, K), 2.5)
  expect_error(.resolve_threshold("other", N, J, K), "Unknown")
  b <- baseline_fit(); f <- b$fit; f$N <- NULL
  expect_error(.dif_inference_itemcov(f, M = 20), "number of persons")
})

test_that("difficulty metric: residual-DIF p-values are computed on the logit scale", {
  # With no item covariates and the estimate center, the logit-scale draws are
  # identical in the two metrics, so the p-values agree up to exact-tie rounding.
  b <- baseline_fit(); M <- 2000
  pd <- .dif_inference_itemcov(b$fit, W = NULL, center = "estimate", M = M, seed = 4, metric = "difficulty")
  pl <- .dif_inference_itemcov(b$fit, W = NULL, center = "estimate", M = M, seed = 4, metric = "logit")
  dp <- abs(pd$residual$p - pl$residual$p)
  expect_lte(max(dp), 10 / M)
  expect_lt(mean(dp > 0), 0.05)
  expect_equal(pd$residual$est, -pl$residual$est / rep(b$fit$alpha.vec, b$K), tolerance = 1e-10)
})

test_that("category-specific impacts (difficulty metric)", {
  b <- baseline_fit(); des <- baseline_design(b$J)
  inf <- .dif_inference_itemcov(b$fit, des$W, des, M = 300, seed = 3)
  ic <- inf$impact_by_category; ce <- inf$category_effects
  expect_equal(nrow(ic), 3 * b$K)
  expect_equal(ic$est, inf$impact$est[match(ic$covariate, inf$impact$covariate)] -
                 ce$est[match(paste(ic$covariate, ic$term), paste(ce$covariate, ce$term))], tolerance = 1e-12)
  expect_equal(as.vector(tapply(ic$est, ic$covariate, mean)[inf$impact$covariate]), inf$impact$est, tolerance = 1e-12)
  expect_null(.dif_inference_itemcov(b$fit, des$W, des, M = 50, seed = 3, metric = "logit")$impact_by_category)
  expect_error(.dif_inference_itemcov(b$fit, list(des$W, des$W, des$W), des, M = 50), "single shared W")
})

test_that("end to end: simulated data with item-covariate DIF (logit metric)", {
  skip_on_cran()
  set.seed(2026)
  Js <- 30; Ns <- 2000
  items <- data.frame(pos = rep(c("noun", "noun", "adj", "verb", "noun", "adj"), 5), aoa = stats::rnorm(Js))
  dsim  <- .build_item_design(items, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"))
  xs    <- cbind(female = stats::rbinom(Ns, 1, 0.5), age = stats::rnorm(Ns))
  lam   <- stats::runif(Js, 0.8, 2); nu0 <- stats::runif(Js, -1.5, 1.5)
  gam   <- cbind(female = c(aoa = 0.3, pos_adj = 0.4, pos_verb = -0.3), age = c(aoa = -0.2, pos_adj = 0, pos_verb = 0.25))
  eps   <- matrix(0, Js, 2); eps[c(4, 11, 20), 1] <- c(0.9, -0.9, 0.9); eps[c(7, 16, 25), 2] <- c(-0.7, 0.7, 0.7)
  alpha <- c(0.4, 0.3); omega <- c(0.2, -0.1)
  sim <- .simulate_itemcov_data(xs, lam, nu0, dsim$W, gam, eps, alpha, omega, seed = 99, metric = "logit")
  for (c in 1:2) {
    S <- which(eps[, c] != 0)
    expect_true(.ml1_condition_check(lam, dsim$W, S, sign(eps[S, c]), metric = "logit", margin = FALSE)$holds)
  }
  fs <- .fit_stage1(sim$dat, xs, start = list(a = rep(1, Js), b = rep(0, Js), g = matrix(0, Js, 2),
                                              beta = c(0, 0), eta = c(0, 0)))
  inf <- .dif_inference_itemcov(fs, W = dsim$W, design = dsim, metric = "logit", M = 5000, seed = 5)
  expect_true(all(abs((inf$impact$est - alpha) / inf$impact$se) < 3))
  expect_true(all(abs((inf$gamma$est - as.vector(gam)) / inf$gamma$se) < 3))
  expect_true(all(inf$residual$flag[as.vector(eps != 0)]))
  expect_lte(sum(inf$residual$flag[as.vector(eps == 0)]), 2)
})

test_that("end to end: difficulty metric (J = 60, N = 2500)", {
  skip_on_cran()
  set.seed(31)
  J <- 60
  items <- data.frame(aoa = stats::rnorm(J, 7, 1.5), pos = sample(c("noun", "verb", "adjective"), J, TRUE, prob = c(.5, .25, .25)))
  des <- .build_item_design(items, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"), standardize = TRUE)
  lam <- stats::runif(J, 0.7, 2.2)
  gam <- cbind(female = c(0.05, -0.12, 0.04), age = c(-0.03, 0.02, 0.06))
  eps <- matrix(0, J, 2); eps[c(5, 12, 20, 33, 41, 50), 1] <- c(0.5, -0.4, 0.6, 0.45, -0.5, 0.4)
  eps[c(8, 17, 26, 44, 55), 2] <- c(-0.3, 0.35, 0.4, -0.3, 0.3)
  alpha <- c(female = 0.1, age = 0.6)
  x <- cbind(female = stats::rbinom(2500, 1, 0.5), age = as.vector(scale(stats::rnorm(2500))))
  sim <- .simulate_itemcov_data(x, lam, stats::rnorm(J, 0, 1.2), des$W, gam, eps, alpha, c(-0.2, 0.1), seed = 7)
  expect_equal(sim$truth$nu, -lam * (des$W %*% gam + eps), tolerance = 1e-12)
  fit <- .fit_stage1(sim$dat, x, list(a = rep(1, J), b = stats::qlogis(colMeans(sim$dat)), g = matrix(0, J, 2),
                                      beta = c(0, 0), eta = c(0, 0)))
  colnames(fit$gamma.mat) <- c("female", "age")
  inf <- .dif_inference_itemcov(fit, des$W, des, M = 3000, seed = 3)
  expect_identical(inf$settings$metric, "difficulty")
  expect_true(all(inf$impact$ci_low <= alpha & alpha <= inf$impact$ci_high))
  true_eff <- unlist(lapply(1:2, function(c) .full_category_effects(stats::setNames(gam[, c], colnames(des$W)), des)))
  ce <- inf$category_effects
  expect_gte(sum(ce$ci_low <= true_eff & true_eff <= ce$ci_high), 6)
  tz <- as.vector(eps) == 0
  expect_lte(sum(inf$residual$flag & tz), 2)
  expect_gt(mean(sign(inf$residual$est[!tz]) == sign(as.vector(eps)[!tz])), 0.8)
})
