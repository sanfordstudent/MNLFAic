# Sampling weights in Stage 1.

weights_setup <- function() {
  set.seed(81)
  J <- 10; N <- 600
  x <- cbind(g1 = stats::rbinom(N, 1, 0.5), z = as.vector(scale(stats::rnorm(N))))
  eps <- matrix(0, J, 2); eps[2, 1] <- 0.6; eps[5, 2] <- -0.4
  sim <- .simulate_itemcov_data(x, stats::runif(J, 0.8, 2), stats::rnorm(J, 0, 0.8), NULL, NULL, eps,
                                c(0.3, 0.4), c(0, 0.1), seed = 82)
  dat <- sim$dat; dat[sample(length(dat), 300)] <- NA
  st <- list(a = rep(1, J), b = stats::qlogis(pmin(pmax(colMeans(dat, na.rm = TRUE), .01), .99)),
             g = matrix(0, J, 2), beta = c(0, 0), eta = c(0, 0))
  list(J = J, N = N, x = x, dat = dat, st = st, w = sample(1:3, N, TRUE))
}

test_that("weights of 1 reproduce the unweighted fit exactly", {
  skip_on_cran()
  s <- weights_setup()
  f0 <- .fit_stage1(s$dat, s$x, s$st, anchor = 1)
  f1 <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, weights = rep(1, s$N))
  for (nm in c("beta", "eta", "alpha.vec", "beta.vec", "gamma.mat", "info", "vcov", "loglik"))
    expect_identical(f1[[nm]], f0[[nm]], label = nm)
})

test_that("integer weights equal duplicated rows; weighted score is zero at the estimate", {
  skip_on_cran()
  s <- weights_setup()
  fw <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, weights = s$w)
  rows <- rep(seq_len(s$N), s$w)
  fd <- .fit_stage1(s$dat[rows, ], s$x[rows, ], s$st, anchor = 1)
  expect_equal(fw$gamma.mat, fd$gamma.mat, tolerance = 1e-6)
  expect_equal(c(fw$beta, fw$eta), c(fd$beta, fd$eta), tolerance = 1e-6)
  expect_equal(fw$alpha.vec, fd$alpha.vec, tolerance = 1e-6)
  expect_equal(fw$loglik, fd$loglik, tolerance = 1e-8)
  expect_equal(fw$info, fd$info, tolerance = 1e-5)
  S  <- .stage1_person_scores(fw, s$dat, s$x)
  sc <- .score_stage1(fw$alpha.vec, fw$beta.vec, fw$gamma.mat, s$x, s$dat, fw$beta, fw$eta, fw$post, fw$nodes, w = s$w)
  expect_equal(colSums(s$w * S), sc, tolerance = 1e-8)
  expect_lt(max(abs(sc) / sqrt(colSums((s$w * S)^2))), 1e-3)
  # the anchor does not change the weighted fit
  fa <- .fit_stage1(s$dat, s$x, s$st, anchor = 4, weights = s$w)
  expect_equal(fa$alpha.vec, fw$alpha.vec, tolerance = 1e-5)
})

test_that("weight input checks", {
  s <- weights_setup()
  expect_error(.fit_stage1(s$dat, s$x, s$st, weights = s$w[-1]), "one non-negative")
  expect_error(.fit_stage1(s$dat, s$x, s$st, weights = replace(s$w, 1, NA)), "one non-negative")
  expect_error(.fit_stage1(s$dat, s$x, s$st, weights = replace(s$w, 1, -1)), "one non-negative")
})

# JK2 replicates with both halves of every zone: zone h = schools 2h - 1, 2h.
jk2_full <- function(w, school) {
  H <- max(school) / 2
  R <- matrix(w, length(w), 2 * H)
  for (h in seq_len(H)) {
    a <- school == 2 * h - 1; b <- school == 2 * h
    R[a, 2 * h - 1] <- 2 * w[a]; R[b, 2 * h - 1] <- 0
    R[a, 2 * h]     <- 0;        R[b, 2 * h]     <- 2 * w[b]
  }
  R
}

test_that("replicate contributions: JK2 pairs collapse exactly; constants", {
  set.seed(5)
  N <- 40; D <- 4; school <- rep(1:8, each = 5); w <- stats::runif(N, 1, 3)
  S <- matrix(stats::rnorm(N * D), N, D)
  A <- matrix(stats::rnorm(D * D), D, D); info <- crossprod(A) + diag(D)
  Rw <- cbind(jk2_full(w, school), w)                   # plus one unused replicate
  rp <- .replicate_contributions(S, Rw - w, info)
  expect_identical(c(rp$n_replicates, rp$n_units, rp$n_paired), c(9L, 4L, 4L))
  d <- t(solve(info) %*% crossprod(S, Rw - w))
  expect_equal(0.5 * crossprod(rp$B), 0.5 * crossprod(d), tolerance = 1e-10)
  expect_equal(.replicate_scale("JK2", 150, 0.5), 0.5)
  expect_equal(.replicate_scale("JK2-half", 75, 0.5), 1)
  expect_equal(.replicate_scale("BRR", 80, 0.5), 1 / 20)
  # wild replicate draws: exact covariance over all sign patterns
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), rp$n_units)))
  Z <- .wild_replicate_draws(rp$B, 0.5, nrow(signs), weights = signs)
  expect_equal(crossprod(Z) / nrow(signs), 0.5 * crossprod(rp$B), tolerance = 1e-12)
})

test_that("replicate covariance: one-step agrees with refitting every replicate", {
  skip_on_cran()
  s <- weights_setup()
  school <- rep(1:20, each = s$N / 20)                  # 10 zones
  Rw <- jk2_full(s$w, school)
  fr <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, weights = s$w, replicates = Rw)
  expect_identical(fr$vcov_type, "replicate")
  expect_identical(c(fr$replicate$n_units, fr$t_df), c(10L, 10L))
  refit <- sapply(seq_len(ncol(Rw)), function(r) {
    f <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, weights = Rw[, r])
    .pack_stage1(f$beta, f$eta, f$alpha.vec, f$beta.vec, f$gamma.mat)
  })
  est <- .pack_stage1(fr$beta, fr$eta, fr$alpha.vec, fr$beta.vec, fr$gamma.mat)
  V_refit <- 0.5 * tcrossprod(refit - est)
  r <- sqrt(diag(fr$vcov) / diag(V_refit))
  expect_lt(abs(stats::median(r) - 1), 0.05)
  expect_gt(min(r), 0.8); expect_lt(max(r), 1.25)
  # weighted CR1 with weights of 1 is the unweighted CR1
  f1 <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, cluster = school, weights = rep(1, s$N))
  f0 <- .fit_stage1(s$dat, s$x, s$st, anchor = 1, cluster = school)
  expect_identical(f1$vcov, f0$vcov)
  expect_error(.fit_stage1(s$dat, s$x, s$st, replicates = Rw), "full-sample")
  expect_error(.fit_stage1(s$dat, s$x, s$st, weights = s$w, cluster = school, replicates = Rw), "not both")
})

test_that("Stage 3 with replicate weights: sandwich and wild draws, method change", {
  skip_on_cran()
  s <- weights_setup()
  school <- rep(1:20, each = s$N / 20)
  fr <- .fit_stage1(s$dat, s$x, s$st, weights = s$w, replicates = jk2_full(s$w, school))
  is <- .dif_inference_itemcov(fr, M = 2000, seed = 3)
  iw <- .dif_inference_itemcov(fr, M = 2000, seed = 3, cluster_method = "wild")
  expect_true(is$settings$vcov_type == "replicate" && is$settings$draw_df == 10 && is$settings$n_units == 10)
  expect_true(iw$settings$wild_weights == "webb" && is.na(iw$settings$draw_df))      # 2^10 < M
  expect_lt(max(abs(iw$variance$se / sqrt(diag(fr$vcov))[.stage1_index(s$J, 2)$eta] - 1)), 0.1)
  fh <- .set_replicate_vcov(fr, "JK2-half")
  expect_equal(fh$vcov, 2 * fr$vcov)
  s2 <- .restore_stage1_fit(.compact_stage1_fit(fr))
  expect_equal(s2$vcov, fr$vcov, tolerance = 1e-10)
  expect_identical(s2$t_df, fr$t_df)
})

test_that("ilsa_replicates(): JK2 and JK2-half construction", {
  w <- c(10, 12, 9, 11, 8, 10)
  R <- ilsa_replicates(w, jkzone = c(1, 1, 1, 2, 2, 2), jkrep = c(0, 1, 1, 0, 1, 0))
  expect_equal(dim(R), c(6, 4))
  expect_equal(R[, 1], c(0, 24, 18, 11, 8, 10))
  expect_equal(R[, 2], c(20, 0, 0, 11, 8, 10))
  expect_equal(R[, 3], c(10, 12, 9, 0, 16, 0))
  expect_equal(R[, 4], c(10, 12, 9, 22, 0, 20))
  Rh <- ilsa_replicates(w, c(1, 1, 1, 2, 2, 2), c(0, 1, 1, 0, 1, 0), method = "JK2-half", n_zones = 3)
  expect_equal(Rh, unname(cbind(R[, 1], R[, 3], w)))
  expect_error(ilsa_replicates(w, c(1, 1, 1, 2, 2, 2), c(0, 2, 1, 0, 1, 0)), "0 or 1")
})

ilsa_data <- function(n = c(320, 480), seed = 1) {
  set.seed(seed)
  J <- 10; N <- sum(n)
  pd <- data.frame(country = factor(rep(c("A", "B"), n)), school = rep(seq_len(N / 20), each = 20))
  pd$wt <- stats::runif(N, 1, 4) * ifelse(pd$country == "A", 50, 5)        # very different country totals
  pd$jkzone <- (pd$school + 1) %/% 2
  pd$jkzone <- pd$jkzone - ifelse(pd$country == "B", min(pd$jkzone[pd$country == "B"]) - 1, 0)   # zones within country
  pd$jkrep <- (pd$school + 1) %% 2
  res <- matrix(0, J, 1); res[3, 1] <- 0.6
  sim <- simulate_mnlfa_ic(pd, ~ country, slopes = stats::runif(J, .8, 2), intercepts = stats::rnorm(J, 0, .7),
                           impact = 0.3, residual_dif = res, seed = seed + 1)
  list(pd = pd, resp = sim$responses, J = J, N = N)
}
ctl_w <- function(...) mnlfa_control(M = 200, n_pilot = 100, seed = 1, ...)

test_that("mnlfa_ic() with weights: scaling, labels, robust covariance, checks", {
  skip_on_cran()
  d <- ilsa_data()
  f0 <- mnlfa_ic(d$resp, d$pd, ~ country, control = ctl_w())
  f1 <- mnlfa_ic(d$resp, d$pd, ~ country, weights = rep(7, d$N), control = ctl_w())   # constant weights
  expect_equal(f1$stage1$gamma.mat, f0$stage1$gamma.mat, tolerance = 1e-10)
  expect_identical(f1$stage1$vcov_type, "cluster")
  expect_true(f1$stage1$cluster_is_person)
  expect_output(print(f1), "robust")
  expect_output(print(f1), "sampling weights")
  fs <- mnlfa_ic(d$resp, d$pd, ~ country, weights = "wt", weight_by = "country", control = ctl_w())
  fh <- mnlfa_ic(d$resp, d$pd, ~ country, weights = "wt", weight_by = "country", control = ctl_w(weight_scale = "house"))
  ws <- tapply(fs$data$weights, d$pd$country, sum); wh <- tapply(fh$data$weights, d$pd$country, sum)
  expect_equal(as.vector(ws), c(400, 400)); expect_equal(as.vector(wh), c(320, 480))
  expect_gt(max(abs(fs$stage1$alpha.vec - fh$stage1$alpha.vec)), 1e-4)   # pooled slopes depend on the scaling
  expect_output(print(fs), "senate scaling within 2 groups")
  expect_error(mnlfa_ic(d$resp, d$pd, ~ country, weights = "wt", control = ctl_w(vcov = "model")), "not valid")
  expect_error(infer(fs, vcov = "model"), "not valid")
  expect_error(infer(fs, weight_scale = "house"), "new fit")
  expect_error(mnlfa_ic(d$resp, d$pd, ~ country, weight_by = "country"), "need `weights`")
  r2 <- d$resp; r2[d$pd$country == "B", 4] <- NA                    # item 4 not administered in B
  expect_error(mnlfa_ic(r2, d$pd, ~ country, control = ctl_w()), "cannot be estimated")
})

test_that("mnlfa_ic() with replicate weights: JK2 by country, method change, wild", {
  skip_on_cran()
  d <- ilsa_data()
  Rw <- ilsa_replicates(d$pd$wt, d$pd$jkzone, d$pd$jkrep)
  fr <- mnlfa_ic(d$resp, d$pd, ~ country, weights = "wt", weight_by = "country", replicate_weights = Rw,
                 control = ctl_w())
  Z <- max(d$pd$jkzone)
  expect_identical(fr$stage1$vcov_type, "replicate")
  expect_identical(fr$inference$settings$vcov_type, "replicate")
  expect_identical(fr$stage1$replicate$n_units, as.integer(Z))
  expect_output(print(fr), sprintf("replicate weights, JK2 (%d zones)", Z), fixed = TRUE)
  sl <- confint(fr, "slopes")
  expect_equal(sl$ci_high - sl$est, stats::qt(.975, Z) * sl$se)
  fh <- infer(fr, replicate_method = "JK2-half")
  expect_equal(confint(fh, "slopes")$se, sqrt(2) * sl$se)
  expect_equal(unname(vcov(fh)), unname(2 * vcov(fr)))
  fw <- infer(fr, cluster_method = "wild")
  expect_output(print(fw), "JK2, wild bootstrap")
  expect_identical(fw$inference$settings$wild_weights, "rademacher")
  # replicate weights by column name; with the same scaling as the full weights
  pd2 <- cbind(d$pd, as.data.frame(Rw)); nm <- paste0("V", seq_len(ncol(Rw)))
  fc <- mnlfa_ic(d$resp, pd2, ~ country, weights = "wt", weight_by = "country", replicate_weights = nm,
                 inference = FALSE, control = ctl_w())
  expect_equal(fc$stage1$vcov, fr$stage1$vcov, tolerance = 1e-10)
  expect_error(mnlfa_ic(d$resp, d$pd, ~ country, cluster = "school", weights = "wt", replicate_weights = Rw), "not both")
})
