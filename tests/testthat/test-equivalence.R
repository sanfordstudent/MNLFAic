# Equivalence with the research code (chen-item-cov-ext/R/), via fixtures.
# Tolerances: deterministic steps (coding, Stage 2, Stage 3 given a fit,
# diagnostics, simulator) must agree to rounding; Stage-1 fits differ only
# through the Gauss-Hermite rule (base R here, fastGHQuad in the research
# code; they agree to about 1e-14), so estimates agree to far below their SEs.

test_that("simulator reproduces the research data exactly", {
  fx <- fixture("fixture_small")$input
  tr <- fx$truth
  sim <- .simulate_itemcov_data(fx$x, tr$lambda, tr$nu0, fx$design$W, tr$gamma, tr$eps, tr$alpha, tr$omega,
                                seed = fx$sim_seed)
  expect_identical(sim$dat, fx$sim_dat)
  expect_equal(sim$theta, fx$sim_theta, tolerance = 1e-14)
  fx2 <- fixture("fixture_missing_cluster")$input
  tr2 <- fx2$truth
  sim2 <- .simulate_itemcov_data(fx2$x, tr2$lambda, tr2$nu0, fx2$design$W, tr2$gamma, tr2$eps, tr2$alpha,
                                 tr2$omega, seed = fx2$sim_seed, cluster = fx2$cluster, icc = fx2$icc)
  expect_identical(sim2$dat, fx2$sim_dat)
})

test_that("simulator does not change the caller's random-number stream", {
  set.seed(5); a <- stats::runif(1)
  set.seed(5); fx <- fixture("fixture_small")$input; tr <- fx$truth
  invisible(.simulate_itemcov_data(fx$x, tr$lambda, tr$nu0, fx$design$W, tr$gamma, tr$eps, tr$alpha,
                                   tr$omega, seed = 1))
  expect_identical(stats::runif(1), a)
})

test_that("item-covariate coding matches", {
  fx <- fixture("fixture_small")$input
  d <- .build_item_design(fx$items, continuous = "aoa", categorical = "pos",
                          omit = list(pos = "noun"), standardize = TRUE)
  expect_equal(d$W, fx$design$W, tolerance = 1e-15)
  expect_identical(d$coding$pos$columns, fx$design$coding$pos$columns)
})

test_that("anchor choice and item-rest correlations match", {
  fx <- fixture("fixture_small")
  expect_identical(.choose_stage1_anchor(fx$input$dat, fx$input$x), fx$anchor_choice)
  expect_equal(.item_rest_correlation(fx$input$dat), fx$item_rest, tolerance = 1e-14)
})

test_that("Stage 1 reproduces the research fits (automatic anchor and anchor 1)", {
  fx <- fixture("fixture_small"); inp <- fx$input
  for (case in c("auto", "anchor1")) {
    ref <- if (case == "auto") fx$fit_auto else fx$fit_anchor1
    f <- suppressWarnings(.fit_stage1(inp$dat, inp$x, start_values(inp$dat, 2),
                                      anchor = if (case == "auto") "auto" else 1))
    se <- sqrt(diag(solve(ref$info)))
    ix <- .stage1_index(ncol(inp$dat), 2)
    expect_identical(f$anchor, ref$anchor)
    expect_equal(f$loglik, ref$loglik, tolerance = 1e-8)
    expect_lt(max(abs(f$alpha.vec - ref$alpha.vec) / se[ix$a]), 1e-4)
    expect_lt(max(abs(f$beta.vec - ref$beta.vec) / se[ix$b]), 1e-4)
    expect_lt(max(abs(c(f$beta, f$eta) - c(ref$beta, ref$eta)) / se[c(ix$beta, ix$eta)]), 1e-4)
    expect_lt(max(abs(f$gamma.mat - ref$gamma.mat)), 1e-4 * max(se))
    expect_equal(sqrt(diag(solve(f$info))), se, tolerance = 1e-4)
  }
})

test_that("per-person scores match and sum to the score", {
  fx <- fixture("fixture_small"); inp <- fx$input
  f <- suppressWarnings(.fit_stage1(inp$dat, inp$x, start_values(inp$dat, 2), anchor = 1))
  S <- .stage1_person_scores(f, inp$dat, inp$x)
  expect_equal(dim(S), dim(fx$person_scores))
  expect_lt(max(abs(S - fx$person_scores)), 1e-4)
  sc <- .score_stage1(f$alpha.vec, f$beta.vec, f$gamma.mat, inp$x, inp$dat, f$beta, f$eta, f$post, f$nodes)
  expect_equal(colSums(S), sc, tolerance = 1e-8)
})

test_that("Stage 2 matches exactly (difficulty, logit, no item covariates)", {
  fx <- fixture("fixture_small"); f <- fx$fit_auto; W <- fx$input$design$W
  for (case in c("difficulty", "logit", "noW")) {
    ref <- fx$stage2[[case]]
    r <- .ml1_itemcov(f$gamma.mat, f$alpha.vec, f$beta, if (case == "noW") NULL else W,
                      if (case == "logit") "logit" else "difficulty")
    expect_equal(r$delta, ref$delta, tolerance = 1e-12)
    expect_equal(r$gamma, ref$gamma, tolerance = 1e-12)
    expect_equal(r$eps, ref$eps, tolerance = 1e-12)
    expect_equal(r$nu_total, ref$nu_total, tolerance = 1e-12)
    expect_equal(r$alpha, ref$alpha, tolerance = 1e-12)
  }
})

test_that("Stage 3 matches exactly given the same fit and seed", {
  fx <- fixture("fixture_small"); f <- as_fit(fx$fit_auto, nrow(fx$input$dat))
  W <- fx$input$design$W; des <- fx$input$design
  args <- list(threshold = list(), estimate = list(center = "estimate"), null = list(center = "null"),
               logit = list(metric = "logit"))
  for (case in names(args)) {
    ref <- fx$stage3[[case]]
    # the research code's default DIF interval shape was equal-tailed
    r <- do.call(.dif_inference_itemcov, c(list(f, W, des, M = 1000, seed = 11, ci_shape_dif = "equal"), args[[case]]))
    for (tab in c("residual", "total", "gamma", "category_effects", "joint", "impact", "impact_by_category"))
      if (!is.null(ref[[tab]])) expect_same_table(r[[tab]], ref[[tab]])
    expect_equal(as.numeric(r$settings$threshold_value), as.numeric(ref$settings$threshold_value))
    expect_equal(r$settings$n_residuals_zeroed, ref$settings$n_residuals_zeroed)
  }
  ref <- fx$stage3$noW
  r <- .dif_inference_itemcov(f, NULL, M = 1000, seed = 11, ci_shape_dif = "equal")
  expect_same_table(r$residual, ref$residual)
  expect_same_table(r$impact, ref$impact)
})

test_that("Stage 3 does not change the caller's random-number stream", {
  fx <- fixture("fixture_small"); f <- as_fit(fx$fit_auto, nrow(fx$input$dat))
  set.seed(9); a <- stats::runif(1)
  set.seed(9); invisible(.dif_inference_itemcov(f, NULL, M = 50, seed = 1))
  expect_identical(stats::runif(1), a)
})

test_that("diagnostics match", {
  fx <- fixture("fixture_small"); inp <- fx$input; W <- inp$design$W
  for (c in 1:2) {
    r <- .ml1_condition_check(inp$truth$lambda, W, inp$truth$eps[, c] != 0, sign(inp$truth$eps[, c]))
    expect_identical(r$holds, fx$diagnostics$condition[[c]]$holds)
    expect_equal(r$rho, fx$diagnostics$condition[[c]]$rho, tolerance = 1e-10)
  }
  tol <- .ml1_tolerance_table(fx$fit_auto$alpha.vec, W, split(seq_len(20), inp$items$pos), n_random = 20, seed = 3)
  expect_equal(tol, fx$diagnostics$tolerance)
  dr <- .ml1_design_report(fx$fit_auto$alpha.vec, W)
  expect_equal(dr$condition_number, fx$diagnostics$design_report$condition_number, tolerance = 1e-12)
  expect_equal(dr$cor_slope_with_W, fx$diagnostics$design_report$cor_slope_with_W, tolerance = 1e-12)
})

test_that("missing responses and clustering: fit, sandwich and Stage 3 match", {
  fx <- fixture("fixture_missing_cluster"); inp <- fx$input
  f <- .fit_stage1(inp$dat, inp$x, start_values(inp$dat, 2), cluster = inp$cluster)
  ref <- fx$fit
  se <- sqrt(diag(solve(ref$info)))
  expect_identical(f$anchor, ref$anchor)
  expect_identical(f$n_clusters, ref$n_clusters)
  expect_lt(max(abs(f$gamma.mat - ref$gamma.mat)), 1e-4 * max(se))
  expect_equal(sqrt(diag(f$vcov)), sqrt(diag(ref$vcov)), tolerance = 1e-3)
  g <- as_fit(ref, nrow(inp$dat))
  for (case in c("cluster", "model", "normal")) {
    r <- .dif_inference_itemcov(g, inp$design$W, inp$design, M = 1000, seed = 22, ci_shape_dif = "equal",
                                vcov = if (case == "model") "model" else "auto",
                                cluster_t = case != "normal")
    rr <- fx[[paste0("stage3_", case)]]
    expect_same_table(r$impact, rr$impact)
    expect_same_table(r$residual, rr$residual)
    expect_identical(r$settings$draw_df, rr$settings$draw_df)
  }
})

test_that("compact storage round-trips (model-based and cluster-robust)", {
  for (nm in c("fixture_small", "fixture_missing_cluster")) {
    fx <- fixture(nm)
    s <- .restore_stage1_fit(fx$compact)
    ref <- if (nm == "fixture_small") fx$fit_auto else fx$fit
    expect_equal(s$info, ref$info, tolerance = 1e-14)
    expect_equal(s$vcov, ref$vcov, tolerance = 1e-10)
    expect_identical(s$vcov_type, ref$vcov_type)
  }
})

test_that("long test (J = 80, N = 1500) matches the research code", {
  skip_on_cran()
  fx <- fixture("fixture_long"); inp <- fx$input
  f <- .fit_stage1(inp$dat, inp$x, start_values(inp$dat, 2))
  ref <- fx$fit
  se <- sqrt(diag(solve(ref$info)))
  expect_identical(f$anchor, ref$anchor)
  expect_equal(f$loglik, ref$loglik, tolerance = 1e-8)
  expect_lt(max(abs(f$gamma.mat - ref$gamma.mat)), 1e-4 * max(se))
  r <- .dif_inference_itemcov(as_fit(ref, nrow(inp$dat)), inp$design$W, inp$design, M = 2000, seed = 33,
                              ci_shape_dif = "equal")
  expect_same_table(r$residual, fx$stage3$residual)
  expect_same_table(r$category_effects, fx$stage3$category_effects)
})
