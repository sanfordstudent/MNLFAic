# User interface: mnlfa_ic(), mnlfa_control(), methods, simulate.

small_data <- function(N = 400, J = 10, seed = 1) {
  set.seed(seed)
  items <- data.frame(item = paste0("w", 1:J), pos = rep(c("noun", "verb", "adj", "noun", "noun"), length.out = J),
                      aoa = round(stats::rnorm(J, 7, 1.5), 1))
  pd <- data.frame(female = stats::rbinom(N, 1, .5), age = stats::rnorm(N), grp = factor(sample(c("a", "b", "c"), N, TRUE)))
  res <- matrix(0, J, 2); res[3, 1] <- 0.8
  sim <- simulate_mnlfa_ic(pd, ~ female + age, items, ~ aoa + pos, slopes = stats::runif(J, .8, 2),
                           intercepts = stats::rnorm(J, 0, .7), impact = c(0.2, 0.5),
                           effects = matrix(c(0, 0.1, -0.1, 0, 0, 0), 3, 2), residual_dif = res, seed = seed + 1)
  list(items = items, pd = pd, resp = sim$responses, sim = sim)
}
ctl_fast <- function(...) mnlfa_control(M = 200, n_pilot = 100, seed = 1, ...)

test_that("mnlfa_ic() reproduces the engine on the same data", {
  skip_on_cran()                    # the engine itself is covered by test-equivalence.R
  fx <- fixture("fixture_small"); inp <- fx$input
  pd <- as.data.frame(inp$x); it <- inp$items
  fit <- mnlfa_ic(inp$dat, pd, ~ female + age, it, ~ aoa + pos,
                  control = mnlfa_control(item_omit = list(pos = "noun"), M = 1000, seed = 11))
  ref <- fx$fit_auto
  # point estimates match the research code
  expect_equal(unname(fit$stage1$gamma.mat), unname(ref$gamma.mat), tolerance = 1e-6)
  expect_equal(unname(fit$coding$design$W), unname(inp$design$W), tolerance = 1e-14)
  expect_equal(fit$inference$impact$est, fx$stage3$threshold$impact$est, tolerance = 1e-6)
  expect_equal(fit$inference$category_effects$est, fx$stage3$threshold$category_effects$est, tolerance = 1e-6)
  expect_equal(fit$inference$residual$est, fx$stage3$threshold$residual$est, tolerance = 1e-6)
  # inference equals the engine on the same fit and settings (bootstrap draws
  # depend on the covariance to rounding, so the comparison uses the same fit)
  eng <- .dif_inference_itemcov(fit$stage1, fit$coding$design$W, fit$coding$design, M = 1000, seed = 11)
  expect_same_table(fit$inference$impact, eng$impact)
  expect_same_table(fit$inference$residual, eng$residual)
  expect_same_table(fit$inference$category_effects, eng$category_effects)
})

test_that("coding: factor person covariates, interactions, item covariates", {
  d <- small_data()
  fit <- mnlfa_ic(d$resp, d$pd, ~ female * age + grp, inference = FALSE, control = ctl_fast())
  expect_identical(fit$data$covariates, c("female", "age", "grpb", "grpc", "female:age"))
  expect_equal(unname(fit$data$x[, "female:age"]), d$pd$female * d$pd$age)
  expect_true(all(fit$data$x[, "grpb"] %in% 0:1))
  fit2 <- mnlfa_ic(d$resp, d$pd, ~ female, d$items, ~ aoa + pos, inference = FALSE, control = ctl_fast())
  W <- fit2$coding$design$W
  expect_identical(colnames(W), c("aoa", "pos_noun", "pos_verb"))     # first level ("adj") implied
  expect_equal(stats::sd(W[, "aoa"]), 1, tolerance = 1e-12)
  expect_null(fit2$inference)
})

test_that("input errors", {
  d <- small_data(); r <- d$resp
  expect_error(mnlfa_ic(r * 2, d$pd, ~ female), "only 0, 1")
  expect_error(mnlfa_ic(r[, 1:2], d$pd, ~ female), "At least 3 items")
  r1 <- r; r1[, 2] <- 1
  expect_error(mnlfa_ic(r1, d$pd, ~ female), "without variation")
  expect_error(mnlfa_ic(r, d$pd[-1, ], ~ female), "rows")
  expect_error(mnlfa_ic(r, d$pd, female ~ age), "one-sided formula")
  expect_error(mnlfa_ic(r, d$pd, ~ female, d$items[-1, ], ~ aoa), "rows")
  it2 <- d$items; it2$item <- rev(it2$item)
  expect_error(mnlfa_ic(r, d$pd, ~ female, it2, ~ aoa), "do not match")
  expect_error(mnlfa_ic(r, d$pd, ~ female, d$items, ~ aoa * pos), "Interactions")
  # numeric row names left over from subsetting are not item names
  it_sub <- d$items[c(10, 1:9), ][order(c(10, 1:9)), ]; it_sub$item <- NULL
  expect_false(is.null(.prepare_items(~ aoa, it_sub, colnames(r))))
  # responses without column names take item names from item_data
  r_un <- unname(r)
  p_un <- .prepare_all(r_un, d$pd, ~ female, d$items, ~ aoa, NULL, ctl_fast(), "fail")
  expect_identical(colnames(p_un$dat), d$items$item)
  expect_error(mnlfa_ic(r, d$pd, ~ female, NULL, ~ aoa), "item_data")
  pd2 <- d$pd; pd2$female[5] <- NA
  expect_error(mnlfa_ic(r, pd2, ~ female), "missing person covariates")
  pd3 <- d$pd; pd3$one <- 1
  expect_error(mnlfa_ic(r, pd3, ~ female + one), "no variation")
  pd4 <- d$pd; pd4$f2 <- 2 * pd4$female
  expect_error(mnlfa_ic(r, pd4, ~ female + f2), "collinear")
  it3 <- d$items; it3$aoa[2] <- NA
  expect_error(mnlfa_ic(r, d$pd, ~ female, it3, ~ aoa), "missing values")
  it4 <- d$items; it4$aoa2 <- 3 * it4$aoa
  expect_error(mnlfa_ic(r, d$pd, ~ female, it4, ~ aoa + aoa2), "collinear")
  expect_error(mnlfa_ic(r[, 1:3], d$pd, ~ female, d$items[1:3, ], ~ aoa + pos), "Too few items")
  expect_error(mnlfa_ic(r, d$pd, ~ female, d$items, ~ pos, control = mnlfa_control(item_omit = list(aoa = 1))), "item_omit")
  expect_error(mnlfa_ic(r, d$pd, ~ female, cluster = 1:3), "one value per person")
  expect_error(mnlfa_ic(r, d$pd, ~ female, control = list()), "mnlfa_control")
})

test_that("missing covariates with na_action = 'omit', persons without responses, warnings", {
  d <- small_data(); r <- d$resp
  pd2 <- d$pd; pd2$female[5] <- NA
  expect_message(f <- mnlfa_ic(r, pd2, ~ female, na_action = "omit", inference = FALSE, control = ctl_fast()), "Dropping 1")
  expect_equal(f$data$N, nrow(r) - 1)
  r2 <- r; r2[7, ] <- NA
  expect_message(f2 <- mnlfa_ic(r2, d$pd, ~ female, inference = FALSE, control = ctl_fast()), "no responses")
  expect_false(7 %in% f2$data$kept)
  # the warnings come from input preparation; tested there, without fitting
  pd5 <- d$pd; pd5$age10 <- pd5$age + 10
  expect_warning(.prepare_all(r, pd5, ~ female + age10, NULL, NULL, NULL, ctl_fast(), "fail"), "far from centered")
  r3 <- r; r3[d$pd$female == 1, 4] <- 1; r3[which(d$pd$female == 1)[1:2], 4] <- 0
  expect_warning(.prepare_all(r3, d$pd, ~ female, NULL, NULL, NULL, ctl_fast(), "fail"), "fewer than 5")
  expect_silent(.prepare_all(r, d$pd, ~ female + age, NULL, NULL, NULL, ctl_fast(), "fail"))
})

test_that("mnlfa_control() validates its arguments", {
  expect_s3_class(mnlfa_control(), "mnlfa_control")
  expect_error(mnlfa_control(M = 0), "M")
  expect_error(mnlfa_control(level = 1.2), "level")
  expect_error(mnlfa_control(fdr = 0), "fdr")
  expect_error(mnlfa_control(threshold = "other"), "threshold")
  expect_error(mnlfa_control(anchor = "first"), "anchor")
  expect_error(mnlfa_control(item_omit = "noun"), "named list")
  expect_error(mnlfa_control(center = "middle"))
  expect_output(print(mnlfa_control()), "settings")
})

test_that("methods: print, summary, coef, confint, vcov, dif_table, diagnose", {
  d <- small_data()
  fit <- mnlfa_ic(d$resp, d$pd, ~ female + age, d$items, ~ aoa + pos, control = ctl_fast())
  expect_output(print(fit), "Impact on the latent mean")
  s <- summary(fit)
  expect_s3_class(s, "summary.mnlfa_ic")
  expect_output(print(s), "Item-covariate effects")
  expect_equal(unname(coef(fit, "mean")), fit$inference$impact$est)
  expect_equal(coef(fit), coef(fit, "mean"))
  expect_equal(unname(coef(fit, "effects")), fit$inference$category_effects$est, tolerance = 1e-12)
  expect_equal(dim(coef(fit, "residual")), c(10, 2))
  expect_length(coef(fit, "slopes"), 10)
  ci <- confint(fit, "mean")
  expect_equal(ci$ci_low, fit$inference$impact$ci_low)
  expect_equal(confint(fit), ci)
  expect_error(confint(fit, "mean", level = 0.9), "infer")
  expect_equal(nrow(confint(fit, "variance", level = 0.9)), 2)
  V <- vcov(fit)
  expect_equal(dim(V), rep(2 * 2 + 2 * 10 + 9 * 2, 2))
  expect_true(all(c("impact:female", "slope:w1", "logvar:age") %in% rownames(V)))
  dt <- dif_table(fit)
  expect_equal(nrow(dt), 20)
  expect_equal(dt$total, dt$residual + dt$explained, tolerance = 1e-12)
  expect_true(all(c("aoa", "pos") %in% names(dt)))
  expect_equal(nrow(dif_table(fit, flagged_only = TRUE)), sum(fit$inference$residual$flag))
  dg <- diagnose(fit, n_random = 10)
  expect_equal(nrow(dg$margin), 2)
  expect_true(all(dg$margin$rho >= 0))
  expect_false(is.null(dg$tolerance))
  expect_true(all(dg$design$leverage_top$item %in% fit$data$items))
  expect_output(print(dg), "Identification margin")
  expect_equal(nrow(diagnose(fit, flags = "offplane", tolerance = FALSE)$margin), 2)
})

test_that("infer() reproduces the stored inference and changes settings", {
  d <- small_data()
  fit <- mnlfa_ic(d$resp, d$pd, ~ female + age, d$items, ~ pos, control = ctl_fast())
  again <- infer(fit)
  expect_identical(again$inference$residual, fit$inference$residual)
  f2 <- infer(fit, center = "estimate")
  expect_identical(f2$inference$settings$center, "estimate")
  expect_error(infer(fit, n_nodes = 31), "new fit")
  f3 <- infer(fit, form = "slope-intercept")
  expect_identical(f3$stage2$metric, "logit")
  f0 <- mnlfa_ic(d$resp, d$pd, ~ female, inference = FALSE, control = ctl_fast())
  expect_error(summary(f0), "No inference")
  expect_output(print(f0), "No inference yet")
  expect_s3_class(infer(f0)$inference$impact, "data.frame")
})

test_that("simulate(): reproducible, keeps missing pattern and RNG stream", {
  d <- small_data(); r <- d$resp; r[1:5, 2] <- NA
  fit <- mnlfa_ic(r, d$pd, ~ female + age, d$items, ~ pos, control = ctl_fast())
  s1 <- simulate(fit, nsim = 2, seed = 4); s2 <- simulate(fit, nsim = 2, seed = 4)
  expect_identical(s1, s2)
  expect_length(s1, 2)
  expect_true(all(is.na(s1[[1]][1:5, 2])))
  expect_identical(colnames(s1[[1]]), colnames(r))
  set.seed(8); u <- stats::runif(1); set.seed(8); invisible(simulate(fit, seed = 1)); expect_identical(stats::runif(1), u)
  tr <- attr(simulate(fit, residual = "none", seed = 1), "truth")
  expect_true(all(tr$residual_dif == 0))
})

test_that("simulate_mnlfa_ic() checks its inputs", {
  pd <- data.frame(g = rep(0:1, 50))
  expect_error(simulate_mnlfa_ic(pd, ~ g, slopes = rep(1, 5), intercepts = rep(0, 4), impact = 0), "same length")
  expect_error(simulate_mnlfa_ic(pd, ~ g, slopes = rep(1, 5), intercepts = rep(0, 5), impact = c(0, 1)), "one value per")
  it <- data.frame(pos = rep(c("a", "b"), length.out = 5))
  expect_error(simulate_mnlfa_ic(pd, ~ g, it, ~ pos, slopes = rep(1, 5), intercepts = rep(0, 5), impact = 0,
                                 effects = matrix(0, 2, 1)), "matrix")
  s <- simulate_mnlfa_ic(pd, ~ g, slopes = rep(1, 5), intercepts = rep(0, 5), impact = 0.5, seed = 1)
  expect_equal(dim(s$responses), c(100, 5))
})

test_that("cluster argument as a column name", {
  skip_on_cran()
  d <- small_data(N = 600)
  d$pd$school <- rep(1:30, each = 20)
  expect_warning(fit <- mnlfa_ic(d$resp, d$pd, ~ female, cluster = "school", control = ctl_fast()), NA)
  expect_identical(fit$stage1$vcov_type, "cluster")
  expect_equal(fit$inference$settings$draw_df, 29)
})
