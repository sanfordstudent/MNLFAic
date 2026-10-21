# Identification condition, tolerance table, design report. Ported from
# chen-item-cov-ext/tests/test_stage2.R, sections 6-8 (logit metric pinned
# where written for it).

toy_pos <- function() .build_item_design(data.frame(
  pos = factor(c(rep("noun", 10), rep("adj", 6), rep("verb", 4)), levels = c("noun", "adj", "verb"))),
  categorical = "pos", omit = list(pos = "noun"))$W

test_that("condition check reproduces Chen et al.'s Proposition 1 cases", {
  a1 <- rep(1, 10)
  chk <- function(...) .ml1_condition_check(a1, NULL, ..., metric = "logit")
  expect_true(chk(support = 9:10, signs = 1)$holds)                 # Chen et al., Fig. 2
  expect_true(chk(support = 7:10, signs = 1)$holds)                 # 4 of 10 same sign
  expect_false(chk(support = 6:10, signs = 1)$holds)                # 5 of 10: fails
  expect_true(chk(support = 5:10, signs = c(1, 1, 1, -1, -1, -1))$holds)  # balanced
  expect_equal(chk(support = 7:10, signs = 1)$rho, 4 / 6, tolerance = 1e-4)
  expect_equal(chk(support = integer(0), signs = 1)$rho, 0)
  expect_error(chk(support = 1:2, signs = c(1, 2)), "signs")
  r <- .ml1_condition_check(a1[1:3], cbind(w = c(1, -1, 0)), support = 1:2, signs = 1)
  expect_false(r$holds); expect_identical(r$rho, Inf)
})

test_that("condition check agrees with exact noise-free recovery (random cases)", {
  n_case <- if (identical(Sys.getenv("NOT_CRAN"), "true")) 300 else 60
  set.seed(4)
  Wpos <- toy_pos()
  agree <- logical(n_case); is_tie <- logical(n_case); tie_says_fail <- TRUE
  for (k in seq_len(n_case)) {
    Jk <- 20; lamk <- stats::runif(Jk, 0.6, 2)
    Wk <- if (k %% 2 == 0) Wpos else cbind(aoa = stats::rnorm(Jk))
    S  <- sample(Jk, sample(2:9, 1)); sg <- sample(c(-1, 1), length(S), replace = TRUE)
    epsk <- numeric(Jk); epsk[S] <- sg * stats::runif(length(S), 0.3, 1.5)
    nu_k <- as.vector(Wk %*% stats::rnorm(ncol(Wk))) + epsk + lamk * stats::rnorm(1)
    rec <- .ml1_itemcov(cbind(nu_k), lamk, W = Wk, metric = "logit")
    recovered <- isTRUE(all.equal(unname(rec$eps[, 1]), epsk, tolerance = 1e-7, scale = 1))
    cc <- .ml1_condition_check(lamk, Wk, S, sg, metric = "logit")
    is_tie[k] <- abs(cc$rho - 1) < 1e-3
    if (is_tie[k]) tie_says_fail <- tie_says_fail && !cc$holds
    agree[k] <- (cc$holds == recovered)
  }
  expect_true(all(agree[!is_tie]))
  expect_true(tie_says_fail)                        # exact ties count as failing
  expect_true(all(which(is_tie) %% 2 == 0))         # ties only with the categorical item covariate
})

test_that("tolerance table", {
  Wpos <- toy_pos()
  tt <- .ml1_tolerance_table(rep(1, 20), Wpos, n_random = 50, metric = "logit",
                             groups = list(noun = 1:10, adj = 11:16, verb = 17:20))
  expect_identical(tt$by_lambda, c(4, 2, 1))
  expect_identical(tt$random_median, c(4, 2, 1))
  set.seed(1); lam <- stats::runif(20, 0.8, 1.6)
  tt2 <- .ml1_tolerance_table(lam, Wpos, n_random = 50, metric = "logit",
                              groups = list(noun = 1:10, adj = 11:16, verb = 17:20))
  verbs <- (17:20)[order(lam[17:20], decreasing = TRUE)]
  k_v <- tt2$by_lambda[3]
  expect_lt(k_v, 3)
  if (k_v > 0) expect_true(.ml1_condition_check(lam, Wpos, verbs[seq_len(k_v)], 1, metric = "logit", margin = FALSE)$holds)
  expect_false(.ml1_condition_check(lam, Wpos, verbs[seq_len(k_v + 1)], 1, metric = "logit", margin = FALSE)$holds)
  set.seed(7); u1 <- stats::runif(1); set.seed(7)
  invisible(.ml1_tolerance_table(lam, Wpos, n_random = 5, groups = list(verb = 17:20)))
  expect_identical(stats::runif(1), u1)             # caller's RNG stream unchanged
})

test_that("design report", {
  set.seed(1); lam <- stats::runif(20, 0.8, 1.6)
  set.seed(2); aoa <- stats::rnorm(20, 6, 1.5)
  rep8 <- .ml1_design_report(lam, cbind(toy_pos(), aoa = aoa - mean(aoa)), metric = "logit")
  expect_length(rep8$cor_slope_with_W, 3)
  expect_true(all(abs(rep8$cor_slope_with_W) <= 1))
  expect_true(all(rep8$leverage_top$leverage <= 1))
  expect_equal(rep8$leverage_mean, 4 / 20, tolerance = 1e-12)
  expect_null(.ml1_design_report(lam)$cor_slope_with_W)
})
