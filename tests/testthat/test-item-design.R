# Ported from chen-item-cov-ext/tests/test_stage2.R, section 2.

items_toy <- function() data.frame(
  pos = factor(c(rep("noun", 10), rep("adj", 6), rep("verb", 4)), levels = c("noun", "adj", "verb")),
  aoa = seq(3, 9, length.out = 20))

test_that("effect coding and centering", {
  items <- items_toy()
  des <- .build_item_design(items, continuous = "aoa", categorical = "pos", omit = list(pos = "noun"))
  E_manual <- cbind(ifelse(items$pos == "adj", 1, ifelse(items$pos == "noun", -1, 0)),
                    ifelse(items$pos == "verb", 1, ifelse(items$pos == "noun", -1, 0)))
  expect_identical(colnames(des$W), c("aoa", "pos_adj", "pos_verb"))
  expect_equal(unname(des$W[, c("pos_adj", "pos_verb")]), E_manual)
  expect_equal(mean(des$W[, "aoa"]), 0, tolerance = 1e-12)
  ds <- .build_item_design(items, continuous = "aoa", standardize = TRUE)
  expect_equal(stats::sd(ds$W[, "aoa"]), 1, tolerance = 1e-12)
})

test_that("implied category effect and sum-to-zero", {
  des <- .build_item_design(items_toy(), continuous = "aoa", categorical = "pos", omit = list(pos = "noun"))
  fe <- .full_category_effects(c(aoa = 0.4, pos_adj = 0.6, pos_verb = 0.2), des)
  expect_equal(fe[["pos_noun"]], -0.8, tolerance = 1e-12)
  expect_equal(sum(fe[c("pos_noun", "pos_adj", "pos_verb")]), 0, tolerance = 1e-12)
  expect_equal(fe[["aoa"]], 0.4)
})

test_that("input errors are caught", {
  items <- items_toy()
  expect_error(.build_item_design(items), "No item covariates")
  expect_error(.build_item_design(items, continuous = "freq"), "not found")
  bad <- items; bad$aoa[3] <- NA
  expect_error(.build_item_design(bad, continuous = "aoa"), "missing values")
  expect_error(.build_item_design(data.frame(pos = rep("noun", 5)), categorical = "pos"), "fewer than 2 levels")
  expect_error(.build_item_design(items, categorical = "pos", omit = list(pos = "pronoun")), "not found")
  expect_error(.build_item_design(data.frame(aoa = rep(1, 5)), continuous = "aoa", standardize = TRUE), "constant")
  expect_error(.build_item_design(as.matrix(items), continuous = "aoa"), "data frame")
})
