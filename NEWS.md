# MNLFAic (development version)

* `confint(fit, "slopes")` gives Wald standard errors and intervals for the
  item slopes, from the Stage-1 covariance (cluster-robust with `cluster =`).
  Like the variance impact, the slopes are not changed by the identifying
  step, so they do not need the bootstrap, and any `level` can be used.
* `vcov()` now returns the covariance that the inference uses: the
  model-based one after `infer(fit, vcov = "model")` on a clustered fit,
  instead of always the cluster-robust one.

* `print()` and `summary()` now label the covariance that the inference
  actually used. Previously, a fit with `cluster =` was labelled
  "cluster-robust" even after `infer(fit, vcov = "model")`, which reruns the
  inference with the model-based covariance. The estimates, intervals and
  p-values were not affected.

# MNLFAic 0.1.1

* First CRAN submission. No changes to the package code: a test that
  reproduces bootstrap draws made on one platform is skipped on CRAN, since
  `MASS::mvrnorm()` gives platform-dependent draws.

# MNLFAic 0.1.0

* First release: anchor-free inference for uniform (i.e., difficulty, or
  optionally, intercept) DIF in moderated nonlinear factor analysis with
  person and item covariates (`mnlfa_ic()`), with summaries, item-level
  tables, simulation from chosen values or from a
  fitted model, diagnostics, stage functions for advanced use, and two vignettes.
* The cluster-robust covariance (`cluster =`) is experimental.
