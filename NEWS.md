# MNLFAic (development version)

* Experimental: sampling weights and replicate-weight variance, for
  international assessments such as TIMSS, PIRLS and PISA. New arguments of
  `mnlfa_ic()`: `weights`, `weight_by` (groups, usually countries, within
  which the weights are rescaled; `mnlfa_control(weight_scale = )` chooses
  senate (default) or house scaling) and `replicate_weights`
  (`mnlfa_control(replicate_method = c("JK2", "JK2-half", "BRR"), fay = )`).
  Stage 1 maximizes the weighted log-likelihood; the replicate covariance
  uses a one-step approximation, so replicates are not refitted. Without
  replicate weights, the covariance is weighted cluster-robust (with
  `cluster =`) or weighted robust. `ilsa_replicates()` builds TIMSS/PIRLS
  jackknife replicate weights from the zone variables. The wild bootstrap
  (`cluster_method = "wild"`) also works with replicate weights.
* `mnlfa_ic()` now stops when an item has no responses among persons with a
  nonzero value of a person covariate (e.g., an item not administered in a
  country), whose DIF cannot be estimated.

* Experimental: score wild cluster bootstrap for clustered data,
  `mnlfa_control(cluster_method = "wild")`, also available without refitting
  through `infer(fit, cluster_method = "wild")`. Each draw combines the
  clusters' Stage-1 score totals with random weights (Rademacher, or Webb's
  six-point weights with few clusters), so the draws have the cluster-robust
  covariance without assuming its shape, need no refits, and are identical
  across platforms for a given seed. With it, the intervals for the variance
  impact and the item slopes also come from the draws. The default is
  unchanged (`"sandwich"`).

* `confint(fit, "slopes")` gives Wald standard errors and intervals for the
  item slopes, from the Stage-1 covariance (cluster-robust with `cluster =`).
  Like the variance impact, the slopes are not changed by the identifying
  step, so they do not need the bootstrap, and any `level` can be used
  (except with the wild cluster bootstrap, below).
* `vcov()` now returns the covariance that the inference uses: the
  model-based one after `infer(fit, vcov = "model")` on a clustered fit,
  instead of always the cluster-robust one.

* `print()` and `summary()` now label the covariance that the inference
  actually used. Previously, a fit with `cluster =` was labeled
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
