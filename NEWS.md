# MNLFAic 0.0.0.9000

* First version: anchor-free inference for uniform (i.e., difficulty, or
  optionally, intercept) DIF in moderated nonlinear factor analysis with
  person and item covariates (`mnlfa_ic()`), with summaries, item-level
  tables, identification diagnostics, simulation from chosen values or from a
  fit, the stage functions for advanced use, and two vignettes.
* The cluster-robust covariance (`cluster =`) is experimental.
* `diagnose()` and `summary()` no longer report an identification margin rho
  computed from estimated residual-DIF flags: in simulations it did not
  indicate when inference was unreliable, and the sparsity assumption cannot
  be tested from the data. `diagnose()` reports the tolerance table and the
  Stage-2 design report, which depend only on the slopes and item covariates;
  its `flags` argument is removed. `ml1_condition_check()` still checks a
  hypothesized pattern of residual DIF.
