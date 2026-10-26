#' MNLFAic: anchor-free DIF inference in moderated nonlinear factor analysis with item covariates (MNLFA-IC)
#'
#' Estimation and inference for uniform (i.e., difficulty or intercept) differential item
#' functioning with respect to an arbitrary number of continuous or
#' categorical person covariates, in moderated nonlinear factor analysis with
#' item covariates. Person covariates affect the mean and variance of the
#' latent trait (impact) and each item's difficulty (DIF), and each item's DIF
#' is split into a part explained by item covariates and item-specific
#' residual DIF. The model is identified without pre-specified anchor items by
#' a sparsity assumption on residual DIF (a generalization of the minimal-L1
#' condition of Chen et al., 2023).
#'
#' The main function is [mnlfa_ic()]. Its results:
#' * [summary.mnlfa_ic()], [dif_table()] and [diagnose()]: the main tables,
#'   item-level results and identification diagnostics;
#' * [coef.mnlfa_ic()] and [confint.mnlfa_ic()]: estimates and confidence
#'   intervals, one type of parameter at a time;
#' * [infer()]: rerunning the bootstrap with other settings;
#' * [simulate_mnlfa_ic()] and [simulate.mnlfa_ic()]: simulating data from
#'   chosen values or from a fit;
#' * the [stage functions][stage-functions]: the three stages and the
#'   identification diagnostics for advanced use.
#'
#' `vignette("MNLFAic")` works through an example, and `vignette("method")`
#' describes the model and the procedure.
#'
#' The package is under development; the user-facing interface is not final.
#'
#' @references Chen, Y., Li, C., Ouyang, J., & Xu, G. (2023). DIF statistical
#'   inference without knowing anchoring items. *Psychometrika, 88*, 1097-1122.
#'   \doi{10.1007/s11336-023-09930-9}
#' @keywords internal
#' @importFrom stats cor cov dnorm median p.adjust pchisq plogis qlogis
#'   quantile rbinom rchisq rnorm sd setNames coef confint vcov simulate
"_PACKAGE"
