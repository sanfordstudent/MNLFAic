#' aDIFic: anchor-free DIF inference with item covariates in MNLFA
#'
#' Estimation and inference for uniform (difficulty) differential item
#' functioning in moderated nonlinear factor analysis with item covariates,
#' identified without anchor items by a sparsity assumption on item-specific
#' residual DIF (a generalization of the minimal-L1 condition of Chen et al.,
#' 2023).
#'
#' The main function is [mnlfa_ic()]; see [summary.mnlfa_ic()], [dif_table()]
#' and [diagnose()] for its results. `vignette("aDIFic")` works through an
#' example, and `vignette("method")` describes the model and the procedure.
#'
#' The package is under development; the user-facing interface is not final.
#'
#' @references Chen, Y., Li, C., Ouyang, J., & Xu, G. (2023). DIF statistical
#'   inference without knowing anchoring items. *Psychometrika, 88*, 1097-1122.
#' @keywords internal
#' @importFrom stats cor cov dnorm median p.adjust pchisq plogis qlogis
#'   quantile rbinom rchisq rnorm sd setNames coef confint vcov simulate
"_PACKAGE"
