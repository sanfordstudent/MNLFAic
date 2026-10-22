# Stage functions for advanced use (thin wrappers around the engine, with
# matrix inputs). Most users need only mnlfa_ic().

#' Stage functions for advanced use
#'
#' The three stages of [mnlfa_ic()] and the identification diagnostics, with
#' matrix inputs. Useful for simulation studies and for re-using fits.
#'
#' * `mnlfa_fit()`: Stage 1, anchored marginal maximum likelihood with all DIF
#'   free except one bookkeeping anchor. Returns the engine fit (slopes
#'   `alpha.vec`, intercepts `beta.vec`, DIF `gamma.mat`, impact `beta`,
#'   log-variance impact `eta`, information `info`, covariance `vcov`, ...).
#' * `ml1_transform()`: Stage 2, the weighted median regression that chooses
#'   the representative satisfying the sparsity assumption.
#' * `ml1_inference()`: Stage 3, the bootstrap.
#' * `ml1_condition_check()`: the identification condition and its margin rho
#'   for a given residual-DIF pattern.
#' * `ml1_tolerance_table()`: how many same-signed residual-DIF items each
#'   group of items can absorb.
#'
#' @param responses N x J matrix of 0/1/NA.
#' @param x N x K numeric matrix of person covariates (no intercept).
#' @param start Optional starting values (see [mnlfa_ic()]).
#' @param cluster Optional cluster identifier.
#' @param control Settings from [mnlfa_control()].
#' @param fit A fit from `mnlfa_fit()`.
#' @param dif J x K matrix of anchored DIF estimates (logits).
#' @param slopes Item slopes (length J).
#' @param impact Optional anchored impact (length K).
#' @param W `NULL`, a J x M item-covariate matrix, or a list of such matrices
#'   (one per person covariate).
#' @param design Optional item design (as stored in an `mnlfa_ic` fit,
#'   `fit$coding$design`) for implied category effects and joint tests.
#' @param form `"difficulty"` or `"slope-intercept"`.
#' @param support Items with residual DIF (logical vector or indices).
#' @param signs Signs (+1/-1) of their residual DIF.
#' @param margin Also compute the margin rho.
#' @param groups Named list of item-index vectors.
#' @param sign Sign of the added residual DIF (+1 or -1).
#' @param n_random,seed Random orders for the tolerance table.
#' @return See Details.
#' @name stage-functions
NULL

#' @rdname stage-functions
#' @export
mnlfa_fit <- function(responses, x, start = NULL, cluster = NULL, control = mnlfa_control()) {
  dat <- .prepare_responses(responses); attr(dat, "named") <- NULL; x <- as.matrix(x)
  if (nrow(x) != nrow(dat)) stop("`x` must have one row per person.", call. = FALSE)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  if (is.null(start)) start <- .default_start(dat, ncol(x))
  fit <- .fit_stage1(dat, x, start, anchor = control$anchor, cluster = cluster,
                     em_tol = control$em_tol, em_stall = control$em_stall, newton_tol = control$newton_tol,
                     max_newton = control$max_newton, n_nodes = control$n_nodes,
                     max_fallback = control$max_fallback, fallback_em_iter = control$fallback_em_iter,
                     verbose = control$verbose)
  colnames(fit$gamma.mat) <- colnames(x)
  fit$N <- nrow(dat)
  fit
}

#' @rdname stage-functions
#' @export
ml1_transform <- function(dif, slopes, impact = NULL, W = NULL, form = c("difficulty", "slope-intercept")) {
  .ml1_itemcov(dif, slopes, impact, W, .form_metric(match.arg(form)))
}

#' @rdname stage-functions
#' @export
ml1_inference <- function(fit, W = NULL, design = NULL, control = mnlfa_control()) {
  .dif_inference_itemcov(fit, W = W, design = design, metric = .form_metric(control$form),
                         center = control$center, M = control$M, sig = 1 - control$level, seed = control$seed,
                         threshold_z = control$threshold, n_pilot = min(control$n_pilot, control$M), fdr = control$fdr,
                         ci_shape_effects = control$ci_effects, ci_shape_dif = control$ci_dif,
                         vcov = control$vcov, cluster_t = control$cluster_t)
}

#' @rdname stage-functions
#' @export
ml1_condition_check <- function(slopes, W = NULL, support, signs = 1, form = c("difficulty", "slope-intercept"),
                                margin = TRUE) {
  .ml1_condition_check(slopes, W, support, signs, .form_metric(match.arg(form)), margin)
}

#' @rdname stage-functions
#' @export
ml1_tolerance_table <- function(slopes, W = NULL, groups, sign = 1, form = c("difficulty", "slope-intercept"),
                                n_random = 200, seed = 1) {
  .ml1_tolerance_table(slopes, W, groups, sign, .form_metric(match.arg(form)), n_random, seed)
}
