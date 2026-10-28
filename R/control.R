#' Settings for mnlfa_ic()
#'
#' The settings for estimation (Stage 1), the identification step
#' (Stage 2) and the bootstrap inference (Stage 3). The defaults are
#' recommended for most analyses.
#'
#' @param form Parameterization of DIF: `"difficulty"` (default; DIF and
#'   item-covariate effects are shifts in item difficulty in latent-SD units,
#'   positive = harder) or `"slope-intercept"` (shifts in the logit intercept,
#'   positive = easier). With item covariates the two define different models
#'   (see Details).
#' @param item_omit Named list giving, for categorical item covariates, the
#'   level whose effect is implied (effect coding), e.g.
#'   `list(domain = "number")`. Default: the first level. The choice changes
#'   only which category's effect is reported as implied; all category effects,
#'   residual DIF and flags are the same.
#' @param item_standardize Standardize continuous item covariates over items
#'   (effects per SD of the item covariate). Default `TRUE`.
#' @param anchor Item whose DIF is fixed at 0 for bookkeeping in Stage 1:
#'   `"auto"` (default) or an item number. Results do not depend on it; it only
#'   affects numerical stability.
#' @param n_nodes Number of adaptive Gauss-Hermite quadrature nodes. This and
#'   the Stage-1 numerical settings that follow (`em_tol` through
#'   `fallback_em_iter`) rarely need to be changed.
#' @param em_tol,em_stall EM tolerance and stall ratio before switching to
#'   Newton steps.
#' @param newton_tol Convergence criterion: every Newton step smaller than this
#'   many standard errors.
#' @param max_newton,max_fallback,fallback_em_iter Newton iterations per round,
#'   and the EM fallback used when Newton fails.
#' @param M Number of bootstrap draws.
#' @param n_pilot Number of draws used to estimate the standard errors that set
#'   the threshold for threshold centering.
#' @param center Centering of the bootstrap: `"threshold"` (default; threshold
#'   centering, recommended), `"estimate"` (centered at the point estimates, the
#'   standard parametric bootstrap; conservative for item-covariate effects in
#'   simulations) or `"null"` (residual DIF centered at 0; intervals for
#'   item-covariate effects were too narrow in simulations).
#' @param threshold Threshold rule for threshold centering, in standard errors:
#'   `"sqrt2loglogN"` (default; \eqn{\sqrt{2 \log \log N}}{sqrt(2 log log N)}),
#'   `"sqrtlogN"` (\eqn{\sqrt{\log N}}{sqrt(log N)}), `"bonferroni"`
#'   (\eqn{\sqrt{2 \log(JK)}}{sqrt(2 log(JK))}), or a number.
#' @param level Confidence level of the intervals.
#' @param fdr Benjamini-Hochberg false discovery rate for flagging residual DIF.
#' @param ci_effects,ci_dif Interval shape, `"symmetric"` (default) or `"equal"`
#'   (equal-tailed): `ci_effects` for impact and item-covariate effects, `ci_dif`
#'   for residual and total DIF. A symmetric interval
#'   is the estimate plus or minus a quantile of the absolute bootstrap
#'   deviations, so it is centered on the estimate; an equal-tailed interval
#'   uses the lower and upper quantiles of the deviations separately (e.g.,
#'   2.5% and 97.5%), so it can be asymmetric around the estimate.
#' @param vcov Covariance used for the bootstrap draws: `"auto"` (cluster-robust
#'   if a cluster was given, otherwise model-based) or `"model"`. Used only with
#'   `cluster`, an experimental feature.
#' @param cluster_t With a cluster-robust covariance, draw from a multivariate
#'   t distribution with G - 1 degrees of freedom (default `TRUE`). Used only
#'   with `cluster`, an experimental feature.
#' @param seed Optional seed for the bootstrap draws. The caller's random-number
#'   stream is left unchanged. A seed reproduces results exactly on the same
#'   machine; across platforms the multivariate normal draws can differ at the
#'   level of bootstrap noise.
#' @param verbose Print progress messages.
#'
#' @details With item covariates, the difficulty and slope-intercept forms are
#'   not reparameterizations of each other: a common effect of an item covariate
#'   is a common shift in latent-SD units in one and a common shift in logits in
#'   the other.
#'
#' @return A list of class `mnlfa_control`.
#' @examples
#' mnlfa_control(M = 2000, seed = 1)
#' @export
mnlfa_control <- function(form = c("difficulty", "slope-intercept"), item_omit = NULL,
                          item_standardize = TRUE, anchor = "auto",
                          n_nodes = 21, em_tol = 0.005, em_stall = 0.5, newton_tol = 1e-4,
                          max_newton = 25, max_fallback = 3, fallback_em_iter = 50,
                          M = 10000, n_pilot = 2000, center = c("threshold", "estimate", "null"),
                          threshold = "sqrt2loglogN", level = 0.95, fdr = 0.05,
                          ci_effects = c("symmetric", "equal"), ci_dif = c("symmetric", "equal"),
                          vcov = c("auto", "model"), cluster_t = TRUE, seed = NULL, verbose = FALSE) {
  form <- match.arg(form); center <- match.arg(center)
  ci_effects <- match.arg(ci_effects); ci_dif <- match.arg(ci_dif); vcov <- match.arg(vcov)
  .check_count <- function(v, nm, min = 1) {
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < min || v != round(v))
      stop("`", nm, "` must be a whole number >= ", min, ".", call. = FALSE)
  }
  .check_prob <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v <= 0 || v >= 1)
      stop("`", nm, "` must be a number strictly between 0 and 1.", call. = FALSE)
  }
  .check_pos <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v <= 0)
      stop("`", nm, "` must be a positive number.", call. = FALSE)
  }
  .check_count(n_nodes, "n_nodes"); .check_count(max_newton, "max_newton")
  .check_count(max_fallback, "max_fallback", 0); .check_count(fallback_em_iter, "fallback_em_iter")
  .check_count(M, "M", 10); .check_count(n_pilot, "n_pilot", 10)
  .check_pos(em_tol, "em_tol"); .check_pos(newton_tol, "newton_tol")
  .check_prob(em_stall, "em_stall"); .check_prob(level, "level"); .check_prob(fdr, "fdr")
  if (!(identical(anchor, "auto") || (is.numeric(anchor) && length(anchor) == 1 && anchor >= 1 && anchor == round(anchor))))
    stop("`anchor` must be \"auto\" or an item number.", call. = FALSE)
  if (!(is.numeric(threshold) && length(threshold) == 1 && threshold > 0) &&
      !(is.character(threshold) && length(threshold) == 1 && threshold %in% c("sqrt2loglogN", "sqrtlogN", "bonferroni")))
    stop("`threshold` must be \"sqrt2loglogN\", \"sqrtlogN\", \"bonferroni\" or a positive number.", call. = FALSE)
  if (!is.null(item_omit) && (!is.list(item_omit) || is.null(names(item_omit))))
    stop("`item_omit` must be a named list, e.g. list(pos = \"noun\").", call. = FALSE)
  if (!is.logical(item_standardize) || length(item_standardize) != 1) stop("`item_standardize` must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(cluster_t) || length(cluster_t) != 1) stop("`cluster_t` must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1) stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) stop("`seed` must be NULL or a number.", call. = FALSE)
  structure(list(form = form, item_omit = item_omit, item_standardize = item_standardize, anchor = anchor,
                 n_nodes = n_nodes, em_tol = em_tol, em_stall = em_stall, newton_tol = newton_tol,
                 max_newton = max_newton, max_fallback = max_fallback, fallback_em_iter = fallback_em_iter,
                 M = M, n_pilot = n_pilot, center = center, threshold = threshold, level = level, fdr = fdr,
                 ci_effects = ci_effects, ci_dif = ci_dif, vcov = vcov, cluster_t = cluster_t, seed = seed,
                 verbose = verbose),
            class = "mnlfa_control")
}

# Engine metric name for a form.
.form_metric <- function(form) if (identical(form, "slope-intercept")) "logit" else "difficulty"

#' @export
print.mnlfa_control <- function(x, ...) {
  cat("mnlfa_ic settings\n")
  utils::str(unclass(x), give.head = FALSE, no.list = TRUE)
  invisible(x)
}
