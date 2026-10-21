#' Anchor-free DIF analysis with item covariates (MNLFA-IC)
#'
#' Fits a two-parameter logistic moderated nonlinear factor analysis in which
#' person covariates affect the latent mean (impact), the latent log-variance,
#' and every item's difficulty (uniform DIF), and in which each item's DIF is
#' decomposed into a part explained by item covariates and an item-specific
#' residual. No anchor items are needed: the model is identified by assuming
#' that residual DIF is sparse (the slope-weighted total residual DIF is
#' minimal, separately for each person covariate), generalizing Chen et al.
#' (2023). Inference uses a threshold-centered parametric bootstrap.
#'
#' The analysis has three stages: (1) marginal maximum likelihood of the model
#' with free DIF for every item and one bookkeeping anchor; (2) a weighted
#' median regression that picks the representative satisfying the sparsity
#' assumption; (3) a parametric bootstrap through Stage 2.
#'
#' @param responses N x J matrix or data frame of 0/1 responses (`NA` =
#'   missing; missing responses are handled by full-information maximum
#'   likelihood under MAR).
#' @param person_data Data frame with one row per person, in the order of
#'   `responses`.
#' @param person One-sided formula of person covariates, e.g. `~ female + age`
#'   or `~ treat * baseline`. Built without an intercept: numeric variables
#'   enter as they are, factors as treatment (0/1) dummies with the first level
#'   as reference. The person with all covariates at 0 is the reference
#'   (latent mean 0, variance 1); center continuous covariates yourself if
#'   needed. The sparsity assumption applies to each column separately, so the
#'   coding of person covariates is part of the identifying assumption.
#' @param item_data Data frame with one row per item, in the order of the
#'   columns of `responses`. If it has an `item` column (or row names), they
#'   must match the column names of `responses`.
#' @param item One-sided formula of item covariates, e.g. `~ aoa + pos`.
#'   Numeric variables are centered at their unweighted mean over items (and
#'   standardized by default); factors and character variables are
#'   effect-coded (see [mnlfa_control()]). `NULL` (default): no item
#'   covariates; the analysis is then Chen et al.'s anchor-free method for any
#'   number of person covariates.
#' @param cluster Optional cluster identifier (vector with one value per
#'   person, or the name of a column of `person_data`). Adds a cluster-robust
#'   covariance for inference. Experimental.
#' @param control Settings from [mnlfa_control()].
#' @param inference Run the bootstrap inference (Stage 3). If `FALSE`, only
#'   point estimates are computed; use [infer()] later.
#' @param start Optional starting values (advanced): a list with `a` (slopes),
#'   `b` (intercepts), `g` (J x K DIF), `beta` (impact), `eta` (log-variance
#'   impact). Default: slopes 1, intercepts from the item proportions correct,
#'   everything else 0.
#' @param na_action What to do with persons who have missing person covariates
#'   (or cluster): `"fail"` (default) or `"omit"`.
#'
#' @return An object of class `mnlfa_ic`; see [summary.mnlfa_ic()],
#'   [dif_table()], [coef.mnlfa_ic()], [diagnose()].
#' @references Chen, Y., Li, C., Ouyang, J., & Xu, G. (2023). DIF statistical
#'   inference without knowing anchoring items. *Psychometrika, 88*, 1097-1122.
#' @examples
#' # Small example (a real analysis would use more persons, items and draws)
#' sim <- simulate_mnlfa_ic(
#'   person_data = data.frame(female = rep(0:1, 150)), person = ~ female,
#'   item_data = data.frame(pos = rep(c("noun", "verb", "adj"), length.out = 9)),
#'   item = ~ pos, slopes = rep(1.2, 9), intercepts = seq(-1, 1, length.out = 9),
#'   impact = 0.3, seed = 1)
#' fit <- mnlfa_ic(sim$responses, sim$person_data, ~ female,
#'                 item_data = sim$item_data, item = ~ pos,
#'                 control = mnlfa_control(M = 200, n_pilot = 100, seed = 1))
#' fit
#' summary(fit)
#' @export
mnlfa_ic <- function(responses, person_data, person, item_data = NULL, item = NULL, cluster = NULL,
                     control = mnlfa_control(), inference = TRUE, start = NULL,
                     na_action = c("fail", "omit")) {
  cl_call <- match.call()
  na_action <- match.arg(na_action)
  if (!inherits(control, "mnlfa_control")) stop("`control` must come from mnlfa_control().", call. = FALSE)
  t0 <- proc.time()[3]
  prep <- .prepare_all(responses, person_data, person, item_data, item, cluster, control, na_action)
  dat <- prep$dat; x <- prep$x; J <- ncol(dat); K <- ncol(x)
  if (is.null(start)) start <- .default_start(dat, K)

  fit <- .fit_stage1(dat, x, start, anchor = control$anchor, cluster = prep$cluster,
                     em_tol = control$em_tol, em_stall = control$em_stall, newton_tol = control$newton_tol,
                     max_newton = control$max_newton, n_nodes = control$n_nodes,
                     max_fallback = control$max_fallback, fallback_em_iter = control$fallback_em_iter,
                     verbose = control$verbose)
  colnames(fit$gamma.mat) <- colnames(x); rownames(fit$gamma.mat) <- colnames(dat)
  names(fit$alpha.vec) <- names(fit$beta.vec) <- colnames(dat)
  names(fit$beta) <- names(fit$eta) <- colnames(x)
  fit$N <- nrow(dat)
  t1 <- proc.time()[3]

  obj <- structure(list(
    call = cl_call,
    data = list(N = nrow(dat), J = J, K = K, items = colnames(dat), covariates = colnames(x),
                responses = dat, x = x, cluster = prep$cluster, kept = prep$kept, N_input = prep$N_input),
    coding = list(person_terms = prep$person_terms, person_levels = prep$person_levels, design = prep$design,
                  form = control$form),
    stage1 = fit, stage2 = NULL, inference = NULL, control = control,
    timing = c(stage1 = unname(t1 - t0))), class = "mnlfa_ic")
  obj$stage2 <- .run_stage2(obj)
  if (inference) obj <- .run_stage3(obj)
  obj
}

.W_of <- function(object) if (is.null(object$coding$design)) NULL else object$coding$design$W

.run_stage2 <- function(object) {
  f <- object$stage1
  .ml1_itemcov(f$gamma.mat, f$alpha.vec, f$beta, .W_of(object), .form_metric(object$control$form))
}

.run_stage3 <- function(object) {
  ctl <- object$control; t0 <- proc.time()[3]
  if (ctl$verbose) message("Stage 3: ", ctl$M, " bootstrap draws")
  object$inference <- .dif_inference_itemcov(
    object$stage1, W = .W_of(object), design = object$coding$design, metric = .form_metric(ctl$form),
    center = ctl$center, M = ctl$M, sig = 1 - ctl$level, seed = ctl$seed, threshold_z = ctl$threshold,
    n_pilot = min(ctl$n_pilot, ctl$M), fdr = ctl$fdr, ci_shape_effects = ctl$ci_effects,
    ci_shape_dif = ctl$ci_dif, N = object$data$N, vcov = ctl$vcov, cluster_t = ctl$cluster_t)
  object$timing["stage3"] <- unname(proc.time()[3] - t0)
  object
}

#' Re-run the inference without refitting
#'
#' Runs Stage 3 (and Stage 2, if the parameterization changes) again with some
#' settings changed, reusing the Stage-1 fit.
#'
#' @param object An `mnlfa_ic` fit.
#' @param ... Arguments of [mnlfa_control()] to change, e.g. `M = 5000`,
#'   `center = "estimate"`, `seed = 2`. Stage-1 settings and item-covariate
#'   coding cannot be changed here.
#' @return The updated `mnlfa_ic` object.
#' @examples
#' \donttest{
#' sim <- simulate_mnlfa_ic(person_data = data.frame(g = rep(0:1, 200)), person = ~ g,
#'                          slopes = rep(1, 10), intercepts = seq(-1, 1, length.out = 10),
#'                          impact = 0.3, seed = 2)
#' fit <- mnlfa_ic(sim$responses, sim$person_data, ~ g, inference = FALSE)
#' fit <- infer(fit, M = 300, n_pilot = 100, seed = 1)
#' }
#' @export
infer <- function(object, ...) {
  if (!inherits(object, "mnlfa_ic")) stop("`object` must be an mnlfa_ic fit.", call. = FALSE)
  changes <- list(...)
  fixed <- c("item_omit", "item_standardize", "anchor", "n_nodes", "em_tol", "em_stall", "newton_tol",
             "max_newton", "max_fallback", "fallback_em_iter")
  if (any(names(changes) %in% fixed))
    stop("These settings need a new fit with mnlfa_ic(): ", paste(intersect(names(changes), fixed), collapse = ", "),
         ".", call. = FALSE)
  old <- object$control
  object$control <- do.call(mnlfa_control, utils::modifyList(unclass(old), changes))
  object$coding$form <- object$control$form
  if (!identical(old$form, object$control$form)) object$stage2 <- .run_stage2(object)
  .run_stage3(object)
}
