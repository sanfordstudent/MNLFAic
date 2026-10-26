#' Simulate data from an MNLFA-IC model
#'
#' Generates 0/1 responses from the model used by [mnlfa_ic()], with person
#' and item covariates coded by the same rules.
#'
#' @param person_data Data frame of person covariates (one row per person).
#' @param person One-sided formula of person covariates (as in [mnlfa_ic()]).
#' @param item_data Data frame of item covariates (one row per item), or
#'   `NULL`.
#' @param item One-sided formula of item covariates, or `NULL`.
#' @param slopes,intercepts Baseline item slopes and intercepts (length J; the
#'   item parameters of the person with all covariates at 0).
#' @param impact Impact on the latent mean (length K, in the order of the
#'   columns built from `person`).
#' @param variance_impact Impact on the latent log-variance (length K;
#'   default 0).
#' @param effects Item-covariate effects on DIF: a matrix with one row per
#'   coded item-covariate column and one column per person covariate (default
#'   0). In the difficulty form these are latent-SD shifts, + = harder. Rows
#'   follow the coded item-covariate columns: the item covariates in formula
#'   order, each categorical one expanded to one column per level except the
#'   omitted one (e.g. `reading`, `domain_data`, `domain_geometry`). If
#'   `effects` has row (or column) names, they must match the coded
#'   item-covariate columns (or the person covariates) and are used to order
#'   the rows (or columns).
#' @param residual_dif J x K matrix of item-specific DIF (default 0), in the
#'   same units.
#' @param form `"difficulty"` (default) or `"slope-intercept"`.
#' @param cluster,icc Optional clustering of persons: with `icc > 0`, a share
#'   `icc` of the latent residual variance is shared within clusters, keeping
#'   each person's marginal distribution unchanged. For use with the
#'   cluster-robust option of [mnlfa_ic()], which is currently an experimental
#'   feature.
#' @param item_omit,item_standardize Coding of item covariates, as in
#'   [mnlfa_control()].
#' @param seed Optional seed; the caller's random-number stream is left
#'   unchanged.
#' @return A list with `responses` (N x J matrix), `theta` (the persons'
#'   latent trait values), `person_data`, `item_data`, the coded person matrix
#'   `x` and item matrix `W`, and `truth`: the generating values, named as the
#'   arguments (`slopes`, `intercepts`, `impact`, `variance_impact`, `effects`,
#'   `residual_dif`, `form`).
#' @examples
#' sim <- simulate_mnlfa_ic(person_data = data.frame(female = rep(0:1, 100)),
#'                          person = ~ female, slopes = rep(1, 8),
#'                          intercepts = seq(-1, 1, length.out = 8), impact = 0.3, seed = 1)
#' dim(sim$responses)
#' @export
simulate_mnlfa_ic <- function(person_data, person, item_data = NULL, item = NULL, slopes, intercepts,
                              impact, variance_impact = NULL, effects = NULL, residual_dif = NULL,
                              form = c("difficulty", "slope-intercept"), cluster = NULL, icc = 0,
                              item_omit = NULL, item_standardize = TRUE, seed = NULL) {
  form <- match.arg(form)
  J <- length(slopes)
  if (length(intercepts) != J) stop("`slopes` and `intercepts` must have the same length.", call. = FALSE)
  if (any(!(slopes > 0))) stop("`slopes` must be positive.", call. = FALSE)
  pp <- .prepare_person(person, person_data, nrow(person_data))
  if (any(pp$na_rows)) stop("Person covariates contain missing values.", call. = FALSE)
  x <- suppressWarnings(.person_matrix(pp, rep(TRUE, nrow(person_data))))
  K <- ncol(x)
  item_names <- paste0("Item", seq_len(J))
  if (!is.null(item_data) && "item" %in% names(item_data)) item_names <- as.character(item_data$item)
  des <- .prepare_items(item, item_data, item_names, item_omit, item_standardize)
  W <- if (is.null(des)) NULL else des$W
  if (length(impact) != K) stop("`impact` must have one value per person covariate (", K, ": ",
                                paste(colnames(x), collapse = ", "), ").", call. = FALSE)
  if (is.null(variance_impact)) variance_impact <- rep(0, K)
  if (length(variance_impact) != K) stop("`variance_impact` must have length ", K, ".", call. = FALSE)
  M <- if (is.null(W)) 0 else ncol(W)
  if (is.null(effects)) effects <- matrix(0, M, K)
  effects <- as.matrix(effects)
  if (M > 0 && !all(dim(effects) == c(M, K)))
    stop("`effects` must be a ", M, " x ", K, " matrix (item-covariate columns: ", paste(colnames(W), collapse = ", "), ").", call. = FALSE)
  if (M > 0 && !is.null(rownames(effects))) {
    if (anyDuplicated(rownames(effects)) || !setequal(rownames(effects), colnames(W)))
      stop("The row names of `effects` must be the coded item-covariate columns: ", paste(colnames(W), collapse = ", "),
           ".", call. = FALSE)
    effects <- effects[colnames(W), , drop = FALSE]
  }
  if (M > 0 && !is.null(colnames(effects))) {
    if (anyDuplicated(colnames(effects)) || !setequal(colnames(effects), colnames(x)))
      stop("The column names of `effects` must be the person covariates: ", paste(colnames(x), collapse = ", "),
           ".", call. = FALSE)
    effects <- effects[, colnames(x), drop = FALSE]
  }
  if (is.null(residual_dif)) residual_dif <- matrix(0, J, K)
  residual_dif <- as.matrix(residual_dif)
  if (!all(dim(residual_dif) == c(J, K))) stop("`residual_dif` must be a ", J, " x ", K, " matrix.", call. = FALSE)
  cl <- .prepare_cluster(cluster, person_data)
  sim <- .simulate_itemcov_data(x, slopes, intercepts, W, if (M > 0) effects else NULL, residual_dif,
                                impact, variance_impact, metric = .form_metric(form), seed = seed,
                                cluster = cl, icc = icc)
  colnames(sim$dat) <- item_names
  truth <- list(slopes = stats::setNames(slopes, item_names), intercepts = stats::setNames(intercepts, item_names),
                impact = stats::setNames(impact, colnames(x)),
                variance_impact = stats::setNames(variance_impact, colnames(x)),
                effects = if (M > 0) `dimnames<-`(effects, list(colnames(W), colnames(x))) else NULL,
                residual_dif = `dimnames<-`(residual_dif, list(item_names, colnames(x))), form = form)
  list(responses = sim$dat, theta = sim$theta, person_data = person_data,
       item_data = if (is.null(item_data)) NULL else item_data, x = x, W = W, truth = truth)
}

#' Simulate data sets from a fitted mnlfa_ic model
#'
#' Generates response data from the fitted model: the estimated slopes,
#' intercepts, impact on the latent mean and log-variance, and item-covariate
#' effects, with residual DIF set as chosen. The persons' covariates and the
#' pattern of missing responses are those of the data. Persons are simulated
#' as independent (no clustering). Analyzing each data set again with
#' [mnlfa_ic()] and comparing the results with the truth shows how the
#' procedure performs in a design like yours.
#'
#' @param object An `mnlfa_ic` fit.
#' @param nsim Number of data sets.
#' @param seed Optional seed; the caller's random-number stream is left
#'   unchanged.
#' @param residual Residual DIF in the simulated truth: `"flagged"` (default;
#'   the estimates for cells flagged by Benjamini-Hochberg, 0 elsewhere),
#'   `"none"`, or `"estimated"` (the estimates everywhere).
#' @param ... Unused.
#' @return A list of `nsim` response matrices, with the simulated truth as
#'   attribute `"truth"`: a list with `slopes`, `intercepts`, `impact` (on the
#'   latent mean), `variance_impact` (on the latent log-variance), `effects`
#'   (item-covariate effects, or `NULL`), `residual_dif` and `form`, named as
#'   the arguments of [simulate_mnlfa_ic()].
#' @export
simulate.mnlfa_ic <- function(object, nsim = 1, seed = NULL, residual = c("flagged", "none", "estimated"), ...) {
  residual <- match.arg(residual)
  d <- object$data; f <- object$stage1; s2 <- object$stage2
  eps <- switch(residual,
                none = matrix(0, d$J, d$K),
                estimated = unname(s2$eps),
                flagged = { .need_inference(object); unname(s2$eps * matrix(object$inference$residual$flag, d$J, d$K)) })
  W <- .W_of(object)
  gam <- if (is.null(W)) NULL else sapply(s2$gamma, identity)
  if (!is.null(gam) && is.null(dim(gam))) gam <- matrix(gam, nrow = 1)
  miss <- is.na(d$responses)
  sims <- .with_seed(seed, lapply(seq_len(nsim), function(i) {
    s <- .simulate_itemcov_data(d$x, unname(f$alpha.vec), unname(f$beta.vec), W, gam, eps,
                                unname(s2$alpha), unname(f$eta), metric = .form_metric(object$control$form))
    r <- s$dat; r[miss] <- NA; colnames(r) <- d$items; r
  }))
  attr(sims, "truth") <- list(slopes = f$alpha.vec, intercepts = f$beta.vec, impact = s2$alpha,
                              variance_impact = f$eta, effects = gam, residual_dif = eps, form = object$control$form)
  sims
}
