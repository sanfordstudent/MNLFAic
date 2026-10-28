# Stage 3: bootstrap inference through the generalized ML1 step (ported from
# chen-item-cov-ext/R/stage3_inference.R).
#
# Chen et al.'s (2023) Algorithm 2, a parametric bootstrap, with Stage 2
# replaced by .ml1_itemcov():
#   1. draw M perturbations Z_m of the anchored Stage-1 parameters from
#      N(0, V), V = I^{-1} (or the cluster-robust covariance; then a
#      multivariate t with G - 1 df by default);
#   2. apply Stage 2 to center + Z_m;
#   3. e_m = G(center + Z_m) - G(center) for every output quantity G;
#   4. intervals, p-values and SEs from the e_m.
# center = "threshold" (default): the estimates, except that residual DIF
#   smaller than threshold_z bootstrap SEs (from a pilot of estimate-centered
#   draws) is moved exactly onto the fitted plane (Fang & Santos's selection
#   estimator of the directional derivative).
# center = "estimate": the estimates (Chen et al.'s paper).
# center = "null": DIF centered at zero for every item (Chen et al.'s code).
# Interval rule: est - (upper, lower) quantiles of e for "estimate"/"threshold";
# est + (lower, upper) for "null". p-value: mean(|e| > |est|).

# Apply Stage 2 to every draw. Returns a list of M-row matrices.
.stage2_on_draws <- function(noise, nu_center, lambda_center, alpha_center, W, metric, J, K, anchor = 1) {
  ix <- .stage1_index(J, K)
  M  <- nrow(noise)
  out <- NULL
  n_nonunique <- 0
  for (m in seq_len(M)) {
    nu_m <- nu_center
    nu_m[-anchor, ] <- nu_m[-anchor, ] + matrix(noise[m, ix$g], J - 1, K, byrow = TRUE)   # anchor fixed
    lam_m <- lambda_center + noise[m, ix$a]
    al_m  <- alpha_center + noise[m, ix$beta]
    r <- .ml1_itemcov(nu_m, lam_m, al_m, W, metric)
    g_flat <- unlist(r$gamma, use.names = TRUE)
    if (is.null(out)) {
      out <- list(eps = matrix(0, M, J * K), nu_total = matrix(0, M, J * K),
                  gamma = matrix(0, M, length(g_flat), dimnames = list(NULL, names(g_flat))),
                  delta = matrix(0, M, K), alpha = matrix(0, M, K))
      if (metric == "difficulty") out$eps_logit <- out$total_logit <- matrix(0, M, J * K)
    }
    out$eps[m, ]      <- as.vector(r$eps)        # column-major: items within covariate
    out$nu_total[m, ] <- as.vector(r$nu_total)
    if (metric == "difficulty") {                # the same on the logit scale, for p-values
      out$eps_logit[m, ]   <- as.vector(-lam_m * r$eps)
      out$total_logit[m, ] <- as.vector(-lam_m * r$nu_total)
    }
    out$gamma[m, ] <- g_flat
    out$delta[m, ] <- r$delta
    out$alpha[m, ] <- r$alpha
    n_nonunique <- n_nonunique + any(r$nonunique)
  }
  out$nonunique_rate <- n_nonunique / M
  out
}

.summarise_draws <- function(est, e, sig, ci_rule, shape = c("equal", "symmetric")) {
  shape <- match.arg(shape)
  lo <- apply(e, 2, stats::quantile, sig / 2, names = FALSE)
  hi <- apply(e, 2, stats::quantile, 1 - sig / 2, names = FALSE)
  qa <- apply(abs(e), 2, stats::quantile, 1 - sig, names = FALSE)    # symmetric (Hong & Li, 2018)
  out <- data.frame(
    est     = est,
    se      = apply(e, 2, stats::sd),
    ci_low  = NA_real_,
    ci_high = NA_real_,
    ci_eq_low   = if (ci_rule == "null") est + lo else est - hi,
    ci_eq_high  = if (ci_rule == "null") est + hi else est - lo,
    ci_sym_low  = est - qa,
    ci_sym_high = est + qa,
    p       = colMeans(abs(e) > rep(1, nrow(e)) %o% abs(est)))
  out$ci_low  <- if (shape == "equal") out$ci_eq_low  else out$ci_sym_low
  out$ci_high <- if (shape == "equal") out$ci_eq_high else out$ci_sym_high
  out
}

# Threshold (in bootstrap SEs) for threshold centering: a number, or a rule.
#   "sqrt2loglogN": sqrt(2 log log N)  (default; grows with N as the Fang-Santos
#                   rate conditions require; ~2.0 at N = 1190)
#   "sqrtlogN"    : sqrt(log N)
#   "bonferroni"  : sqrt(2 log(J K))
.resolve_threshold <- function(threshold_z, N, J, K) {
  if (is.numeric(threshold_z)) return(threshold_z)
  switch(threshold_z,
         sqrt2loglogN = sqrt(2 * log(log(N))),
         sqrtlogN     = sqrt(log(N)),
         bonferroni   = sqrt(2 * log(J * K)),
         stop("Unknown threshold rule: ", threshold_z, call. = FALSE))
}

# fit       : Stage-1 fit (.fit_stage1()); needs beta, alpha.vec, gamma.mat, info
# W         : NULL, J x M item-covariate matrix, or list of length C
# design    : optional .build_item_design() output (shared W): adds implied
#             omitted-category effects, joint tests per categorical variable,
#             and (difficulty metric) category-specific impacts
# M, sig    : number of draws, interval level 1 - sig
# seed      : if not NULL, the draws use this seed and the caller's RNG state
#             is restored afterwards
# vcov      : "auto" (fit$vcov when the fit is cluster-robust, else I^{-1}) or
#             "model" (always I^{-1})
# cluster_t : with a cluster-robust vcov, multivariate t draws with G - 1 df
# fdr       : Benjamini-Hochberg level for flagging residual DIF (all J x C tests)
# ci_shape_effects : interval shape for impact and item-covariate effects
# ci_shape_dif     : interval shape for residual and total DIF
.dif_inference_itemcov <- function(fit, W = NULL, design = NULL,
                                   metric = c("difficulty", "logit"),
                                   center = c("threshold", "estimate", "null"),
                                   M = 10000, sig = 0.05, seed = NULL,
                                   threshold_z = "sqrt2loglogN", n_pilot = 2000, fdr = 0.05,
                                   ci_shape_effects = c("symmetric", "equal"),
                                   ci_shape_dif = c("symmetric", "equal"),
                                   N = NULL, keep_draws = FALSE, vcov = c("auto", "model"),
                                   cluster_t = TRUE) {
  vcov <- match.arg(vcov)
  metric <- match.arg(metric)
  center <- match.arg(center)
  ci_shape_effects <- match.arg(ci_shape_effects)
  ci_shape_dif     <- match.arg(ci_shape_dif)
  J <- length(fit$alpha.vec); K <- ncol(as.matrix(fit$gamma.mat))
  if (is.null(fit$info)) stop("The Stage-1 fit has no information matrix.", call. = FALSE)
  if (!is.null(design) && !is.matrix(W)) stop("`design` can only be used with a single shared W.", call. = FALSE)
  cov_names <- if (!is.null(colnames(fit$gamma.mat))) colnames(fit$gamma.mat) else paste0("x", 1:K)

  anchor <- if (!is.null(fit$anchor)) fit$anchor else 1L
  use_cluster <- vcov == "auto" && identical(fit$vcov_type, "cluster")
  varian <- if (use_cluster) fit$vcov else solve(fit$info)
  vcov_type <- if (use_cluster) "cluster" else "model"
  est <- .ml1_itemcov(fit$gamma.mat, fit$alpha.vec, fit$beta, W, metric)

  draw_df <- NA
  noise <- .with_seed(seed, {
    nz <- MASS::mvrnorm(n = M, mu = rep(0, nrow(varian)), Sigma = varian)
    if (use_cluster && cluster_t) {              # multivariate t, G - 1 df (after the normals)
      draw_df <- fit$n_clusters - 1
      nz <- nz / sqrt(stats::rchisq(M, draw_df) / draw_df)
    }
    nz
  })

  # ── center ──────────────────────────────────────────────────────────────
  nu_c <- fit$gamma.mat; alpha_c <- fit$beta
  n_zeroed <- NA_integer_; z_thr <- NA_real_
  if (center == "null") {
    nu_c <- matrix(0, J, K)
  } else if (center == "threshold") {
    pilot <- .stage2_on_draws(noise[seq_len(min(n_pilot, M)), , drop = FALSE],
                              fit$gamma.mat, fit$alpha.vec, fit$beta, W, metric, J, K, anchor)
    eps_hat <- as.vector(est$eps)
    e_pilot <- pilot$eps - rep(1, nrow(pilot$eps)) %o% eps_hat
    se_eps  <- apply(e_pilot, 2, stats::sd)
    if (is.null(N)) N <- if (!is.null(fit$post)) nrow(fit$post) else if (!is.null(fit$N)) fit$N else NA
    if (is.character(threshold_z) && threshold_z %in% c("sqrt2loglogN", "sqrtlogN") && is.na(N))
      stop("threshold_z = '", threshold_z, "' needs the number of persons: pass N = ...", call. = FALSE)
    z_thr <- .resolve_threshold(threshold_z, N, J, K)
    small <- abs(eps_hat) < z_thr * se_eps
    # move small residuals onto the fitted plane (on the logit scale: a
    # difficulty residual e is -lambda * e in logits); Stage 2 is unchanged
    to_logit <- if (metric == "difficulty") -rep(fit$alpha.vec, K) else 1
    nu_c <- fit$gamma.mat - matrix(ifelse(small, eps_hat * to_logit, 0), J, K)
    n_zeroed <- sum(small)
  }
  # Stage 2 at the center is written down, not solved: the centers are
  # degenerate LAD problems on which the simplex solver can cycle forever.
  g_center <- est
  if (center == "null") {
    g_center$delta[] <- 0
    g_center$gamma   <- lapply(est$gamma, function(g) 0 * g)
    g_center$eps[]   <- 0
    g_center$nu_total[] <- 0
    g_center$alpha   <- alpha_c
  } else if (center == "threshold") {
    moved <- matrix(ifelse(small, eps_hat, 0), J, K)
    g_center$eps      <- est$eps - moved
    g_center$nu_total <- est$nu_total - moved
  }

  # ── draws ───────────────────────────────────────────────────────────────
  d <- .stage2_on_draws(noise, nu_c, fit$alpha.vec, alpha_c, W, metric, J, K, anchor)
  dev <- function(draws, at) draws - rep(1, M) %o% at
  e_eps   <- dev(d$eps, as.vector(g_center$eps))
  e_nu    <- dev(d$nu_total, as.vector(g_center$nu_total))
  e_gamma <- dev(d$gamma, unlist(g_center$gamma))
  e_alpha <- dev(d$alpha, g_center$alpha)
  ci_rule <- if (center == "null") "null" else "paper"

  # Difficulty metric: estimates and intervals in difficulty units, but
  # "DIF = 0" is tested on the logit scale (the same null hypothesis; dividing
  # by each draw's slope makes difficulty-scale draws heavy-tailed and costs
  # the small p-values BH relies on).
  lam_cell <- rep(fit$alpha.vec, K)
  p_logit <- function(draws_logit, center_vals, est_vals)
    colMeans(abs(dev(draws_logit, -lam_cell * as.vector(center_vals))) >
             rep(1, M) %o% abs(-lam_cell * as.vector(est_vals)))
  item_df <- function(est_mat, e, p = NULL) {
    s <- .summarise_draws(as.vector(est_mat), e, sig, ci_rule, ci_shape_dif)
    if (!is.null(p)) s$p <- p
    cbind(data.frame(item = rep(seq_len(J), K), covariate = rep(cov_names, each = J)), s)
  }
  logit_p <- metric == "difficulty"
  residual <- item_df(est$eps, e_eps, if (logit_p) p_logit(d$eps_logit, g_center$eps, est$eps))
  residual$p_bh <- stats::p.adjust(residual$p, "BH")
  residual$flag <- residual$p_bh < fdr
  total <- item_df(est$nu_total, e_nu, if (logit_p) p_logit(d$total_logit, g_center$nu_total, est$nu_total))
  impact <- cbind(data.frame(covariate = cov_names),
                  .summarise_draws(est$alpha, e_alpha, sig, ci_rule, ci_shape_effects))

  gamma <- NULL; category_effects <- NULL; joint <- NULL; impact_by_category <- NULL
  g_est <- unlist(est$gamma)
  if (length(g_est) > 0) {
    gamma <- cbind(data.frame(term = names(g_est)), .summarise_draws(g_est, e_gamma, sig, ci_rule, ci_shape_effects))
    if (!is.null(design)) {
      # linear map from one covariate's gamma to all effects, incl. implied categories
      cols <- colnames(W)
      Tm <- do.call(cbind, lapply(cols, function(v)
        .full_category_effects(stats::setNames(as.numeric(cols == v), cols), design)))
      colnames(Tm) <- cols
      ce <- list(); jt <- list(); ic <- list()
      cat_terms <- unlist(lapply(names(design$coding), function(v) {
        cd <- design$coding[[v]]
        if (cd$type == "categorical") paste0(v, "_", cd$levels) else NULL }))
      for (c in seq_len(K)) {
        idx <- (c - 1) * length(cols) + seq_along(cols)
        est_full <- as.vector(Tm %*% g_est[idx])
        e_full   <- e_gamma[, idx, drop = FALSE] %*% t(Tm)
        ce[[c]]  <- cbind(data.frame(covariate = cov_names[c], term = rownames(Tm)),
                          .summarise_draws(est_full, e_full, sig, ci_rule, ci_shape_effects))
        # Impact for an item with no residual DIF in each category (difficulty
        # metric): alpha_c - gamma_c,level, latent-SD units. Effect coding makes
        # the unweighted mean over a variable's levels equal alpha_c.
        if (metric == "difficulty" && length(cat_terms) > 0) {
          k <- match(cat_terms, rownames(Tm))
          ic[[c]] <- cbind(data.frame(covariate = cov_names[c], term = cat_terms),
                           .summarise_draws(est$alpha[c] - est_full[k], e_alpha[, c] - e_full[, k, drop = FALSE],
                                            sig, ci_rule, ci_shape_effects))
        }
        for (v in names(design$coding)) {          # joint Wald test per categorical variable
          cd <- design$coding[[v]]
          if (cd$type != "categorical") next
          jdx <- idx[match(cd$columns, cols)]
          Vi  <- solve(stats::cov(e_gamma[, jdx, drop = FALSE]))
          stat   <- drop(t(g_est[jdx]) %*% Vi %*% g_est[jdx])
          stat_m <- rowSums((e_gamma[, jdx, drop = FALSE] %*% Vi) * e_gamma[, jdx, drop = FALSE])
          jt[[length(jt) + 1]] <- data.frame(covariate = cov_names[c], variable = v, df = length(jdx),
                                             stat = stat, p_boot = mean(stat_m >= stat),
                                             p_chisq = stats::pchisq(stat, length(jdx), lower.tail = FALSE))
        }
      }
      category_effects <- do.call(rbind, ce)
      joint <- if (length(jt) > 0) do.call(rbind, jt) else NULL
      impact_by_category <- if (length(ic) > 0) do.call(rbind, ic) else NULL
    }
  }

  out <- list(
    residual = residual, total = total, gamma = gamma,
    category_effects = category_effects, joint = joint, impact = impact,
    impact_by_category = impact_by_category, point = est,
    # cells the threshold centering kept off the fitted plane (J x K logical;
    # NULL for other centers): a less conservative flag set for diagnostics
    off_plane = if (center == "threshold") matrix(!small, J, K, dimnames = list(NULL, cov_names)) else NULL,
    settings = list(center = center, metric = metric, M = M, sig = sig, seed = seed,
                    threshold_z = if (center == "threshold") threshold_z else NA,
                    threshold_value = z_thr,
                    ci_shape_effects = ci_shape_effects, ci_shape_dif = ci_shape_dif,
                    n_residuals_zeroed = n_zeroed, fdr = fdr,
                    nonunique_rate = d$nonunique_rate, vcov_type = vcov_type, draw_df = draw_df,
                    n_clusters = if (vcov_type == "cluster") fit$n_clusters else NA))
  if (keep_draws) out$draws <- list(e_eps = e_eps, e_nu = e_nu, e_gamma = e_gamma,
                                    e_alpha = e_alpha, noise = noise)
  out
}
