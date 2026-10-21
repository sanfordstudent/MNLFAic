# Stage 1: anchored maximum likelihood (ported from
# chen-item-cov-ext/R/stage1_fit.R; adaptive quadrature only).
#
# EM converges very slowly along the impact/DIF "slide", so its own stopping
# rule stops far from the optimum. The fit therefore
#   1. runs EM (adaptive quadrature) with a loose tolerance to get close;
#   2. takes Newton steps theta <- theta + I^{-1} S with the analytic score
#      (Fisher's identity) and the observed information (Louis' identity),
#      halving each step until the log-likelihood does not decrease;
#   3. stops when the last step moved every parameter by less than newton_tol
#      standard errors.
# When Newton fails (singular information, failed line search) or uses up
# max_newton iterations, up to fallback_em_iter EM iterations are run and
# Newton resumes, at most max_fallback times; after that the fit is returned
# with converged = FALSE and a warning.

.pack_stage1 <- function(beta, eta, a, b, g) c(beta, eta, a, b, as.vector(t(g[-1, , drop = FALSE])))

.unpack_stage1 <- function(p, J, K) {
  g <- matrix(0, J, K)
  g[-1, ] <- matrix(p[(2 * K + 2 * J + 1):length(p)], J - 1, K, byrow = TRUE)
  list(beta = p[1:K], eta = p[(K + 1):(2 * K)],
       a = p[(2 * K + 1):(2 * K + J)], b = p[(2 * K + J + 1):(2 * K + 2 * J)], g = g)
}

# Fit with item 1 as the anchor (its DIF fixed at 0).
# dat, x: N x J responses (0/1/NA), N x K person covariates
# start : list(a, b, g, beta, eta) (g: J x K, row 1 = 0)
# Returns the estimates (beta, eta, alpha.vec = slopes, beta.vec = intercepts,
# gamma.mat = J x K DIF), post (N x n posterior weights), nodes (N x n), info
# (observed information, anchor removed), loglik, convergence.
.fit_stage1_anchor1 <- function(dat, x, start, em_tol = 0.005, em_stall = 0.5, newton_tol = 1e-4,
                                max_newton = 25, n_nodes = 21, max_fallback = 3,
                                fallback_em_iter = 50, verbose = FALSE) {
  t_start <- proc.time()[3]
  J <- ncol(dat); K <- ncol(x)
  start_g <- if (is.null(dim(start$g))) matrix(start$g, J, K) else start$g

  em <- .em_aq(start$a, start$b, start_g, x, dat, start$beta, start$eta, n_nodes = n_nodes,
               tol = em_tol, stall_ratio = em_stall, verbose = verbose)
  em_iter <- em$iterations
  last_mode <- em$mode
  e_step <- function(u) {
    es <- .estep_aq(u$a, u$b, u$g, x, dat, u$beta, u$eta, n_nodes, start = last_mode)
    last_mode <<- es$mode
    es
  }
  info_at <- function(u, es) .info_stage1(u$a, u$b, u$g, x, dat, u$beta, u$eta, es$q, es$Theta)
  if (verbose) message(sprintf("EM: %d iterations, %.0f s", em_iter, proc.time()[3] - t_start))

  em_more <- function(u) .em_aq(u$a, u$b, u$g, x, dat, u$beta, u$eta, n_nodes = n_nodes,
                                tol = em_tol / 10, max_iter = fallback_em_iter)
  can_fallback <- function() n_fallback < max_fallback
  fallback <- function(reason, it) {
    n_fallback <<- n_fallback + 1
    em2 <- em_more(u)
    p   <<- .pack_stage1(em2$beta, em2$eta, em2$alpha.vec, em2$beta.vec, em2$gamma.mat)
    u   <<- .unpack_stage1(p, J, K)
    es  <<- e_step(u); ll <<- es$loglik
    ll_trace <<- c(ll_trace, ll)
    if (verbose) message(sprintf("Newton %2d: %s; EM fallback %d (%d iterations): loglik %.6f",
                                 it, reason, n_fallback, em2$iterations, ll))
  }

  p  <- .pack_stage1(em$beta, em$eta, em$alpha.vec, em$beta.vec, em$gamma.mat)
  u  <- .unpack_stage1(p, J, K)
  es <- e_step(u)
  ll <- es$loglik
  ll_trace  <- ll
  converged <- FALSE
  max_ratio <- NA_real_
  n_fallback <- 0
  stop_reason <- "max_newton"

  repeat {
    for (it in seq_len(max_newton)) {
      score <- .score_stage1(u$a, u$b, u$g, x, dat, u$beta, u$eta, es$q, es$Theta)
      info  <- info_at(u, es)
      step  <- tryCatch(solve(info, score), error = function(e) NULL)
      failed <- if (is.null(step)) "singular information" else NA_character_
      if (is.na(failed)) {
        se <- sqrt(pmax(diag(solve(info)), .Machine$double.eps))
        # Step halving: never accept a real decrease in log-likelihood. The
        # tolerance allows for quadrature noise (the person-specific nodes move
        # with beta/eta); at convergence it is far below 1e-6.
        frac <- 1
        repeat {
          p_new  <- p + frac * step
          u_new  <- .unpack_stage1(p_new, J, K)
          es_new <- e_step(u_new)
          if (is.finite(es_new$loglik) && es_new$loglik >= ll - 1e-6) break
          frac <- frac / 2
          if (frac < 1e-6) { failed <- "line search failed"; break }
        }
      }
      if (!is.na(failed)) {
        if (!can_fallback()) { stop_reason <- failed; break }
        fallback(failed, it)
        next
      }
      max_ratio <- max(abs(frac * step) / se)
      p <- p_new; u <- u_new; es <- es_new; ll <- es$loglik
      ll_trace <- c(ll_trace, ll)
      if (verbose) message(sprintf("Newton %2d: loglik %.6f  step fraction %.3g  max|step|/SE %.2e  (%.0f s)",
                                   it, ll, frac, max_ratio, proc.time()[3] - t_start))
      if (max_ratio < newton_tol) { converged <- TRUE; stop_reason <- "converged"; break }
    }
    if (converged || stop_reason != "max_newton" || !can_fallback()) break
    fallback("max_newton reached", max_newton)
  }
  if (!converged) warning("Stage 1: Newton did not converge (", stop_reason, ").", call. = FALSE)

  info  <- info_at(u, es)
  info  <- (info + t(info)) / 2           # exactly symmetric
  score <- .score_stage1(u$a, u$b, u$g, x, dat, u$beta, u$eta, es$q, es$Theta)
  list(beta = u$beta, eta = u$eta, alpha.vec = u$a, beta.vec = u$b, gamma.mat = u$g,
       post = es$q, nodes = es$Theta, info = info, loglik = ll,
       convergence = list(converged = converged, stop_reason = stop_reason, em_fallbacks = n_fallback,
                          quadrature = "adaptive", n_nodes = n_nodes, em_iterations = em_iter,
                          newton_iterations = length(ll_trace) - 1, max_step_over_se = max_ratio,
                          max_abs_score = max(abs(score)), loglik_trace = ll_trace,
                          seconds = unname(proc.time()[3] - t_start)))
}

# Automatic Stage-1 anchor. The anchor is bookkeeping (every Stage-2/3 result
# is invariant to it) but matters numerically: a poorly determined anchor
# spreads its uncertainty to every item, making fits ill-conditioned.
# Rule: among items whose item-rest correlation is at least the median, the
# one whose smallest count of correct / incorrect responses, overall and within
# each group of every binary person covariate, is largest.
.choose_stage1_anchor <- function(dat, x) {
  cnt <- function(rows) {
    d <- dat[rows, , drop = FALSE]
    pmin(colSums(d == 1, na.rm = TRUE), colSums(d == 0, na.rm = TRUE))
  }
  score <- cnt(rep(TRUE, nrow(dat)))
  for (k in seq_len(ncol(x))) {
    v <- x[, k]
    if (all(v %in% c(0, 1))) score <- pmin(score, cnt(v == 1), cnt(v == 0))
  }
  r_rest <- .item_rest_correlation(dat)
  eligible <- !is.na(r_rest) & r_rest >= stats::median(r_rest, na.rm = TRUE)
  if (!any(eligible)) eligible[] <- TRUE
  which.max(ifelse(eligible, score, -Inf))
}

# Correlation of each item with the person's proportion correct on the other
# items (observed responses only); NA for an item without variation.
.item_rest_correlation <- function(dat) {
  vapply(seq_len(ncol(dat)), function(j) {
    rest <- rowMeans(dat[, -j, drop = FALSE], na.rm = TRUE)
    ok <- !is.na(dat[, j]) & is.finite(rest)
    if (length(unique(dat[ok, j])) < 2 || length(unique(rest[ok])) < 2) return(NA_real_)
    stats::cor(dat[ok, j], rest[ok])
  }, numeric(1))
}

# Stage-1 fit with any anchor, results in the original item order.
#   anchor : "auto" (default; .choose_stage1_anchor()) or an item index.
#   cluster: optional, one value per person (no NA). Adds the CR1 cluster-
#            robust covariance vcov = I^{-1} [G/(G-1) sum_g S_g S_g'] I^{-1}
#            from the cluster sums S_g of per-person scores; estimates are
#            unchanged.
# Items are reordered so the anchor comes first, fitted, and every output is
# mapped back: alpha.vec, beta.vec, gamma.mat (row `anchor` is the fixed 0);
# info (and vcov) laid out as beta (K), eta (K), a (J), b (J), then the DIF of
# the non-anchor items in original order (row-major). Adds anchor, vcov,
# vcov_type ("model"/"cluster"), and with a cluster n_clusters and
# cluster_scores (G x D).
.fit_stage1 <- function(dat, x, start, anchor = "auto", cluster = NULL, ...) {
  J <- ncol(dat); K <- ncol(x)
  if (!is.null(cluster)) {
    if (length(cluster) != nrow(dat)) stop("`cluster` must have one value per person (row of the data).", call. = FALSE)
    if (anyNA(cluster)) stop("`cluster` contains missing values.", call. = FALSE)
  }
  a <- if (identical(anchor, "auto")) .choose_stage1_anchor(dat, x) else as.integer(anchor)
  ord <- c(a, setdiff(seq_len(J), a)); inv <- order(ord)        # original item j sits at position inv[j]
  dat_o <- dat[, ord, drop = FALSE]
  if (a == 1L) {
    fit <- .fit_stage1_anchor1(dat, x, start, ...)
    P <- seq_len(nrow(fit$info))
  } else {
    g <- if (is.null(dim(start$g))) matrix(start$g, J, K) else start$g
    st <- start; st$a <- start$a[ord]; st$b <- start$b[ord]
    st$g <- g[ord, , drop = FALSE]; st$g[1, ] <- 0
    fit <- .fit_stage1_anchor1(dat_o, x, st, ...)
    ix_a <- 2 * K + seq_len(J); ix_b <- 2 * K + J + seq_len(J)
    g_pos <- function(r) 2 * K + 2 * J + (r - 2) * K + seq_len(K)   # DIF block of reordered item r >= 2
    P <- c(seq_len(2 * K), ix_a[inv], ix_b[inv],
           unlist(lapply(setdiff(seq_len(J), a), function(j) g_pos(inv[j]))))
  }
  S <- if (!is.null(cluster)) .stage1_person_scores(fit, dat_o, x) else NULL
  if (a != 1L) {
    fit$alpha.vec <- fit$alpha.vec[inv]; fit$beta.vec <- fit$beta.vec[inv]
    fit$gamma.mat <- fit$gamma.mat[inv, , drop = FALSE]
    fit$info <- fit$info[P, P]
  }
  fit$anchor <- a
  .add_stage1_vcov(fit, if (is.null(S)) NULL else S[, P, drop = FALSE], cluster)
}

# Attach vcov (and, with a cluster, the CR1 sandwich ingredients) to a fit.
.add_stage1_vcov <- function(fit, S = NULL, cluster = NULL) {
  Vm <- solve(fit$info); Vm <- (Vm + t(Vm)) / 2
  if (is.null(cluster)) {
    fit$vcov <- Vm; fit$vcov_type <- "model"
    return(fit)
  }
  Sg <- rowsum(S, as.character(cluster), reorder = FALSE)
  G  <- nrow(Sg)
  if (G < 20) warning("Only ", G, " clusters; cluster-robust standard errors may be too small.", call. = FALSE)
  fit$cluster_scores <- Sg
  fit$n_clusters <- G
  fit$vcov <- .stage1_cluster_vcov(fit$info, Sg)
  fit$vcov_type <- "cluster"
  fit
}

# CR1 sandwich from the information and the cluster score sums.
.stage1_cluster_vcov <- function(info, cluster_scores) {
  G <- nrow(cluster_scores)
  Vm <- solve(info)
  V  <- Vm %*% (G / (G - 1) * crossprod(cluster_scores)) %*% Vm
  (V + t(V)) / 2
}

# The Stage-1 covariance for inference: fit$vcov when present, else I^{-1};
# type = "model" forces I^{-1}.
.stage1_vcov <- function(fit, type = c("auto", "model")) {
  type <- match.arg(type)
  if (type == "model" || is.null(fit$vcov)) { V <- solve(fit$info); return((V + t(V)) / 2) }
  fit$vcov
}
