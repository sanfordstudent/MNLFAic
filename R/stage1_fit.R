# Stage 1: anchored maximum likelihood with adaptive quadrature.
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

# Item types. A 1PL item shares one slope with the other 1PL items: a linear
# restriction theta = T phi of the 2PL parameters, with T a 0/1 matrix copying
# the common slope to every 1PL item. Score and information of phi are T's and
# T'IT, so the 2PL score and Louis information are reused unchanged; the
# model-based covariance of theta is T (T'IT)^{-1} T' (.inv_info()). T = NULL
# when every item is 2PL (the unrestricted code path, unchanged).
.restriction_matrix <- function(J, K, onepl) {
  if (is.null(onepl) || !any(onepl)) return(NULL)
  D <- 2 * K + 2 * J + (J - 1) * K
  col <- seq_len(D)
  s1 <- 2 * K + which(onepl)
  col[s1] <- s1[1]
  keep <- sort(unique(col))
  Tm <- matrix(0, D, length(keep))
  Tm[cbind(seq_len(D), match(col, keep))] <- 1
  Tm
}

# T of a fit (original item order).
.fit_T <- function(fit) .restriction_matrix(length(fit$alpha.vec), ncol(as.matrix(fit$gamma.mat)), fit$onepl)

# Model-based covariance from the information: I^{-1}, or T (T'IT)^{-1} T'.
.inv_info <- function(info, Tm = NULL) {
  if (is.null(Tm)) return(solve(info))
  Tm %*% solve(crossprod(Tm, info %*% Tm), t(Tm))
}

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
# w     : optional sampling weights (one per person); the fit maximizes the
#         weighted (pseudo-) log-likelihood
# onepl : optional logical (one per item, fitted order): 1PL items, which share
#         one slope
# Returns the estimates (beta, eta, alpha.vec = slopes, beta.vec = intercepts,
# gamma.mat = J x K DIF), post (N x n posterior weights), nodes (N x n), info
# (observed information, anchor removed), loglik, convergence.
.fit_stage1_anchor1 <- function(dat, x, start, em_tol = 0.005, em_stall = 0.5, newton_tol = 1e-4,
                                max_newton = 25, n_nodes = 21, max_fallback = 3,
                                fallback_em_iter = 50, verbose = FALSE, w = NULL, onepl = NULL) {
  t_start <- proc.time()[3]
  J <- ncol(dat); K <- ncol(x)
  start_g <- if (is.null(dim(start$g))) matrix(start$g, J, K) else start$g
  Tm <- .restriction_matrix(J, K, onepl)
  if (!is.null(Tm)) start$a[onepl] <- mean(start$a[onepl])        # start on the restriction

  em <- .em_aq(start$a, start$b, start_g, x, dat, start$beta, start$eta, n_nodes = n_nodes,
               tol = em_tol, stall_ratio = em_stall, verbose = verbose, w = w, onepl = onepl)
  em_iter <- em$iterations
  last_mode <- em$mode
  e_step <- function(u) {
    es <- .estep_aq(u$a, u$b, u$g, x, dat, u$beta, u$eta, n_nodes, start = last_mode, w = w)
    last_mode <<- es$mode
    es
  }
  info_at <- function(u, es) .info_stage1(u$a, u$b, u$g, x, dat, u$beta, u$eta, es$q, es$Theta, w = w)
  if (verbose) message(sprintf("EM: %d iterations, %.0f s", em_iter, proc.time()[3] - t_start))

  em_more <- function(u) .em_aq(u$a, u$b, u$g, x, dat, u$beta, u$eta, n_nodes = n_nodes,
                                tol = em_tol / 10, max_iter = fallback_em_iter, w = w, onepl = onepl)
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
      score <- .score_stage1(u$a, u$b, u$g, x, dat, u$beta, u$eta, es$q, es$Theta, w = w)
      info  <- info_at(u, es)
      step  <- tryCatch(if (is.null(Tm)) solve(info, score) else
                          as.vector(Tm %*% solve(crossprod(Tm, info %*% Tm), crossprod(Tm, score))),
                        error = function(e) NULL)
      failed <- if (is.null(step)) "singular information" else NA_character_
      if (is.na(failed)) {
        se <- sqrt(pmax(diag(.inv_info(info, Tm)), .Machine$double.eps))
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
  score <- .score_stage1(u$a, u$b, u$g, x, dat, u$beta, u$eta, es$q, es$Theta, w = w)
  if (!is.null(Tm)) score <- as.vector(crossprod(Tm, score))        # restricted score
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
#   weights: optional sampling weights (one per person, >= 0), used as given
#            (scaling is done by the caller). The estimates maximize the
#            weighted log-likelihood; the cluster sums use weighted scores.
#   replicates: optional N x R replicate weights (needs `weights`, scaled the
#            same way; not with `cluster`). Adds the replicate-weight
#            covariance by one-step linearization (.replicate_contributions());
#            replicate_method and fay set its constant (.set_replicate_vcov()).
#   onepl  : optional logical, one per item (original order): 1PL items, which
#            share one slope (see .restriction_matrix()). Stored in the fit.
# Items are reordered so the anchor comes first, fitted, and every output is
# mapped back: alpha.vec, beta.vec, gamma.mat (row `anchor` is the fixed 0);
# info (and vcov) laid out as beta (K), eta (K), a (J), b (J), then the DIF of
# the non-anchor items in original order (row-major). Adds anchor, vcov,
# vcov_type ("model"/"cluster"), and with a cluster n_clusters and
# cluster_scores (G x D).
.fit_stage1 <- function(dat, x, start, anchor = "auto", cluster = NULL, weights = NULL, replicates = NULL,
                        replicate_method = "JK2", fay = 0.5, onepl = NULL, ...) {
  J <- ncol(dat); K <- ncol(x)
  if (!is.null(onepl) && (length(onepl) != J || anyNA(onepl))) stop("`onepl` must have one value per item.", call. = FALSE)
  if (!is.null(onepl) && !any(onepl)) onepl <- NULL
  if (!is.null(weights) && (length(weights) != nrow(dat) || anyNA(weights) || any(weights < 0)))
    stop("`weights` must have one non-negative, non-missing value per person.", call. = FALSE)
  if (!is.null(replicates)) {
    if (is.null(weights)) stop("Replicate weights need the full-sample `weights`.", call. = FALSE)
    if (!is.null(cluster)) stop("Use either `cluster` or replicate weights, not both.", call. = FALSE)
    replicates <- as.matrix(replicates)
    if (nrow(replicates) != nrow(dat) || anyNA(replicates) || any(replicates < 0))
      stop("Replicate weights must be a matrix with one row per person and no missing or negative values.",
           call. = FALSE)
  }
  if (!is.null(cluster)) {
    if (length(cluster) != nrow(dat)) stop("`cluster` must have one value per person (row of the data).", call. = FALSE)
    if (anyNA(cluster)) stop("`cluster` contains missing values.", call. = FALSE)
  }
  a <- if (identical(anchor, "auto")) .choose_stage1_anchor(dat, x) else as.integer(anchor)
  ord <- c(a, setdiff(seq_len(J), a)); inv <- order(ord)        # original item j sits at position inv[j]
  dat_o <- dat[, ord, drop = FALSE]
  if (a == 1L) {
    fit <- .fit_stage1_anchor1(dat, x, start, w = weights, onepl = onepl, ...)
    P <- seq_len(nrow(fit$info))
  } else {
    g <- if (is.null(dim(start$g))) matrix(start$g, J, K) else start$g
    st <- start; st$a <- start$a[ord]; st$b <- start$b[ord]
    st$g <- g[ord, , drop = FALSE]; st$g[1, ] <- 0
    fit <- .fit_stage1_anchor1(dat_o, x, st, w = weights, onepl = onepl[ord], ...)
    ix_a <- 2 * K + seq_len(J); ix_b <- 2 * K + J + seq_len(J)
    g_pos <- function(r) 2 * K + 2 * J + (r - 2) * K + seq_len(K)   # DIF block of reordered item r >= 2
    P <- c(seq_len(2 * K), ix_a[inv], ix_b[inv],
           unlist(lapply(setdiff(seq_len(J), a), function(j) g_pos(inv[j]))))
  }
  S <- if (!is.null(cluster) || !is.null(replicates)) .stage1_person_scores(fit, dat_o, x) else NULL
  if (!is.null(S) && !is.null(weights) && is.null(replicates)) S <- weights * S
  fit$weights <- weights
  if (a != 1L) {
    fit$alpha.vec <- fit$alpha.vec[inv]; fit$beta.vec <- fit$beta.vec[inv]
    fit$gamma.mat <- fit$gamma.mat[inv, , drop = FALSE]
    fit$info <- fit$info[P, P]
  }
  fit$anchor <- a
  fit$onepl <- onepl
  if (!is.null(replicates)) {
    fit$replicate <- .replicate_contributions(S[, P, drop = FALSE], replicates - weights, fit$info, .fit_T(fit))
    return(.set_replicate_vcov(fit, replicate_method, fay))
  }
  .add_stage1_vcov(fit, if (is.null(S)) NULL else S[, P, drop = FALSE], cluster)
}

# Attach vcov (and, with a cluster, the CR1 sandwich ingredients) to a fit.
.add_stage1_vcov <- function(fit, S = NULL, cluster = NULL) {
  Tm <- .fit_T(fit)
  Vm <- .inv_info(fit$info, Tm); Vm <- (Vm + t(Vm)) / 2
  if (is.null(cluster)) {
    fit$vcov <- Vm; fit$vcov_type <- "model"
    return(fit)
  }
  Sg <- rowsum(S, as.character(cluster), reorder = FALSE)
  G  <- nrow(Sg)
  if (G < 20) warning("Only ", G, " clusters; cluster-robust standard errors may be too small.", call. = FALSE)
  fit$cluster_scores <- Sg
  fit$n_clusters <- G
  fit$t_df <- G - 1
  fit$vcov <- .stage1_cluster_vcov(fit$info, Sg, Tm)
  fit$vcov_type <- "cluster"
  fit
}

# Replicate-weight covariance by one-step linearization. Replicate r's estimates
# are approximated by one Newton step from the full-sample estimates, where the
# full-sample weighted score is zero:
#   d_r = I_w^{-1} sum_i (w_ir - w_i) s_i,
# and V = c sum_r d_r d_r' with the replication constant c. Replicates whose
# weight changes are exact negatives of each other (JK2 with both halves of
# every zone) have d_s = -d_r; each such pair is kept once with multiplicity 2,
# so the wild bootstrap gets one weight per zone. Replicates with no weight
# change (unused zones) are dropped.
# S: N x D unweighted per-person scores; delta: N x R; info: the weighted
# information. Returns B (rows sqrt(multiplicity) d_r, so V = c B'B),
# n_replicates (R), n_units (rows of B) and n_paired.
.replicate_contributions <- function(S, delta, info, Tm = NULL) {
  R <- ncol(delta)
  Vm <- .inv_info(info, Tm); Vm <- (Vm + t(Vm)) / 2
  d <- t(Vm %*% crossprod(S, delta))
  nrm <- sqrt(colSums(delta^2))
  C <- crossprod(delta)
  keep <- nrm > 0; mult <- rep(1, R)
  for (r in which(keep)) {
    if (!keep[r] || mult[r] == 2) next
    tol <- 1e-10 * nrm[r] * nrm
    s <- which(keep & seq_len(R) > r & mult == 1 & abs(C[r, ] + nrm[r] * nrm) <= tol & abs(nrm - nrm[r]) <= 1e-10 * nrm[r])
    if (length(s) > 0) { keep[s[1]] <- FALSE; mult[r] <- 2 }
  }
  list(B = sqrt(mult[keep]) * d[keep, , drop = FALSE], n_replicates = R, n_units = sum(keep),
       n_paired = sum(mult[keep] == 2))
}

# Replication constant and degrees of freedom.
#   JK2      : both halves of every zone (IEA, 2015 onward): c = 1/2, df = zones
#   JK2-half : one replicate per zone (earlier IEA cycles): c = 1, df = zones
#   BRR      : Fay's method with factor k (PISA: k = 0.5, R = 80): c = 1 / (R (1 - k)^2), df = R
.replicate_scale <- function(method, R, fay) {
  switch(method, JK2 = 0.5, `JK2-half` = 1, BRR = 1 / (R * (1 - fay)^2),
         stop("Unknown replicate method: ", method, call. = FALSE))
}

.set_replicate_vcov <- function(fit, method = "JK2", fay = 0.5) {
  rp <- fit$replicate
  c0 <- .replicate_scale(method, rp$n_replicates, fay)
  V <- c0 * crossprod(rp$B)
  fit$vcov <- (V + t(V)) / 2
  fit$vcov_type <- "replicate"
  fit$rep_method <- method; fit$rep_fay <- fay; fit$rep_scale <- c0
  fit$t_df <- if (method == "BRR") rp$n_replicates else rp$n_units
  fit
}

# CR1 sandwich from the information and the cluster score sums.
.stage1_cluster_vcov <- function(info, cluster_scores, Tm = NULL) {
  G <- nrow(cluster_scores)
  Vm <- .inv_info(info, Tm)
  V  <- Vm %*% (G / (G - 1) * crossprod(cluster_scores)) %*% Vm
  (V + t(V)) / 2
}

# The Stage-1 covariance for inference: fit$vcov when present, else I^{-1};
# type = "model" forces I^{-1}.
.stage1_vcov <- function(fit, type = c("auto", "model")) {
  type <- match.arg(type)
  if (type == "model" || is.null(fit$vcov)) { V <- .inv_info(fit$info, .fit_T(fit)); return((V + t(V)) / 2) }
  fit$vcov
}
