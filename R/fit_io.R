# Compact storage of Stage-1 fits (ported from chen-item-cov-ext/R/fit_io.R),
# so Stage 2/3 can be re-run without refitting. The information matrix is
# symmetric, so only its upper triangle is kept; for a cluster-robust fit the
# G x D cluster score sums are kept and the sandwich is recomputed on restore.

.compact_stage1_fit <- function(fit, extra = list()) {
  info <- (fit$info + t(fit$info)) / 2
  c(list(beta = fit$beta, eta = fit$eta, alpha.vec = fit$alpha.vec,
         beta.vec = fit$beta.vec, gamma.mat = fit$gamma.mat,
         info_upper = info[upper.tri(info, diag = TRUE)], D = nrow(info),
         N = if (!is.null(fit$post)) nrow(fit$post) else fit$N,
         loglik = fit$loglik, convergence = fit$convergence,
         anchor = if (!is.null(fit$anchor)) fit$anchor else 1L,
         vcov_type = if (!is.null(fit$vcov_type)) fit$vcov_type else "model",
         cluster_scores = fit$cluster_scores, n_clusters = fit$n_clusters), extra)
}

.restore_stage1_fit <- function(s) {
  info <- matrix(0, s$D, s$D)
  info[upper.tri(info, diag = TRUE)] <- s$info_upper
  s$info <- info + t(info) - diag(diag(info))
  s$info_upper <- NULL
  if (identical(s$vcov_type, "cluster")) {
    s$vcov <- .stage1_cluster_vcov(s$info, s$cluster_scores)
  } else {
    Vm <- solve(s$info); s$vcov <- (Vm + t(Vm)) / 2; s$vcov_type <- "model"
  }
  s
}
