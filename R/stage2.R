# Stage 2: generalized minimal-L1 (ML1) transformation with item covariates
# (ported from chen-item-cov-ext/R/stage2_ml1_itemcov.R).
#
# For each person covariate c, a median (LAD) regression without intercept of
# the anchored Stage-1 DIF estimates on the Stage-2 design:
#   logit metric      : design [lambda, W]
#   difficulty metric : design [lambda, lambda * W], i.e. a slope-weighted
#                       median regression of difficulty DIF on [1, W]
# giving delta_c (the slide), gamma_c (item-covariate effects) and the
# residual DIF. Impact alpha_c = alpha_tilde_c + delta_c.
#
# metric = "difficulty" (default): gamma, residual and total DIF are shifts in
#   item difficulty, latent-SD units, + = harder:
#   nu_jc = -lambda_j (gamma_c' w_j + eps_jc).
# metric = "logit": item covariates shift the logit intercept,
#   nu_jc = gamma_c' w_j + eps_jc; logits, + = easier.
# With no item covariates both metrics reduce to Chen et al.'s (2023) ML1
# step and give the same impact.

# Normalize W into a list with one entry (matrix or NULL) per person
# covariate. W may be NULL, a J x M matrix shared by all covariates, or a list
# of length C.
.as_W_list <- function(W, C, J) {
  if (is.null(W)) return(vector("list", C))
  if (is.matrix(W)) W <- rep(list(W), C)
  if (!is.list(W) || length(W) != C) stop("W must be NULL, a J x M matrix, or a list of length C.", call. = FALSE)
  for (c in seq_len(C)) {
    Wc <- W[[c]]
    if (!is.null(Wc)) {
      if (!is.matrix(Wc) || nrow(Wc) != J) stop("W[[", c, "]] must be a matrix with J rows.", call. = FALSE)
      if (is.null(colnames(Wc))) colnames(W[[c]]) <- paste0("w", seq_len(ncol(Wc)))
    }
  }
  W
}

# Stage-2 design matrix for one person covariate.
.ml1_design <- function(lambda, Wc, metric = c("difficulty", "logit")) {
  metric <- match.arg(metric)
  if (is.null(Wc) || ncol(Wc) == 0) {
    X <- matrix(lambda, ncol = 1)
  } else if (metric == "logit") {
    X <- cbind(lambda, Wc)
  } else {
    X <- cbind(lambda, lambda * Wc)
  }
  colnames(X) <- c("delta", if (!is.null(Wc)) colnames(Wc))
  X
}

# One LAD fit (Barrodale-Roberts simplex, as rq(y ~ -1 + X, tau = 0.5));
# records whether quantreg flagged the solution as possibly nonunique.
# Guard: an all-zero response has the exact solution 0; the simplex solver can
# cycle forever on it, so it is not called.
.lad_fit <- function(X, y) {
  if (all(y == 0)) {
    coef <- stats::setNames(numeric(ncol(X)), colnames(X))
    return(list(coef = coef, resid = numeric(length(y)), nonunique = FALSE))
  }
  nonunique <- FALSE
  fit <- withCallingHandlers(
    quantreg::rq.fit(X, y, tau = 0.5, method = "br"),
    warning = function(w) {
      if (grepl("nonunique", conditionMessage(w))) {
        nonunique <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  coef <- fit$coefficients
  names(coef) <- colnames(X)
  list(coef = coef, resid = as.vector(y - X %*% coef), nonunique = nonunique)
}

# nu_hat     : J x C anchored Stage-1 DIF estimates (logits)
# lambda_hat : J-vector of Stage-1 slopes
# alpha_hat  : optional C-vector of anchored Stage-1 impact
# W          : NULL, J x M matrix, or list of length C
# Returns list(delta, gamma, eps, nu_total, alpha, nonunique, metric); eps and
# nu_total (total DIF) are J x C in the metric's units.
.ml1_itemcov <- function(nu_hat, lambda_hat, alpha_hat = NULL, W = NULL,
                         metric = c("difficulty", "logit")) {
  metric <- match.arg(metric)
  nu_hat <- as.matrix(nu_hat)
  J <- nrow(nu_hat); C <- ncol(nu_hat)
  if (length(lambda_hat) != J) stop("lambda_hat must have length J.", call. = FALSE)
  W_list <- .as_W_list(W, C, J)
  cov_names <- if (!is.null(colnames(nu_hat))) colnames(nu_hat) else paste0("x", seq_len(C))

  delta     <- stats::setNames(numeric(C), cov_names)
  gamma     <- stats::setNames(vector("list", C), cov_names)
  eps       <- matrix(0, J, C, dimnames = list(rownames(nu_hat), cov_names))
  nonunique <- stats::setNames(logical(C), cov_names)
  for (c in seq_len(C)) {
    X   <- .ml1_design(lambda_hat, W_list[[c]], metric)
    fit <- .lad_fit(X, nu_hat[, c])
    delta[c]     <- fit$coef[1]
    gamma[[c]]   <- fit$coef[-1]
    eps[, c]     <- fit$resid
    nonunique[c] <- fit$nonunique
  }
  nu_total <- nu_hat - outer(lambda_hat, delta)
  dimnames(nu_total) <- dimnames(eps)
  if (metric == "difficulty") {             # report on the difficulty scale, + = harder
    gamma    <- lapply(gamma, function(g) -g)
    eps      <- -eps / lambda_hat
    nu_total <- -nu_total / lambda_hat
  }
  list(delta = delta, gamma = gamma, eps = eps, nu_total = nu_total,
       alpha = if (!is.null(alpha_hat)) alpha_hat + delta else NULL,
       nonunique = nonunique, metric = metric)
}
