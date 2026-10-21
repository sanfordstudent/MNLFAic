# Diagnostics for the generalized ML1 identification condition (ported from
# chen-item-cov-ext/R/diagnostics.R).
#
# The condition, for one person covariate: with S = items with residual DIF,
# s_j their signs, and x_j item j's row of the Stage-2 design,
#     sum_{j in S} s_j x_j'u  <  sum_{j not in S} |x_j'u|   for every u != 0.
# With b = sum_{j in S} s_j x_j, the margin
#     rho = max_{u != 0} b'u / sum_{j not in S} |x_j'u|
# must be < 1. Computation: the condition holds exactly when a median
# regression on constructed data (the DIF-free items with response 0, plus one
# pseudo-item with row t * b and a huge response) returns coefficients of
# exactly zero; rho is found by bisection over t.

# Strict inequality: "holds" is tested at a slightly inflated multiplier, so
# ties (rho = 1) and margins within 1e-4 of 1 count as failing (the safe side).
.strict <- 1 + 1e-4

# Interior-point solver ("fn"): the simplex method can cycle forever on these
# degenerate problems. The decision is scale-free: if the condition holds at t,
# the solution is u = 0; if not, it moves until t * b'u ~ big.
.ml1_holds_at <- function(Xc, b, t) {
  if (all(b == 0)) return(TRUE)
  big <- 1e3 * (1 + max(abs(Xc)))
  fit <- suppressWarnings(quantreg::rq.fit(rbind(Xc, t * b), c(rep(0, nrow(Xc)), big),
                                           tau = 0.5, method = "fn"))
  t * abs(sum(b * fit$coefficients)) < 1e-3 * big
}

.ml1_rho <- function(Xc, b, rel_tol = 1e-4) {
  if (all(b == 0)) return(0)
  holds1 <- .ml1_holds_at(Xc, b, 1)
  if (holds1) {                            # bracket [t_lo, t_hi]: holds at t_lo, fails at t_hi
    t_lo <- 1; t_hi <- 2
    while (.ml1_holds_at(Xc, b, t_hi)) {
      t_lo <- t_hi; t_hi <- 2 * t_hi
      if (t_hi > 1e12) return(0)
    }
  } else {
    t_hi <- 1; t_lo <- 0.5
    while (!.ml1_holds_at(Xc, b, t_lo)) {
      t_hi <- t_lo; t_lo <- t_lo / 2
      if (t_lo < 1e-12) return(Inf)
    }
  }
  while ((t_hi - t_lo) / t_lo > rel_tol) {
    t_mid <- sqrt(t_lo * t_hi)
    if (.ml1_holds_at(Xc, b, t_mid)) t_lo <- t_mid else t_hi <- t_mid
  }
  1 / sqrt(t_lo * t_hi)
}

# Check the condition for one person covariate.
#   lambda  : J-vector of slopes
#   W       : NULL or J x M item-covariate matrix (as used in Stage 2)
#   support : logical J-vector (TRUE = residual DIF) or item indices
#   signs   : J-vector (or length-|S| vector) of +1/-1 for the residual-DIF items
#   margin  : also compute rho (about 20-40 extra LAD fits)
# Returns list(holds, rho, reason).
.ml1_condition_check <- function(lambda, W = NULL, support, signs = 1,
                                 metric = c("difficulty", "logit"), margin = TRUE) {
  metric <- match.arg(metric)
  X <- .ml1_design(lambda, W, metric)
  J <- nrow(X); p <- ncol(X)
  S <- if (is.logical(support)) which(support) else as.integer(support)
  s <- if (length(signs) == J) signs[S] else rep_len(signs, length(S))
  if (any(!(s %in% c(-1, 1)))) stop("signs must be +1 or -1.", call. = FALSE)
  Xc <- X[setdiff(seq_len(J), S), , drop = FALSE]
  if (nrow(Xc) == 0 || qr(Xc)$rank < p)
    return(list(holds = FALSE, rho = Inf,
                reason = "items without residual DIF do not span the Stage-2 design"))
  b <- colSums(s * X[S, , drop = FALSE])
  holds <- .ml1_holds_at(Xc, b, .strict)
  rho   <- if (margin) .ml1_rho(Xc, b) else NA_real_
  list(holds = holds, rho = rho,
       reason = if (holds) "condition holds" else "residual-DIF items outweigh the rest in some direction")
}

# How many same-signed residual-DIF items each group of items can absorb before
# the condition fails, with no residual DIF elsewhere.
#   groups   : named list of item-index vectors
#   sign     : +1 or -1 (all added items share it; the hardest case)
#   n_random : number of random orders for the "typical" figure
# by_lambda: items added in decreasing slope order (exact worst case for a
# category with equal slopes; otherwise a heuristic); random: median and 5th
# percentile over random orders. Each order is scanned by binary search.
# seed: the random orders use this seed; the caller's RNG state is restored.
.ml1_tolerance_table <- function(lambda, W = NULL, groups, sign = 1,
                                 metric = c("difficulty", "logit"), n_random = 200, seed = 1) {
  metric <- match.arg(metric)
  X <- .ml1_design(lambda, W, metric)
  J <- nrow(X)
  holds_with <- function(S) {
    if (length(S) == 0) return(TRUE)
    Xc <- X[setdiff(seq_len(J), S), , drop = FALSE]
    if (qr(Xc)$rank < ncol(X)) return(FALSE)
    .ml1_holds_at(Xc, colSums(sign * X[S, , drop = FALSE]), .strict)
  }
  tolerated <- function(ord) {
    lo <- 0; hi <- length(ord)             # largest k with holds(first k items)
    if (holds_with(ord)) return(hi)
    while (hi - lo > 1) {
      mid <- (lo + hi) %/% 2
      if (holds_with(ord[seq_len(mid)])) lo <- mid else hi <- mid
    }
    lo
  }
  .with_seed(seed, {
    rows <- lapply(names(groups), function(g) {
      idx  <- groups[[g]]
      rand <- replicate(n_random, tolerated(idx[sample.int(length(idx))]))
      data.frame(group = g, n_items = length(idx),
                 by_lambda = tolerated(idx[order(lambda[idx], decreasing = TRUE)]),
                 random_median = stats::median(rand),
                 random_p05 = unname(stats::quantile(rand, 0.05)))
    })
    do.call(rbind, rows)
  })
}

# Summary of the Stage-2 design: collinearity of slopes with item covariates
# (which makes impact hard to separate from item-covariate DIF) and leverage.
.ml1_design_report <- function(lambda, W = NULL, metric = c("difficulty", "logit"), n_top = 5) {
  metric <- match.arg(metric)
  X <- .ml1_design(lambda, W, metric)
  Xs <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  H  <- Xs %*% solve(crossprod(Xs), t(Xs))
  lev <- diag(H)
  list(J = nrow(X), p = ncol(X),
       cor_slope_with_W = if (ncol(X) > 1) drop(stats::cor(X[, 1], X[, -1, drop = FALSE])) else NULL,
       condition_number = kappa(Xs, exact = TRUE),
       leverage_mean = ncol(X) / nrow(X),
       leverage_top = utils::head(data.frame(item = seq_along(lev), leverage = lev)[order(-lev), ], n_top))
}
