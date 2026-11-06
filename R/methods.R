# Methods and helpers for mnlfa_ic objects.

.units <- function(object) {
  if (identical(object$control$form, "slope-intercept")) "logits; + = easier for higher values of the person covariate"
  else "latent SD; + = harder for higher values of the person covariate"
}

.need_inference <- function(object) {
  if (is.null(object$inference)) stop("No inference yet: run infer(object) or mnlfa_ic(..., inference = TRUE).", call. = FALSE)
}

# Variance impact (log-variance coefficients) with Wald intervals from Stage 1;
# t with G - 1 df when the covariance is cluster-robust.
.variance_table <- function(object, level = object$control$level) {
  f <- object$stage1; K <- object$data$K
  V <- if (identical(object$control$vcov, "model")) .stage1_vcov(f, "model") else f$vcov
  se <- sqrt(diag(V)[(K + 1):(2 * K)])
  cluster <- identical(f$vcov_type, "cluster") && !identical(object$control$vcov, "model")
  crit <- if (cluster) stats::qt(1 - (1 - level) / 2, f$n_clusters - 1) else stats::qnorm(1 - (1 - level) / 2)
  p <- if (cluster) 2 * stats::pt(-abs(f$eta / se), f$n_clusters - 1) else 2 * stats::pnorm(-abs(f$eta / se))
  data.frame(covariate = object$data$covariates, est = unname(f$eta), se = se,
             ci_low = unname(f$eta) - crit * se, ci_high = unname(f$eta) + crit * se, p = p, row.names = NULL)
}

#' @export
print.mnlfa_ic <- function(x, digits = 3, ...) {
  d <- x$data; f <- x$stage1
  cat("MNLFA-IC: anchor-free DIF analysis\n")
  cat(sprintf("  %d persons, %d items; person covariates: %s\n", d$N, d$J, paste(d$covariates, collapse = ", ")))
  des <- x$coding$design
  cat("  item covariates:", if (is.null(des)) "none" else paste(names(des$coding), collapse = ", "), "\n")
  cat(sprintf("  form: %s; Stage 1 %s (%d Newton steps, %d EM fallbacks)%s\n", x$control$form,
              if (f$convergence$converged) "converged" else "DID NOT CONVERGE",
              f$convergence$newton_iterations, f$convergence$em_fallbacks,
              if (identical(f$vcov_type, "cluster")) sprintf("; cluster-robust (%d clusters)", f$n_clusters) else ""))
  cat("\nImpact on the latent mean (latent SD per unit):\n")
  if (is.null(x$inference)) {
    print(round(x$stage2$alpha, digits))
    cat("\nNo inference yet: run infer(fit).\n")
  } else {
    im <- x$inference$impact
    tab <- data.frame(estimate = im$est, ci_low = im$ci_low, ci_high = im$ci_high, p = im$p, row.names = im$covariate)
    print(round(tab, digits))
    fl <- tapply(x$inference$residual$flag, factor(x$inference$residual$covariate, levels = d$covariates), sum)
    cat(sprintf("\nResidual DIF flagged (BH %g%%): %s\n", 100 * x$control$fdr,
                paste(sprintf("%s %d", names(fl), fl), collapse = ", ")))
  }
  cat("Use summary(), dif_table(), diagnose() for more.\n")
  invisible(x)
}

#' Summarize an mnlfa_ic fit
#'
#' The main results of an analysis: impact on the latent mean and
#' log-variance, item-covariate effects on DIF, joint tests for categorical
#' item covariates, impact for an item of each category, and flagged residual
#' DIF. See `vignette("MNLFAic")` for how to read each table.
#'
#' @param object An `mnlfa_ic` fit.
#' @param ... Unused.
#' @return An object of class `summary.mnlfa_ic` with tables `impact`,
#'   `variance`, `effects` (item-covariate effects, including implied
#'   categories), `joint` (joint tests per categorical item covariate),
#'   `category_impact` (impact for an item of each category with no residual
#'   DIF) and `flagged` (residual DIF flagged by Benjamini-Hochberg).
#' @example man-roxygen/example-fit.R
#' @examples
#' s <- summary(fit)
#' s
#' s$effects
#' @export
summary.mnlfa_ic <- function(object, ...) {
  .need_inference(object)
  inf <- object$inference
  sel <- function(df, cols) if (is.null(df)) NULL else df[, intersect(cols, names(df)), drop = FALSE]
  cols <- c("covariate", "term", "est", "se", "ci_low", "ci_high", "p")
  effects <- if (!is.null(inf$category_effects)) inf$category_effects else
             if (!is.null(inf$gamma)) cbind(covariate = rep(object$data$covariates, each = nrow(inf$gamma) / object$data$K), inf$gamma) else NULL
  flagged <- dif_table(object, flagged_only = TRUE)
  structure(list(call = object$call, data = object$data[c("N", "J", "K", "items", "covariates")],
                 form = object$control$form, units = .units(object), level = object$control$level,
                 fdr = object$control$fdr, convergence = object$stage1$convergence,
                 vcov_type = object$stage1$vcov_type,
                 impact = sel(inf$impact, cols), variance = .variance_table(object),
                 effects = sel(effects, cols), joint = inf$joint, category_impact = sel(inf$impact_by_category, cols),
                 flagged = flagged),
            class = "summary.mnlfa_ic")
}

#' @export
print.summary.mnlfa_ic <- function(x, digits = 3, ...) {
  rnd <- function(df) { if (is.null(df)) return(NULL); num <- vapply(df, is.numeric, logical(1)); df[num] <- lapply(df[num], round, digits); df }
  cat("MNLFA-IC: anchor-free DIF analysis\n")
  cat(sprintf("%d persons, %d items, %d person covariate(s); form: %s; %g%% intervals%s\n", x$data$N, x$data$J,
              x$data$K, x$form, 100 * x$level, if (identical(x$vcov_type, "cluster")) "; cluster-robust" else ""))
  if (!x$convergence$converged) cat("WARNING: Stage 1 did not converge (", x$convergence$stop_reason, ").\n", sep = "")
  cat("\nImpact on the latent mean (latent SD per unit of the covariate):\n"); print(rnd(x$impact), row.names = FALSE)
  cat("\nImpact on the latent log-variance (Wald):\n"); print(rnd(x$variance), row.names = FALSE)
  if (!is.null(x$effects)) {
    cat("\nItem-covariate effects on DIF (", x$units, "):\n", sep = ""); print(rnd(x$effects), row.names = FALSE)
  }
  if (!is.null(x$joint)) { cat("\nJoint tests for categorical item covariates:\n"); print(rnd(x$joint), row.names = FALSE) }
  if (!is.null(x$category_impact)) {
    cat("\nImpact for an item of each category with no residual DIF (latent SD):\n")
    print(rnd(x$category_impact), row.names = FALSE)
  }
  cat(sprintf("\nResidual DIF flagged (Benjamini-Hochberg, FDR %g%%; %s):\n", 100 * x$fdr, x$units))
  if (nrow(x$flagged) == 0) cat("  none\n") else
    print(rnd(x$flagged[, c("item", "covariate", "est", "ci_low", "ci_high", "p_bh", "total")]), row.names = FALSE)
  invisible(x)
}

#' Coefficients of an mnlfa_ic fit
#'
#' Point estimates from an `mnlfa_ic` fit, one type of parameter at a time.
#'
#' @param object An `mnlfa_ic` fit.
#' @param type `"mean"` (default; impact on the latent
#'   mean), `"variance"` (impact on the latent log-variance), `"effects"`
#'   (item-covariate effects on DIF, with implied categories), `"residual"` or
#'   `"total"` (residual and total DIF, items in rows and person covariates in
#'   columns), `"slopes"` or `"intercepts"` (baseline item parameters, for the
#'   reference person with all covariates at 0).
#' @param ... Unused.
#' @return A named vector, or a matrix for `"residual"` and `"total"`.
#' @example man-roxygen/example-fit.R
#' @examples
#' coef(fit)              # impact on the latent mean
#' coef(fit, "effects")   # item-covariate effects on DIF
#' coef(fit, "residual")  # residual DIF, items x person covariates
#' @export
coef.mnlfa_ic <- function(object, type = c("mean", "variance", "effects", "residual", "total", "slopes", "intercepts"), ...) {
  type <- match.arg(type)
  s2 <- object$stage2; f <- object$stage1
  switch(type,
    mean = s2$alpha,
    variance = f$eta,
    effects = {
      des <- object$coding$design
      if (is.null(des)) return(numeric(0))
      out <- unlist(lapply(object$data$covariates, function(c) {
        e <- .full_category_effects(s2$gamma[[c]], des); stats::setNames(e, paste(c, names(e), sep = ":")) }))
      out
    },
    residual = s2$eps,
    total = s2$nu_total,
    slopes = f$alpha.vec,
    intercepts = f$beta.vec)
}

#' Confidence intervals for an mnlfa_ic fit
#'
#' Confidence intervals from an `mnlfa_ic` fit, one type of parameter at a
#' time.
#'
#' @param object An `mnlfa_ic` fit with inference.
#' @param parm `"mean"` (default; impact on the latent mean),
#'   `"variance"` (impact on the latent log-variance), `"effects"`
#'   (item-covariate effects on DIF, with implied categories), `"residual"` or
#'   `"total"` (residual and total DIF, one row per item and person
#'   covariate).
#' @param level Confidence level. The bootstrap intervals are computed when the
#'   model is fitted, at the level set by `mnlfa_control(level = )` (default
#'   0.95), so `level` must equal that value. For intervals at a different
#'   level, rerun the bootstrap with `infer(object, level = )`, which reuses the
#'   Stage-1 fit. The Wald intervals for `"variance"` can be computed at any
#'   level. As a result, `level` can only be changed here for `"variance"`.
#' @param ... Unused.
#' @return A data frame with the covariate (and the term or item, where
#'   relevant), the estimate `est`, and the confidence limits `ci_low` and
#'   `ci_high`.
#' @example man-roxygen/example-fit.R
#' @examples
#' confint(fit)                            # impact on the latent mean
#' confint(fit, "effects")
#' confint(fit, "variance", level = 0.9)   # Wald intervals, any level
#' @export
confint.mnlfa_ic <- function(object, parm = c("mean", "variance", "effects", "residual", "total"),
                             level = object$control$level, ...) {
  parm <- match.arg(parm)
  if (parm == "variance") return(.variance_table(object, level)[, c("covariate", "est", "ci_low", "ci_high")])
  .need_inference(object)
  if (!isTRUE(all.equal(level, object$control$level)))
    stop("The inference was run at level ", object$control$level, "; use infer(object, level = ", level, ").", call. = FALSE)
  inf <- object$inference
  out <- switch(parm,
    mean = inf$impact[, c("covariate", "est", "ci_low", "ci_high")],
    effects = if (!is.null(inf$category_effects)) inf$category_effects[, c("covariate", "term", "est", "ci_low", "ci_high")] else
              if (!is.null(inf$gamma)) inf$gamma[, c("term", "est", "ci_low", "ci_high")] else NULL,
    residual = dif_table(object)[, c("item", "covariate", "est", "ci_low", "ci_high")],
    total = dif_table(object, type = "total")[, c("item", "covariate", "est", "ci_low", "ci_high")])
  if (!is.null(out)) rownames(out) <- NULL
  out
}

#' @export
vcov.mnlfa_ic <- function(object, ...) {
  f <- object$stage1; d <- object$data
  V <- f$vcov
  others <- setdiff(seq_len(d$J), f$anchor)
  nm <- c(paste0("impact:", d$covariates), paste0("logvar:", d$covariates), paste0("slope:", d$items),
          paste0("intercept:", d$items),
          as.vector(t(outer(d$items[others], d$covariates, function(i, c) paste0("dif:", i, ":", c)))))
  dimnames(V) <- list(nm, nm)
  V
}

#' Item-level DIF table
#'
#' One row per item and person covariate, with residual (item-specific) DIF or
#' total DIF, its confidence interval and tests, and the part explained by the
#' item covariates.
#'
#' @param object An `mnlfa_ic` fit with inference.
#' @param type `"residual"` (default) or `"total"` DIF for the estimates,
#'   confidence intervals and p-values. The `flag` column always refers to
#'   residual DIF.
#' @param flagged_only Only rows flagged by Benjamini-Hochberg (residual DIF).
#' @return A data frame with one row per item and person covariate and columns:
#'   * `item`, `covariate`: the item and the person covariate;
#'   * `est`, `se`, `ci_low`, `ci_high`, `p`: estimate, bootstrap standard
#'     error, confidence interval and p-value of the DIF chosen by `type`;
#'   * `p_bh` (residual DIF only): Benjamini-Hochberg-adjusted p-value;
#'   * `flag`: residual DIF flagged at the false discovery rate set in
#'     [mnlfa_control()];
#'   * `total`, `explained`, `residual`: point estimates of total DIF and its
#'     two parts (total = explained + residual);
#'   * the item covariates of each item (if any).
#' @example man-roxygen/example-fit.R
#' @examples
#' dif_table(fit)
#' dif_table(fit, type = "total")
#' dif_table(fit, flagged_only = TRUE)
#' @export
dif_table <- function(object, type = c("residual", "total"), flagged_only = FALSE) {
  type <- match.arg(type)
  .need_inference(object)
  inf <- object$inference; d <- object$data
  r <- inf$residual; t <- inf$total
  src <- if (type == "residual") r else t
  out <- data.frame(item = d$items[src$item], covariate = src$covariate, est = src$est, se = src$se,
                    ci_low = src$ci_low, ci_high = src$ci_high, p = src$p, stringsAsFactors = FALSE)
  if (type == "residual") out$p_bh <- r$p_bh
  out$flag <- r$flag
  out$total <- t$est; out$residual <- r$est; out$explained <- t$est - r$est
  des <- object$coding$design
  if (!is.null(des)) out <- cbind(out, des$data[src$item, , drop = FALSE], row.names = NULL)
  if (flagged_only) out <- out[out$flag, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Identification diagnostics
#'
#' Two descriptions of how much residual DIF the identification can absorb,
#' computed from the estimated slopes and the item covariates: a tolerance
#' table (how many same-signed residual-DIF items each group of items could
#' contain before identification fails) and a report on the Stage-2 design.
#'
#' Neither uses the estimated residual DIF. The sparsity assumption concerns
#' the true, unknown pattern of residual DIF, and cannot be tested from the
#' data. To check whether a hypothesized pattern of residual DIF would be
#' identified, and by what margin, use [ml1_condition_check()].
#'
#' @param object An `mnlfa_ic` fit.
#' @param groups Named list of item-index vectors for the tolerance table.
#'   Default: the levels of each categorical item covariate and the lowest and
#'   highest 10% of each continuous one.
#' @param tolerance Compute the tolerance table.
#' @param n_random Number of random orders in the tolerance table.
#' @return A list of class `mnlfa_ic_diagnostics` with
#'   * `tolerance` (or `NULL` without item covariates or with
#'     `tolerance = FALSE`): a data frame with, per group of items, the number
#'     of items (`n_items`) and the number of same-signed residual-DIF items
#'     the group could absorb when they are added in order of decreasing slope
#'     (`by_lambda`, typically the worst case) and over random orders (median
#'     `random_median` and 5th percentile `random_p05`);
#'   * `design`: a list with the correlations of the slopes with the item
#'     covariates, the condition number and the highest-leverage items.
#' @example man-roxygen/example-fit.R
#' @examples
#' dg <- diagnose(fit, n_random = 50)
#' dg
#' dg$tolerance
#' @export
diagnose <- function(object, groups = NULL, tolerance = TRUE, n_random = 200) {
  d <- object$data; f <- object$stage1
  W <- .W_of(object); metric <- .form_metric(object$control$form)
  tol <- NULL
  if (tolerance && !is.null(W)) {
    if (is.null(groups)) {
      groups <- list(); des <- object$coding$design
      for (v in names(des$coding)) {
        cd <- des$coding[[v]]; xv <- des$data[[v]]
        if (cd$type == "categorical") {
          for (lv in cd$levels) groups[[paste0(v, "=", lv)]] <- which(as.character(xv) == lv)
        } else {
          q <- stats::quantile(xv, c(0.1, 0.9))
          groups[[paste0(v, " lowest 10%")]] <- which(xv <= q[1])
          groups[[paste0(v, " highest 10%")]] <- which(xv >= q[2])
        }
      }
    }
    tol <- .ml1_tolerance_table(f$alpha.vec, W, groups, sign = 1, metric = metric, n_random = n_random)
  }
  design <- .ml1_design_report(f$alpha.vec, W, metric)
  design$leverage_top$item <- d$items[design$leverage_top$item]
  rownames(design$leverage_top) <- NULL
  structure(list(tolerance = tol, design = design), class = "mnlfa_ic_diagnostics")
}

#' @export
print.mnlfa_ic_diagnostics <- function(x, digits = 3, ...) {
  if (!is.null(x$tolerance)) {
    cat("Tolerance: same-signed residual-DIF items each group could absorb before identification fails\n",
        "(by_lambda: highest-slope items first; random_median / random_p05: over random orders):\n", sep = "")
    print(x$tolerance, row.names = FALSE)
    cat("\n")
  }
  ds <- x$design
  cat(sprintf("Stage-2 design: %d items, %d columns; condition number %.2f\n", ds$J, ds$p, ds$condition_number))
  if (!is.null(ds$cor_slope_with_W)) {
    cat("Correlation of slopes with item covariates:\n")
    print(round(ds$cor_slope_with_W, digits))
  }
  cat(sprintf("Highest-leverage items (mean leverage %.3f):\n", ds$leverage_mean))
  lt <- ds$leverage_top; lt$leverage <- round(lt$leverage, digits)
  print(lt, row.names = FALSE)
  invisible(x)
}
