# Replicate weights for IEA studies (TIMSS, PIRLS).

#' Jackknife replicate weights for TIMSS and PIRLS
#'
#' Builds the jackknife (JK2) replicate weights of IEA studies from the
#' sampling-zone variables, for use as `replicate_weights` in [mnlfa_ic()].
#' IEA data files provide the zone (`JKZONE`) and the replicate indicator
#' (`JKREP`) instead of the replicate weights themselves. Part of the
#' experimental support for sampling weights.
#'
#' Schools are paired into jackknife zones. With `method = "JK2"` (IEA studies
#' from 2015 onward), each zone gives two replicates: in the first, the weights
#' of students with `jkrep == 1` in that zone are doubled and those with
#' `jkrep == 0` are set to 0; in the second, the reverse. Students in other
#' zones keep their weights. With `method = "JK2-half"` (earlier cycles), each
#' zone gives only the first of these. Use the same `method` in
#' [mnlfa_control()], which sets the variance factor (1/2 and 1).
#'
#' When several countries are pooled, zones are numbered within each country,
#' so replicate r perturbs zone r in every country at once, as in the IEA's
#' pooled analyses.
#'
#' @param weights Total student weights (e.g., `TOTWGT`), or senate or house
#'   weights; [mnlfa_ic()] rescales them in any case.
#' @param jkzone Jackknife zone of each student (e.g., `JKZONE`), integers from
#'   1.
#' @param jkrep Replicate indicator of each student (e.g., `JKREP`), 0 or 1.
#' @param method `"JK2"` (default; both halves of every zone) or `"JK2-half"`.
#' @param n_zones Number of zones. Default: the largest value of `jkzone`
#'   (zones without students give replicates equal to the full sample, which
#'   add nothing to the variance).
#' @return A matrix with one row per student and one column per replicate:
#'   2 x `n_zones` columns for `"JK2"` (zone h in columns 2h - 1 and 2h) and
#'   `n_zones` for `"JK2-half"`.
#' @references Foy, P., & LaRoche, S. (2020). Estimating standard errors in the
#'   TIMSS 2019 results. In M. O. Martin, M. von Davier, & I. V. S. Mullis
#'   (Eds.), *Methods and Procedures: TIMSS 2019 Technical Report*
#'   (chap. 14). TIMSS & PIRLS International Study Center, Boston College.
#' @examples
#' w <- c(10, 12, 9, 11, 8, 10)
#' ilsa_replicates(w, jkzone = c(1, 1, 1, 2, 2, 2), jkrep = c(0, 1, 1, 0, 1, 0))
#' @export
ilsa_replicates <- function(weights, jkzone, jkrep, method = c("JK2", "JK2-half"), n_zones = NULL) {
  method <- match.arg(method)
  N <- length(weights)
  if (length(jkzone) != N || length(jkrep) != N)
    stop("`weights`, `jkzone` and `jkrep` must have the same length.", call. = FALSE)
  if (anyNA(weights) || anyNA(jkzone) || anyNA(jkrep)) stop("Missing values in `weights`, `jkzone` or `jkrep`.", call. = FALSE)
  if (!all(jkrep %in% c(0, 1))) stop("`jkrep` must be 0 or 1.", call. = FALSE)
  if (!is.numeric(jkzone) || any(jkzone < 1) || any(jkzone != round(jkzone)))
    stop("`jkzone` must contain whole numbers from 1.", call. = FALSE)
  if (is.null(n_zones)) n_zones <- max(jkzone)
  if (max(jkzone) > n_zones) stop("`jkzone` exceeds `n_zones`.", call. = FALSE)
  per_zone <- if (method == "JK2") 2 else 1
  R <- matrix(weights, N, per_zone * n_zones)
  for (h in unique(jkzone)) {
    z <- jkzone == h
    col <- per_zone * (h - 1) + 1
    R[z, col] <- 2 * weights[z] * jkrep[z]
    if (method == "JK2") R[z, col + 1] <- 2 * weights[z] * (1 - jkrep[z])
  }
  R
}
