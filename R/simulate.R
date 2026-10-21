# Data generation from the model (ported from chen-item-cov-ext/R/simulate.R).
#
#   logit P(y_ij = 1 | theta_i) = lambda_j theta_i + nu0_j + sum_c nu_jc x_ci
#   theta_i ~ N(x_i' alpha, exp(x_i' omega))
#   nu_jc = -lambda_j (w_j' gamma_c + eps_jc)   (difficulty metric, default:
#           gamma, eps in latent-SD units, + = harder)
#        or w_j' gamma_c + eps_jc               (logit metric: logits, + = easier)

# x            : N x K person covariates
# lambda, nu0  : J-vectors
# W            : J x M item-covariate matrix (or NULL)
# gamma        : M x K matrix of item-covariate effects (column c = gamma_c)
# eps          : J x K residual DIF
# alpha, omega : K-vectors (impact on mean and log-variance)
# seed         : if not NULL, the data use this seed and the caller's RNG state
#                is restored afterwards
# cluster, icc : optional clustering of persons. With a cluster vector and
#                icc > 0, theta's residual is split into a cluster part and a
#                person part, theta_i = mu_i + sigma_i (sqrt(1 - icc) e_i +
#                sqrt(icc) u_c), so each person's marginal distribution is the
#                model's and only within-cluster dependence is added.
.simulate_itemcov_data <- function(x, lambda, nu0, W = NULL, gamma = NULL, eps, alpha, omega,
                                   metric = c("difficulty", "logit"), seed = NULL,
                                   cluster = NULL, icc = 0) {
  metric <- match.arg(metric)
  .with_seed(seed, {
    N <- nrow(x); J <- length(lambda); K <- ncol(x)
    explained <- if (is.null(W)) matrix(0, J, K) else W %*% gamma
    nu <- if (metric == "difficulty") -lambda * (explained + eps) else explained + eps
    if (is.null(cluster) || icc == 0) {
      theta <- stats::rnorm(N, as.vector(x %*% alpha), sqrt(exp(as.vector(x %*% omega))))
    } else {
      stopifnot(length(cluster) == N, icc > 0, icc < 1)
      e  <- stats::rnorm(N)
      cf <- factor(cluster); u <- stats::rnorm(nlevels(cf))[as.integer(cf)]
      theta <- as.vector(x %*% alpha) + sqrt(exp(as.vector(x %*% omega))) * (sqrt(1 - icc) * e + sqrt(icc) * u)
    }
    lp  <- outer(theta, lambda) + rep(1, N) %*% t(nu0) + x %*% t(nu)
    dat <- matrix(stats::rbinom(N * J, 1, stats::plogis(lp)), N, J)
    list(dat = dat, x = x, theta = theta,
         truth = list(lambda = lambda, nu0 = nu0, nu = nu, W = W, gamma = gamma,
                      eps = eps, alpha = alpha, omega = omega, metric = metric))
  })
}
