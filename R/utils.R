# Small internal helpers.

# Evaluate `expr` with the random-number generator seeded by `seed`, then
# restore the caller's RNG state, so functions with a `seed` argument do not
# change the global stream (CRAN policy). With seed = NULL, `expr` uses the
# current stream as is.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  expr
}

# Positions of the parameter blocks in the packed (anchor-removed) Stage-1
# vector: alpha (K), omega (K), slopes (J), intercepts (J), then the DIF of the
# non-anchor items, row-major ((J - 1) K).
.stage1_index <- function(J, K) {
  list(beta = 1:K, eta = (K + 1):(2 * K), a = (2 * K + 1):(2 * K + J),
       b = (2 * K + J + 1):(2 * K + 2 * J),
       g = (2 * K + 2 * J + 1):(2 * K + 2 * J + (J - 1) * K))
}
