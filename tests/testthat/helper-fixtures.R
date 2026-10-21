# Fixtures: outputs of the research code (chen-item-cov-ext/R/) on fixed
# inputs, made by package/dev/make_fixtures.R in the project folder. The
# equivalence tests compare the package against them.
fixture <- function(name) {
  path <- testthat::test_path("fixtures", paste0(name, ".rds"))
  if (!file.exists(path)) testthat::skip(paste("fixture not available:", name))
  readRDS(path)
}

start_values <- function(dat, K) {
  J <- ncol(dat)
  list(a = rep(1, J), b = stats::qlogis(pmin(pmax(colMeans(dat, na.rm = TRUE), .01), .99)),
       g = matrix(0, J, K), beta = rep(0, K), eta = rep(0, K))
}

# Research-code fit as stored in a fixture, made usable by the package's Stage 3.
as_fit <- function(f, N) { f$N <- N; f }

# Compare two data frames of inference output column by column.
expect_same_table <- function(a, b, tol = 1e-10) {
  testthat::expect_equal(names(a), names(b))
  for (v in names(a)) {
    if (is.numeric(a[[v]])) testthat::expect_equal(a[[v]], b[[v]], tolerance = tol, label = v)
    else testthat::expect_identical(as.character(a[[v]]), as.character(b[[v]]), label = v)
  }
}
