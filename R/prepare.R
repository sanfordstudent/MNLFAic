# Input checks and preparation for mnlfa_ic() and simulate_mnlfa_ic().

# Responses: N x J numeric 0/1/NA matrix with item names.
.prepare_responses <- function(responses) {
  if (!is.matrix(responses) && !is.data.frame(responses))
    stop("`responses` must be a matrix or data frame of 0/1 responses.", call. = FALSE)
  nm <- colnames(responses)
  dat <- as.matrix(responses)
  if (is.logical(dat)) storage.mode(dat) <- "double"
  if (!is.numeric(dat)) stop("`responses` must contain only 0, 1 (or TRUE/FALSE) and NA.", call. = FALSE)
  storage.mode(dat) <- "double"
  if (any(!(dat[!is.na(dat)] %in% c(0, 1)))) stop("`responses` must contain only 0, 1 and NA.", call. = FALSE)
  if (ncol(dat) < 3) stop("At least 3 items are needed.", call. = FALSE)
  named <- !is.null(nm)
  if (!named) nm <- paste0("Item", seq_len(ncol(dat)))
  colnames(dat) <- nm
  rownames(dat) <- NULL
  attr(dat, "named") <- named
  dat
}

# Item identifiers in item_data: an `item` column, or row names unless they are
# just numbers (e.g. left over from subsetting a data frame). NULL if none.
.item_ids <- function(item_data) {
  if (is.null(item_data)) return(NULL)
  if ("item" %in% names(item_data)) return(as.character(item_data$item))
  rn <- rownames(item_data)
  if (!is.null(rn) && !all(grepl("^[0-9]+$", rn))) return(rn)
  NULL
}

# Person covariates from a one-sided formula: model.matrix() without the
# intercept (factors as treatment dummies; interactions as products).
.prepare_person <- function(person, person_data, N) {
  if (!inherits(person, "formula") || length(person) != 2)
    stop("`person` must be a one-sided formula, e.g. ~ female + age.", call. = FALSE)
  if (!is.data.frame(person_data)) stop("`person_data` must be a data frame.", call. = FALSE)
  if (nrow(person_data) != N)
    stop("`person_data` has ", nrow(person_data), " rows but `responses` has ", N, ".", call. = FALSE)
  tt <- stats::terms(person, data = person_data)
  if (length(attr(tt, "term.labels")) == 0) stop("`person` must contain at least one covariate.", call. = FALSE)
  attr(tt, "intercept") <- 1L
  mf <- stats::model.frame(tt, data = person_data, na.action = stats::na.pass)
  na_rows <- !stats::complete.cases(mf)
  list(tt = tt, mf = mf, na_rows = na_rows)
}

.person_matrix <- function(prep, keep) {
  mf <- prep$mf[keep, , drop = FALSE]
  X <- stats::model.matrix(prep$tt, mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL
  rownames(X) <- NULL
  sds <- apply(X, 2, stats::sd)
  if (any(!(sds > 0)))
    stop("Person covariate(s) with no variation: ", paste(colnames(X)[!(sds > 0)], collapse = ", "), ".", call. = FALSE)
  qx <- qr(X)
  if (qx$rank < ncol(X))
    stop("Person covariates are collinear: ", paste(colnames(X)[qx$pivot[(qx$rank + 1):ncol(X)]], collapse = ", "),
         " can be written as a combination of the others.", call. = FALSE)
  binary <- apply(X, 2, function(v) all(v %in% c(0, 1)))
  off <- !binary & abs(colMeans(X)) > sds
  if (any(off))
    warning("Person covariate(s) far from centered: ", paste(colnames(X)[off], collapse = ", "),
            ". Impact and baseline item parameters refer to the person with all covariates at 0; ",
            "consider centering them.", call. = FALSE)
  X
}

# Item covariates from a one-sided formula: numeric -> continuous (centered,
# optionally standardized); factor/character/logical -> effect-coded.
.prepare_items <- function(item, item_data, item_names, item_omit = NULL, standardize = TRUE) {
  J <- length(item_names)
  if (is.null(item)) return(NULL)
  if (!inherits(item, "formula") || length(item) != 2)
    stop("`item` must be a one-sided formula, e.g. ~ aoa + pos.", call. = FALSE)
  if (is.null(item_data)) stop("`item` was given but `item_data` is missing.", call. = FALSE)
  if (!is.data.frame(item_data)) stop("`item_data` must be a data frame.", call. = FALSE)
  if (nrow(item_data) != J)
    stop("`item_data` has ", nrow(item_data), " rows but there are ", J, " items.", call. = FALSE)
  ids <- .item_ids(item_data)
  if (!is.null(ids) && !identical(ids, item_names))
    stop("The items in `item_data` (its `item` column or row names) do not match the column names of ",
         "`responses` in the same order.", call. = FALSE)
  tt <- stats::terms(item, data = item_data)
  if (any(attr(tt, "order") > 1))
    stop("Interactions in the `item` formula are not supported; create interaction terms manually as columns of `item_data`.", call. = FALSE)
  if (length(attr(tt, "term.labels")) == 0) return(NULL)
  mf <- stats::model.frame(tt, data = item_data, na.action = stats::na.pass)
  if (anyNA(mf)) stop("Item covariates contain missing values.", call. = FALSE)
  vars <- names(mf)
  is_cat <- vapply(mf, function(v) is.factor(v) || is.character(v) || is.logical(v), logical(1))
  items_df <- as.data.frame(lapply(mf, function(v) if (is.logical(v)) factor(v) else v), check.names = FALSE)
  if (!is.null(item_omit)) {
    bad <- setdiff(names(item_omit), vars[is_cat])
    if (length(bad)) stop("`item_omit` refers to non-categorical or unknown item covariates: ",
                          paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  des <- .build_item_design(items_df, continuous = vars[!is_cat], categorical = vars[is_cat],
                            omit = if (is.null(item_omit)) list() else item_omit, standardize = standardize)
  rownames(des$W) <- item_names
  if (J < ncol(des$W) + 2)
    stop("Too few items (", J, ") for ", ncol(des$W), " item-covariate columns.", call. = FALSE)
  if (qr(cbind(1, des$W))$rank < ncol(des$W) + 1)
    stop("Item covariates are collinear (with each other or with a constant).", call. = FALSE)
  des$data <- items_df
  des
}

.prepare_cluster <- function(cluster, person_data) {
  if (is.null(cluster)) return(NULL)
  if (is.character(cluster) && length(cluster) == 1 && cluster %in% names(person_data)) return(person_data[[cluster]])
  cluster
}

# Items whose smallest correct/incorrect count, overall or within a group of a
# binary person covariate, is below `min_count`.
.sparse_items <- function(dat, x, min_count = 5) {
  cnt <- function(rows) {
    d <- dat[rows, , drop = FALSE]
    pmin(colSums(d == 1, na.rm = TRUE), colSums(d == 0, na.rm = TRUE))
  }
  m <- cnt(rep(TRUE, nrow(dat)))
  for (k in seq_len(ncol(x))) {
    v <- x[, k]
    if (all(v %in% c(0, 1))) m <- pmin(m, cnt(v == 1), cnt(v == 0))
  }
  which(m < min_count)
}

.default_start <- function(dat, K) {
  J <- ncol(dat)
  list(a = rep(1, J), b = stats::qlogis(pmin(pmax(colMeans(dat, na.rm = TRUE), .01), .99)),
       g = matrix(0, J, K), beta = rep(0, K), eta = rep(0, K))
}

# Everything mnlfa_ic() needs, checked.
.prepare_all <- function(responses, person_data, person, item_data, item, cluster, control, na_action) {
  dat <- .prepare_responses(responses)
  N0 <- nrow(dat)
  # without column names in `responses`, take item names from item_data
  ids <- if (!is.null(item)) .item_ids(item_data) else NULL
  if (!attr(dat, "named") && !is.null(ids) && length(ids) == ncol(dat)) colnames(dat) <- ids
  attr(dat, "named") <- NULL
  pp <- .prepare_person(person, person_data, N0)
  cl <- .prepare_cluster(cluster, person_data)
  if (!is.null(cl) && length(cl) != N0)
    stop("`cluster` must have one value per person (", N0, ").", call. = FALSE)
  drop_na <- pp$na_rows | (if (is.null(cl)) FALSE else is.na(cl))
  if (any(drop_na)) {
    if (na_action == "fail")
      stop(sum(drop_na), " person(s) have missing person covariates or cluster. Remove them, or use na_action = \"omit\".",
           call. = FALSE)
    message("Dropping ", sum(drop_na), " person(s) with missing person covariates or cluster.")
  }
  no_resp <- rowSums(!is.na(dat)) == 0
  if (any(no_resp & !drop_na)) message("Dropping ", sum(no_resp & !drop_na), " person(s) with no responses.")
  keep <- !drop_na & !no_resp
  dat <- dat[keep, , drop = FALSE]
  x <- .person_matrix(pp, keep)
  if (!is.null(cl)) cl <- cl[keep]
  p_item <- colMeans(dat, na.rm = TRUE)
  flat <- !is.finite(p_item) | p_item %in% c(0, 1)
  if (any(flat)) stop("Item(s) without variation in the responses: ", paste(colnames(dat)[flat], collapse = ", "), ".", call. = FALSE)
  sp <- .sparse_items(dat, x)
  if (length(sp))
    warning("Item(s) with fewer than 5 correct or incorrect responses overall or within a group of a binary ",
            "person covariate: ", paste(colnames(dat)[sp], collapse = ", "),
            ". Their DIF may be poorly determined; consider screening near-ceiling or near-floor items.", call. = FALSE)
  des <- .prepare_items(item, item_data, colnames(dat), control$item_omit, control$item_standardize)
  if (is.numeric(control$anchor) && control$anchor > ncol(dat)) stop("`anchor` exceeds the number of items.", call. = FALSE)
  list(dat = dat, x = x, cluster = cl, design = des, kept = which(keep), person_terms = pp$tt,
       person_levels = stats::.getXlevels(pp$tt, pp$mf), N_input = N0)
}
