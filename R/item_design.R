# Item-covariate design matrix W (ported from chen-item-cov-ext/R/item_design.R).
#
# Coding conventions (part of the identifying assumption, not only of
# interpretation):
#   - continuous item covariates are centered at the UNWEIGHTED mean over
#     items (optionally also divided by their SD over items);
#   - categorical item covariates are EFFECT-coded (unweighted): one column
#     per non-omitted level, coded 1 for that level, -1 for the omitted level,
#     0 otherwise. The omitted level's effect is implied as minus the sum of
#     the other levels' effects (see .full_category_effects()).

# items       : data frame with one row per item (J rows, in item order)
# continuous  : names of continuous columns
# categorical : names of categorical columns (character or factor)
# omit        : named list giving the omitted level per categorical column
#               (default: the first level)
# standardize : also divide continuous columns by their SD over items
# Returns list(W = J x M matrix with column names, coding = list per variable)
.build_item_design <- function(items, continuous = character(0), categorical = character(0),
                               omit = list(), standardize = FALSE) {
  if (!is.data.frame(items)) stop("`items` must be a data frame.", call. = FALSE)
  vars <- c(continuous, categorical)
  if (length(vars) == 0) stop("No item covariates requested.", call. = FALSE)
  missing_cols <- setdiff(vars, names(items))
  if (length(missing_cols) > 0)
    stop("Item covariate columns not found: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (anyNA(items[vars])) stop("Item covariates contain missing values.", call. = FALSE)

  cols <- list(); coding <- list()
  for (v in continuous) {
    xv  <- as.numeric(items[[v]])
    ctr <- mean(xv)
    scl <- if (standardize) stats::sd(xv) else 1
    if (standardize && !(scl > 0)) stop("Item covariate '", v, "' is constant.", call. = FALSE)
    cols[[v]]   <- (xv - ctr) / scl
    coding[[v]] <- list(type = "continuous", center = ctr, scale = scl, columns = v)
  }
  for (v in categorical) {
    fv   <- factor(items[[v]])
    levs <- levels(fv)
    if (length(levs) < 2) stop("Categorical item covariate '", v, "' has fewer than 2 levels.", call. = FALSE)
    om <- if (!is.null(omit[[v]])) as.character(omit[[v]]) else levs[1]
    if (!(om %in% levs)) stop("Omitted level '", om, "' not found in '", v, "'.", call. = FALSE)
    coded_levs <- setdiff(levs, om)
    col_names  <- paste0(v, "_", coded_levs)
    for (k in seq_along(coded_levs))
      cols[[col_names[k]]] <- ifelse(fv == coded_levs[k], 1, ifelse(fv == om, -1, 0))
    coding[[v]] <- list(type = "categorical", levels = levs, omitted = om,
                        coded_levels = coded_levs, columns = col_names, counts = table(fv))
  }
  W <- do.call(cbind, cols)
  colnames(W) <- names(cols)
  rownames(W) <- rownames(items)
  list(W = W, coding = coding)
}

# Expand one person covariate's coefficients into the full set of effects,
# including the implied omitted-category effect of each categorical variable.
#   gamma  : named numeric vector (names = columns of W)
#   design : output of .build_item_design()
.full_category_effects <- function(gamma, design) {
  out <- numeric(0)
  for (v in names(design$coding)) {
    cd <- design$coding[[v]]
    if (cd$type == "continuous") {
      out[v] <- gamma[[v]]
    } else {
      eff <- gamma[cd$columns]
      names(eff) <- paste0(v, "_", cd$coded_levels)
      out <- c(out, eff)
      out[paste0(v, "_", cd$omitted)] <- -sum(eff)
    }
  }
  out
}
