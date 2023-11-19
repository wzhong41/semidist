#' Switch the representation of a categorical object
#'
#' @description
#' Categorical data with n observations and R levels can typically be represented as two forms in R:
#' a factor with length n, or an n by K indicator matrix with elements being 0 or 1.
#' This function is to switch the form of a categorical object from one to the another.
#'
#' @param obj an object representing categorical data,
#' either a factor or an indicator matrix with each row representing an observation.
#'
#' @return categorical object in the another form.
#' @export
switch_cat_repr <- function(obj) {
  if (any(class(obj) == "factor")) {
    Y <- vapply(
      levels(obj),
      function(r) {
        obj == r
      },
      rep(1, length(obj))
    )
    return(Y)
  } else if (any(class(obj) == "matrix")) {
    R <- ncol(obj)
    if (is.null(colnames(obj))) {
      lev <- 1:R
    } else {
      lev <- colnames(obj)
    }
    ind <- obj > 0
    y <- factor(rep(NA, nrow(obj)), levels = lev)
    for (r in 1:R) {
      y[ind[, r]] <- lev[[r]]
    }
    return(y)
  } else {
    stop('Class of `obj` should be either "factor" or "matrix"!')
  }
}

# Estimate the trace of the covariance matrix and its square (R implementation for verified)
tr_estimate_R_impl <- function(X) {
  n <- nrow(X)

  XXt <- tcrossprod(X)
  Y1 <- sum(diag(crossprod(X))) / n
  Y2 <- (sum(XXt^2) - sum(diag(XXt^2))) / (n * (n - 1))
  Y3 <- (sum(XXt) - Y1 * n) / (n * (n - 1))

  XXt2 <- XXt %*% XXt
  vec_diag_XXt <- rep(diag(XXt), each = n - 1)
  vec_offdiag_XXt <- as.vector(XXt)[!as.logical(diag(n))]
  Y4 <-
    (sum(XXt2) - sum(diag(XXt2)) - 2 * sum(vec_diag_XXt * vec_offdiag_XXt)) / (n * (n - 1) * (n - 2))
  Y5 <- ((n * (n - 1) * Y3)^2 - 2 * n * (n - 1) * Y2 - 4 * n * (n - 1) * (n - 2) * Y4) /
    (n * (n - 1) * (n - 2) * (n - 3))

  Tn1 <- Y1 - Y3
  Tn2 <- Y2 - 2 * Y4 + Y5

  list(
    tr_S = Tn1,
    tr_S2 = Tn2
  )
}
