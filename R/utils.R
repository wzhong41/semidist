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


#' Print Method for Independence Tests Between Categorical and Continuous Variables
#' 
#' @description
#' Printing object of class `"indtest"`, by simple [print] method.
#' 
#' @param x `"indtest"` class object.
#' @param digits minimal number of *significant* digits.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' # Man-made functionally dependent data --------------------------------------
#' n <- 30; R <- 3
#' x <- rep(0, n)
#' x[1:10] <- 0.3; x[11:20] <- 0.2; x[21:30] <- -0.1
#' y <- factor(rep(1:3, each = 10))
#' test <- mv_test(x, y)
#' print(test)
#' test_asym <- mv_test(x, y, test_type = "asym")
#' print(test_asym)
#' 
#' @returns None
#' @export
print.indtest <- function(x, digits = getOption("digits"), ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("Data: ", x$name_data, ",\t", "Sample size = ", x$n, "\n", sep = "")
  if (x$method == "MV Independence Test (Asymptotic Test)") {
    crit_vals <- x$crit_vals
    cat("Test statistic = ", x$stat, "\n")
    cat("Asymptotic critical value: \t")
    for (i in 1:length(crit_vals)) {
      cat(names(crit_vals)[i], " -- ", crit_vals[i], "\t", sep = "")
    }
    cat("\n")
  } else {
    fp <- format.pval(x$pvalue, digits = max(1L, digits - 3L))
    out_pvalue <- paste("p-value", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
    cat("Test statistic = ", x$stat, ",\t", out_pvalue, "\n", sep = "")
  }
  cat("Alternative hypothesis: Two random variables are not independent", sep = "")

  invisible(x)
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
