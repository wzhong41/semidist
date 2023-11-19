#' Mean Variance (MV) statistics
#' @description Compute the statistics of mean variance (MV) index, which can
#'   measure the dependence between a univariate continuous variable and a
#'   categorical variable. See Cui, Li and Zhong (2015); Cui and Zhong (2019) 
#'   for details.
#'
#' @param x Data of univariate continuous variables, which should be a vector of
#'   length \eqn{n}.
#' @param y Data of categorical variables, which should be a factor of length
#'   \eqn{n}.
#' @param return_mat A boolean. If `FALSE` (the default), only the calculated
#'   statistic is returned. If `TRUE`, also return the matrix of the indicator
#'   for x <= x_i, which is useful for the permutation test.
#'
#' @return The value of the corresponding sample statistic.
#'
#'   If the argument `return_mat` of `mv()` is set as `TRUE`, a list with
#'   elements
#'   * `mv`: the MV index statistic;
#'   * `mat_x`: the matrices of the distances of the indicator for x <= x_i;
#'
#'   will be returned.
#' 
#' @seealso 
#' * [mv_test()] for implementing independence test via MV index;
#' * [mv_sis()] for implementing feature screening via MV index.
#' 
#' @examples
#' x <- mtcars[, "mpg"]
#' y <- factor(mtcars[, "am"])
#' print(mv(x, y))
#'
#' # Man-made independent data -------------------------------------------------
#' n <- 30; R <- 5; prob <- rep(1/R, R)
#' x <- rnorm(n)
#' y <- factor(sample(1:R, size = n, replace = TRUE, prob = prob), levels = 1:R)
#' print(mv(x, y))
#'
#' # Man-made functionally dependent data --------------------------------------
#' n <- 30; R <- 3
#' x <- rep(0, n)
#' x[1:10] <- 0.3; x[11:20] <- 0.2; x[21:30] <- -0.1
#' y <- factor(rep(1:3, each = 10))
#' print(mv(x, y))
#' 
#' @export
mv <- function(x, y, return_mat = FALSE) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  if (is.matrix(x) && ncol(x) > 1) {
    stop("`x` should be a numeric vector")
  } else {
    x <- drop(x)
  }
  
  n <- length(x)
  R <- nlevels(y)
  Y <- switch_cat_repr(y)
  n_all <- colSums(Y)
  
  # Column i stands for the indicator for x <= x[i]
  mat_x <- vapply(
    x, 
    function(x_i) x <= x_i, 
    rep(1, n)
  )
  
  # Marginal CDF of x (evaluated at each x_i) 
  F_0 <- colMeans(mat_x)
  # F_0 <- rank(x, ties.method = "max") / n # Equivalent!
  
  # Conditional CDF of x given y = r (evaluated at each x_i) 
  # Column r stands for the F_r
  F_all <- crossprod(mat_x, Y) %*% diag(1 / n_all)
  
  mv <- sum((F_all - F_0)^2 %*% diag(n_all/n)) / n
  
  if (return_mat) {
    return(
      list(
        mv = mv,
        mat_x = mat_x
      )
    )
  } else {
    return(mv)
  }
}


#' MV independence test (via Permutation Test)
#'
#' @description Implement the MV independence test via permutation test.
#'
#' @param x Data of univariate continuous variables, which should be a vector of
#'   length \eqn{n}.
#' @param y Data of categorical variables, which should be a factor of length
#'   \eqn{n}.
#' @param num_perm The number of replications in permutation test.
#'
#' @returns A list with class `"indtest"` containing the following components
#' * `method`: name of the test;
#' * `name_data`: names of the `x` and `y`;
#' * `n`: sample size of the data;
#' * `num_perm`: number of replications in permutation test;
#' * `stat`: test statistic;
#' * `pvalue`: computed p-value.
#'
#' @examples
#' x <- mtcars[, "mpg"]
#' y <- factor(mtcars[, "am"])
#' test <- mv_test(x, y)
#' print(test)
#'
#' # Man-made independent data -------------------------------------------------
#' n <- 30; R <- 5; prob <- rep(1/R, R)
#' x <- rnorm(n)
#' y <- factor(sample(1:R, size = n, replace = TRUE, prob = prob), levels = 1:R)
#' test <- mv_test(x, y)
#' print(test)
#'
#' # Man-made functionally dependent data --------------------------------------
#' n <- 30; R <- 3
#' x <- rep(0, n)
#' x[1:10] <- 0.3; x[11:20] <- 0.2; x[21:30] <- -0.1
#' y <- factor(rep(1:3, each = 10))
#' test <- mv_test(x, y)
#' print(test)
#'
#' @export
mv_test <- function(x, y, num_perm = 10000) {
  mv_obj <- mv(x, y, return_mat = TRUE)
  
  n <- length(y)
  Y <- switch_cat_repr(y)
  n_all <- colSums(Y)
  name_data <- paste(
    deparse(substitute(x)), "and", deparse(substitute(y))
  )
  
  mv_n <- mv_obj$mv
  mat_x <- mv_obj$mat_x
  F_0 <- colMeans(mat_x)
  
  name_method <- paste0(
    "MV Independence Test (Permutation Test with K = ",
    num_perm, ")"
  )
  num_rej <- 0
  for (perm in 1:num_perm) {
    F_all <- crossprod(mat_x, Y[sample(1:n), ]) %*% diag(1 / n_all)
    
    mv_perm <- sum((F_all - F_0)^2 %*% diag(n_all/n)) / n
    
    if (mv_perm >= mv_n) {
      num_rej <- num_rej + 1
    }
  }
  pvalue <- (num_rej + 1) / (num_perm + 1)
  
  indtest <- list(
    method = name_method,
    name_data = name_data,
    n = n,
    num_perm = num_perm,
    stat = n * mv_n,
    pvalue = pvalue
  )
  class(indtest) <- "indtest"
  
  indtest
}

#' Feature screening via MV Index
#'
#' @description Implement the feature screening for the classification problem
#'   via MV index.
#'
#' @param X Data of multivariate covariates, which should be an
#'   \eqn{n}-by-\eqn{p} matrix.
#' @param y Data of categorical response, which should be a factor of length
#'   \eqn{n}.
#' @param d An integer specifying how many features should be kept after
#'   screening. Defaults to `NULL`. If `NULL`, then it will be set as \eqn{[n /
#'   log(n)]}, where \eqn{[x]} denotes the integer part of x.
#'
#' @param parallel A boolean indicating whether to calculate parallelly via
#'   `furrr::future_map`. Defaults to `FALSE`.
#'
#' @returns A list of the objects about the implemented feature screening:
#' * `measurement`: sample MV index calculated for each single covariate;
#' * `selected`: indicies or names (if avaiable as colnames of `X`) of
#'   covariates that are selected after feature screening;
#' * `ordering`: order of the calculated measurements of each single covariate.
#'   The first one is the largest, and the last is the smallest.
#'
#' @examples
#' X <- mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")]
#' y <- factor(mtcars[, "am"])
#'
#' mv_sis(X, y, d = 4)
#'
#' @export
mv_sis <- function(X, y, d = NULL, parallel = FALSE) {
  p <- ncol(X)
  n <- nrow(X)
  
  if (is.null(d)) {
    d <- floor(n / log(n))
  }
  
  if (is.null(colnames(X))) {
    names_X <- 1:p
    colnames(X) <- names_X
  } else {
    names_X <- colnames(X)
  }
  
  omega <- furrr::future_map_dbl(names_X, ~ mv(X[, .x], y))
  ordering <- order(omega, decreasing = TRUE)
  ordering_names <- names_X[ordering]
  D_hat <- ordering_names[seq(d)]
  
  list(
    measurement = omega,
    selected = D_hat,
    ordering = ordering_names
  )
}
