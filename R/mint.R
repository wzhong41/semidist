MIsemi <- function(X, Y, ks, perm_term = NULL, dist_x = NULL, n_x = NULL) {
  n <- nrow(X)
  R <- ncol(Y)
  
  kmax <- max(ks)
  
  rho_class <- matrix(nrow = n, ncol = kmax)
  for (r in 1:R) {
    ind_r <- as.logical(Y[, r])
    if (nrow(X[ind_r, , drop = FALSE]) == 0) {
      next
    } else {
      rho_class[ind_r, ] <- knn.dist(X[ind_r, , drop = FALSE], k = kmax)
    }
  }
  
  if (is.null(perm_term)) {
    freq <- colSums(Y)
    dist_x <- as.matrix(dist(X))
    n_x <- drop(Y %*% freq)
  }
  
  # for (k in ks) {
  #   for (i in 1:n) {
  #     m[i, as.character(k)] <- sum(dist_x[i, -i] <= rho_class[i, k])
  #   }
  # }
  m <- matrix(nrow = n, ncol = length(ks), dimnames = list(NULL, ks))
  for (i in 1:n) {
    m[i, ] <- colSums(outer(dist_x[i, -i], rho_class[i, ks], "<="))
  }
  
  if (is.null(perm_term)) {
    MI <- digamma(n) - mean(digamma(n_x)) + digamma(ks) - colMeans(digamma(m))
  } else {
    MI <- perm_term - colMeans(digamma(m))
  }
  names(MI) <- ks
  
  MI
}

#' Mutual information independence test (categorical-continuous case)
#'
#' @description Implement the mutual information independence test (MINT)
#'   (Berrett and Samworth, 2019), but with some modification in estimating the
#'   mutual informaion (MI) between a categorical random variable and a
#'   continuous variable. The modification is based on the idea of Ross (2014).
#'
#'   `MINTsemiperm()` implements the permutation independence test via
#'   mutual information, but the parameter `k` should be pre-specified.
#'   
#'   `MINTsemiauto()` automatically selects an appropriate `k` based on a 
#'   data-driven procedure, and conducts `MINTsemiperm()` with the `k` chosen.
#'
#' @param X Data of multivariate continuous variables, which should be an
#'   \eqn{n}-by-\eqn{p} matrix, or, a vector of length \eqn{n} (for univariate
#'   variable).
#' @param y Data of categorical variables, which should be a factor of length
#'   \eqn{n}.
#' @param k Number of nearest neighbor. See References for details.
#' @param kmax Maximum `k` in the automatic search for optimal `k`.
#' @param B,B1,B2 Number of permutations to use. Defaults to 1000.
#' 
#' @returns A list with class `"indtest"` containing the following components
#' * `method`: name of the test;
#' * `name_data`: names of the `X` and `y`;
#' * `n`: sample size of the data;
#' * `num_perm`: number of replications in permutation test;
#' * `stat`: test statistic;
#' * `pvalue`: computed p-value.
#' 
#' For `MINTsemiauto()`, the list also contains
#' * `kmax`: maximum `k` in the automatic search for optimal `k`;
#' * `kopt`: optimal `k` chosen.
#'
#' @examples
#' X <- mtcars[, c("mpg", "disp", "drat", "wt")]
#' y <- factor(mtcars[, "am"])
#' 
#' MINTsemiperm(X, y, 5)
#' MINTsemiauto(X, y, kmax = 32)
#'
#' @references 
#' 1. Berrett, Thomas B., and Richard J. Samworth. "Nonparametric
#'   independence testing via mutual information." *Biometrika* 106, no. 3
#'   (2019): 547-566. 
#' 1. Ross, Brian C. "Mutual information between discrete and
#'   continuous data sets." *PloS one* 9, no. 2 (2014): e87357.
#' @export
MINTsemiperm <- function(X, y, k, B = 1000) {
  n <- nrow(X)
  Y <- switch_cat_repr(y)
  name_data <- paste(
    deparse(substitute(X)), "and", deparse(substitute(y))
  )
  
  dist_x <- as.matrix(dist(X))
  freq <- colSums(Y)
  n_x <- drop(Y %*% freq)
  perm_term <- digamma(n) - mean(digamma(n_x)) + digamma(k)
  
  name_method <- paste0(
    "MINT Independence Test (Permutation Test with B = ",
    B, " and k = ", k, ")"
  )
  H <- MIsemi(X, Y, k, perm_term, dist_x, n_x)
  
  Hp <- rep(0, B)
  for (b in 1:B) {
    Yperm <- Y[sample(n), ]
    Hp[b] <- MIsemi(X, Yperm, k, perm_term, dist_x, n_x)
  }
  pvalue <- (1 + sum(Hp >= H)) / (B + 1)
  
  indtest <- list(
    method = name_method,
    name_data = name_data,
    n = n,
    num_perm = B,
    stat = H,
    pvalue = pvalue
  )
  class(indtest) <- "indtest"
  
  indtest
}

#' @rdname MINTsemiperm
#' @export
MINTsemiauto <- function(X, y, kmax, B1 = 1000, B2 = 1000) {
  if (!is.matrix(X) && is.numeric(X)) {
    X <- as.matrix(X, ncol = 1)
  }
  n <- nrow(X)
  Y <- switch_cat_repr(y)
  
  dist_x <- as.matrix(dist(X))
  freq <- colSums(Y)
  kmax <- min(min(freq)-1, kmax)
  n_x <- drop(Y %*% freq)
  ks <- 1:kmax
  perm_term <- digamma(n) - mean(digamma(n_x)) + digamma(ks)
  
  Hp <- matrix(rep(0, B1 * kmax), ncol = kmax)
  for (b in 1:B1) {
    Y1 <- Y[sample(n), ]
    Y2 <- Y[sample(n), ]
    
    # t_start <- proc.time()
    Hp[b, ] <- MIsemi(X, Y1, ks, perm_term, dist_x, n_x) - MIsemi(X, Y2, ks, perm_term, dist_x, n_x)
    # t <- proc.time() - t_start
  }
  MSE <- (1 / B1) * colSums(Hp^2)
  kopt <- which.min(MSE)
  
  
  indtest <- MINTsemiperm(X, y, kopt, B = B2)
  
  indtest$kopt <- kopt
  indtest$kmax <- kmax
  
  indtest
}