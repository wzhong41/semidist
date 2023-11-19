#' Semi-distance covariance and correlation statistics
#'
#' @description Compute the statistics (or sample estimates) of semi-distance
#'   covariance and correlation. The semi-distance correlation is a standardized
#'   version of semi-distance covariance, and it can measure the dependence
#'   between a *multivariate* continuous variable and a categorical variable.
#'   See Details for the definition of semi-distance covariance and
#'   semi-distance correlation.
#'
#' @param X Data of multivariate continuous variables, which should be an
#'   \eqn{n}-by-\eqn{p} matrix, or, a vector of length \eqn{n} (for univariate
#'   variable).
#' @param y Data of categorical variables, which should be a factor of length
#'   \eqn{n}.
#' @param type Type of statistic: `"V"` (the default) or `"U"`. See Details.
#' @param return_mat A boolean. If `FALSE` (the default), only the calculated
#'   statistic is returned. If `TRUE`, also return the matrix of the distances
#'   of X and the divergences of y, which is useful for the permutation test.
#'
#' @details For \eqn{\bm{X} \in \mathbb{R}^{p}} and \eqn{Y \in \{1, 2, \cdots,
#'   R\}}, the (population-level) semi-distance covariance is defined as
#'   \deqn{\mathrm{SDcov}(\bm{X}, Y) =
#'   \mathrm{E}\left[\|\bm{X}-\widetilde{\bm{X}}\|\left(1-\sum_{r=1}^R
#'   I(Y=r,\widetilde{Y}=r)/p_r\right)\right],} where \eqn{p_r = P(Y = r)} and
#'   \eqn{(\widetilde{\bm{X}}, \widetilde{Y})} is an iid copy of \eqn{(\bm{X},
#'   Y)}.
#'   The (population-level) semi-distance correlation is defined as
#'   \deqn{\mathrm{SDcor}(\bm{X}, Y) = \dfrac{\mathrm{SDcov}(\bm{X},
#'   Y)}{\mathrm{dvar}(\bm{X})\sqrt{R-1}},} where \eqn{\mathrm{dvar}(\bm{X})} is
#'   the distance variance (Szekely, Rizzo, and Bakirov 2007) of \eqn{\bm{X}}.
#'
#'   With \eqn{n} observations \eqn{\{(\bm{X}_i, Y_i)\}_{i=1}^{n}}, `sdcov()`
#'   and `sdcor()` can compute the sample estimates for the semi-distance
#'   covariance and correlation.
#'
#'   If `type = "V"`, the semi-distance covariance statistic is computed as a
#'   V-statistic, which takes a very similar form as the energy-based statistic
#'   with double centering, and is always non-negative. Specifically,
#'   \deqn{\text{SDcov}_n(\bm{X}, y) = \frac{1}{n^2} \sum_{k=1}^{n}
#'   \sum_{l=1}^{n} A_{kl} B_{kl},}
#'   where \deqn{A_{kl} = a_{kl} - \bar{a}_{k.} - \bar{a}_{.l} + \bar{a}_{..}}
#'   is the double centering (Szekely, Rizzo, and Bakirov 2007) of
#'   \eqn{a_{kl} = \| \bm{X}_k - \bm{X}_l \|,} and \deqn{B_{kl} =
#'   1 - \sum_{r=1}^{R} I(Y_k = r) I(Y_l = r) / \hat{p}_r} with \eqn{\hat{p}_r =
#'   n_r / n = n^{-1}\sum_{i=1}^{n} I(Y_i = r)}.
#'   The semi-distance correlation statistic is \deqn{\text{SDcor}_n(\bm{X}, y)
#'   = \dfrac{\text{SDcov}_n(\bm{X}, y)}{\text{dvar}_n(\bm{X})\sqrt{R - 1}},}
#'   where \eqn{\text{dvar}_n(\bm{X})} is the V-statistic of distance variance
#'   of \eqn{\bm{X}}.
#'
#'   If `type = "U"`, then the semi-distance covariance statistic is computed as
#'   an ``estimated U-statistic'', which is utilized in the independence test
#'   statistic and is not necessarily non-negative. Specifically,
#'   \deqn{\widetilde{\text{SDcov}}_n(\bm{X}, y) = \frac{1}{n(n-1)}
#'   \sum_{i \ne j} \| \bm{X}_i - \bm{X}_j \| \left(1 - \sum_{r=1}^{R}
#'   I(Y_i = r) I(Y_j = r) / \tilde{p}_r\right),}
#'   where \eqn{\tilde{p}_r = (n_r-1) / (n-1) = (n-1)^{-1}(\sum_{i=1}^{n} I(Y_i
#'   = r) - 1)}. Note that the test statistic of the semi-distance independence
#'   test is \deqn{T_n = n \cdot \widetilde{\text{SDcov}}_n(\bm{X}, y).}
#'
#' @returns The value of the corresponding sample statistic.
#'
#'   If the argument `return_mat` of `sdcov()` is set as `TRUE`, a list with
#'   elements
#'   * `sdcov`: the semi-distance covariance statistic;
#'   * `mat_x, mat_y`: the matrices of the distances of X and the divergences
#'   of y, respectively;
#'
#'   will be returned.
#'
#' @seealso
#'   * [sd_test()] for implementing independence test via semi-distance
#' covariance;
#'   * [sd_sis()] for implementing groupwise feature screening via
#' semi-distance correlation.
#'
#' @examples
#' X <- mtcars[, c("mpg", "disp", "drat", "wt")]
#' y <- factor(mtcars[, "am"])
#' print(sdcov(X, y))
#' print(sdcor(X, y))
#'
#' # Man-made independent data -------------------------------------------------
#' n <- 30; R <- 5; p <- 3; prob <- rep(1/R, R)
#' X <- matrix(rnorm(n*p), n, p)
#' y <- factor(sample(1:R, size = n, replace = TRUE, prob = prob), levels = 1:R)
#' print(sdcov(X, y))
#' print(sdcor(X, y))
#'
#' # Man-made functionally dependent data --------------------------------------
#' n <- 30; R <- 3; p <- 3
#' X <- matrix(0, n, p)
#' X[1:10, 1] <- 1; X[11:20, 2] <- 1; X[21:30, 3] <- 1
#' y <- factor(rep(1:3, each = 10))
#' print(sdcov(X, y))
#' print(sdcor(X, y))
#'
#' @export
sdcov <- function(X, y, type = "V", return_mat = FALSE) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  } else if (is.numeric(X)) {
    X <- as.matrix(X, ncol = 1) # Univariate case
  } else if (!is.matrix(X)) {
    stop("`X` should be a matrix or a numeric vector")
  }
  
  if (!is.factor(y)) {
    stop("`y` should be a factor")
  } else {
    Y <- switch_cat_repr(y) # Convert to indicator matrix representation
  }
  
  n <- nrow(Y)
  n_r <- colSums(Y)
  
  if (type == "U") {
    Y <- Y[, which(n_r > 1), drop = FALSE] # Filter: `r` such that n_r > 1
    prob <- (colSums(Y) - 1) / (n - 1)
    
    denom <- n * (n - 1)
    mat_x <- as.matrix(dist(X))
    mat_y <- 1 - Y %*% diag(1 / prob) %*% t(Y)
  } else if (type == "V") {
    prob <- colSums(Y) / n
    
    denom <- n^2
    mat_x <- as.matrix(dist(X))
    # mat_x <- energy::Dcenter(X) # Equivalent!
    mat_y <- 1 - Y %*% diag(1 / prob) %*% t(Y)
  } else {
    stop('`type` should be either "U" or "V"')
  }
  
  sdcov <- sum(mat_x * mat_y) / denom
  
  if (return_mat) {
    return(
      list(
        sdcov = sdcov,
        mat_x = mat_x,
        mat_y = mat_y
      )
    )
  } else {
    return(sdcov)
  }
}

#' @rdname sdcov
#' @export
sdcor <- function(X, y) {
  R <- length(levels(y))
  
  dvar_X <- energy::dcov(X, X)
  sdcor <- sdcov(X, y) / (dvar_X * sqrt(R - 1))
  
  sdcor
}

#' Semi-distance independence test
#'
#' @description Implement the semi-distance independence test via permutation
#'   test, or via the asymptotic approximation when the dimensionality of
#'   continuous variables \eqn{p} is high.
#'
#' @inheritParams sdcov
#' @param test_type Type of the test:
#'   * `"perm"` (the default): Implement the test via permutation test;
#'   * `"asym"`: Implement the test via the asymptotic approximation when the
#'   dimension of continuous variables \eqn{p} is high.
#'
#'   See the Reference for details.
#' @param num_perm The number of replications in permutation test. Defaults to
#'   10000. See Details and Reference.
#'
#' @details The semi-distance independence test statistic is \deqn{T_n = n \cdot
#'   \widetilde{\text{SDcov}}_n(X, y),} where the
#'   \eqn{\widetilde{\text{SDcov}}_n(X, y)} can be computed by `sdcov(X, y, type
#'   = "U")`.
#'
#'   For the permutation test (`test_type = "perm"`), totally \eqn{K}
#'   replications of permutation will be conducted, and the argument `num_perm`
#'   specifies the \eqn{K} here. The p-value of permutation test is computed by
#'   \deqn{\text{p-value} = (\sum_{k=1}^K I(T^{\ast (k)}_{n} \ge T_{n}) + 1) /
#'   (K + 1),} where \eqn{T_{n}} is the semi-distance test statistic and
#'   \eqn{T^{\ast (k)}_{n}} is the test statistic with \eqn{k}-th permutation
#'   sample.
#'
#'   When the dimension of the continuous variables is high, the asymptotic
#'   approximation approach can be applied (`test_type = "asym"`), which is
#'   computationally faster since no permutation is needed.
#'
#' @returns A list with class `"indtest"` containing the following components
#' * `method`: name of the test;
#' * `name_data`: names of the `X` and `y`;
#' * `n`: sample size of the data;
#' * `test_type`: type of the test;
#' * `num_perm`: number of replications in permutation test, if
#'   `test_type = "perm"`;
#' * `stat`: test statistic;
#' * `pvalue`: computed p-value.
#'
#' @seealso [sdcov()] for computing the statistic of semi-distance covariance.
#'
#' @examples
#' X <- mtcars[, c("mpg", "disp", "drat", "wt")]
#' y <- factor(mtcars[, "am"])
#' test <- sd_test(X, y)
#' print(test)
#'
#' # Man-made independent data -------------------------------------------------
#' n <- 30; R <- 5; p <- 3; prob <- rep(1/R, R)
#' X <- matrix(rnorm(n*p), n, p)
#' y <- factor(sample(1:R, size = n, replace = TRUE, prob = prob), levels = 1:R)
#' test <- sd_test(X, y)
#' print(test)
#'
#' # Man-made functionally dependent data --------------------------------------
#' n <- 30; R <- 3; p <- 3
#' X <- matrix(0, n, p)
#' X[1:10, 1] <- 1; X[11:20, 2] <- 1; X[21:30, 3] <- 1
#' y <- factor(rep(1:3, each = 10))
#' test <- sd_test(X, y)
#' print(test)
#'
#' #' Man-made high-dimensionally independent data -----------------------------
#' n <- 30; R <- 3; p <- 100
#' X <- matrix(rnorm(n*p), n, p)
#' y <- factor(rep(1:3, each = 10))
#' test <- sd_test(X, y)
#' print(test)
#'
#' test <- sd_test(X, y, test_type = "asym")
#' print(test)
#'
#' # Man-made high-dimensionally dependent data --------------------------------
#' n <- 30; R <- 3; p <- 100
#' X <- matrix(0, n, p)
#' X[1:10, 1] <- 1; X[11:20, 2] <- 1; X[21:30, 3] <- 1
#' y <- factor(rep(1:3, each = 10))
#' test <- sd_test(X, y)
#' print(test)
#'
#' test <- sd_test(X, y, test_type = "asym")
#' print(test)
#'
#' @export
sd_test <- function(X, y, test_type = "perm", num_perm = 10000) {
  obj <- sdcov(X, y, type = "U", return_mat = TRUE)
  
  n <- length(y)
  name_data <- paste(
    deparse(substitute(X)), "and", deparse(substitute(y))
  )
  
  sdcov_n <- obj$sdcov
  mat_x <- obj$mat_x
  mat_y <- obj$mat_y
  
  stat <- n * sdcov_n
  
  if (test_type == "perm") {
    name_method <- paste0(
      "Semi-Distance Independence Test (Permutation Test with K = ",
      num_perm, ")"
    )
    
    num_rej <- 0
    for (perm in 1:num_perm) {
      ind_perm <- sample(1:n)
      mat_y_perm <- mat_y[ind_perm, ind_perm]
      sdcov_star <- sum(mat_x * mat_y_perm)/ (n * (n - 1))
      if (sdcov_star >= sdcov_n) {
        num_rej <- num_rej + 1
      }
    }
    pvalue <- (num_rej + 1) / (num_perm + 1)
  } else if (test_type == "asym") {
    name_method <- "Semi-distance Independence (Asymptotic) Test"
    
    tr_obj <- tr_estimate(X)
    tr_S <- tr_obj$tr_S
    tr_S2 <- tr_obj$tr_S2
    
    n_r <- colSums(switch_cat_repr(y))
    d <- tr_S2 / tr_S * (sum(n_r / (n_r - 1)) - n / (n - 1))
    pvalue <- pnorm(stat, sd = sqrt(d), lower.tail = FALSE) # one-sided
  } else {
    stop('`test_type` should be either "perm" or "asym"!')
  }
  
  indtest <- list(
    method = name_method,
    name_data = name_data,
    n = n,
    test_type = test_type,
    stat = n * sdcov_n,
    pvalue = pvalue
  )
  if (test_type == "perm") {
    indtest$num_perm = num_perm
  }
  class(indtest) <- "indtest"
  
  indtest
}

#' Feature screening via semi-distance correlation
#'
#' @description Implement the (grouped) feature screening for the classification
#'   problem via semi-distance correlation.
#'
#' @param X Data of multivariate covariates, which should be an
#'   \eqn{n}-by-\eqn{p} matrix.
#' @param y Data of categorical response, which should be a factor of length
#'   \eqn{n}.
#' @param group_info A list specifying the group information, with elements
#'   being sets of indicies of covariates in a same group. For example,
#'   `list(c(1, 2, 3), c(4, 5))` specifies that covariates 1, 2, 3 are in a
#'   group and covariates 4, 5 are in another group.
#'
#'   Defaults to `NULL`. If `NULL`, then it will be set as `list(1, 2, ..., p)`,
#'   that is, treat each single covariate as a group.
#'
#'   If `X` has colnames, then the colnames can be used to specified the
#'   `group_info`. For example, `list(c("a", "b"), c("c", "d"))`.
#'
#'   The names of the list can help recoginize the group. For example,
#'   `list(grp_ab = c("a", "b"), grp_cd = c("c", "d"))`. If names of the list
#'   are not specified, `c("Grp 1", "Grp 2", ..., "Grp J")` will be applied.
#'
#' @param d An integer specifying *at least* how many (single) features should
#'   be kept after screening. For example, if `group_info = list(c(1, 2), c(3,
#'   4))` and `d = 3`, then all features 1, 2, 3, 4 must be selected since it
#'   should guarantee at least 3 features are kept.
#'
#'   Defaults to `NULL`. If `NULL`, then it will be set as \eqn{[n / log(n)]},
#'   where \eqn{[x]} denotes the integer part of x.
#'
#' @param parallel A boolean indicating whether to calculate parallelly via
#'   `furrr::future_map`. Defaults to `FALSE`.
#'
#' @returns A list of the objects about the implemented feature screening:
#' * `group_info`: group information;
#' * `measurement`: sample semi-distance correlations calculated for the groups
#'   specified in `group_info`;
#' * `selected`: indicies/names of (single) covariates that are selected after
#'   feature screening;
#' * `ordering`: order of the calculated measurements of the groups specified in
#'   `group_info`. The first one is the largest, and the last is the smallest.
#'
#' @seealso [sdcor()] for calculating the sample semi-distance correlation.
#'
#' @examples
#' X <- mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")]
#' y <- factor(mtcars[, "am"])
#'
#' sd_sis(X, y, d = 4)
#'
#' # Suppose we have prior information for the group structure as
#' # ("mpg", "drat"), ("disp", "hp") and ("wt", "qsec")
#' group_info <- list(
#'   mpg_drat = c("mpg", "drat"),
#'   disp_hp = c("disp", "hp"),
#'   wt_qsec = c("wt", "qsec")
#' )
#' sd_sis(X, y, group_info, d = 4)
#'
#' @export
sd_sis <- function(X, y, group_info = NULL, d = NULL, parallel = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(d)) {
    d <- floor(n / log(n))
  }
  
  if (is.null(colnames(X))) {
    names_X <- 1:p
    colnames(X) <- names_X
  } else {
    names_X <- colnames(X)
  }
  
  if (is.null(group_info)) {
    group_info <- split(names_X, 1:p)
    names(group_info) <- paste("Grp", names_X)
  } else if (is.null(names(group_info))) {
    names(group_info) <- paste("Grp", 1:length(group_info))
  }
  
  J <- length(group_info)
  
  if (parallel) {
    omega <- furrr::future_map_dbl(group_info, ~ sdcor(X[, .x, drop = FALSE], y))
  } else {
    omega <- map_dbl(group_info, ~ sdcor(X[, .x], y))
  }
  
  ordering <- order(omega, decreasing = TRUE)
  ordering_names <- names(group_info)[ordering]
  
  if (length(group_info[[ordering[1]]]) >= d) {
    D_hat <- group_info[[ordering[1]]]
  } else {
    D_hat <- c()
    for (j in 1:J) {
      D_hat <- unique(append(D_hat, group_info[[ordering[j]]]))
      if (length(D_hat) >= d) {
        break
      } 
    }
  }
  
  list(
    group_info = group_info,
    measurement = omega,
    selected = D_hat,
    ordering = ordering_names
  )
}
