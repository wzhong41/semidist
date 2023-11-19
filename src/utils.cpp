#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Estimate the trace of the covariance matrix and its square
//' 
//' @description 
//' For a design matrix \eqn{\mathbf{X}}, estimate the trace of its covariance matrix \eqn{\Sigma = \mathrm{cov}(\mathbf{X})},
//' and the square of covariance matrix \eqn{\Sigma^2}.
//'
//' @param X The design matrix.
//' @return A list with elements:
//' * `tr_S`: estimate for trace of \eqn{\Sigma};
//' * `tr_S2`: estimate for trace of \eqn{\Sigma^2}.
// [[Rcpp::export]]
List tr_estimate(const arma::mat& X) {
  int n = X.n_rows;
  
  arma::mat XXt = X * X.t();
  double Y1 = arma::accu(arma::diagvec(XXt)) / n;
  double Y3 = (arma::accu(XXt) - Y1 * n) / (n * (n - 1));
  double Y2 = (arma::accu(arma::pow(XXt, 2)) - arma::accu(arma::pow(arma::diagvec(XXt), 2))) / (n * (n - 1));
  
  arma::mat XXt2 = XXt * XXt;
  arma::vec diag_XXt = arma::vectorise(arma::repelem(arma::diagvec(XXt), n - 1, 1));
  arma::vec offdiag_XXt = arma::vectorise(XXt(arma::find(arma::eye<arma::umat>(n, n) == 0)));
  double Y4 = (arma::accu(XXt2) - arma::accu(arma::diagvec(XXt2)) - 2 * arma::accu(diag_XXt % offdiag_XXt)) /
    (n * (n - 1) * (n - 2));
  double Y5 = ((n * (n - 1) * Y3) * (n * (n - 1) * Y3) - 2 * n * (n - 1) * Y2 - 4 * n * (n - 1) * (n - 2) * Y4) /
    (n * (n - 1) * (n - 2) * (n - 3));
  
  double Tn1 = Y1 - Y3;
  double Tn2 = Y2 - 2 * Y4 + Y5;
  
  return List::create(
    Named("tr_S") = Tn1,
    Named("tr_S2") = Tn2
  );
}