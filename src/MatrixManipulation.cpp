#include <RcppArmadillo.h>
#include <cmath>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Function that inverts a non-singular matrix
//' @param a matrix to be inverted
//' @returns the inverse matrix
// [[Rcpp::export]]
mat Transp(const mat &x)
{
  return trans(x);
}

//' Function that transposes
//' @param a matrix to be transposed
//' @returns the transposed matrix
// [[Rcpp::export]]
mat Inverse(const mat &x)
{
  return inv(x);
}
