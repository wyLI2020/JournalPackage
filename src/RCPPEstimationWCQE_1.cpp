#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double rhotau(double tau, double u) {
  // rhotau() - \rho_\tau(u) = (\tau - I(u<0)) u
  double result = 0.0;
  if(u < 0) {
    result = (tau - 1.0) * u;
  } else {
    result = tau * u;
  }
  return result;
}
double Psitau(double tau, double u) {
  // Psitau() - \Psi_\tau(u) = \tau - I(u<0)
  double result = 0.0;
  if(u < 0) {
    result = tau - 1.0;
  } else {
    result = tau;
  }
  return result;
}
// [[Rcpp::export]]
double wqrLF(arma::vec varphi, int N, arma::vec y, arma::mat xwith1, double tau, arma::vec weights) {
  // wqrLF() - loss function of quantile
  arma::vec eleV(N);
  for (int i = 0; i < N; i++) {
    double u = y(i) - dot(xwith1.row(i), varphi);
    eleV(i) = weights(i) * rhotau(tau, u);
  }
  double result = sum(eleV);
  return result;
}
// [[Rcpp::export]]
arma::vec gr_wqrLF(arma::vec varphi, int N, arma::vec y, arma::mat xwith1, double tau, arma::vec weights) {
  // gr_wqrLF() - gradient of loss function of quantile
  arma::mat eleM(varphi.n_elem,N);
  arma::mat xwith1T = xwith1.t(); // x^{\prime}
  for (int i = 0; i < N; i++) {
    double u = y(i) - dot(xwith1T.col(i), varphi);
    eleM.col(i) = weights(i) * Psitau(tau, u) * xwith1T.col(i);
  }
  arma::vec result = sum(eleM, 1);
  return result;
}
