#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec f_betac_tilde(arma::mat varphitautildeM) {
  // f_betac_tilde - \tilde{\beta}_{c} 
  arma::mat varphitautildeM_abs = abs(varphitautildeM);
  // \tilde{Q}_{\eta}(\tau_{k}), k = 1,...,K
  arma::rowvec Qetatautilde = varphitautildeM.row(0);
  arma::rowvec Qetatautilde_abs = abs(Qetatautilde);
  // \tilde{\beta}_{c} 
  arma::vec sum_varphitautildeM_abs = sum(varphitautildeM_abs, 1);
  double sum_Qetatautilde_abs = sum(Qetatautilde_abs);
  arma::vec betactilde = sum_varphitautildeM_abs / sum_Qetatautilde_abs;
  return betactilde;
}
// [[Rcpp::export]]
arma::vec f_weights_rq(int N, arma::mat Xawith1, arma::vec betactilde) {
  // weights in WCQE: 
  arma::vec weights_rq(N);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    weights_rq(i) = 1 / dot(xawith1i, betactilde);
  }
  return weights_rq;
}
double rhotau(double tau, double u) {
  // rhotau() - \rho_\tau(u)
  double result = 0.0;
  if(u < 0) {
    result = (tau - 1) * u;
  } else {
    result = tau * u;
  }
  return result;
}
// [[Rcpp::export]]
double wqrLF(arma::vec varphi, int N, arma::vec y, arma::mat x, double tau, arma::vec weights) {
  // wqrLF() - loss function of quantile
  arma::vec elev(N);
  for (int i = 0; i < N; i++) {
    double u = y(i) - dot(x.row(i), varphi);
    elev(i) = weights(i) * rhotau(tau, u);
  }
  double result = sum(elev);
  return result;
}
  


