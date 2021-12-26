#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat f_H_est(int K, arma::vec tauK, arma::vec XitauKestV) {
  // f_H_est - estimation of matrix H
  arma::mat Hest(K,K); Hest.fill(0.0);
  for (int k1 = 0; k1 < K; k1++) {
    for (int k2 = 0; k2 < K; k2++) {
      double tauk1 = tauK(k1); double tauk2 = tauK(k2);
      double Xiest_tauk1 = XitauKestV(k1);
      double Xiest_tauk2 = XitauKestV(k2);
      arma::vec tauk1tauk2(2); tauk1tauk2(0) = tauk1; tauk1tauk2(1) = tauk2;
      double Stauk1tauk2 = min(tauk1tauk2) - tauk1*tauk2;
      Hest(k1,k2) = (1.0/Xiest_tauk1) * Stauk1tauk2 * (1.0/Xiest_tauk2);
    }
  }
  return Hest;
}
// [[Rcpp::export]]
arma::vec f_WQAE(arma::mat varphitauKhatM, arma::vec qest, arma::mat Hest) {
  // f_WQAE - \hat{\beta}(\hat{\pi}_{opt,K}) = \sum_{k=1}^{K} \pi_{k} \hat{\varphi}(\tau_{k})
  arma::vec pioptKhat = Hest.i()*qest / dot(qest, Hest.i()*qest); // \hat{\pi}_{opt,K}
  // \hat{\beta}(\hat{\pi}_{opt,K})
  arma::vec betahat = varphitauKhatM * pioptKhat;
  return betahat;
}
// [[Rcpp::export]]
arma::vec f_ASD_WQAE(int N, arma::vec qest, arma::mat Hest, arma::mat Dhat) {
  // f_ASD_WQAE - estimation for ASD of WQAE
  arma::mat CovM = 1.0/dot(qest, Hest.i()*qest) * Dhat.i() / (N*1.0); // covariance matrix of \hat{\beta}
  // ASD
  arma::vec ASDsquarebetahat = CovM.diag(0);
  arma::vec ASDbetahat = sqrt(ASDsquarebetahat);
  return ASDbetahat;
}
