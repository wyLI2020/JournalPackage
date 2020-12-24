#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double IQR_rcpp(arma::vec x){
  // IQR - interquartile range of the x values
  Function f("IQR");
  NumericVector result = f(Named("x") = x);
  return result(0);
}
double dnorm_rcpp(double x){
  // dnorm - Gaussian kernel function
  Function f("dnorm");
  NumericVector result = f(Named("x") = x);
  return result(0);
}
arma::vec quantile_rcpp(arma::vec x, arma::vec tauseq){
  // quantile - sample quantile
  Function f("quantile");
  NumericVector resultr = f(Named("x") = x, Named("probs") = tauseq);
  arma::vec result(resultr.begin(), resultr.size(), false);
  return result;
}
arma::vec f_eta_tilde(int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec betactilde) {
  // f_eta_tilde - residuals, \tilde{\eta}_{i} = \hat{\vartheta}_{i} / (x_{ai}^\prime \tilde{\beta}_{c})
  arma::vec etatilde(N);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    etatilde(i) = varthetahat(i) / dot(xawith1i, betactilde);
  }
  return etatilde;
}
arma::vec f_fetaQetatau_tilde(arma::vec tauseq, int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec betactilde, arma::vec etatilde) { // Gaussian kernel
  // f_fetaQetatau_tilde - estimation of the density of \eta at Q_{\eta}(\tau), \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau))
  int taunum = tauseq.n_elem;
  arma::vec QetataucheckV = quantile_rcpp(etatilde, tauseq); // \check{Q}_{\eta}(\tau), that is \tau-th sample quantile of residuals \tilde{\eta}_{i}'s
  // bandwidth h
  double s = stddev(etatilde); // standard deviation of the residuals
  double Rhat = IQR_rcpp(etatilde); // interquartile of the residuals
  arma::vec vecFORmin(2); vecFORmin(0) = s; vecFORmin(1) = Rhat/1.34;
  double h = 0.9 * pow((N*1.0), -0.2) * min(vecFORmin);
  // \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau)) = (Nh)^{-1} \sum K((\check{Q}_{\eta}(\tau) - \tilde{\eta}_{i})/h)
  arma::vec fetaQetatautildeV(taunum); // the j-th element is \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau_{j}))
  for (int j = 0; j < taunum; j++) {
    arma::vec Ki(N); // the i-th element is K((\check{Q}_{\eta}(\tau_{j}) - \tilde{\eta}_{i})/h)
    for (int i = 0; i < N; i++) {
      Ki(i) = dnorm_rcpp((QetataucheckV(j) - etatilde(i)) / h);
    }
    double fetatilde = sum(Ki) / (N*h); // \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau_{j}))
    fetaQetatautildeV(j) = fetatilde;
  }
  return fetaQetatautildeV;
}
// [[Rcpp::export]]
arma::vec f_beta_hat(arma::vec tauK, int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec betactilde, arma::mat varphitauKhatM) {
  // f_beta_hat - \hat{\beta}(\hat{\pi}_{opt,K})
  int K = tauK.n_elem;
  arma::vec etatilde = f_eta_tilde(N, Xawith1, varthetahat, betactilde); // residual, \tilde{\eta}
  arma::vec fetaQetatautildeV = f_fetaQetatau_tilde(tauK, N, Xawith1, varthetahat, betactilde, etatilde); // \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau))
  // \hat{\pi}_{opt,K}
  arma::rowvec QetatauKhat = varphitauKhatM.row(0);
  arma::vec qhat(K); for (int k = 0; k < K; k++) {qhat(k) = QetatauKhat(k);} // \hat{q}
  arma::mat Htilde(K,K); Htilde.fill(0.0); // \tilde{H}
  for (int k1 = 0; k1 < K; k1++) {
    for (int k2 = 0; k2 < K; k2++) {
      double tauk1 = tauK(k1); double tauk2 = tauK(k2);
      double Xitildetauk1 = fetaQetatautildeV(k1);
      double Xitildetauk2 = fetaQetatautildeV(k2);
      arma::vec tauk1tauk2(2); tauk1tauk2(0) = tauk1; tauk1tauk2(1) = tauk2;
      double Stauk1tauk2 = min(tauk1tauk2) - tauk1*tauk2;
      Htilde(k1,k2) = (1.0/Xitildetauk1) * Stauk1tauk2 * (1.0/Xitildetauk2);
    }
  }
  arma::vec pioptKhat = Htilde.i()*qhat / dot(qhat, Htilde.i()*qhat);
  // \hat{\beta}(\hat{\pi}_{opt,K})
  arma::mat piKvarphiK(varphitauKhatM.n_rows, varphitauKhatM.n_cols); piKvarphiK.fill(0.0); // the k-th column is \hat{\pi}_{k,opt,K} \hat{\varphi}(\tau_{k})
  for (int k = 0; k < K; k++) {
    double pik = pioptKhat(k);
    arma::vec varphihattauk = varphitauKhatM.col(k);
    piKvarphiK.col(k) = pik * varphihattauk;
  }
  arma::vec betahat = sum(piKvarphiK, 1);
  return betahat;
}
