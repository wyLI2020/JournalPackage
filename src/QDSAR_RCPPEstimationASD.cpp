#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double IQR_rcpp2(arma::vec x){
  // IQR_rcpp2 - IQR() in R, interquartile range of the x values
  Function f("IQR");
  NumericVector result = f(Named("x") = x);
  return result(0);
}
double dnorm_rcpp2(double x){
  // dnorm_rcpp2 - dnorm() in R, Gaussian kernel function
  Function f("dnorm");
  NumericVector result = f(Named("x") = x);
  return result(0);
}
arma::vec quantile_rcpp2(arma::vec x, arma::vec tauseq){
  // quantile_rcpp2 - quantile(), sample quantile
  Function f("quantile");
  NumericVector resultr = f(Named("x") = x, Named("probs") = tauseq);
  arma::vec result(resultr.begin(), resultr.size(), false);
  return result;
}
arma::vec f_eta_est(int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec beta_est) {
  // f_eta_est - residuals, \hat{\eta}_{i} = \hat{\vartheta}_{i} / (x_{ai}^\prime \hat{\beta}(\hat{\pi}_{opt,K})), or \tilde{\eta}_{i} = \hat{\vartheta}_{i} / (x_{ai}^\prime \tilde{\beta}_{c})
  arma::vec eta_est(N);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    eta_est(i) = varthetahat(i) / dot(xawith1i, beta_est);
  }
  return eta_est;
}
arma::vec f_fetaQetatau_est(arma::vec tauseq, int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec beta_est, arma::vec eta_est) { // Gaussian kernel
  // f_fetaQetatau_est - estimation of the density of \eta at Q_{\eta}(\tau), \hat{f}_{\eta}(\check{Q}_{\eta}(\tau)) or \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau))
  int taunum = tauseq.n_elem;
  arma::vec QetataucheckV = quantile_rcpp2(eta_est, tauseq); // \check{Q}_{\eta}(\tau), that is \tau-th sample quantile of residuals \hat{\eta}_{i}'s or \tilde{\eta}_{i}'s
  // bandwidth h
  double s = stddev(eta_est); // standard deviation of the residuals
  double Rhat = IQR_rcpp2(eta_est); // interquartile of the residuals
  arma::vec vecFORmin(2); vecFORmin(0) = s; vecFORmin(1) = Rhat/1.34;
  double h = 0.9 * pow((N*1.0), -0.2) * min(vecFORmin);
  // \hat{f}_{\eta}(\check{Q}_{\eta}(\tau)) = (Nh)^{-1} \sum K((\check{Q}_{\eta}(\tau) - \hat{\eta}_{i})/h), or \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau)) = (Nh)^{-1} \sum K((\check{Q}_{\eta}(\tau) - \tilde{\eta}_{i})/h)
  arma::vec fetaQetatau_estV(taunum); // the j-th element is \hat{f}_{\eta}(\check{Q}_{\eta}(\tau_{j})) or \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau_{j}))
  for (int j = 0; j < taunum; j++) {
    arma::vec Ki(N); // the i-th element is K((\check{Q}_{\eta}(\tau_{j}) - \hat{\eta}_{i})/h) or K((\check{Q}_{\eta}(\tau_{j}) - \tilde{\eta}_{i})/h)
    for (int i = 0; i < N; i++) {
      Ki(i) = dnorm_rcpp2((QetataucheckV(j) - eta_est(i)) / h);
    }
    double feta_est = sum(Ki) / (N*h); // \hat{f}_{\eta}(\check{Q}_{\eta}(\tau_{j})) or \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau_{j}))
    fetaQetatau_estV(j) = feta_est;
  }
  return fetaQetatau_estV;
}
// [[Rcpp::export]]
arma::mat f_D_est(int N, arma::mat Xawith1, arma::vec beta_est) {
  // f_D_est - estimation of D_{N}, \hat{D}_{N} or \tilde{D}_{N}
  int dimvarphi = Xawith1.n_cols;
  arma::cube DC_est(dimvarphi,dimvarphi,N); DC_est.fill(0.0);
  arma::mat xawith1iM(1,dimvarphi);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    xawith1iM.row(0) = Xawith1.row(i);
    DC_est.slice(i) = (xawith1iM.t()*xawith1iM) / (dot(xawith1i, beta_est)*dot(xawith1i, beta_est));
  }
  arma::mat D_est = mean(DC_est, 2);
  return D_est;
}
// [[Rcpp::export]]
arma::mat f_ASD_varphi(arma::vec tauseq, int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec betactilde, arma::mat varphitauhatM, arma::mat Dtilde) {
  // f_ASD_varphi - estimation for ASD of \hat{\varphi}(\tau), with \tau in tauseq
  int dimvarphi = Xawith1.n_cols;
  int taunum = tauseq.n_elem;
  arma::vec etatilde = f_eta_est(N, Xawith1, varthetahat, betactilde); // residual, \tilde{\eta}
  arma::vec fetaQetatautildeV = f_fetaQetatau_est(tauseq, N, Xawith1, varthetahat, betactilde, etatilde); // \tilde{f}_{\eta}(\check{Q}_{\eta}(\tau))
  // ASD
  arma::mat ASDM(dimvarphi,taunum); // the j-th column is the estimation of ASD of \hat{\varphi}(\tau_{j})
  for (int j = 0; j < taunum; j++) {
    double tau = tauseq(j);
    double Xihattau = fetaQetatautildeV(j); // \tilde{\Xi}(\tau_{j})
    arma::mat CovM = tau*(1-tau) * 1.0/(Xihattau*Xihattau) * Dtilde.i() / (N*1.0); // covariance matrix of \hat{\varphi}(\tau_{j})
    arma::vec ASDvarphitauhat(dimvarphi);
    for (int l = 0; l < dimvarphi; l++) {
      ASDvarphitauhat(l) = sqrt(CovM(l,l));
    }
    ASDM.col(j) = ASDvarphitauhat;
  }
  return ASDM;
}
// [[Rcpp::export]]
arma::vec f_ASD_beta(arma::vec tauK, int N, int T, arma::mat Xawith1, arma::vec varthetahat, arma::vec betahat, arma::mat varphitauKhatM, arma::mat Dhat) {
  // f_ASD_beta - estimation for ASD of \hat{beta}
  int dimbeta = Xawith1.n_cols;
  int K = tauK.n_elem;
  arma::vec etahat = f_eta_est(N, Xawith1, varthetahat, betahat); // residual, \hat{\eta}
  arma::vec fetaQetatauKhatV = f_fetaQetatau_est(tauK, N, Xawith1, varthetahat, betahat, etahat); // \hat{f}_{\eta}(\check{Q}_{\eta}(\tau))
  // \hat{H}
  arma::rowvec QetatauKhat = varphitauKhatM.row(0);
  arma::vec qhat(K); for (int k = 0; k < K; k++) {qhat(k) = QetatauKhat(k);} // \hat{q}
  arma::mat Hhat(K,K); Hhat.fill(0.0);
  for (int k1 = 0; k1 < K; k1++) {
    for (int k2 = 0; k2 < K; k2++) {
      double tauk1 = tauK(k1); double tauk2 = tauK(k2);
      double Xihattauk1 = fetaQetatauKhatV(k1);
      double Xihattauk2 = fetaQetatauKhatV(k2);
      arma::vec tauk1tauk2(2); tauk1tauk2(0) = tauk1; tauk1tauk2(1) = tauk2;
      double Stauk1tauk2 = min(tauk1tauk2) - tauk1*tauk2;
      Hhat(k1,k2) = (1.0/Xihattauk1) * Stauk1tauk2 * (1.0/Xihattauk2);
    }
  }
  // covariance matrix of \hat{\beta}
  arma::mat CovM = 1.0/dot(qhat, Hhat.i()*qhat) * Dhat.i() / (N*1.0);
  // ASD
  arma::vec ASDbetahat(dimbeta);
  for (int l = 0; l < dimbeta; l++) {
    ASDbetahat(l) = sqrt(CovM(l,l));
  }
  return ASDbetahat;
}