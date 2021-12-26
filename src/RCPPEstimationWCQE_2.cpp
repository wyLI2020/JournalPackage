#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat f_WCQE_tausM(int p, int N, arma::vec y, arma::mat x, arma::mat xwith1, int tausNUM, arma::vec tausV, arma::vec weightsV) {
  // WCQE
  Environment OPT = Environment::namespace_env("QuantileDSAR");
  Function f = OPT["Rf_WCQE"];
  // Function f("Rf_WCQE");
  arma::mat WCQEM(xwith1.n_cols,tausNUM); WCQEM.fill(0.0);
  for (int k = 0; k < tausNUM; k++) {
    double tauk = tausV(k);
    NumericVector result_rcpp = f(p, N, y, x, xwith1, tauk, weightsV);
    arma::vec result_arma(result_rcpp.begin(), result_rcpp.size(), false);
    WCQEM.col(k) = result_arma;
  }
  return WCQEM;
}
// [[Rcpp::export]]
arma::vec f_betac_tilde(arma::mat varphitautildeM) {
  // f_betac_tilde - \tilde{\beta}_{c} = \frac {\sum_{k}|\tilde{\varphi}(\tau_{k})|} {\sum_{k}|\tilde{Q}_{\eta}(\tau_{k})|}
  arma::mat varphitautildeM_abs = abs(varphitautildeM);
  // \tilde{Q}_{\eta}(\tau_{k}), k = 1,...,K
  arma::rowvec Qetatautilde = varphitautildeM.row(0);
  arma::rowvec Qetatautilde_abs = abs(Qetatautilde);
  // \tilde{\beta}_{c} = \frac{\sum_{k}|\tilde{\varphi}(\tau_{k})|}{\sum_{k}|\tilde{Q}_{\eta}(\tau_{k})|}
  arma::vec sum_varphitautildeM_abs = sum(varphitautildeM_abs, 1);
  double sum_Qetatautilde_abs = sum(Qetatautilde_abs);
  arma::vec betactilde = sum_varphitautildeM_abs / sum_Qetatautilde_abs;
  return betactilde;
}
// [[Rcpp::export]]
arma::vec f_weights_rq(int N, arma::mat Xawith1, arma::vec betactilde) {
  // weights in WCQE: (x_{a1}^\prime \tilde{\beta}_{c}, ..., x_{aN}^\prime \tilde{\beta}_{c})^\prime
  arma::vec weights_rq(N);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    weights_rq(i) = 1 / dot(xawith1i, betactilde);
  }
  return weights_rq;
}
// [[Rcpp::export]]
arma::vec f_eta_est(int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec betaest) {
  // f_eta_est - residuals, \eta_{i} = \vartheta_{i} / (x_{ai}^\prime \beta)
  arma::vec etaest(N);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    etaest(i) = varthetahat(i) / dot(xawith1i, betaest);
  }
  return etaest;
}
double IQR_rcpp(arma::vec x){
  // IQR_rcpp - IQR() in R, interquartile range of the x values
  Function f("IQR");
  NumericVector result = f(Named("x") = x);
  return result(0);
}
double dnorm_rcpp(double x){
  // dnorm_rcpp - dnorm() in R, Gaussian kernel function
  Function f("dnorm");
  NumericVector result = f(Named("x") = x);
  return result(0);
}
// [[Rcpp::export]]
arma::vec f_feta_tildeV(arma::vec xV, int N, arma::mat Xawith1, arma::vec varthetahat, arma::vec betactilde) { // Gaussian kernel
  // f_feta_tildeV - estimations of the density of \eta at points x's in xV, \tilde{f}_{\eta}(x) = (Nh)^{-1} \sum_{i=1}^{N} K((x - \tilde{\eta}_{i})/h) with \tilde{\eta}_{i} = \vartheta_{i} / (x_{ai}^{\prime} \tilde{\beta}_{c})
  arma::vec etatilde = f_eta_est(N, Xawith1, varthetahat, betactilde); // \tilde{\eta}
  // bandwidth h
  double s = stddev(etatilde); // standard deviation of the residuals
  double Rhat = IQR_rcpp(etatilde); // interquartile of the residuals
  arma::vec vecFORmin(2); vecFORmin(0) = s; vecFORmin(1) = Rhat/1.34;
  double h = 0.9 * pow((N*1.0), -0.2) * min(vecFORmin);
  // \hat{f}_{\eta}(x_{j}) with x_{j}'s in xV
  int eleNUM = xV.n_elem;
  arma::vec fetatildeV(eleNUM); // the j-th element is \tilde{f}_{\eta}(x_{j})
  for (int j = 0; j < eleNUM; j++) {
    double xj = xV(j);
    arma::vec KV(N); // the i-th element is K((x_{j} - \tilde{\eta}_{i})/h)
    for (int i = 0; i < N; i++) {
      KV(i) = dnorm_rcpp((xj - etatilde(i)) / h);
    }
    double fetatilde_xj = sum(KV) / (N*h); // \tilde{f}_{\eta}(x)
    fetatildeV(j) = fetatilde_xj;
  }
  return fetatildeV;
}
// [[Rcpp::export]]
arma::mat f_D_est(int N, arma::mat Xawith1, arma::vec betaest) {
  // f_D_est - estimation of D_{N}, \hat{D}_{N} or \tilde{D}_{N}
  int dimvarphi = Xawith1.n_cols;
  arma::cube DC_est(dimvarphi,dimvarphi,N); DC_est.fill(0.0);
  arma::mat xawith1iM(1,dimvarphi);
  for (int i = 0; i < N; i++) {
    arma::rowvec xawith1i = Xawith1.row(i);
    xawith1iM.row(0) = Xawith1.row(i);
    DC_est.slice(i) = (xawith1iM.t()*xawith1iM) / (dot(xawith1i, betaest)*dot(xawith1i, betaest));
  }
  arma::mat Dest = mean(DC_est, 2);
  return Dest;
}
// [[Rcpp::export]]
arma::mat f_ASD_WCQE(arma::vec tauseq, arma::vec XitauseqtildeV, arma::mat Dtilde, int N, int p) {
  // f_ASD_WCQE - estimation for ASD of WCQE, with \tau in tauseq
  int dimvarphi = p + 1;
  int taunum = tauseq.n_elem;
  // ASD
  arma::mat ASDM(dimvarphi,taunum); // the j-th column is the estimation of ASD of \hat{\varphi}(\tau_{j})
  for (int j = 0; j < taunum; j++) {
    double tauj = tauseq(j);
    double Xitilde_tauj = XitauseqtildeV(j); // \tilde{\Xi}(\tau_{j})
    arma::mat CovMj = tauj*(1-tauj) * 1.0/(Xitilde_tauj*Xitilde_tauj) * Dtilde.i() / (N*1.0); // covariance matrix of \hat{\varphi}(\tau_{j})
    arma::vec ASDsquarevarphihat_tauj = CovMj.diag(0);
    arma::vec ASDvarphihat_tauj = sqrt(ASDsquarevarphihat_tauj);
    ASDM.col(j) = ASDvarphihat_tauj;
  }
  return ASDM;
}
