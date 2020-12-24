#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double RCPPlnlOFdelta(arma::vec delta, int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1) {
  // RCPPlnlOFdelta() - \ln l(delta)
  int dimbeta = Xawith1.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat rowvec1_T(1, T); rowvec1_T.fill(1.0);
  arma::mat J_T = vec1_T * rowvec1_T;
  arma::mat I_NT(N*T,N*T); I_NT.eye(N*T,N*T);
  arma::mat I_T(T,T); I_T.eye(T,T);
  // \lambda and \beta^\star
  double lambda = delta(0);
  arma::vec betastar = delta.tail(dimbeta);
  // \Omega(\beta^\star)
  arma::mat A(N,N); A.fill(0.0); // 1/\sigma^{2}_{\epsilon} Cov(\vartheta)
  for(int i = 0; i < N; i++) {
    A(i,i) = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar);
  }
  arma::mat Omega = kron(J_T, A) + I_NT;
  // S(\lambda)
  arma::mat B = I_N - lambda * W;
  arma::mat S = kron(I_T, B);
  // \hat{phi}(\delta)
  arma::mat left = Ztilde.t() * Omega.i() * Ztilde;
  arma::vec right = Ztilde.t() * Omega.i() * S * Y;
  arma::vec phihat = left.i() * right;
  // \hat{\sigma_\epsilon^2}(\delta)
  arma::vec Vhat = S * Y - Ztilde * phihat;
  double sigmaepsilon2hat = dot(Vhat, Omega.i()*Vhat) / (N*T*1.0);
  // \ln l(\delta)
  double valueOmega, signOmega; log_det(valueOmega, signOmega, Omega);
  // double detOmega = exp(valueOmega)*signOmega;
  double valueS, signS; log_det(valueS, signS, S);
  // double detS = exp(valueS)*signS;
  // double lnldelta = - N*T/2.0*(log(2*3.14)+1) - N*T/2.0*log(sigmaepsilon2hat) - 1.0/2.0*log(detOmega) + log(detS);
  double lnldelta = - N*T/2.0*(log(2*3.14)+1) - N*T/2.0*log(sigmaepsilon2hat) - 1.0/2.0*valueOmega + valueS;
  return lnldelta;
}
// // [[Rcpp::export]]
// double rhotau(double tau, double u) {
//   // rhotau() - \rho_\tau(u)
//   double result = 0.0;
//   if(u < 0) {
//     result = (tau - 1) * u;
//   } else {
//     result = tau * u;
//   }
//   return result;
// }
