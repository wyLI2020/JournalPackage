#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec f_chi_hat(arma::vec deltahat, int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1) {
  // f_chi_hat() - \hat{chi}(\hat{delta}): \hat{phi}(\hat{delta}) and \hat{\sigma}^{2}_{\epsilon}(\hat{delta})
  int dimbeta = Xawith1.n_cols;
  int dimphi = Ztilde.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat rowvec1_T(1, T); rowvec1_T.fill(1.0);
  arma::mat J_T = vec1_T * rowvec1_T;
  arma::mat I_NT(N*T,N*T); I_NT.eye(N*T,N*T);
  arma::mat I_T(T,T); I_T.eye(T,T);
  // \lambda and \beta^{\star}
  double lambdahat = deltahat(0);
  arma::vec betastarhat = deltahat.tail(dimbeta);
  // \Omega(\beta^{\star})
  arma::mat Ahat(N,N); Ahat.fill(0.0); // 1/\sigma^{2}_{\epsilon} Cov(\vartheta)
  for(int i = 0; i < N; i++) {
    Ahat(i,i) = dot(Xawith1.row(i), betastarhat) * dot(Xawith1.row(i), betastarhat);
  }
  arma::mat Omegahat = kron(J_T, Ahat) + I_NT;
  // S(\lambda)
  arma::mat Bhat = I_N - lambdahat * W;
  arma::mat Shat = kron(I_T, Bhat);
  // \hat{phi}(\delta)
  arma::mat lefthat = Ztilde.t() * Omegahat.i() * Ztilde;
  arma::vec righthat = Ztilde.t() * Omegahat.i() * Shat * Y;
  arma::vec phihat = lefthat.i() * righthat;
  // \hat{\sigma}^{2}_{\epsilon}(\delta)
  arma::vec Vhat = Shat * Y - Ztilde * phihat;
  double sigmaepsilon2hat = dot(Vhat, Omegahat.i()*Vhat) / (N*T*1.0);
  // \hat{\chi}(\delta)
  arma::vec chihat(dimphi+1);
  chihat.head(dimphi) = phihat;
  chihat(dimphi) = sigmaepsilon2hat;
  return chihat;
}
// [[Rcpp::export]]
arma::vec f_zeta_hat(int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec deltahat) {
  // f_zeta_hat() - QMLE \hat{zeta} = (\hat{\phi}^\prime, \hat{\sigma}^{2}_{\epsilon}, \hat{\delta}^\prime)^\prime = (\hat{\phi}^\prime, \hat^{2}_{\sigma_\epsilon}, \hat{lambda}, \hat{\beta}^{\star \prime})^\prime
  int dimphi = Ztilde.n_cols;
  int dimbeta = Xawith1.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat vec1_N(N, 1); vec1_N.fill(1.0);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat vec1_NT(N*T, 1); vec1_NT.fill(1.0);
  // \hat{chi}, \chi = (\phi^\prime, \sigma^{2}_\epsilon)^\prime
  arma::vec chihat = f_chi_hat(deltahat, N, T, W, Y, Ztilde, Xawith1);
  // \hat{zeta}, \zeta = (\phi^\prime, \sigma^{2}_\epsilon, \lambda, \beta^\star)^\prime
  arma::vec zetahat(dimphi + dimbeta + 2);
  zetahat.head(dimphi+1) = chihat;
  zetahat.tail(dimbeta+1) = deltahat;
  return zetahat;
}

