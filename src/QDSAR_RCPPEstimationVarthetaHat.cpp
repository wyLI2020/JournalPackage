#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec f_vartheta_hat(arma::vec zetahat, int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde) {
  // f_vartheta_hat - \hat{\vartheta} = \hat{\theta} - Xwith1 \hat{\psi}
  int dimphi = Ztilde.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::vec phihat = zetahat.head(dimphi);
  double lambdahat = zetahat(dimphi+1);
  arma::mat Bhat = I_N - lambdahat * W;
  arma::mat varthetahat_mat(N, T); // matrix of which the t-th column is the estimation of \vartheta in time t
  for (int t = 0; t < T; t++) {
    arma::vec yt = Y.subvec(t*N, (t+1)*N-1);
    arma::mat Ztildet = Ztilde.submat(t*N, 0, (t+1)*N-1, Ztilde.n_cols-1);
    varthetahat_mat.col(t) = Bhat * yt - Ztildet * phihat;
  }
  arma::vec varthetahat = mean(varthetahat_mat, 1); // \hat{\vartheta}
  return varthetahat;
}