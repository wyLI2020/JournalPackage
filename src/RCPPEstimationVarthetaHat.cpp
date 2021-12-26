#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat f3_asmat(arma::vec vec1, int nrow, int ncol) {
  // Fill matrix with elements of vector
  arma::mat vec1_mat(nrow*ncol, 1); vec1_mat.col(0) = vec1;
  vec1_mat.reshape(nrow, ncol);
  return vec1_mat;
}
arma::vec f3_asvec(arma::mat mat1) {
  // Matrix straighten
  int nrow = mat1.n_rows;
  int ncol = mat1.n_cols;
  mat1.reshape(nrow*ncol, 1);
  return mat1.col(0);
}
arma::vec f3_SVec(arma::vec someVec, int N, int T, arma::mat B) {
  // f3_SVec() - S(\lambda) Vec, with S(\lambda) = I_{T} \otimes B_{N}
  arma::mat Vec_mat = f3_asmat(someVec, N, T);
  arma::vec Svec = f3_asvec(B * Vec_mat);
  return Svec;
}
// [[Rcpp::export]]
arma::vec f_V_hat(arma::vec zetahat, int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde) {
  // f_V_hat() - \hat{V}_{NT}, the estimate of V_{NT} = \iota \otimes \vartheta + \varepsilon
  int dimphi = Ztilde.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::vec phihat = zetahat.head(dimphi);
  double lambdahat = zetahat(dimphi+1);
  arma::mat Bhat = I_N - lambdahat * W; // B(\hat{\lambda}) = I_{N} - \hat{\lambda}
  arma::vec ShatY = f3_SVec(Y, N, T, Bhat); // S(\hat{\lambda}) Y
  arma::vec Vhat = ShatY - Ztilde*phihat; // \hat{V}_{NT} = S(\hat{\lambda}) Y - \tilde{Z} \hat{\phi}
  return Vhat;
}
// [[Rcpp::export]]
arma::vec f_vartheta_hat(arma::vec Vhat, int N, int T) {
  // f_vartheta_hat() - \hat{\vartheta} with \vartheta = \theta - Xwith1 \psi
  arma::mat vhat_mat = f3_asmat(Vhat, N, T);
  arma::vec varthetahat = mean(vhat_mat, 1); // \hat{\vartheta}
  return varthetahat;
}
