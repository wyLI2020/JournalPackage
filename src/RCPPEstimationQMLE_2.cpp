#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat f2_asmat(arma::vec vec1, int nrow, int ncol) {
  // Fill matrix with elements of vector
  arma::mat vec1_mat(nrow*ncol, 1); vec1_mat.col(0) = vec1;
  vec1_mat.reshape(nrow, ncol);
  return vec1_mat;
}
arma::vec f2_asvec(arma::mat mat1) {
  // Matrix straighten
  int nrow = mat1.n_rows;
  int ncol = mat1.n_cols;
  mat1.reshape(nrow*ncol, 1);
  return mat1.col(0);
}
arma::vec f2_SVec(arma::vec someVec, int N, int T, arma::mat B) {
  // f2_SVec() - S(\lambda) Vec, with S(\lambda) = I_{T} \otimes B_{N}
  arma::mat Vec_mat = f2_asmat(someVec, N, T);
  arma::vec Svec = f2_asvec(B * Vec_mat);
  return Svec;
}
arma::vec f2_IotimesWVec(arma::vec someVec, int N, int T, arma::mat W) {
  // f_IotimesWVec - (I_{T} \otimes W_{N}) Vec
  arma::mat Vec_mat = f2_asmat(someVec, N, T);
  arma::vec result = f2_asvec(W * Vec_mat);
  return result;
}
arma::vec f2_OmegainverseVec(arma::vec someVec, int N, int T, arma::mat IplusTAinverse) {
  // f2_OmegainverseVec() - \Omega(\beta^{\star})^{-1} Vec, with \Omega(\beta^{\star})^{-1} = \frac{1}{T} J_{T} \otimes (I_{N} + T A_{N})^{-1} + I_{NT} - \frac{1}{T} J_{T} \otimes I_{N}
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  // \Omega(\beta^{\star})^{-1} Vec = \frac{1}{T} J_{T} \otimes (I_{N} + T A_{N})^{-1} Vec + Vec - \frac{1}{T} J_{T} \otimes I_{N} Vec
  arma::mat Vec_mat = f2_asmat(someVec, N, T);
  arma::vec OmegainverseVec_1 = 1.0/(T*1.0) * kron(vec1_T, sum(IplusTAinverse*Vec_mat, 1)); // \frac{1}{T} J_{T} \otimes (I_{N} + T A_{N})^{-1} Vec
  arma::vec OmegainverseVec_2 = someVec; // Vec
  arma::vec OmegainverseVec_3 = - 1.0/(T*1.0) * kron(vec1_T, sum(Vec_mat, 1)); // - \frac{1}{T} J_{T} \otimes I_{N} Vec
  arma::vec OmegainverseVec = OmegainverseVec_1 + OmegainverseVec_2 + OmegainverseVec_3;
  return OmegainverseVec;
}
arma::vec f2_DO1betastarOmegainverseVec_l(int l, arma::vec someVec, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f2_DO1betastarOmegainverseVec_l() - \frac{\partial}{\partial\beta^{\star}_{l}}\Omega \Omega^{-1} Vec, with the left matrix = J_{T} \otimes [\diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})} (I_{NT} + T A_{N})^{-1}]
  arma::vec diagVec(N); // diagonal of matrix [\diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})} (I_{NT} + T A_{N})^{-1}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d1aii_l = 2.0 * Xawith1(i, l) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})}
    diagVec(i) = d1aii_l * (1.0 / (1.0 + 1.0 * T * aii));
  }
  arma::mat RightM(N,N); RightM.fill(0.0); RightM.diag() = diagVec; // original matirx = J_{T} \otimes RightM
  arma::mat Vec_mat = f2_asmat(someVec, N, T);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::vec result = kron(vec1_T, sum(RightM*Vec_mat, 1));
  return result;
}
arma::vec f_OmegainverseDO1betastarOmegainverseVec_l(int l, arma::vec someVec, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f_OmegainverseDO1betastarOmegainverseVec_l() - \Omega^{-1} \frac{\partial}{\partial\beta^{\star}_{l}}\Omega \Omega^{-1} Vec, with the left matirx = J_{T} \otimes [(I_{NT} + T A_{N})^{-1} \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})} (I_{NT} + T A_{N})^{-1}]
  arma::vec diagVec(N); // diagonal of matrix [(I_{NT} + T A_{N})^{-1} \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})} (I_{NT} + T A_{N})^{-1}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d1aii_l = 2.0 * Xawith1(i, l) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})}
    diagVec(i) = (1.0 / (1.0 + 1.0 * T * aii)) * d1aii_l * (1.0 / (1.0 + 1.0 * T * aii));
  }
  arma::mat RightM(N,N); RightM.fill(0.0); RightM.diag() = diagVec; // original matirx = J_{T} \otimes RightM
  arma::mat Vec_mat = f2_asmat(someVec, N, T);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::vec result = kron(vec1_T, sum(RightM*Vec_mat, 1));
  return result;
}
arma::vec f_DO2betastarOmegainverseVec_l2l1(int l1, int l2, arma::vec someVec, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f_DO2betastarOmegainverseVec_l2l1() - \frac{\partial^2}{\partial\beta^{\star}_{l2}\partial\beta^{\star}_{l1}}\Omega \Omega^{-1} Vec, with the left matrix = J_{T} \otimes [\diag{2 x_{al2i} x_{al1i}} (I_{NT} + T A_{N})^{-1}]
  arma::vec diagVec(N); // diagonal of matrix [\diag{2 x_{al2i} x_{al1i}} (I_{NT} + T A_{N})^{-1}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d2aii_l2l1 = 2.0 * Xawith1(i, l2) * Xawith1(i, l1); // element (i,i) in matrix \diag{2 x_{al2i} x_{al1i}}
    diagVec(i) = d2aii_l2l1 * (1.0 / (1.0 + 1.0 * T * aii));
  }
  arma::mat RightM(N,N); RightM.fill(0.0); RightM.diag() = diagVec; // original matirx = J_{T} \otimes RightM
  arma::mat Vec_mat = f2_asmat(someVec, N, T);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::vec result = kron(vec1_T, sum(RightM*Vec_mat, 1));
  return result;
}
double f_trOmegainverseDO1betastarOmegainverseDO1betastar_l2l1(int l1, int l2, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f_trOmegainverseDO1betastarOmegainverseDO1betastar_l2l1() - tr[\Omega^{-1} \frac{\partial}{\partial\beta^{\star}_{l2}}\Omega \Omega^{-1} \frac{\partial}{\partial\beta^{\star}_{l1}}\Omega], with the matrix = T J_{T} \otimes [(I_{NT} + T A_{N})^{-1} \diag{2 x_{al2i} x_{ai}^{\prime} \beta^{\star}} (I_{NT} + T A_{N})^{-1} \diag{2 x_{al1i} x_{ai}^{\prime} \beta^{\star}}]
  arma::vec diagVec(N); // diagonal of matrix [(I_{NT} + T A_{N})^{-1} \diag{2 x_{al2i} x_{ai}^{\prime} \beta^{\star}} (I_{NT} + T A_{N})^{-1} \diag{2 x_{al1i} x_{ai}^{\prime} \beta^{\star}}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d1aii_l1 = 2.0 * Xawith1(i, l1) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix \diag{2 x_{al1i} (x_{ai}^{\prime} \beta^{\star})}
    double d1aii_l2 = 2.0 * Xawith1(i, l2) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix \diag{2 x_{al2i} (x_{ai}^{\prime} \beta^{\star})}
    diagVec(i) = (1.0 / (1.0 + 1.0 * T * aii)) * d1aii_l2 * (1.0 / (1.0 + 1.0 * T * aii)) * d1aii_l1;
  }
  double result = 1.0 * T * T * sum(diagVec); // original matrix = T J_{T} \otimes diag{diagVec}
  return result;
}
double f_trOmegainverseDO2betastar_l2l1(int l1, int l2, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f_trOmegainverseDO2betastar_l2l1() - tr[\Omega^{-1} \frac{\partial^2}{\partial\beta^{\star}_{l2}\partial\beta^{\star}_{l1}}\Omega], with the matrix = J_{T} \otimes [(I_{NT} + T A_{N})^{-1} \diag{2 x_{al2i} x_{al1i}}]
  arma::vec diagVec(N); // diagonal of matrix [(I_{NT} + T A_{N})^{-1} \diag{2 x_{al2i} x_{al1i}}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d2aii_l2l1 = 2.0 * Xawith1(i, l2) * Xawith1(i, l1); // element (i,i) in matrix \diag{2 x_{al2i} x_{al1i}}
    diagVec(i) = (1.0 / (1.0 + 1.0 * T * aii)) * d2aii_l2l1;
  }
  double result = 1.0 * T * sum(diagVec); // original matrix = J_{T} \otimes diag{diagVec}
  return result;
}
arma::vec Cf_score_b(int b, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat, arma::vec Vhat) {
  // Cf_score_b() - the bootstrapped score function
  Environment OPT = Environment::namespace_env("QuantileDSAR");
  Function f = OPT["Rf_score_b"];
  // Function f("Rf_score_b");
  NumericVector result_rcpp = f(b, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat);
  arma::vec result_arma(result_rcpp.begin(), result_rcpp.size(), false);
  return result_arma;
}
arma::mat f_boots_VCmatrix(int Bnum, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat, arma::vec Vhat) {
  // f_boots_VCmatrix() - in b-th bootstrap, firstly obtain T matched booststrap samples {\hat{r}_{1}^{b}, ..., \hat{r}_{T}^{b}}, then calculate the estimation of VC matrix, which follows (Su and Yang, 2015)
  int dimzeta = zetahat.n_elem;
  arma::mat scoreBnum(dimzeta, Bnum); scoreBnum.fill(0.0); // the b-th column is scoreb
  arma::cube scorescoreTBnum(dimzeta, dimzeta, Bnum); scorescoreTBnum.fill(0.0); // the b-th slice is (scoreb %*% score^\prime), where scoreb is the estimation of score function in b-th bootstrap
  for (int b = 0; b < Bnum; b++) {
    arma::vec scoreb = Cf_score_b(b, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat);
    arma::mat scorebM(dimzeta, 1); scorebM.col(0) = scoreb;
    scoreBnum.col(b) = scoreb;
    scorescoreTBnum.slice(b) = scorebM * scorebM.t();
  }
  // calculate the estimation of VC matrix, that is Sigma^{\star}_{N} + \Sigma_{N} = E(\frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta} \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta^{\prime}})
  arma::mat sec1 = mean(scorescoreTBnum, 2); // \frac{1}{b} \sum\limits_{b=1}^{Bnum} (scoreb %*% scoreb^{\prime})
  arma::mat sec2(dimzeta, 1); // \frac{1}{B} \sum\limits_{b=1}^{Bnum} scoreb
  sec2.col(0) = mean(scoreBnum, 1);
  arma::mat VCmatrixhat = sec1 - sec2 * sec2.t();
  return VCmatrixhat;
}
arma::mat f_SigmaPlusSigmastarhat(int Bnum, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat, arma::vec Vhat) {
  // f_SigmaPlusSigmastarhat() - estimation of Sigma^{\star}_{N} + \Sigma_{N} = \frac{1}{NT} VC matrix
  arma::mat VCmatrixhat = f_boots_VCmatrix(Bnum, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat); // estimation of VC matrix = E(\frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta} \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta^{\prime}})
  arma::mat SigmaPlusSigmastarhat = VCmatrixhat / (N*T*1.0); 
  return SigmaPlusSigmastarhat;
}
arma::mat f_Sigmahat(int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat, arma::vec Vhat) {
  // f_Sigmahat() - estimation of matrix \Sigma_{N}, \hat{\Sigma}_{N} = - \frac{1}{NT} \frac{\partial^2}{\partial zeta \partial zeta^\prime}\lnL_{NT}(\zeta_{0})
  int dimzeta = Ztilde.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  double sigmaepsilon2hat = zetahat(p+q+2);
  double lambdahat = zetahat(p+q+3);
  arma::vec betastarhat = zetahat.tail(p+1);
  arma::mat IplusTAhatinverse(N,N); IplusTAhatinverse.fill(0.0); // (I_{N} + T \hat{A}_{N})^{-1}
  for(int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastarhat) * dot(Xawith1.row(i), betastarhat); // A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    IplusTAhatinverse(i,i) = 1.0 / (1.0 + 1.0 * T * aii);
  }
  arma::mat Bhat = I_N - lambdahat * W; // B(\hat{\lambda}) = I_{N} - \hat{\lambda}
  arma::vec OmegainverseVhat = f2_OmegainverseVec(Vhat, N, T, IplusTAhatinverse); // \Omega^{-1}(\hat{\beta}^{\star}) \hat{V}
  arma::mat OmegainverseZtildehat(N*T, Ztilde.n_cols); OmegainverseZtildehat.fill(0.0); // \Omega^{-1}(\hat{\beta}^{\star}) \tide{Z}
  for(int h = 0; h < dimzeta; h++) {
    arma::vec ztilde_h = Ztilde.col(h);
    arma::vec Omegainverseztildehat_h = f2_OmegainverseVec(ztilde_h, N, T, IplusTAhatinverse); 
    OmegainverseZtildehat.col(h) = Omegainverseztildehat_h;
  }
  arma::vec IotimesWY = f2_IotimesWVec(Y, N, T, W); // (I_{T} \otimes W_{N}) Y
  arma::vec OmegainverseIotimesWYhat = f2_OmegainverseVec(IotimesWY, N, T, IplusTAhatinverse); // \Omega^{-1}(\hat{\beta}^{\star}) (I_{T} \otimes W_{N}) Y
  // Second partial derivatives
  arma::mat D2_phiphi = - 1.0/sigmaepsilon2hat * Ztilde.t()*OmegainverseZtildehat;
  arma::mat D2_sigmaphi = - 1.0/(sigmaepsilon2hat*sigmaepsilon2hat) * (OmegainverseZtildehat.t()*Vhat).t();
  double D2_sigmasigma = - 1.0/(sigmaepsilon2hat*sigmaepsilon2hat*sigmaepsilon2hat) * dot(Vhat, OmegainverseVhat) + 1.0*N*T/(2.0*sigmaepsilon2hat*sigmaepsilon2hat);
  arma::mat D2_lambdaphi = - 1.0/sigmaepsilon2hat * (OmegainverseZtildehat.t()*IotimesWY).t();
  double D2_lambdasigma = - 1.0/(sigmaepsilon2hat*sigmaepsilon2hat) * dot(IotimesWY, OmegainverseVhat);
  double D2_lambdalambda = - 1.0*T * trace(Bhat.i()*W*Bhat.i()*W) - 1.0/sigmaepsilon2hat * dot(IotimesWY, OmegainverseIotimesWYhat);
  arma::mat D2_betastarphi(1+p,p+q+2); D2_betastarphi.fill(0.0);
  arma::mat D2_betastarsigma(1+p,1); D2_betastarsigma.fill(0.0);
  arma::mat D2_betastarlambda(1+p,1); D2_betastarlambda.fill(0.0);
  for (int l = 0; l < (1+p); l++) {
    arma::vec DO1betastarOmegainverseVhat_l = f2_DO1betastarOmegainverseVec_l(l, Vhat, N, T, betastarhat, Xawith1); // \frac{\partial}{\partial \beta^{\star}_{l}}\Omega(\hat{\beta}^{\star}) \Omega^{-1}(\hat{\beta}^{\star}) \hat{V}
    arma::mat DO1betastarOmegainverseZtildehat_l(N*T, Ztilde.n_cols); DO1betastarOmegainverseZtildehat_l.fill(0.0); // \frac{\partial}{\partial \beta^{\star}_{l}}\Omega(\hat{\beta}^{\star}) \Omega^{-1}(\hat{\beta}^{\star}) \tilde{Z}
    for(int h = 0; h < dimzeta; h++) {
      arma::vec ztilde_h = Ztilde.col(h);
      arma::vec DO1betastarOmegainverseztildehat_lh = f2_DO1betastarOmegainverseVec_l(l, ztilde_h, N, T, betastarhat, Xawith1);
      DO1betastarOmegainverseZtildehat_l.col(h) = DO1betastarOmegainverseztildehat_lh;
    }
    D2_betastarphi.row(l) = - 1.0/sigmaepsilon2hat * (DO1betastarOmegainverseZtildehat_l.t()*OmegainverseVhat).t();
    D2_betastarsigma(l,0) = - 1.0/(2.0*sigmaepsilon2hat*sigmaepsilon2hat) * dot(OmegainverseVhat, DO1betastarOmegainverseVhat_l);
    D2_betastarlambda(l,0) = - 1.0/sigmaepsilon2hat * dot(OmegainverseIotimesWYhat, DO1betastarOmegainverseVhat_l);
  }
  arma::mat D2_betastarbetastar(1+p,1+p); D2_betastarbetastar.fill(0.0);
  for (int l2 = 0; l2 < (1+p); l2++) {
    for (int l1 = 0; l1 < (1+p); l1++) {
      arma::vec OmegainverseDO1betastarOmegainverseVhat_l1 = f_OmegainverseDO1betastarOmegainverseVec_l(l1, Vhat, N, T, betastarhat, Xawith1); // \Omega^{-1}(\hat{\beta}^{\star}) \frac{\partial}{\partial \beta^{\star}_{l1}}\Omega(\hat{\beta}^{\star}) \Omega^{-1}(\hat{\beta}^{\star}) \hat{V}
      arma::vec DO1betastarOmegainverseVhat_l2 = f2_DO1betastarOmegainverseVec_l(l2, Vhat, N, T, betastarhat, Xawith1); // \frac{\partial}{\partial \beta^{\star}_{l2}}\Omega(\hat{\beta}^{\star}) \Omega^{-1}(\hat{\beta}^{\star}) \hat{V}
      arma::vec DO2betastarOmegainverseVhat_l2l1 = f_DO2betastarOmegainverseVec_l2l1(l1, l2, Vhat, N, T, betastarhat, Xawith1); // \frac{\partial^2}{\partial \beta^{\star}_{l2} \partial \beta^{\star}_{l1}}\Omega(\hat{\beta}^{\star}) \Omega^{-1}(\hat{\beta}^{\star}) \hat{V}
      double trOmegainverseDO1betastarOmegainverseDO1betastarhat_l2l1 = f_trOmegainverseDO1betastarOmegainverseDO1betastar_l2l1(l1, l2, N, T, betastarhat, Xawith1); // tr[\Omega^{-1}(\hat{\beta}^{\star}) \frac{\partial}{\partial \beta^{\star}_{l2}}\Omega(\hat{\beta}^{\star}) \Omega^{-1}(\hat{\beta}^{\star}) \frac{\partial}{\partial \beta^{\star}_{l1}}\Omega(\hat{\beta}^{\star})]
      double trOmegainverseDO2betastarhat_l2l1 = f_trOmegainverseDO2betastar_l2l1(l1, l2, N, T, betastarhat, Xawith1); // tr[\Omega^{-1}(\hat{\beta}^{\star}) \frac{\partial^2}{\partial \beta^{\star}_{l2} \partial \beta^{\star}_{l1}}\Omega(\hat{\beta}^{\star})]
      D2_betastarbetastar(l2,l1) = - 1.0/2.0 * (- trOmegainverseDO1betastarOmegainverseDO1betastarhat_l2l1 + trOmegainverseDO2betastarhat_l2l1) + 1.0/(2.0*sigmaepsilon2hat) * (-2*dot(DO1betastarOmegainverseVhat_l2, OmegainverseDO1betastarOmegainverseVhat_l1) + dot(OmegainverseVhat, DO2betastarOmegainverseVhat_l2l1));
    }
  }
  arma::mat D2(2*p+q+5, 2*p+q+5); D2.fill(0.0);
  D2.submat(0,0,p+q+1,p+q+1) = D2_phiphi;
  D2.submat(p+q+2,0,p+q+2,p+q+1) = D2_sigmaphi; D2(p+q+2,p+q+2) = D2_sigmasigma;
  D2.submat(p+q+3,0,p+q+3,p+q+1) = D2_lambdaphi; D2(p+q+3,p+q+2) = D2_lambdasigma; D2(p+q+3,p+q+3) = D2_lambdalambda;
  D2.submat(p+q+4,0,2*p+q+4,p+q+1) = D2_betastarphi; D2.submat(p+q+4,p+q+2,2*p+q+4,p+q+2) = D2_betastarsigma; D2.submat(p+q+4,p+q+3,2*p+q+4,p+q+3) = D2_betastarlambda; D2.submat(p+q+4,p+q+4,2*p+q+4,2*p+q+4) = D2_betastarbetastar;
  D2.submat(0,p+q+2,p+q+1,p+q+2) = D2_sigmaphi.t();
  D2.submat(0,p+q+3,p+q+1,p+q+3) = D2_lambdaphi.t(); D2(p+q+2,p+q+3) = D2_lambdasigma;
  D2.submat(0,p+q+4,p+q+1,2*p+q+4) = D2_betastarphi.t(); D2.submat(p+q+2,p+q+4,p+q+2,2*p+q+4) = D2_betastarsigma.t(); D2.submat(p+q+3,p+q+4,p+q+3,2*p+q+4) = D2_betastarlambda.t();
  // \hat{\Sigma}_{N} = - \frac{1}{NT} \frac{\partial^2 \lnL_{NT}(\zeta_{0})}{\partial zeta \partial zeta^\prime} in \Sigma_{N}
  arma::mat Sigmahat = - D2 / (N*T*1.0);
  return Sigmahat;
}
// [[Rcpp::export]]
arma::mat f_CovM_QMLE(int Bnum, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat, arma::vec Vhat) {
  // f_CovM_QMLE() - estimation for the covariance matrix of QMLE \hat{\zeta}
  arma::mat SigmaPlusSigmastarhat = f_SigmaPlusSigmastarhat(Bnum, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat);
  arma::mat Sigmahat = f_Sigmahat(N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat);
  // estimation of the covariance matrix
  arma::mat CovMzetahat = Sigmahat.i() * SigmaPlusSigmastarhat * Sigmahat.i() / (N*T*1.0);
  return CovMzetahat;
}
// [[Rcpp::export]]
arma::mat f_CovM_secbetacheck(int p, int q, arma::vec zetahat, arma::mat CovMzetahat) {
  // f_CovM_secbetacheck() - estimation for the covariance matrix of c(\check{\beta}_{1}, ..., \check{\beta}_{p}) in \check{\beta}
  int dimzeta = zetahat.n_elem;
  arma::vec betastarhat = zetahat.tail(1+p); // \hat{\beta}^{\star}
  double betastar0hat = betastarhat(0); // \hat{\beta}^{\star}_{0}
  arma::mat derivamat(dimzeta, p); derivamat.fill(0.0); // derivative matrix in delta method
  for (int l = 0; l < p; l++) {
    arma::mat derivavec_l(dimzeta, 1); derivavec_l.fill(0.0); // the l-th column of derivative matrix, that is derivative vector \frac{\partial}{\partial \zeta}h(\zeta) of transformation function h(\zeta) = \beta^{\star}_{l} / \beta^{\star}_{0}
    double betastarlhat = betastarhat(l+1); // \hat{\beta}^{\star}_{l}
    derivavec_l(4+p+q, 0) = - betastarlhat / (betastar0hat*betastar0hat); // \frac{\partial}{\partial \beta^{\star}_{0}} h(\zeta)
    derivavec_l(4+p+q+(l+1), 0) = 1.0 / betastar0hat; // \frac{\partial}{\partial \beta^{\star}_{l}} h(\zeta)
    derivamat.col(l) = derivavec_l.col(0);
  }
  arma::mat CovMsecbetacheck = derivamat.t() * CovMzetahat * derivamat;
  return CovMsecbetacheck;
}
