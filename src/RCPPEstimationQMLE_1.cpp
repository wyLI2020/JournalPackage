#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat f1_asmat(arma::vec vec1, int nrow, int ncol) {
  // Fill matrix with elements of vector
  arma::mat vec1_mat(nrow*ncol, 1); vec1_mat.col(0) = vec1;
  vec1_mat.reshape(nrow, ncol);
  return vec1_mat;
}
arma::vec f1_asvec(arma::mat mat1) {
  // Matrix straighten
  int nrow = mat1.n_rows;
  int ncol = mat1.n_cols;
  mat1.reshape(nrow*ncol, 1);
  return mat1.col(0);
}
arma::vec f1_SVec(arma::vec someVec, int N, int T, arma::mat B) {
  // f1_SVec() - S(\lambda) Vec, with S(\lambda) = I_{T} \otimes B_{N}
  arma::mat Vec_mat = f1_asmat(someVec, N, T);
  arma::vec Svec = f1_asvec(B * Vec_mat);
  return Svec;
}
arma::vec f1_OmegainverseVec(arma::vec someVec, int N, int T, arma::mat IplusTAinverse) {
  // f1_OmegainverseVec() - \Omega(\beta^{\star})^{-1} Vec, with \Omega(\beta^{\star})^{-1} = \frac{1}{T} J_{T} \otimes (I_{N} + T A_{N})^{-1} + I_{NT} - \frac{1}{T} J_{T} \otimes I_{N}
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  // \Omega(\beta^{\star})^{-1} Vec = \frac{1}{T} J_{T} \otimes (I_{N} + T A_{N})^{-1} Vec + Vec - \frac{1}{T} J_{T} \otimes I_{N} Vec
  arma::mat Vec_mat = f1_asmat(someVec, N, T);
  arma::vec OmegainverseVec_1 = 1.0/(T*1.0) * kron(vec1_T, sum(IplusTAinverse*Vec_mat, 1)); // \frac{1}{T} J_{T} \otimes (I_{N} + T A_{N})^{-1} Vec
  arma::vec OmegainverseVec_2 = someVec; // Vec
  arma::vec OmegainverseVec_3 = - 1.0/(T*1.0) * kron(vec1_T, sum(Vec_mat, 1)); // - \frac{1}{T} J_{T} \otimes I_{N} Vec
  arma::vec OmegainverseVec = OmegainverseVec_1 + OmegainverseVec_2 + OmegainverseVec_3;
  return OmegainverseVec;
}
arma::vec f_chi_delta(arma::vec delta, int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1) {
  // f_chi_delta() - the function \hat{chi}(\delta): \hat{phi}(\delta) and \hat{\sigma}^{2}_{\varepsilon}(\delta)
  int dimbeta = Xawith1.n_cols;
  int dimphi = Ztilde.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  // \lambda and \beta^{\star}
  double lambda = delta(0);
  arma::vec betastar = delta.tail(dimbeta);
  // (I_{N} + T A_{N}(\beta^{\star}))^{-1} in \Omega(\beta^{\star})^{-1}, A_{N} = 1/\sigma^{2}_{\varepsilon} Cov(\vartheta)
  arma::mat IplusTAinverse(N,N); IplusTAinverse.fill(0.0); // (I_{N} + T A_{N})^{-1}
  for(int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    IplusTAinverse(i,i) = 1.0 / (1.0 + 1.0 * T * aii);
  }
  // B(\lambda)
  arma::mat B = I_N - lambda * W;
  // \hat{phi}(\delta)
  arma::mat OmegainverseZtilde(N*T,dimphi); OmegainverseZtilde.fill(0.0); // \Omega^{-1}(\beta^{\star}) \tilde{Z}
  for(int l = 0; l < dimphi; l++) {
    arma::vec ztilel = Ztilde.col(l); // the l-th column of \tilde{Z}
    OmegainverseZtilde.col(l) = f1_OmegainverseVec(ztilel, N, T, IplusTAinverse); 
  }
  arma::vec SY = f1_SVec(Y, N, T, B); // S(\lambda) Y = (I_{T} \otimes B_{N}) Y
  arma::vec OmegainverseSY = f1_OmegainverseVec(SY, N, T, IplusTAinverse); // \Omega^{-1}(\beta^{\star}) S(\lambda) Y
  arma::mat left = Ztilde.t() * OmegainverseZtilde;
  arma::vec right = Ztilde.t() * OmegainverseSY;
  arma::vec phi = left.i() * right;
  // \hat{\sigma}^{2}_{\varepsilon}(\delta)
  arma::vec V = SY - Ztilde * phi;
  arma::vec OmegainverseV = f1_OmegainverseVec(V, N, T, IplusTAinverse);
  double sigmaepsilon2 = dot(V, OmegainverseV) / (N*T*1.0);
  // \hat{\chi}(\delta)
  arma::vec chi(dimphi+1);
  chi.head(dimphi) = phi;
  chi(dimphi) = sigmaepsilon2;
  return chi;
}
// [[Rcpp::export]]
double RCPPlnlOFdelta(arma::vec delta, int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1) {
  // RCPPlnlOFdelta() - \ln l(delta)
  int dimbeta = Xawith1.n_cols; 
  int dimphi = Ztilde.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat rowvec1_T(1, T); rowvec1_T.fill(1.0);
  arma::mat J_T = vec1_T * rowvec1_T;
  arma::mat I_T(T,T); I_T.eye(T,T);
  // \lambda and \beta^\star
  double lambda = delta(0);
  arma::vec betastar = delta.tail(dimbeta);
  // \hat{\sigma}^{2}_{\varepsilon}(\delta)
  arma::vec chi = f_chi_delta(delta, N, T, W, Y, Ztilde, Xawith1);
  double sigmaepsilon2 = chi(dimphi);
  // A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
  arma::vec Adiag(N); // diagonal of A_{N}(\beta^{\star})
  for(int i = 0; i < N; i++) {
    Adiag(i) = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar);
  }
  // B(\lambda)
  arma::mat B = I_N - lambda * W;
  // \ln l(\delta)
  double logdetOmega = sum(log(1.0 + 1.0 * T * Adiag));  // |\Omega| = \Pi_{i=1}^{N} (1 + T a_{ii}) with a_{ii} denotes the element of A_{N} at i-th row ang i-th column, and then ln |\Omega| = \sum_{i=1}^{N} ln(1 + T a_{ii})
  double logdetB, signdetB; log_det(logdetB, signdetB, B);
  double logdetS = 1.0 * T * logdetB; // |S_{NT}| = |B_{N}|^{T}, and then ln |S_{NT}| = T ln |B_{N}|
  double lnldelta = - N*T/2.0*(log(2*3.14)+1) - N*T/2.0*log(sigmaepsilon2) - 1.0/2.0*logdetOmega + logdetS;
  return lnldelta;
}
// [[Rcpp::export]]
arma::vec f_zeta_hat(int N, int T, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec deltahat) {
  // f_zeta_hat() - QMLE \hat{zeta} = (\hat{\phi}^\prime, \hat{\sigma}^{2}_{\varepsilon}, \hat{\delta}^\prime)^\prime = (\hat{\phi}^\prime, \hat_{\sigma}^{2}_{\varepsilon}, \hat{lambda}, \hat{\beta}^{\star \prime})^\prime
  int dimphi = Ztilde.n_cols;
  int dimbeta = Xawith1.n_cols;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat vec1_N(N, 1); vec1_N.fill(1.0);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat vec1_NT(N*T, 1); vec1_NT.fill(1.0);
  // \hat{chi}, \chi = (\phi^\prime, \sigma^{2}_\varepsilon)^\prime
  arma::vec chihat = f_chi_delta(deltahat, N, T, W, Y, Ztilde, Xawith1);
  // \hat{zeta}, \zeta = (\phi^\prime, \sigma^{2}_\varepsilon, \lambda, \beta^\star)^\prime
  arma::vec zetahat(dimphi + dimbeta + 2);
  zetahat.head(dimphi+1) = chihat;
  zetahat.tail(dimbeta+1) = deltahat;
  return zetahat;
}
arma::vec f1_IotimesWVec(arma::vec someVec, int N, int T, arma::mat W) {
  // f_IotimesWVec - (I_{T} \otimes W_{N}) Vec
  arma::mat Vec_mat = f1_asmat(someVec, N, T);
  arma::vec result = f1_asvec(W * Vec_mat);
  return result;
}
arma::vec f1_DO1betastarOmegainverseVec_l(int l, arma::vec someVec, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f1_DO1betastarOmegainverseVec_l() - \frac{\partial}{\partial\beta^{\star}_{l}}\Omega \Omega^{-1} Vec, with the left matrix = J_{T} \otimes [\diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})} (I_{NT} + T A_{N})^{-1}]
  arma::vec diagVec(N); // diagonal of matrix [\diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})} (I_{NT} + T A_{N})^{-1}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d1aii_l = 2.0 * Xawith1(i, l) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})}
    diagVec(i) = d1aii_l * (1.0 / (1.0 + 1.0 * T * aii));
  }
  arma::mat RightM(N,N); RightM.fill(0.0); RightM.diag() = diagVec; // original matirx = J_{T} \otimes RightM
  arma::mat Vec_mat = f1_asmat(someVec, N, T);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::vec result = kron(vec1_T, sum(RightM*Vec_mat, 1));
  return result;
}
double f_trOmegainverseDO1betastar_l(int l, int N, int T, arma::vec betastar, arma::mat Xawith1) {
  // f_trOmegainverseDO1betastar_l() - tr[\Omega^{-1} \frac{\partial}{\partial\beta^{\star}_{l}}\Omega], the matrix = J_{T} \otimes [(I_{NT} + T A_{N})^{-1} \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})}]
  arma::vec diagVec(N); // diagonal of matrix [(I_{NT} + T A_{N})^{-1} \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})}]
  for (int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    double d1aii_l = 2.0 * Xawith1(i, l) * dot(Xawith1.row(i), betastar); // element (i,i) in matrix \diag{2 x_{ali} (x_{ai}^{\prime} \beta^{\star})}
    diagVec(i) = (1.0 / (1.0 + 1.0 * T * aii)) * d1aii_l;
  }
  double result = 1.0 * T * sum(diagVec); // original matrix = J_{T} \otimes diag{diagVec}
  return result;
}
arma::vec f_score(int N, int T, int p, int q, double sigmaepsilon2, arma::vec betastar, arma::vec Y, arma::vec V, arma::mat W, arma::mat B, arma::mat Ztilde, arma::mat Xawith1) {
  // f_score() - score function, that is \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta}
  arma::mat IplusTAinverse(N,N); IplusTAinverse.fill(0.0); // (I_{N} + T A_{N})^{-1}
  for(int i = 0; i < N; i++) {
    double aii = dot(Xawith1.row(i), betastar) * dot(Xawith1.row(i), betastar); // A_{N}(\beta^\star) = \diag{(x_{a1}^{\prime} \beta^{\star})^{2}, ..., (x_{aN}^{\prime} \beta^{\star})^{2}}
    IplusTAinverse(i,i) = 1.0 / (1.0 + 1.0 * T * aii);
  }
  arma::vec OmegainverseV = f1_OmegainverseVec(V, N, T, IplusTAinverse); // \Omega^{-1} V
  arma::vec IotimesWY = f1_IotimesWVec(Y, N, T, W); // (I_{T} \otimes W_{N}) Y
  arma::vec D1_phi = 1.0/sigmaepsilon2 * Ztilde.t() * OmegainverseV;
  double D1_sigma = 1.0/(2.0*sigmaepsilon2*sigmaepsilon2) * dot(V, OmegainverseV) - 1.0*N*T/(2.0*sigmaepsilon2);
  double D1_lambda = - 1.0*T * trace(B.i()*W) + 1.0/sigmaepsilon2 * dot(IotimesWY, OmegainverseV);
  arma::vec D1_betastar(1+p);
  for (int l = 0; l < (1+p); l++) {
    double trOmegainverseDO1betastar_l = f_trOmegainverseDO1betastar_l(l, N, T, betastar, Xawith1); // tr[\Omega^{-1} \frac{\partial}{\partial \beta^{\star}_{l}}\Omega]
    arma::vec DO1betastarOmegainverseV = f1_DO1betastarOmegainverseVec_l(l, V, N, T, betastar, Xawith1); // \frac{\partial}{\partial \beta^{\star}_{l}}\Omega \Omega^{-1} V
    D1_betastar(l) = - 1.0/2.0 * trOmegainverseDO1betastar_l + 1.0/(2.0*sigmaepsilon2) * dot(OmegainverseV, DO1betastarOmegainverseV);
  }
  arma::mat D1(2*p+q+5,1);
  D1.submat(0,0,p+q+1,0) = D1_phi;
  D1(p+q+2,0) = D1_sigma;
  D1(p+q+3,0) = D1_lambda;
  D1.submat(p+q+4,0,2*p+q+4,0) = D1_betastar;
  // score function, that is \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta}
  arma::vec score = D1.col(0);
  return score;
}
// [[Rcpp::export]]
arma::vec f_score_b(arma::vec indicatorb, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat, arma::vec Vhat) {
  // f_score_b() - the estimation of score function in b-th bootstrap
  arma::vec Ylag1 = Ztilde.col(0);
  arma::vec y0 = Ylag1.head(N); // y_{0}
  arma::mat Z = Ztilde.cols(1, q);
  arma::mat Xwith1 = Ztilde.submat(0, q+1, N-1, (Ztilde.n_cols-1));
  arma::vec phihat = zetahat.head(p+q+2);
  double alphahat = phihat(0);
  arma::vec gammahat = phihat.subvec(1, q);
  arma::vec psihat = phihat.tail(1+p);
  double sigmaepsilon2hat = zetahat(p+q+2);
  double lambdahat = zetahat(p+q+3);
  arma::vec betastarhat = zetahat.tail(p+1);
  // (1) Matrix (\hat{v}_{1}, ..., \hat{v}_{T}) and matrix (\hat{r}_{1}, ..., \hat{r}_{T}), where \hat{v}_{t} is the estimation of v_{t} and we sample the mathched boostrap samples from {\hat{r}_{1}, ..., \hat{r}_{T}}
  arma::mat vhat_mat = f1_asmat(Vhat, N, T); // the t-th column is \hat{v}_{t}, t = 1, ..., T, v_{it} = x_{i}^\prime \beta \eta_{i} + \varepsilon_{it}, t = 1, ..., T
  arma::mat Ahat(N,N); Ahat.fill(0.0); // A_{N}(\hat(\beta)^{\star}) = 1/\sigma^{2}_{\varepsilon} Cov(\vartheta)
  for(int i = 0; i < N; i++) {
    arma::rowvec xai = Xawith1.row(i);
    Ahat(i,i) = dot(xai, betastarhat) * dot(xai, betastarhat);
  }
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat AhatplusI = Ahat + I_N; // A_{N}(\hat(\beta)^{\star}) + I_{N} = 1/\sigma^{2}_{\varepsilon} Cov(\vartheta + \varepsilon_{t})
  arma::mat AhatplusI_sqrt = sqrt(AhatplusI); // (A_{N}(\hat(\beta)^{\star}) + I_{N})^{1/2} and A_{N}(\hat(\beta)^{\star}) + I_{N} is a diagonal matrix
  arma::mat rhat_mat(N, T); rhat_mat.fill(0.0); // the t-th column is \hat{r}_{t}, t = 1, ..., T, \hat{r}_{t} = (A_{N}(\hat(\beta)^{\star}) + I_{N})^{-1/2} \hat{v}_{t}
  for (int t = 0; t < T; t++) {
    arma::vec rhatt = AhatplusI_sqrt.i() * vhat_mat.col(t);
    rhat_mat.col(t) = rhatt - mean(rhatt);
  }
  // (2,3,4) Bootstrap for Bnum estimations of VC matrix, that is Sigma^{\star}_{N} + \Sigma_{N} = E(\frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta} \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta^{\prime}})
  arma::mat Bhat = I_N - lambdahat * W;
  // T matched booststrap samples {\hat{r}_{1}^{b}, ..., \hat{r}_{T}^{b}} in b-th bootstrap
  arma::mat rhatb_mat(N, T); rhatb_mat.fill(0.0);
  for (int i = 0; i < N; i++) {
    rhatb_mat.row(i) = rhat_mat.row(indicatorb(i)-1);
  }
  // calculate the \hat{V}^{b}_{NT}, \hat{Y}^{b}, \hat{Y}_{-1}^{b}
  arma::mat vhatb_mat(N, T); vhatb_mat.fill(0.0); // the t-th column is \hat{v}^{b}_{t} = (A_{N}(\hat(\beta)^{\star}) + I_{N}) * \hat{r}^{b}_{t}, t = 1, ..., T
  arma::mat yhatb_mat(N, T+1); yhatb_mat.fill(0.0); // the t-th column is \hat{y}^{b}_{t}
  yhatb_mat.col(0) = y0;
  for (int t = 0; t < T; t++) { // t = 1, ..., T
    vhatb_mat.col(t) = AhatplusI_sqrt * rhatb_mat.col(t); // \hat{v}^{b}_{t}
    yhatb_mat.col(t+1) = Bhat.i() * (alphahat * yhatb_mat.col(t) + Z.rows((N*t), (N*t+N-1)) * gammahat + Xwith1 * psihat + vhatb_mat.col(t)); // \hat{y}^{b}_{t}
  }
  arma::vec Vhatb = vectorise(vhatb_mat, 0);
  arma::vec Ylag1hatb = vectorise(yhatb_mat.cols(0, T-1), 0);
  arma::vec Yhatb = vectorise(yhatb_mat.cols(1, T), 0);
  arma::mat Ztildehatb = Ztilde; Ztildehatb.col(0) = Ylag1hatb;
  // calculate the esimation of score function
  arma::vec scoreb = f_score(N, T, p, q, sigmaepsilon2hat, betastarhat, Yhatb, Vhatb, W, Bhat, Ztildehatb, Xawith1);
  return scoreb;
}  
