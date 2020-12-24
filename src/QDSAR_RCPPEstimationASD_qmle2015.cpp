#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec f_score (int N, int T, int p, int q, double sigmaepsilon2, arma::vec betastar, arma::vec Y, arma::vec V, arma::mat I_T, arma::mat J_T, arma::mat W, arma::mat S, arma::mat Omega, arma::mat Ztilde, arma::mat Xawith1) {
  // f_score - score function, that is \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta}
  arma::vec D1_phi = 1.0/sigmaepsilon2 * Ztilde.t()*Omega.i()*V;
  double D1_sigma = 1.0/(2.0*sigmaepsilon2*sigmaepsilon2) * dot(V, Omega.i()*V) - N*T/(2.0*sigmaepsilon2);
  double D1_lambda = - trace(S.i()*kron(I_T, W)) + 1.0/sigmaepsilon2 * dot(Y, kron(I_T, W).t()*Omega.i()*V);
  arma::vec D1_betastar(1+p);
  for (int l = 0; l < (1+p); l++) {
    arma::mat D1lA(N,N); D1lA.fill(0.0);
    for (int i = 0; i < N; i++) {
      D1lA(i,i) = 2.0*Xawith1(i,l)*dot(Xawith1.row(i), betastar);
    }
    arma::mat D1lOmega = kron(J_T, D1lA);
    D1_betastar(l) = - 1.0/2.0 * trace(Omega.i()*D1lOmega) + 1.0/(2.0*sigmaepsilon2) * dot(V, Omega.i()*D1lOmega*Omega.i()*V);
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
arma::vec f_score_b (int b, arma::vec indicatorb, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat) {
  // f_score_b - the estimation of score function in b-th bootstrap
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::vec Ylag1 = Ztilde.col(0);
  arma::mat Z = Ztilde.cols(1, q);
  arma::mat Xwith1 = Ztilde.submat(0, q+1, N-1, (Ztilde.n_cols-1));
  arma::mat y(N, T+1); y.fill(0.0); // the t-th column is y_{t}, t = 0, ..., T
  for (int i = 0; i < N; i++) {
    for (int t = 0; t < T; t++) {
      y(i, t+1) = Y(N * t + i);
    }
  }
  y.col(0) = Ylag1.head(N);
  arma::vec phihat = zetahat.head(p+q+2);
  double alphahat = phihat(0);
  arma::vec gammahat = phihat.subvec(1, q);
  arma::vec psihat = phihat.tail(1+p);
  double sigmaepsilon2hat = zetahat(p+q+2);
  double lambdahat = zetahat(p+q+3);
  arma::vec betastarhat = zetahat.tail(p+1);
  // (1) Matrix (\hat{v}_{1}, ..., \hat{v}_{T}) and matrix (\hat{r}_{1}, ..., \hat{r}_{T}), where \hat{v}_{t} is the estimation of v_{t} and we sample the boostrap samples from {\hat{r}_{1}, ..., \hat{r}_{T}}
  arma::mat Ahat(N,N); Ahat.fill(0.0); // 1/\sigma^{2}_{\epsilon} Cov(\vartheta)
  for(int i = 0; i < N; i++) {
    arma::rowvec xai = Xawith1.row(i);
    Ahat(i,i) = dot(xai, betastarhat) * dot(xai, betastarhat);
  }
  arma::mat variance_matrix = Ahat + I_N; // (A_{N}(\hat(\beta)^{\star}) + I_{N})
  // arma::mat vhat(N, T); vhat.fill(0.0); the t-th column is \hat{v}_{t}, t = 1, ..., T, v_{it} = x_{i}^\prime \beta \eta_{i} + \epsilon_{it}, t = 1, ..., T
  arma::mat rhat(N, T); rhat.fill(0.0); // the t-th column is \hat{r}_{t}, t = 1, ..., T, \hat{r}_{t} = (A_{N}(\hat(\beta)^{\star}) + I_{N})^{-1} \hat{v}_{t}
  for (int t = 0; t < T; t++) {
    arma::vec yt = y.col(t+1);
    arma::vec ytminus1 = y.col(t);
    arma::mat Zt = Z.rows((N*t), (N*t+N-1));
    arma::vec vhatt = yt - alphahat * ytminus1 - lambdahat * W * yt - Zt * gammahat - Xwith1 * psihat;
    arma::vec rhatt = variance_matrix.i() * vhatt;
    // vhat.col(t) = vhatt;
    rhat.col(t) = rhatt - mean(rhatt);
  }
  // (2,3,4) Bootstrap for Bnum estimations of VC matrix, that is Sigma^{\star}_{N} + \Sigma_{N} = E(\frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta} \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta})
  arma::mat I_T(T,T); I_T.eye(T,T);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat rowvec1_T(1, T); rowvec1_T.fill(1.0);
  arma::mat J_T = vec1_T * rowvec1_T;
  arma::mat Bhat = I_N - lambdahat * W;
  arma::mat Shat = kron(I_T, Bhat);
  arma::mat I_NT(N*T,N*T); I_NT.eye(N*T,N*T);
  arma::mat Omegahat = kron(J_T, Ahat) + I_NT;
  // T matched booststrap samples {\hat{r}_{1}^{b}, ..., \hat{r}_{T}^{b}} in b-th bootstrap
  arma::mat rhatb(N, T); rhatb.fill(0.0);
  for (int i = 0; i < N; i++) {
    rhatb.row(i) = rhat.row(indicatorb(i)-1);
  }
  // calculate the \hat{V}^{b}_{NT}, \hat{Y}^{b}, \hat{Y}_{-1}^{b}
  arma::mat vhatb(N, T); vhatb.fill(0.0); // the t-th column is \hat{v}^{b}_{t} = (A_{N}(\hat(\beta)^{\star}) + I_{N}) * \hat{r}^{b}_{t}, t = 1, ..., T
  vhatb.col(0) = variance_matrix * rhatb.col(0); // \hat{v}^{b}_{1}
  arma::mat yhatb(N, T+1); yhatb.fill(0.0); // the t-th column is \hat{y}^{b}_{t}
  yhatb.col(0) = y.col(0); // y_{0}
  yhatb.col(1) = Bhat.i() * (alphahat * yhatb.col(0) + Z.rows(0, N-1) * gammahat + Xwith1 * psihat + vhatb.col(0)); // \hat{y}^{b}_{1}
  for (int t = 1; t < T; t++) { // t = 2, ..., T
    vhatb.col(t) = variance_matrix * rhatb.col(t); // \hat{v}^{b}_{t}
    yhatb.col(t+1) = Bhat.i() * (alphahat * yhatb.col(t) + Z.rows((N*t), (N*t+N-1)) * gammahat + Xwith1 * psihat + vhatb.col(t)); // \hat{y}^{b}_{t}
  }
  arma::vec Vhatb = vectorise(vhatb, 0);
  arma::vec Ylag1hatb = vectorise(yhatb.cols(0, T-1), 0);
  arma::vec Yhatb = vectorise(yhatb.cols(1, T), 0);
  arma::mat Ztildehatb = Ztilde; Ztildehatb.col(0) = Ylag1hatb;
  // calculate the esimation of score function
  arma::vec scoreb = f_score (N, T, p, q, sigmaepsilon2hat, betastarhat, Yhatb, Vhatb, I_T, J_T, W, Shat, Omegahat, Ztildehatb, Xawith1);
  return scoreb;
}  

arma::mat f_boots_VCmatrix (int Bnum, arma::mat scoreBnum, int dimzeta) {
  // f_boots_VCmatrix - in b-th bootstrap, firstly obtain T matched booststrap samples {\hat{r}_{1}^{b}, ..., \hat{r}_{T}^{b}}, then calculate the estimation of VC matrix, which follows (Su and Yang, 2015)
  arma::cube scorescoreTBnum(dimzeta, dimzeta, Bnum); scorescoreTBnum.fill(0.0); // the b-th slice is (scoreb %*% score^\prime), where scoreb is the estimation of score function in b-th bootstrap
  for (int b = 0; b < Bnum; b++) {
    arma::vec scoreb = scoreBnum.col(b);
    arma::mat scorebM(dimzeta, 1); scorebM.col(0) = scoreb;
    scorescoreTBnum.slice(b) = scorebM * scorebM.t();
  }
  // calculate the estimation of VC matrix, that is Sigma^{\star}_{N} + \Sigma_{N} = E(\frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta} \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta})
  arma::mat sec1 = mean(scorescoreTBnum, 2); // \frac{1}{b} \sum\limits_{b=1}^{Bnum} (scoreb %*% scoreb^{\prime})
  arma::mat sec2(dimzeta, 1); // \frac{1}{B} \sum\limits_{b=1}^{Bnum} scoreb
  sec2.col(0) = mean(scoreBnum, 1);
  arma::mat VCmatrixhat = sec1 - sec2 * sec2.t();
  return VCmatrixhat;
}

// [[Rcpp::export]]
arma::vec f_ASD_zeta (int Bnum, arma::mat scoreBnum, int N, int T, int p, int q, arma::mat W, arma::vec Y, arma::mat Ztilde, arma::mat Xawith1, arma::vec zetahat) {
  // f_ASD_zeta - estimation for ASD of \hat{\zeta} with the normal distribution assumption
  int dimzeta = zetahat.n_elem;
  arma::mat I_N(N,N); I_N.eye(N,N);
  arma::mat I_T(T,T); I_T.eye(T,T);
  arma::mat I_NT(N*T,N*T); I_NT.eye(N*T,N*T);
  arma::mat vec1_T(T, 1); vec1_T.fill(1.0);
  arma::mat rowvec1_T(1, T); rowvec1_T.fill(1.0);
  arma::mat J_T = vec1_T * rowvec1_T;
  arma::vec phihat = zetahat.head(p+q+2);
  double sigmaepsilon2hat = zetahat(p+q+2);
  double lambdahat = zetahat(p+q+3);
  arma::vec betastarhat = zetahat.tail(p+1);
  arma::mat Bhat = I_N - lambdahat * W;
  arma::mat Shat = kron(I_T, Bhat);
  arma::vec Vhat = Shat*Y - Ztilde*phihat;
  arma::mat Ahat(N,N); Ahat.fill(0.0); // 1/\sigma^{2}_{\epsilon} Cov(\vartheta)
  for(int i = 0; i < N; i++) {
    Ahat(i,i) = dot(Xawith1.row(i), betastarhat) * dot(Xawith1.row(i), betastarhat);
  }
  arma::mat Omegahat = kron(J_T, Ahat) + I_NT;
  // Second partial derivatives
  arma::mat D2_phiphi = - 1.0/sigmaepsilon2hat * Ztilde.t()*Omegahat.i()*Ztilde;
  arma::mat D2_sigmaphi = - 1.0/(sigmaepsilon2hat*sigmaepsilon2hat) * (Ztilde.t()*Omegahat.i()*Vhat).t();
  double D2_sigmasigma = - 1.0/(sigmaepsilon2hat*sigmaepsilon2hat*sigmaepsilon2hat) * dot(Vhat, Omegahat.i()*Vhat) + N*T/(2.0*sigmaepsilon2hat*sigmaepsilon2hat);
  arma::mat D2_lambdaphi = - 1.0/sigmaepsilon2hat * (Ztilde.t()*Omegahat.i()*kron(I_T,W)*Y).t();
  double D2_lambdasigma = - 1.0/(sigmaepsilon2hat*sigmaepsilon2hat) * dot(Y, kron(I_T,W).t()*Omegahat.i()*Vhat);
  double D2_lambdalambda = - trace(Shat.i()*kron(I_T, W)*Shat.i()*kron(I_T, W)) - 1.0/sigmaepsilon2hat * dot(Y, kron(I_T, W).t()*Omegahat.i()*kron(I_T, W)*Y);
  arma::mat D2_betastarphi(1+p,p+q+2); D2_betastarphi.fill(0.0);
  arma::mat D2_betastarsigma(1+p,1); D2_betastarsigma.fill(0.0);
  arma::mat D2_betastarlambda(1+p,1); D2_betastarlambda.fill(0.0);
  for (int l = 0; l < (1+p); l++) {
    arma::mat D1lA(N,N); D1lA.fill(0.0);
    for (int i = 0; i < N; i++) {
      D1lA(i,i) = 2.0*Xawith1(i,l)*dot(Xawith1.row(i), betastarhat);
    }
    arma::mat D1lOmega = kron(J_T, D1lA);
    D2_betastarphi.row(l) = - 1.0/sigmaepsilon2hat * (Ztilde.t()*Omegahat.i()*D1lOmega*Omegahat.i()*Vhat).t();
    D2_betastarsigma(l,0) = - 1.0/(2.0*sigmaepsilon2hat*sigmaepsilon2hat) * dot(Vhat, Omegahat.i()*D1lOmega*Omegahat.i()*Vhat);
    D2_betastarlambda(l,0) = - 1.0/sigmaepsilon2hat * dot(Y, kron(I_T, W).t()*Omegahat.i()*D1lOmega*Omegahat.i()*Vhat);
  }
  arma::mat D2_betastarbetastar(1+p,1+p); D2_betastarbetastar.fill(0.0);
  for (int l1 = 0; l1 < (1+p); l1++) {
    for (int l2 = 0; l2 < (1+p); l2++) {
      arma::mat D1l1A(N,N); D1l1A.fill(0.0); arma::mat D1l2A(N,N); D1l2A.fill(0.0);
      arma::mat D2l1l2A(N,N); D2l1l2A.fill(0.0);
      for (int i = 0; i < N; i++) {
        D1l1A(i,i) = 2.0*Xawith1(i,l1)*dot(Xawith1.row(i), betastarhat);
        D1l2A(i,i) = 2.0*Xawith1(i,l2)*dot(Xawith1.row(i), betastarhat);
        D2l1l2A(i,i) = 2.0*Xawith1(i,l1)*Xawith1(i,l2);
      }
      arma::mat D1l1Omega = kron(J_T, D1l1A); arma::mat D1l2Omega = kron(J_T, D1l2A);
      arma::mat D2l1l2Omega = kron(J_T, D2l1l2A);
      D2_betastarbetastar(l1,l2) = - 1.0/2.0 * trace(-Omegahat.i()*D1l1Omega*Omegahat.i()*D1l2Omega + Omegahat.i()*D2l1l2Omega) + 1.0/(2.0*sigmaepsilon2hat) * (-2*dot(Vhat, Omegahat.i()*D1l1Omega*Omegahat.i()*D1l2Omega*Omegahat.i()*Vhat) + dot(Vhat, Omegahat.i()*D2l1l2Omega*Omegahat.i()*Vhat));
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
  // \hat{\Sigma}_{N} + \hat{\Sigma}^{\star}_{N} = \frac{1}{NT} VCmatrixhat, where \Sigma_{N} + Sigma^{\star}_{N} = \frac{1}{NT} E(\frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta} \frac{\partial \lnL_{NT}(\zeta_{0})}{\partial zeta})
  arma::mat VCmatrixhat = f_boots_VCmatrix(Bnum, scoreBnum, dimzeta);
  arma::mat SigmaPlusSigmastarhat = VCmatrixhat / (N*T*1.0);
  // estimation of the covariance matrix
  arma::mat CovM = Sigmahat.i() * SigmaPlusSigmastarhat * Sigmahat.i() / (N*T*1.0);
  // ASD
  arma::vec ASDzetahat(dimzeta);
  for (int l = 0; l < dimzeta; l++) {
    ASDzetahat(l) = sqrt(CovM(l,l));
  }
  return ASDzetahat;
}