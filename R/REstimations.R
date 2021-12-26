# library(Rcpp)
# library(RcppArmadillo)
# library(quantreg)
# library(nloptr)
# Rcpp::sourceCpp('RCPPEstimationQMLE_1.cpp')
Rf_score_b <- function(b, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat) {
  indicatorb <- sample(c(1:N), N, replace = TRUE)
  scoreb <- f_score_b(indicatorb, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat)
  return(scoreb)
}
# Rcpp::sourceCpp('RCPPEstimationQMLE_2.cpp')
RlnlOFdelta_Opposite <- function(delta, N, T, W, Y, Ztilde, Xawith1) {
  # RlnlOFdelta_Opposite() ,i.e., -\ln l(delta)
  lnldelta <- RCPPlnlOFdelta(delta, N, T, W, Y, Ztilde, Xawith1)
  lnldelta_Opposite <- -lnldelta
  return(lnldelta_Opposite)
}
# Rcpp::sourceCpp('RCPPEstimationVarthetaHat.cpp')
# Rcpp::sourceCpp('RCPPEstimationWCQE_1.cpp')
Rf_WCQE_withc <- function(p, N, y, xwith1, tau, weights) { # obtain WCQE with constrains 
  fn <- function(varphi) { # objective function
    return(wqrLF(varphi, N, y, xwith1, tau, weights))
  }
  gr <- function(varphi) { # gradient of the objective function
    return(gr_wqrLF(varphi, N, y, xwith1, tau, weights))
  }
  hin <- function(varphi) { # the inequalty constraints, hin(varphi) >= 0
    if (p==1) { # hin(varphi), which is applied to the case p=1
      return(c(varphi[1] * varphi[2]))
    } else if (p==2) { # hin(varphi), which is applied to the case p=2
      return(c(varphi[1] * varphi[2], varphi[1] * varphi[3]))
    }
  }
  varphi_ini <- as.matrix(rep(0,p+1))
  return(auglag(varphi_ini, fn, hin = hin, localsolver = c("COBYLA"))$par)
}
Rf_WCQE <- function(p, N, y, x, xwith1, tau, weights) {
  wcqe <- rq(y~x, tau = tau, weights = weights)$coefficients
  if(all(wcqe >= 0) | all(wcqe <= 0)) {
    result <- wcqe
  } else {
    wcqe_c <- Rf_WCQE_withc(p, N, y, xwith1, tau, weights)
    result <- wcqe_c
  }
  return(result)
}
# Rcpp::sourceCpp('RCPPEstimationWCQE_2.cpp')
# Rcpp::sourceCpp('RCPPEstimationWQAE.cpp')

Rf_ZetaVarphiBetaASD_hat <- function(tauseq, N, T, W, Y, Ylag1, Z, X, K, tauK, Bnum) {
  # prepare of matrix
  q <- ncol(Z)
  p <- ncol(X)
  Xa <- abs(X)
  Xwith1 <- cbind(rep(1,N), X)
  Xawith1 <- cbind(rep(1,N), Xa)
  ZtildeSEC1 <- cbind(Ylag1, Z)
  ZtildeSEC2 <- kronecker(rep(1,T), Xwith1)
  Ztilde <- cbind(ZtildeSEC1, ZtildeSEC2)
  # QMLE \hat{\zeta}
  deltaInitial <- c(0.5,rep(1,1+p)); deltaLower <- c(-0.999999,rep(0.000001,1+p)); deltaUpper <- c(0.999999,rep(Inf,1+p)) # initial value of \delta in optim()
  deltahat <- optim(par = deltaInitial, fn = RlnlOFdelta_Opposite, lower = deltaLower, upper = deltaUpper, method = "L-BFGS-B", N=N, T=T, W=W, Y=Y, Ztilde=Ztilde, Xawith1=Xawith1)$par
  zetahat <- f_zeta_hat(N, T, W, Y, Ztilde, Xawith1, deltahat) # QMLE
  ## ASD of QMLE
  Vhat <- f_V_hat(zetahat, N, T, W, Y, Ztilde) # \hat{V}_{NT}
  CovMzetahat <- f_CovM_QMLE(Bnum, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat, Vhat)
  ASDzetahat <- sqrt(diag(CovMzetahat))
  # \check{\beta}
  dimzeta <- length(zetahat)
  betastarhat <- zetahat[(dimzeta-p):dimzeta]
  betacheck <- betastarhat / betastarhat[1] # \check{\beta}
  ## ASD of \check{\beta}
  CovMsecbetacheck <- f_CovM_secbetacheck(p, q, zetahat, CovMzetahat)
  ASDsecbetacheck <- sqrt(diag(CovMsecbetacheck))
  ASDbetacheck <- c(NA, ASDsecbetacheck)
  # \hat{\vartheta}
  varthetahat <- f_vartheta_hat(Vhat, N, T)
  varthetahat <- varthetahat - mean(varthetahat) # zero mean \hat{\vartheta}
  # WCQE 
  ## \tilde{\beta}_{c}
  varphitauKtildeM <- f_WCQE_tausM(p, N, varthetahat, Xa, Xawith1, K, tauK, rep(1,N)) # unweighted conditional quantile estimation, \tilde{\varphi}(\tau_{k}), k = 1,...,K
  betactilde <- f_betac_tilde(varphitauKtildeM) # \tilde{\beta}_{c}
  ## \hat{\varphi}(\tau) with \tau's in tauseq
  weights_rq <- f_weights_rq(N, Xawith1, betactilde) # weights in WCQE
  varphitauseqhatM <- f_WCQE_tausM(p, N, varthetahat, Xa, Xawith1, length(tauseq), tauseq, weights_rq) # WCQE
  ## ASD of WCQE
  etatilde <- f_eta_est(N, Xawith1, varthetahat, betactilde) # \tilde{\eta}
  QetatauseqbarV <- quantile(etatilde, probs = tauseq) # \bar{Q}_{\eta}(\tau) with \tau's in tauseq
  XitauseqbarV <- f_feta_tildeV(QetatauseqbarV, N, Xawith1, varthetahat, betactilde) # \bar{\Xi}(\tau) with \tau's in tauseq, \bar{\Xi}(\tau) = \tilde{f}_{\eta}(\bar{Q}_{\eta}(\tau))
  Dtilde <- f_D_est(N, Xawith1, betactilde) # \tilde_{D}_{N} = \frac{1}{N} \sum_{i=1}^{N} (x_{ai} x_{ai}^{\prime}) / (x_{ai}^\prime \tilde{\beta}_{c})^{2}
  ASDvarphitauseqhatM <- f_ASD_WCQE(tauseq, XitauseqbarV, Dtilde, N, p)
  # WQAE 
  varphitauKhatM <- f_WCQE_tausM(p, N, varthetahat, Xa, Xawith1, K, tauK, weights_rq) # WCQE \hat{\varphi}(\tau_{k}), k = 1,...,K
  QetatauKhatV <- varphitauKhatM[1,] # \hat{Q}_{\eta}(\tau) with \tau's in tauK
  qhat <- QetatauKhatV # \hat{q}
  XitauKhatV <- f_feta_tildeV(QetatauKhatV, N, Xawith1, varthetahat, betactilde) # \hat{\Xi}(\tau) with \tau's in tauK, \hat{\Xi}(\tau) = \tilde{f}_{\eta}(\hat{Q}_{\eta}(\tau))
  Hhat <- f_H_est(K, tauK, XitauKhatV) # \hat{H}
  betahat <- f_WQAE(varphitauKhatM, qhat, Hhat) # WQAE \hat{\beta}(\hat{\pi}_{opt,K})
  ## ASD of WQAE
  ASDbetahat <- f_ASD_WQAE(N, qhat, Hhat, Dtilde)
  # vector of the estimators
  result_est <- c(zetahat, betacheck[-1], betactilde[-1], betahat[-1], as.vector(varphitauseqhatM))
  result_asd <- c(ASDzetahat, ASDbetacheck[-1], rep(NA, p), ASDbetahat[-1], as.vector(ASDvarphitauseqhatM))
  result <- c(result_est, result_asd)
  namealphahat <- "alphahat"; namegammahat <- paste0("gamma", 1:q, "hat"); namepsihat <- paste0("psi", 0:p, "hat")
  namephihat <- c(namealphahat, namegammahat, namepsihat); namesigmavarepsilon2hat <- "sigmavarepsilon2hat"
  namelambdahat <- "lambdahat"; namebetastarhat <- paste0("betastar", 0:p, "hat")
  namezetahat <- c(namephihat, namesigmavarepsilon2hat, namelambdahat, namebetastarhat)
  namebetacheck <- paste0("beta", 0:p, "check")
  namebetactilde <- paste0("beta", 0:p, "ctilde")
  namebetahat <- paste0("beta", 0:p, "hat")
  namevarphitauseqhatV <- paste0("varphi", rep((0:p), length(tauseq)), "(", rep(tauseq, each=1+p), ")hat")
  nameASDzetahat <- paste0("ASD", namezetahat)
  nameASDbetacheck <- paste0("ASD", namebetacheck)
  nameASDbetactilde <- paste0("ASD", namebetactilde)
  nameASDbetahat <- paste0("ASD", namebetahat)
  nameASDvarphitauseqhatV <- paste0("ASD", namevarphitauseqhatV)
  name_est <- c(namezetahat, namebetacheck[-1], namebetactilde[-1], namebetahat[-1], namevarphitauseqhatV)
  name_asd <- c(nameASDzetahat, nameASDbetacheck[-1], nameASDbetactilde[-1], nameASDbetahat[-1], nameASDvarphitauseqhatV)
  names(result) <- c(name_est, name_asd)
  return(result)
}
