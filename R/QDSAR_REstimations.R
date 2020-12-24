RlnlOFdelta_Opposite <- function(delta, N, T, W, Y, Ztilde, Xawith1) {
  # RlnlOFdelta_Opposite() ,i.e., -\ln l(delta)
  lnldelta <- RCPPlnlOFdelta(delta, N, T, W, Y, Ztilde, Xawith1)
  lnldelta_Opposite <- -lnldelta
  return(lnldelta_Opposite)
}

Rf_ZetaVarphiBeta_hat <- function(tauseq = c(0.25, 0.75), K, N, T, W, Y, Ylag1, Z, X, Bnum, clsnum) {
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
  deltaInitial <- c(0.1,rep(1,1+p)); deltaLower <- c(-1,rep(0,1+p)); deltaUpper <- c(1,rep(Inf,1+p))
  deltahat <- optim(par = deltaInitial, fn = RlnlOFdelta_Opposite, lower = deltaLower, upper = deltaUpper, method = "L-BFGS-B", N=N, T=T, W=W, Y=Y, Ztilde=Ztilde, Xawith1=Xawith1)$par
  zetahat <- f_zeta_hat(N, T, W, Y, Ztilde, Xawith1, deltahat) # QMLE
  namealpha <- "alpha"; namegamma <- paste0("gamma", 1:q); namepsi <- paste0("psi", 0:p)
  namephi <- c(namealpha, namegamma, namepsi); namesigmavarepsilon2 <- "sigmavarepsilon2"
  namelambda <- "lambda"; namebetastar <- paste0("betastar", 0:p)
  namezeta <- c(namephi, namesigmavarepsilon2, namelambda, namebetastar)
  rownames(zetahat) <- namezeta; colnames(zetahat) <- "QMLE"
  # \hat{\vartheta}
  varthetahat <- f_vartheta_hat(zetahat, N, T, W, Y, Ztilde) 
  # WCQE \hat{\varphi}(\tau_{j}) with \tau's in tauseq
  tauK <- round(c(1:K)/(K+1), 2) # (\tau_{1}, ..., \tau_{K})^\prime
  varphitauKtildeM <- rq(varthetahat~Xa, tau = tauK)$coefficients # unweighted conditional quantile estimation, \tilde{\varphi}(\tau_{k}), k = 1,...,K
  betactilde <- f_betac_tilde(varphitauKtildeM) # \tilde{\beta}_{c}
  weights_rq <- f_weights_rq(N, Xawith1, betactilde) # weights in WCQE
  varphitauhatM <- matrix(nrow = 1+p, ncol = 2)
  initialQ1 <- rep(-0.1, 1+p); upQ1 <- rep(0, 1+p)
  initialQ3 <- rep(0.1, 1+p); lowQ3 <- rep(0, 1+p)
  varphitauhatM[,1] <- optim(par = initialQ1, fn = wqrLF, upper = upQ1, method = "L-BFGS-B", N=N, y=varthetahat, x=Xawith1, tau=0.25, weights=weights_rq)$par
  varphitauhatM[,2] <- optim(par = initialQ3, fn = wqrLF, lower = lowQ3, method = "L-BFGS-B", N=N, y=varthetahat, x=Xawith1, tau=0.75, weights=weights_rq)$par
  rownames(varphitauhatM) <- paste0("varphi", 0:p, "(tau)"); colnames(varphitauhatM) <- paste0("tau=", tauseq)
  # WQAE \hat{\beta}(\hat{\pi}_{opt,K})
  varphitauKhatM <- rq(varthetahat~Xa, tau = tauK, weights = weights_rq)$coefficients # WCQE \hat{\varphi}(\tau_{k}), k = 1,...,K
  betahat <- f_beta_hat(tauK, N, Xawith1, varthetahat, betactilde, varphitauKhatM) # WQAE
  rownames(betahat) <- paste0("beta", 0:p); colnames(betahat) <- "WQAE"
  # ASD
  Dtilde <- f_D_est(N, Xawith1, betactilde)
  ASDvarphitauhatM <- f_ASD_varphi(tauseq, N, Xawith1, varthetahat, betactilde, varphitauhatM, Dtilde)
  rownames(ASDvarphitauhatM) <- paste0("varphi", 0:p, "(tau)"); colnames(ASDvarphitauhatM) <- paste0("tau=", tauseq)
  Dhat <- f_D_est(N, Xawith1, betahat)
  ASDbetahat <- f_ASD_beta(tauK, N, T, Xawith1, varthetahat, betahat, varphitauKhatM, Dhat)
  rownames(ASDbetahat) <- paste0("beta", 0:p); colnames(ASDbetahat) <- "ASD"
  ## ASD of QMLE
  set.seed(12); indicatorBnum <- matrix(sample(c(1:N), N*Bnum, replace = TRUE), N, Bnum)
  SEQ <- 1:Bnum
  parFun <- function(b){ # the function in parallel
    indicatorb <- indicatorBnum[,b]
    onerep <- f_score_b(b, indicatorb, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat)
    return(onerep)
  }
  cls <- makeCluster(clsnum, type="FORK")
  scoreBnum <- parSapply(cls, SEQ, parFun) # Carry out the task parFUNCTION parallely
  stopCluster(cls) # stop workers
  ASDzetahat <- f_ASD_zeta(Bnum, scoreBnum, N, T, p, q, W, Y, Ztilde, Xawith1, zetahat)
  rownames(ASDzetahat) <- namezeta; colnames(ASDzetahat) <- "ASD"
  # # \hat{\eta}
  # etahat <- varthetahat / (Xawith1%*%betahat)
  # etatilde <- varthetahat / (Xawith1%*%betactilde)
  # List of the estimators
  result <- list(zetahat=zetahat, betactilde=betactilde, varphitauhat=varphitauhatM, betahat=betahat, ASD_zeta=ASDzetahat, ASD_varphitau=ASDvarphitauhatM, ASDbeta=ASDbetahat)
  # result <- list(zetahat=zetahat, betactilde=betactilde, varphitauhat=varphitauhatM, betahat=betahat, ASD_zeta=ASDzetahat, ASD_varphitau=ASDvarphitauhatM, ASDbeta=ASDbetahat, etahat=etahat, etatilde=etatilde)
  return(result)
}
