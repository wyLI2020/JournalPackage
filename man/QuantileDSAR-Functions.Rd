\name{QuantileDSAR-Functions}
\alias{Rf_ZetaVarphiBetaASD_hat}
\alias{Rf_score_b}
\alias{RlnlOFdelta_Opposite}
\alias{Rf_WCQE_withc}
\alias{Rf_WCQE}
\alias{RCPPlnlOFdelta}
\alias{f_zeta_hat}
\alias{f_score_b}
\alias{f_CovM_QMLE}
\alias{f_CovM_secbetacheck}
\alias{f_V_hat}
\alias{f_vartheta_hat}
\alias{wqrLF}
\alias{gr_wqrLF}
\alias{f_WCQE_tausM}
\alias{f_betac_tilde}
\alias{f_weights_rq}
\alias{f_eta_est}
\alias{f_feta_tildeV}
\alias{f_D_est}
\alias{f_ASD_WCQE}
\alias{f_H_est}
\alias{f_WQAE}
\alias{f_ASD_WQAE}

\title{The estimation function in QuantileDSAR package}
\description{
  To calculate the estimators QMLE, WCQE and WQAE, as well as their ASD in the DSAR model.
}
\usage{
  Rf_ZetaVarphiBetaASD_hat(tauseq, N, T, W, Y, Ylag1, Z, X, K, tauK, Bnum)
}
\arguments{
  \item{tauseq}{vector, the quantiles to be estimated, default = "c(0.25, 0.75)"}
  \item{N}{integer, the number of spatial units}
  \item{T}{integer, the number of time periods}
  \item{W}{(\emph{N}, \emph{N}) matrix, the spatial weights matrix}
  \item{Y}{\emph{NT}-dimensional vector, response}
  \item{Ylag1}{\emph{NT}-dimensional vector, the first lag of response}
  \item{Z}{(\emph{NT}, q) matrix, time-varying regressors}
  \item{X}{(\emph{N}, p) matrix, time-invariant regressors}
  \item{K}{integer, the number of quantile levels in WCQE and WQAE}
  \item{tauK}{\emph{K}-dimensional vector, quantile levels in WCQE and WQAE}
  \item{Bnum}{integer, the number of bootstrap samples}
}
\value{

  \code{Rf_ZetaVarphiBetaASD_hat} returns a vector.

}
\details{
  With the preproccess in following "Examples", we can get the estimators of the DSAR model by the code as follows: 
  
  \code{Rf_ZetaVarphiBetaASD_hat(tauseq=c(0.25,0.75), N, T, W, Y, Ylag1, Z, X, K=9, tauK=c(1:9)/10, Bnum=500)}
}
\examples{
  data(W)
  data(AQIweather_df)
  data(Economy_df)
  N <- nrow(Economy_df)
  NTplus1 <- nrow(AQIweather_df)
  T <- NTplus1 / N - 1
  AQI <- AQIweather_df$AQI; nAQI <- AQI / sd(AQI)
  TEM <- AQIweather_df$TEMweek[-c(1:N)]; nTEM <- TEM / sd(TEM)
  PRE <- AQIweather_df$PREweek[-c(1:N)]; nPRE <- PRE / sd(PRE)
  WIN <- AQIweather_df$WINweek[-c(1:N)]; nWIN <- WIN / sd(WIN)
  GRP <- Economy_df$GRP; nGRP <- GRP / sd(GRP)
  Industry2nd <- Economy_df$Industry2nd; nIndustry2nd <- Industry2nd / sd(Industry2nd)
  Y <- nAQI[(N+1):NTplus1]
  Ylag1 <- nAQI[1:(NTplus1-N)]
  Z <- matrix(NA, nrow = N*T, ncol = 3)
  Z[,1] <- nTEM - mean(nTEM); Z[,2] <- nPRE - mean(nPRE); Z[,3] <- nWIN - mean(nWIN)
  X <- matrix(NA, nrow = N, ncol = 2)
  X[,1] <- nGRP; X[,2] <- nIndustry2nd
  set.seed(6)
  Rf_ZetaVarphiBetaASD_hat(tauseq=c(0.25,0.75), N, T, W, Y, Ylag1, Z, X, K=9, tauK=c(1:9)/10, Bnum=500)
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

