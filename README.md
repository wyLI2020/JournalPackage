# QuantileDSAR

R package "QuantileDSAR" is to implement the hybrid two-stage estimation procedure for the DSAR model with random effect.

## Installation

```R
#install.packages("Rcpp")
#install.packages("RcppArmadillo")
#install.packages("quantreg")
#install.packages("parallel")
library(Rcpp)
library(RcppArmadillo)
library(quantreg)
library(parallel)
system("R CMD INSTALL --build QuantileDSAR_1.0.tar.gz")
library(QuantileDSAR)
```

## Usage

```R
Rf_ZetaVarphiBeta_hat(tauseq = c(0.25, 0.75), K, N, T, W, Y, Ylag1, Z, X, Bnum, clsnum)
```

- **tauseq**:  vector, the quantiles to be estimated, default = "c(0.25, 0.75)"

- **K**:  integer, the number of quantile levels for $\widetilde{\beta}_{c}$ in WCQE
- **N**:  integer, the number of spatial units
- **T**:  integer, the number of time periods
- **W**:  $N \times N$ matrix, the spatial weights matrix
- **Y**:  $NT$-dimensional vector, response
- **Ylag1**:  $NT$-dimensional vector, the first lag of response
- **Z**:  $NT \times q$ matrix, time-varying regressors
- **X**:  $N \times p$ matrix, time-invariant regressors
- **Bnum**:  integer, the number of bootstrap samples
- **clsnum**:  the number of copies of R running in parallel for bootstrap

## Example

```R
library(QuantileDSAR)
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
Rf_ZetaVarphiBeta_hat(tauseq = c(0.25, 0.75), K=9, N, T, W, Y, Ylag1, Z, X, Bnum=500, clsnum=6)
```

## Development

This R package is developed by Wenyu Li. 

