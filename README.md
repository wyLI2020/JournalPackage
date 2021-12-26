# QuantileDSAR

R package "QuantileDSAR" is to implement the hybrid two-stage estimation procedure for the DSAR model with random effect.

## Installation

```R
#install.packages("devtools")
library(devtools)
install_github("wyLI2020/QuantileDSAR")
```

## Usage

```R
Rf_ZetaVarphiBetaASD_hat(tauseq, N, T, W, Y, Ylag1, Z, X, K, tauK, Bnum)
```

- **tauseq**:  vector, the quantiles to be estimated
- **N**:  integer, the number of spatial units
- **T**:  integer, the number of time periods
- **W**:  (N, N) matrix, the spatial weights matrix
- **Y**:  NT​-dimensional vector, response
- **Ylag1**:  NT​-dimensional vector, the first lag of response
- **Z**:  (NT, q​) matrix, time-varying regressors
- **X**:  (N, p) matrix, time-invariant regressors
- **K**:  integer, the number of quantile levels for $\widetilde{\beta}_{c}$ in WCQE and WQAE
- **tauK**:  K​-dimensional vector, quantile levels for $\widetilde{\beta}_{c}$ in WCQE and WQAE
- **Bnum**:  integer, the number of bootstrap samples

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

set.seed(6)
Rf_ZetaVarphiBetaASD_hat(tauseq=c(0.25,0.75), N, T, W, Y, Ylag1, Z, X, K=9, tauK=c(1:9)/10, Bnum=500)
```



