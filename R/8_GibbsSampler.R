#' @description Function to draw the coefficients for the AR(1) Probit model governing regime
#' switching probabilities from a inverse gamma distribution
#' @param regimeVec vector with ones and zeros of the regime draw
#' @return vector with the draws for the coefficients

Beta_fctn <- function(regimeVec) {
  # Set up the priors so that that p = q = .975 (maybe change prior variance)
  meanVec_0 <- c(-qnorm(.975), -qnorm(.025) + qnorm(.975))
  sdMat_0_inv <- Inverse(diag(rep(.5, 2)))
  X <- cbind(rep(1, length(regimeVec) - 1), lag(regimeVec, n = 1)[-1])
  XtX <- Transp(X) %*% X
  Y <- head(regimeVec, -1)
  meanVec_1 <- Inverse(sdMat_0_inv + XtX) %*% (sdMat_0_inv %*% meanVec_0 + Transp(X) %*% Y)
  sdMat_1 <- Inverse(sdMat_0_inv + XtX)
  draw <- rmvnorm(1, mean = meanVec_1, sigma = sdMat_1)
  return(c("beta_0" = draw[, 1], "beta_1" = draw[, 2]))
}


#' @description Function to draw the innovation standard deviation of the measurement equation
#' @param data vector with the noisy measurement
#' @param stateVec matrix of the state vector draw
#' @param regimeVec vector with ones and zeros of the regime draw
#' @param paramList list with the additional parameters
#' @param seasPos position of the cyclical component in the state vector
#' @return draw of the standard deviation

SdEpsilon_fctn <- function(data, stateVec, regimeVec = NULL, paramList, seasPos = 3) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  data <- data[-c(1:10)]
  # Set the priors
  alpha_0 <- 0
  beta_0 <- 0
  # Constrain the parameter space
  lowerLimit <- 1e-5
  upperLimit <- Inf
  # Compute the SSR from the regression of the data on the drift and seasonal term (maybe wider dist better)
  mu <- stateVec[1, ]
  gamma <- stateVec[seasPos, ]
  v_t <- data - mu - gamma
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + length(data))
  beta_1 <- beta_0 + SSR
  # Draw the new value
  epsilon_var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  epsilon_sd_candidate <- sqrt(epsilon_var)
  if (epsilon_sd_candidate < lowerLimit | epsilon_sd_candidate > upperLimit) {
    epsilon_sd_draw <- paramList$epsilon
  } else {
    epsilon_sd_draw <- epsilon_sd_candidate
  }
  return(epsilon_sd_draw)
}


#' @description Function to draw the AR covariates of an AR(2) cycle component
#' @param stateVec matrix of the state vector draw
#' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
#' @param paramList list with the additional parameters
#' @param cyclePos position of the cycle component in the state vector
#' @return draw of the standard deviation

ARPhi_fctn <- function(stateVec, regimeHeteroscVec, paramList, cyclePos) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  regimeHeteroscVec <- regimeHeteroscVec[-c(1:10)]
  # Set the priors
  mean_0 <- c(0, 0)
  sd_0_inv <- Inverse(diag(2))
  # Get necessary values
  phiBackup <- paramList$phi
  Dimens <- NROW(stateVec)
  nPeriods <- NCOL(stateVec)
  eta_0 <- paramList$eta$eta_0
  eta_1 <- paramList$eta$eta_1
  sdVec <- ifelse(regimeHeteroscVec == 1, eta_1, eta_0)
  cycle <- stateVec[cyclePos, -c(1:2)] / sdVec[-c(1:2)]
  cycle_tl1 <- stateVec[cyclePos, -c(1, nPeriods)] / sdVec[-c(1, nPeriods)]
  cycle_tl2 <- stateVec[cyclePos, -c(nPeriods - 1, nPeriods)] / sdVec[-c(nPeriods - 1, nPeriods)]
  Xt <- rbind(cycle_tl1, cycle_tl2)
  X <- Transp(Xt)
  XtX <- Xt %*% X
  XtY <- Xt %*% matrix(cycle, nc = 1)
  # Adjust the priors
  sd_1 <- Inverse(sd_0_inv + XtX)
  mean_1 <- sd_1 %*% (sd_0_inv %*% mean_0 + XtY)
  draw <- rmvnorm(1, mean = mean_1, sigma = sd_1)
  if (StatTest_fctn(draw) == F) draw <- phiBackup
  return(draw)
}


#' @description Function to checks if the coefficients of an AR(2) process are within the unit circle
#' @param paramVec vector with the two coefficients
#' @return TRUE if the coefficients are within the unit circle. Else FALSE

StatTest_fctn <- function(paramVec) {
  if (is.null(dim(paramVec))) paramVec <- matrix(paramVec, 1)
  k <- nrow(paramVec)
  p <- ncol(paramVec) / k
  paramVec <- rbind(paramVec, cbind(diag((p - 1) * k), matrix(0, (p - 1) * k, k)))
  eigv <- eigen(paramVec, only.values = TRUE)$values
  return(all(abs(eigv) < 1))
}


#' @description Function to draw the regime dependent sd of innovations governing the cyclical component
#' @param stateVec matrix of the state vector draw
#' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
#' @param paramList list with the additional parameters
#' @param cyclePos position of the cyclical component in the state vector
#' @return draw of the standard deviation

SdEtaSwitch_fctn <- function(stateVec, regimeHeteroscVec, paramList, cyclePos) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  regimeHeteroscVec <- regimeHeteroscVec[-c(1:10)]
  # Set the priors
  alpha_0 <- 0
  beta_0 <- 0
  alpha_2 <- 0
  beta_2 <- 0
  # Define variables
  nPeriods <- NCOL(stateVec)
  phiVec <- paramList$phi
  cycle <- stateVec[cyclePos, -1]
  cycleAR <- stateVec[cyclePos:(cyclePos + 1), -nPeriods]
  Heterosc_S_t <- regimeHeteroscVec[-1]
  # Get multiplicative specification of the two variances
  eta_0 <- paramList$eta$eta_0
  eta_1 <- paramList$eta$eta_1
  x <- (eta_1^2 / eta_0^2 - eta_0^2)
  correcFactor <- (1 + x * Heterosc_S_t)^(-.5)
  # Compute the SSR from the regression of the trend on its lag and drift terms, corrected for S_t = 1
  resid <- cycle - c(matrix(phiVec, nr = 1) %*% cycleAR)
  v_t_1 <- resid * correcFactor
  SSR_1 <- sum(v_t_1^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + nPeriods)
  beta_1 <- beta_0 + SSR_1
  # Draw the new value xi_0
  eta_0_var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  eta_0_draw <- sqrt(eta_0_var)
  # Compute the SSR from the regression of the trend on its lag and drift terms, corrected for (heterosc.) S_t = 0
  v_t_2 <- Heterosc_S_t * resid / eta_0_draw
  SSR_2 <- sum(v_t_2^2, na.rm = T)
  # Adjust the priors
  alpha_3 <- (alpha_2 + sum(Heterosc_S_t))
  beta_3 <- beta_2 + SSR_2
  # Draw the new value xi_0
  x_draw <- rinvgamma(1, alpha_3 / 2, beta_3 / 2)
  eta_1_draw <- sqrt(eta_0_var * x_draw)
  return(c(eta_0_draw, eta_1_draw))
}


#' @description Function to draw the regime independent sd of innovations governing the cyclical component
#' @param stateVec matrix of the state vector draw
#' @param paramList list with the additional parameters
#' @param cyclePos position of the cyclical component in the state vector
#' @return draw of the standard deviation

SdEtaConst_fctn <- function(stateVec, paramList, cyclePos) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  # Set the priors
  alpha_0 <- 0
  beta_0 <- 0
  # Constrain the parameter space
  lowerLimit <- 1e-5
  upperLimit <- Inf
  # Define variables
  nPeriods <- NCOL(stateVec)
  phiVec <- paramList$phi
  cycle <- stateVec[cyclePos, -1]
  cycleAR <- stateVec[cyclePos:(cyclePos + 1), -nPeriods]
  # Compute the SSR from the regression of the cyclical term on its lags
  v_t <- cycle - c(matrix(phiVec, nr = 1) %*% cycleAR)
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + nPeriods)
  beta_1 <- beta_0 + SSR
  # Draw the new value
  eta_var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  eta_sd_canditate <- sqrt(eta_var)
  if (eta_sd_canditate < lowerLimit | eta_sd_canditate > upperLimit) {
    eta_sd_draw <- paramList$eta
  } else {
    eta_sd_draw <- eta_sd_canditate
  }
  return(eta_sd_draw)
}

#' @description Function to draw the sd of innovations governing the cyclical component
#' @param stateVec matrix of the state vector draw
#' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
#' @param paramList list with the additional parameters
#' @param cyclePos position of the cyclical component in the state vector
#' @return draw of the standard deviation

SdEta_fctn <- function(stateVec, regimeHeteroscVec = NULL, paramList, cyclePos) {
  if (is.null(regimeHeteroscVec)) {
    eta_vec <- SdEtaConst_fctn(stateVec = stateVec, paramList = paramList, cyclePos = 3)
  } else {
    eta_vec <- SdEtaSwitch_fctn(stateVec = stateVec, regimeHeteroscVec = regimeHeteroscVec, paramList = paramList, cyclePos = 3)
  }
  return(eta_vec)
}


#' @description Function to draw the innovation standard deviation of the trend component
#' @param stateVec matrix of the state vector draw
#' @param regimeDriftVec vector with regime realizations governing the trend drift. Vector takes either values 1 or 0
#' @param regimeHeteroscVec vector with regime realizations governing the trend heteroscedasticity. Vector takes either values 1 or 0
#' if NULL, the regimeDriftVec is considered to specify the heteroscedasticity
#' @param paramList list with the additional parameters
#' @return vector of the drawn standard deviation(s)

SdXi_fctn <- function(stateVec, regimeDriftVec, regimeHeteroscVec = NULL, paramList) {
  if (is.null(regimeHeteroscVec)) {
    xiVec <- SdXiConst_fctn(stateVec = stateVec, regimeVec = regimeDriftVec, paramList = paramList)
  } else {
    xiVec <- SdXiSwitch_fctn(stateVec = stateVec, regimeDriftVec = regimeDriftVec, regimeHeteroscVec = regimeHeteroscVec, paramList = paramList)
  }
  return(xiVec)
}
