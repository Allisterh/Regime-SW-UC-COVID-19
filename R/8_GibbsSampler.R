#' @description Function to draw the coefficients for the AR(1) Probit model governing regime
#' switching probabilities from a inverse gamma distribution
#' @param regimeVec vector with ones and zeros of the regime draw
#' @return vector with the draws for the coefficients

Beta_fctn <- function(regimeVec) {
  # Set up the priors so that that p = q = .975 (maybe change prior variance)
  meanVec_0 <- c(-qnorm(.975), -qnorm(.025) + qnorm(.975))
  sdMat_0 <- Inverse(diag(rep(.5, 2)))
  X <- cbind(rep(1, length(regimeVec) - 1), lag(regimeVec, n = 1)[-1])
  XtX <- Transp(X) %*% X
  Y <- head(regimeVec, -1)
  meanVec_1 <- Inverse(sdMat_0 + XtX) %*% (sdMat_0 %*% meanVec_0 + Transp(X) %*% Y)
  sdMat_1 <- Inverse(sdMat_0 + XtX)
  draw <- rmvnorm(1, mean = meanVec_1, sigma = sdMat_1)
  return(c("beta_0" = draw[, 1], "beta_1" = draw[, 2]))
}


#' @description Function to draw the innovation standard deviation of the measurement equation
#' @param data vector with the observations
#' @param stateVec matrix of the state vector draw
#' @param regimeVec vector with ones and zeros of the regime draw
#' @return draw of the standard deviation

SdEpsilon_fctn <- function(data, stateVec, regimeVec = NULL) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  data <- data[-c(1:10)]
  # Set the priors
  alpha_0 <- 0
  beta_0 <- 0
  # Compute the SSR from the regression of the data on the drift and seasonal term (maybe wider dist better)
  mu <- stateVec[1, ]
  gamma <- stateVec[3, ]
  v_t <- data - mu - gamma
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + length(data))
  beta_1 <- beta_0 + SSR
  # Draw the new value
  var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  sd <- sqrt(var)
  return(sd)
}


#' @description Function to draw the innovation standard deviation of the trend component
#' @param data vector with the observations
#' @param stateVec matrix of the state vector draw
#' @param regimeVec vector with ones and zeros of the regime draw
#' @param paramList list with the additional parameters
#' @return draw of the standard deviation

SdXi_fctn <- function(data, stateVec, regimeVec, paramList) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  regimeVec <- regimeVec[-c(1:10)]
  # Set the priors
  alpha_0 <- 1000
  # alpha_0 <- 0
  beta_0 <- 0
  # Compute the SSR from the regression of the trend on its lag and drift terms
  nu_1 <- paramList$nu_1
  mu_t <- stateVec[1, -NCOL(stateVec)]
  nu_0 <- stateVec[2, -NCOL(stateVec)]
  mu_lt <- stateVec[1, -1]
  S_t <- regimeVec[-length(regimeVec)]
  v_t <- mu_t - mu_lt - nu_0 - S_t * nu_1
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + length(data))
  beta_1 <- beta_0 + SSR
  # Draw the new value
  var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  sd <- sqrt(var)
  return(sd)
}


#' @description Function to draw the innovation standard deviation of the trend component in case of regime induced heterocedasticity
#' @param data vector with the observations
#' @param stateVec matrix of the state vector draw
#' @param regimeVec vector with ones and zeros of the regime draw
#' @param paramList list with the additional parameters
#' @return vector of the draw of the two standard deviations

SdXiSwitch_fctn <- function(data, stateVec, regimeDriftVec, regimeHeteroscVec, paramList) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  regimeDriftVec <- regimeDriftVec[-c(1:10)]
  regimeHeteroscVec <- regimeHeteroscVec[-c(1:10)]
  # Set the priors
  alpha_0 <- 0
  beta_0 <- 0
  alpha_2 <- 0
  beta_2 <- 0
  # Define variables
  nu_1 <- paramList$nu_1
  mu_t <- stateVec[1, -NCOL(stateVec)]
  nu_0 <- stateVec[2, -NCOL(stateVec)]
  mu_lt <- stateVec[1, -1]
  Drift_S_t <- regimeDriftVec[-length(regimeDriftVec)]
  Heterosc_S_t <- regimeHeteroscVec[-length(regimeHeteroscVec)]
  # Get multiplicative specification of the two variances
  xi_0 <- paramList$xi$xi_0
  xi_1 <- paramList$xi$xi_1
  x <- (xi_1^2 / xi_0^2 - xi_0^2)
  correcFactor <- (1 + x * Heterosc_S_t)^(-.5)
  # Compute the SSR from the regression of the trend on its lag and drift terms, corrected for S_t = 1
  v_t_1 <- (mu_t - mu_lt - nu_0 - Drift_S_t * nu_1) * correcFactor
  SSR_1 <- sum(v_t_1^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + length(data))
  beta_1 <- beta_0 + SSR_1
  # Draw the new value xi_0
  xi_0_var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  xi_0_draw <- sqrt(xi_0_var)
  # Compute the SSR from the regression of the trend on its lag and drift terms, corrected for S_t = 0
  v_t_2 <- Heterosc_S_t * (mu_t - mu_lt - nu_0 - Drift_S_t * nu_1) / xi_0_draw
  SSR_2 <- sum(v_t_2^2, na.rm = T)
  # Adjust the priors
  alpha_3 <- (alpha_2 + sum(S_t))
  beta_3 <- beta_2 + SSR_2
  # Draw the new value xi_0
  x_draw <- rinvgamma(1, alpha_3 / 2, beta_3 / 2)
  xi_1_draw <- sqrt(xi_0_var * x_draw)

  return(c(xi_0_draw, xi_1_draw))
}


#' @description Function to draw the innovation standard deviation of the seasonal component
#' @param stateVec matrix of the state vector draw
#' @param seasPos position of the seasonal component in the state vector
#' @return draw of the standard deviation

SdOmega_fctn <- function(stateVec, seasPos = 3) {
  # discard first 10 observations because of the high state vector variation due to the filter initialization
  stateVec <- stateVec[, -c(1:10)]
  dimens <- NROW(stateVec)
  # Set the priors so that SdOmega .005 (maybe wider dist better)
  # alpha_0 <- 1e5
  alpha_0 <- 0
  beta_0 <- 0
  # Compute the SSR based on the seasonal unit root (regression of the seasonal term on its determ dummies and the unit root)
  gamma_regress <- stateVec[c(seasPos:8, dimens), -NCOL(stateVec)]
  beta <- matrix(c(rep(-1, 6), 1), nr = 1)
  gamma <- stateVec[seasPos, -1]
  omega_t <- gamma - c(beta %*% gamma_regress)
  SSR <- sum(omega_t[-length(omega_t)]^2, na.rm = T)
  # Adjust the priors
  alpha_1 <- (alpha_0 + NCOL(stateVec))
  beta_1 <- beta_0 + SSR
  # Draw the new value
  var <- rinvgamma(1, alpha_1 / 2, beta_1 / 2)
  sd <- sqrt(var)
  return(sd)
}
