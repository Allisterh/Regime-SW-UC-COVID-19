#' @description Function that simulates 2 infection regimes
#' @param nPeriods length of the simulated process
#' @param heteroscInd boolean. If TRUE, different transition probabilities are used
#' @return vector of simulated regimes

simRegime_fctn <- function(nPeriods, heteroscInd = FALSE) {
  # Set transition probabilities
  q <- ifelse(heteroscInd == TRUE, .97, .8)
  p <- ifelse(heteroscInd == TRUE, .99, .9)
  # Simulate the series
  regime <- rep(0, nPeriods)
  regime[1] <- rbinom(1, 1, .5)
  for (i in 2:nPeriods) {
    prob <- ifelse(regime[i - 1] == 1, p, 1 - q)
    regime[i] <- rbinom(1, 1, prob)
  }
  return(regime)
}


#' @description Function that simulates the trend process
#' @param nPeriods length of the simulated process
#' @param regimeVec vector of 1s and 0s governing the trend drift
#' @param regimeHeteroscVec vector of 1s and 0s governing the trend heteroscedasticity. If NULL,
#' the regimeDrift vector is used.
#' @return vector of the simulated trend

simTrend_fctn <- function(nPeriods, regimeVec, heteroscInd = FALSE, regimeHeteroscVec = NULL) {
  if (is.null(regimeHeteroscVec)) regimeHeteroscVec <- regimeVec
  # Set the drift
  nu_0 <- .04
  nu_1 <- -.06
  # Simulate the innovations
  sdXi_0 <- .05
  sdXi_1 <- ifelse(heteroscInd == FALSE, .05, .1)
  xi_0 <- rnorm(nPeriods, sd = sdXi_0)
  xi_1 <- rnorm(nPeriods, sd = sdXi_1)
  xi <- ifelse(regimeHeteroscVec == 0, xi_0, xi_1)
  # Construct the trend. First period is set to approximate the COVID-19 infections
  trend <- rep(0, nPeriods)
  trend[1] <- 10
  for (i in 2:nPeriods) {
    trend[i] <- trend[i - 1] + xi[i] + nu_0 + nu_1 * regimeVec[i]
  }
  return(trend)
}


#' @description Function that simulates the seasonal unit root
#' @param nPeriods length of the simulated process
#' @return vector of the simulated seasonality

simURSeas_fctn <- function(nPeriods) {
  # Simulate the innovations
  sdOmega <- .05
  omega <- rnorm(nPeriods + 7, sd = sdOmega)
  # Recursively construct the seasonal unit root
  seas <- rep(0, nPeriods + 7)
  x <- rep(0, nPeriods + 7)
  for (i in 8:(nPeriods + 7)) {
    x[i] <- x[i - 7] + omega[i]
    seas[i] <- -sum(
      seas[i - 1], seas[i - 2], seas[i - 3], seas[i - 4],
      seas[i - 5], seas[i - 6]
    ) + x[i - 1]
  }
  seas <- tail(seas, -7)
  return(seas)
}


#' @description Function that simulates the deterministic seasonality
#' @param nPeriods length of the simulated process
#' @return vector of the simulated seasonality

simDetermSeas_fctn <- function(nPeriods) {
  seas <- rep_len(c(0.03, -0.05, -0.1, 0.05, 0.01, 0.04, 0.02), nPeriods)
  return(seas)
}


#' @description Function that simulates the seasonal cyclical component
#' @param nPeriods length of the simulated process
#' @return vector of the simulated cycle

simCycle_fctn <- function(nPeriods, regimeHeteroscVec = NULL) {
  if (is.null(regimeHeteroscVec)) regimeHeteroscVec <- rep(0, nPeriods + 3)
  # Simulate the innovations
  sdEta_0 <- .5
  sdEta_1 <- .1
  phi_1 <- .5
  phi_2 <- -.2
  eta_0 <- rnorm(nPeriods + 3, sd = sdEta_0)
  eta_1 <- rnorm(nPeriods + 3, sd = sdEta_1)
  eta <- ifelse(regimeHeteroscVec == 0, eta_0, eta_1)
  # Recursively construct the cyclical term
  cycle <- rep(0, nPeriods + 3)
  for (i in 3:(nPeriods + 3)) {
    cycle[i] <- phi_1 * cycle[i - 1] + phi_2 * cycle[i - 2] + eta[i]
  }
  cycle <- tail(cycle, -3)
  return(cycle)
}
