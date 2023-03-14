#' @description Function to constrain a parameter vector to the relevant parameter space for the model specification with
#' regime induced AR cycle heterogeneity
#' @param par vector of unconstrained parameters
#' @returns a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:3] <- exp(par[1:3])
  # Defines regimes by nu_1 < 0
  constrPar[4] <- -exp(par[4])
  # Assures stationary AR(2) cycle
  phi_1 <- par[5] / (1 + abs(par[5]))
  phi_2 <- (1 - abs(phi_1)) * par[6] / (1 + abs(par[6])) + abs(phi_1) - phi_1^2
  constrPar[5] <- 2 * phi_1
  constrPar[6] <- -(phi_1^2 + phi_2)
  # Constrains probabilities to >=0 & <=1
  constrPar[7:8] <- 1 / (1 + exp(par[7:8]))
  return(constrPar)
}


#' @description Function that sets up a list with additional parameters for the filter
#' @param par vector of parameters
#' @param constrainPar logical. If yes, the parameter constraining function for the numerical optimization
#' is applied
#' @return a list with the additional parameters

ParList_fctn <- function(par, constrainPar) {
  if (constrainPar == TRUE) {
    # Reverse parameter transformation
    par <- ParConstrain_fctn(par)
  }
  # Load input (provided constraining function is already reversed)
  xi <- par[1] # Sd of innovations to trend component
  omega <- par[2] # Sd of innovations to the seasonal component
  eta <- par[3] # Sd of innovations to the cycle component in regime 0
  nu_1 <- par[4] # Additional drift for down-turning regime 1
  phi_1 <- par[5] # AR covariate of the cycle
  phi_2 <- par[6] # AR covariate of the cycle
  q <- par[7] # Probability of staying in regime 0
  p <- par[8] # Probability of staying in regime 1

  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- paramList$xi$xi_1 <- xi
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- omega
  paramList$epsilon <- 0
  paramList$nu_1 <- nu_1
  paramList$Probs$q <- q
  paramList$Probs$p <- p
  paramList$eta$eta_0 <- paramList$eta$eta_1 <- eta
  paramList$phi <- c(phi_1, phi_2)

  return(paramList)
}


#' @description Function that sets up a list with system matrices for the filter
#' @param paramList list with the additional parameters
#' @return a list with the system matrices

SystemMat_fctn <- function(paramList) {
  # Transition matrix
  Tt_trend <- Tt_trend_fctn()
  Tt_cycle <- Tt_AR2cycle_fctn(phiVec = paramList$phi)
  Tt_seasDeterm <- Tt_seasDeterm_fctn()
  Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  Tt <- as.matrix(bdiag(Tt_trend, Tt_cycle, Tt_seasDeterm, Tt_seasUr))
  # allows the seasonal unit root to load on gamma_t+1
  Tt[5, 11] <- 1
  Dimens <- NCOL(Tt)
  # Measurement matrix
  Z <- matrix(c(1, 0, 1, 0, 1, rep(0, Dimens - 5)), 1, Dimens)
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, 2), 1, rep(0, Dimens - 3)),
    c(rep(0, 10), 1, rep(0, Dimens - 11))
  )
  # Create the output list
  systemList <- list()
  systemList$Tt <- Tt
  systemList$Z <- Z
  systemList$R <- R

  return(systemList)
}


#' @description Function that constructs the model specific random parameter table for the model specification with regime induced cycle
#' heterogeneity
#' @param nRandom number of random parameter vectors
#' @return tibble with the random parameter vectors

ThetaRand_fctn <- function(nRandom) {
  thetaRand <- tibble(
    xi = log(runif(nRandom, 0.01, 0.25)),
    omega = log(runif(nRandom, 0, 0.3)),
    eta = log(runif(nRandom, 0.01, 0.3)),
    nu_1 = log(runif(nRandom, 0.001, 0.3)),
    phi_1 = runif(nRandom, -7, 7),
    phi_2 = runif(nRandom, -7, 7),
    q = runif(nRandom, -5, 5),
    p = runif(nRandom, -5, 5)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the model specification with regime induced cycle
#' heterogeneity
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi = sigma[xi],
    omega = sigma[omega],
    eta = sigma[eta],
    nu_0 = nu[0],
    nu_1 = nu[1],
    phi_1 = Phi[1],
    phi_2 = Phi[2],
    p = italic(p),
    q = italic(q)
  )
  return(labelMap)
}


#' @description Function that draws the additional parameters as part of the Gibbs sampler
#' @param data vector for which the likelihood values are to be computed
#' @param stateVec array with draws of the state vector
#' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
#' @param paramList list with the additional parameters
#' @return tibble with the random parameter vectors

AdditionalParam_fctn <- function(data, stateVec, regimeVec, paramList) {
  # xi
  paramList$xi$xi_0 <- paramList$xi$xi_1 <- SdXi_fctn(stateVec = stateVec, regimeVec = regimeVec, paramList = paramList)
  # nu_1
  paramList$nu_1 <- Nu_1_fctn(stateVec = stateVec, regimeDriftVec = regimeVec, paramList = paramList)
  # omega
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- SdOmega_fctn(stateVec = stateVec, paramList = paramList, seasPos = 5)
  # p and q
  transProbs <- TransProbs_fctn(regimeVec = regimeVec, paramList = paramList)
  paramList$Probs$q <- transProbs[1]
  paramList$Probs$p <- transProbs[2]
  # eta
  eta <- SdEta_fctn(stateVec = stateVec, paramList = paramList, cyclePos = 3)
  paramList$eta$eta_0 <- eta
  paramList$eta$eta_1 <- eta
  # phi
  phi <- ARPhi_fctn(stateVec = stateVec, regimeHeteroscVec = regimeVec, paramList = paramList, cyclePos = 3)
  paramList$phi <- phi

  return(paramList)
}


#' @description Produces a vector of additional parameters from a list
#' @param paramList standardized list of additional parameters
#' @return vector with relevant parameters

outputParam_fctn <- function(paramList) {
  return(c(
    "xi" = paramList$xi$xi_0, "omega" = paramList$omega$omega_0, "eta" = paramList$eta$eta, "nu_1" = paramList$nu_1,
    "phi_1" = paramList$phi[1], "phi_2" = paramList$phi[2], "q" = paramList$Probs$q, "p" = paramList$Probs$p
  ))
}


#' @description Simulates data according to model specifications
#' @param nPeriods Number of time periods to be simulated
#' @return a vector with the simulated series

simData_fctn <- function(nPeriods) {
  # Construct the series
  regime <- simRegime_fctn(nPeriods = nPeriods)
  trend <- simTrend_fcnt(nPeriods = nPeriods, regimeVec = regime)
  seas <- simURSeas_fctn(nPeriods = nPeriods)
  cycle <- simCycle_fctn(nPeriods)
  y <- trend + seas + cycle
  return(y)
}


#' @description Returns initial values for the Gibbs sampler
#' @return a vector with initial values

paramIni_fctn <- function() {
  return(c(
    "xi" = 1e5, "omega" = 1e5, "epsilon" = 1e5,
    "nu_1" - .5, "phi_1" = 0, "phi_2" = 0, "q" = .95, "p" = .95
  ))
}
