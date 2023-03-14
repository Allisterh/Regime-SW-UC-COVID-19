#' @description Function to constrain a parameter vector to the relevant parameter space for the original model specification
#' @param par vector of unconstrained parameters
#' @return a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:3] <- exp(par[1:3])
  # Defines regimes by nu_1 < 0
  constrPar[4] <- -exp(par[4])
  # Constrains probabilities to >=0 & <=1
  constrPar[5:8] <- 1 / (1 + exp(par[5:8]))
  return(constrPar)
}


#' @description Function that sets up a list with system matrices for the filter
#' @param paramList list with the additional parameters
#' @return a list with the system matrices

SystemMat_fctn <- function(paramList = NULL) {
  # Transition matrix
  Tt_trend <- Tt_trend_fctn()
  Tt_seasDeterm <- Tt_seasDeterm_fctn()
  Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm))
  Dimens <- NCOL(Tt)
  # Measurement matrix
  Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), nr = 1, nc = Dimens)
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- matrix(c(1, rep(0, Dimens - 1), rep(0, 2 * Dimens)), nc = 3, byrow = T)
  # Create the output list
  systemList <- list()
  systemList$Tt <- Tt
  systemList$Z <- Z
  systemList$R <- R

  return(systemList)
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
  xi_0 <- par[1] # Sd of innovations to the trend component in heterosc. regime 0
  xi_1 <- par[2] # Sd of innovations to the trend component in heterosc. regime 1
  epsilon <- par[3] # Sd of measurement equation residual
  nu_1 <- par[4] # Additional drift for down-turning drift regime 1
  q <- par[5] # Probability of staying in drift regime 0
  p <- par[6] # Probability of staying in drift regime 1
  r <- par[7] # Probability of staying in heterosc. regime 0
  s <- par[8] # Probability of staying in heterosc. regime 1
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- xi_0
  paramList$xi$xi_1 <- xi_1
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- 0
  paramList$epsilon <- epsilon
  paramList$nu_1 <- nu_1
  paramList$Probs$q <- q
  paramList$Probs$p <- p
  paramList$Probs$r <- r
  paramList$Probs$s <- s
  paramList$eta$eta_0 <- paramList$eta$eta_1 <- 0

  return(paramList)
}


#' @description Function that constructs the model specific random parameter table for the original model specification
#' @param nRandom number of random parameter vectors
#' @return tibble with the random parameter vectors

ThetaRand_fctn <- function(nRandom) {
  thetaRand <- tibble(
    xi_0 = log(runif(nRandom, 0.01, 0.25)),
    xi_1 = log(runif(nRandom, 0.01, 0.25)),
    epsilon = log(runif(nRandom, 0.01, 0.3)),
    nu_1 = log(runif(nRandom, 0.001, 0.3)),
    q = runif(nRandom, -5, 5),
    p = runif(nRandom, -5, 5),
    r = runif(nRandom, -5, 5),
    s = runif(nRandom, -5, 5)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the original model specification
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi_0 = sigma[xi[0]],
    xi_1 = sigma[xi[1]],
    epsilon = sigma[epsilon],
    nu_0 = nu[0],
    nu_1 = nu[1],
    p = italic(p),
    q = italic(q),
    r = italic(r),
    s = italic(s)
  )
  return(labelMap)
}


#' @description Function that draws the additional parameters as part of the Gibbs sampler
#' @param data vector for which the likelihood values are to be computed
#' @param stateVec array with draws of the state vector
#' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
#' @param paramList list with the additional parameters
#' @return tibble with the random parameter vectors

AdditionalParam_fctn <- function(data, stateVec, regimeDriftVec, regimeHeteroscVec, paramList) {
  # epsilon
  paramList$epsilon <- SdEpsilon_fctn(data = data, stateVec = stateVec, paramList = paramList, seasPos = 3)
  # xi
  xiVec <- SdXiSwitch_fctn(
    stateVec = stateVec, regimeDriftVec = regimeDriftVec, regimeHeteroscVec = regimeHeteroscVec,
    paramList = paramList
  )
  paramList$xi$xi_0 <- xiVec[1]
  paramList$xi$xi_1 <- xiVec[2]
  # nu_1
  paramList$nu_1 <- Nu_1_fctn(stateVec = stateVec, regimeDriftVec = regimeDriftVec, paramList = paramList, regimeHeteroscVec = regimeHeteroscVec)
  # p and q
  transProbsDrift <- TransProbs_fctn(
    regimeVec = regimeDriftVec, paramList = paramList, heteroscRegime = FALSE
  )
  paramList$Probs$q <- transProbsDrift[1]
  paramList$Probs$p <- transProbsDrift[2]
  # r and s
  transProbsHeterosc <- TransProbs_fctn(
    regimeVec = regimeHeteroscVec, paramList = paramList, heteroscRegime = TRUE
  )
  paramList$Probs$r <- transProbsHeterosc[1]
  paramList$Probs$s <- transProbsHeterosc[2]

  return(paramList)
}


#' @description Produces a vector of additional parameters from a list
#' @param paramList standardized list of additional parameters
#' @return vector with relevant parameters

outputParam_fctn <- function(paramList) {
  return(c(
    "xi_0" = paramList$xi$xi_0, "xi_1" = paramList$xi$xi_1, "epsilon" = paramList$epsilon,
    "nu_1" = paramList$nu_1, "q" = paramList$Probs$q, "p" = paramList$Probs$p, "r" = paramList$Probs$r, "s" = paramList$Probs$s
  ))
}


#' @description Simulates data according to model specifications
#' @param epsilonInd If TRUE, adds noise to the measurement equation
#' @return a vector with the simulated series

simData_fctn <- function(epsilonInd = TRUE) {
  # Create the noise term
  sdEpsilon <- .15
  epsilon <- rnorm(nTotal, sd = sdEpsilon)
  # Construct the series
  regimeDrift <- simRegime_fctn(nPeriods = nPeriods, heteroscInd = FALSE)
  regimeHeterosc <- simRegime_fctn(nPeriods = nPeriods, heteroscInd = TRUE)
  trend <- simTrend_fctn(nPeriods = nPeriods, regimeVec = regimeDrift, heteroscInd = TRUE, regimeHeteroscVec = regimeHeterosc)
  seas <- simDetermSeas_fctn(nPeriods = nPeriods)
  y <- trend + seas + epsilon
  return(y)
}


#' @description Returns initial values for the Gibbs sampler
#' @return a vector with initial values

paramIni_fctn <- function() {
  return(c(
    "xi_0" = 1e4, "xi_1" = 1e4, "epsilon" = 1e4,
    "nu_1" = -.5, "q" = .95, "p" = .95, "r" = .5, "s" = .5
  ))
}
