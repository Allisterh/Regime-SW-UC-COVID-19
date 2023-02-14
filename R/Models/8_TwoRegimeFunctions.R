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
  omega <- par[3] # Sd of innovations to seasonal component
  epsilon <- par[4] # Sd of measurement equation residual
  nu_1 <- par[5] # Additional drift for down-turning drift regime 1
  q <- par[6] # Probability of staying in drift regime 0
  p <- par[7] # Probability of staying in drift regime 1
  r <- par[8] # Probability of staying in heterosc. regime 0
  s <- par[9] # Probability of staying in heterosc. regime 1
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- xi_0
  paramList$xi$xi_1 <- xi_1
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- omega
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
    omega = log(runif(nRandom, 0, 0.3)),
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
    xi_0 = sigma[xi_0],
    xi_1 = sigma[xi_1],
    omega = sigma[omega],
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

Additionalparam_fctn <- function(data, stateVec, regimeDriftVec, regimeHeteroscVec, paramList) {
  # epsilon
  paramList$epsilon <- SdEpsilon_fctn(data = data, stateVec = stateVec)
  # xi
  xiVec <- SdXiSwitch_fctn(
    data = data, stateVec = stateVec, regimeDriftVec = regimeDriftVec, regimeHeteroscVec = regimeHeteroscVec,
    paramList = paramList
  )
  paramList$xi$xi_0 <- xiVec[1]
  paramList$xi$xi_1 <- xiVec[2]
  # nu_1
  paramList$nu_1 <- Nu_1_TwoRegime_fctn(stateVec = stateVec, regimeDriftVec = regimeDriftVec, regimeHeteroscVec = regimeHeteroscVec, paramList = paramList)
  # omega
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- SdOmega_fctn(stateVec = stateVec, seasPos = 3)
  # p and q
  transProbsDrift <- TransProbsRegime_fctn(regimeDriftVec = regimeDriftVec, regimeHeteroscVec = regimeHeteroscVec, paramList = paramList)
  paramList$Probs$q <- transProbsDrift[1]
  paramList$Probs$p <- transProbsDrift[2]
  transProbsHeterosc <- TransProbsRegime_fctn(
    regimeDriftVec = regimeDriftVec, regimeHeteroscVec = regimeHeteroscVec, paramList = paramList,
    secondRegime = TRUE
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
    "xi_0" = paramList$xi$xi_0, "xi_1" = paramList$xi$xi_1, "omega" = paramList$omega$omega_0, "epsilon" = paramList$epsilon,
    "nu_1" = paramList$nu_1, "q" = paramList$Probs$q, "p" = paramList$Probs$p, "r" = paramList$Probs$r, "s" = paramList$Probs$s
  ))
}
