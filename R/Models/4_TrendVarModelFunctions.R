#' @description Function to constrain a parameter vector to the relevant parameter space for the model specification with
#' regime induced trend heterogeneity
#' @param par vector of unconstrained parameters
#' @returns a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:4] <- exp(par[1:4])
  # Defines regimes by nu_1 < 0
  constrPar[5] <- -exp(par[5])
  # Constrains probabilities to >=0 & <=1
  constrPar[6:7] <- 1 / (1 + exp(par[6:7]))
  return(constrPar)
}


#' @description Function that reverses the ParConstrain_fctn to correct for the SE calculation for the model specification with
#' regime induced trend heterogeneity
#' @param par vector of constrained parameters
#' @returns a diagonal matrix with the correction factors

InvParConstrain_fctn <- function(par) {
  correction <- rep(1, length(par))
  # Assures positive sd of system innovations
  correction[1:5] <- diag(numDeriv::jacobian(exp, par[1:5]))
  # Defines regimes by nu_1 < 0
  correction[5] <- -correction[5]
  # Constrains probabilities to >=0 & <=1
  correction[6:7] <- diag(numDeriv::jacobian(function(x) 1 / (1 + exp(x)), par[6:7]))
  correctionMat <- diag(correction)
  return(correctionMat)
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
  xi_0 <- par[1] # Sd of innovations to trend component in regime 0
  xi_1 <- par[2] # Sd of innovations to trend component in regime 1
  omega <- par[3] # Sd of innovations to seasonal
  epsilon <- par[4] # Sd of measurement equation residual
  nu_1 <- par[5] # Additional drift for down-turning regime 1
  q <- par[6] # Probability of staying in regime 0
  p <- par[7] # Probability of staying in regime 1
  
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- xi_0
  paramList$xi$xi_1 <- xi_1
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- omega
  paramList$epsilon <- epsilon
  paramList$nu_1 <- nu_1
  paramList$Probs$q <- q
  paramList$Probs$p <- p
  paramList$eta$eta_0 <- paramList$eta$eta_1 <- 0
  
  return(paramList)
}


#' @description Function that sets up a list with system matrices for the filter
#' @param paramList list with the additional parameters
#' @return a list with the system matrices

SystemMat_fctn <- function(paramList = NULL) {
  # Transition matrix
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_seasUr))
  # allows the seasonal unit root to load on gamma_t+1
  Tt[3, 9] <- 1
  Dimens <- NCOL(Tt)
  # Measurement matrix
  Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    rep(0, Dimens),
    c(rep(0, 8), 1, rep(0, Dimens - 9))
  )
  # Create the output list
  systemList <- list()
  systemList$Tt <- Tt
  systemList$Z <- Z
  systemList$R <- R
  
  return(systemList)
}


#' @description Function that sets up the Kim filter
#' @param par vector of parameters
#' @param constrainPar logical. If yes, the parameter constraining function for the numerical optimization
#' is applied
#' @return a list with the additional parameters and system matrices

ModelSetup_fctn <- function(par, constrainPar) {
  #---------------------------------------------------------------------------------------#
  # Initialize filter
  #---------------------------------------------------------------------------------------#
  
  if (constrainPar == TRUE){
    # Reverse parameter transformation
    par <- ParConstrain_fctn(par)
  }
  # Load input (provided constraining function is already reversed)
  xi_0 <- par[1] # Sd of innovations to trend component in regime 0
  xi_1 <- par[2] # Sd of innovations to trend component in regime 1
  omega <- par[3] # Sd of innovations to seasonal
  epsilon <- par[4] # Sd of measurement equation residual
  nu_1 <- par[5] # Additional drift for down-turning regime 1
  q <- par[6] # Probability of staying in regime 0
  p <- par[7] # Probability of staying in regime 1
  
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- xi_0
  paramList$xi$xi_1 <- xi_1
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- omega
  paramList$epsilon <- epsilon
  paramList$nu_1 <- nu_1
  paramList$Probs$q <- q
  paramList$Probs$p <- p
  paramList$eta$eta_0 <- paramList$eta$eta_1 <- 0
  
  # Transition matrix
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_seasUr))
  # allows the seasonal unit root to load on gamma_t+1
  Tt[3, 9] <- 1
  Dimens <- NCOL(Tt)
  # Measurement matrix
  Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    rep(0, Dimens),
    c(rep(0, 8), 1, rep(0, Dimens - 9))
  )
  # Create the output list
  systemList <- list()
  systemList$Tt <- Tt
  systemList$Z <- Z
  systemList$R <- R
  
  return(list("paramList" = paramList, "systemList" = systemList))
}


#' @description Function that constructs the model specific random parameter table for the model specification with regime 
#' induced trend heterogeneity
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
    p = runif(nRandom, -5, 5)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the model specification with regime induced 
#' trend heterogeneity
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi_0 = sigma[xi][0],
    xi_1 = sigma[xi][1],
    omega = sigma[omega],
    epsilon = sigma[epsilon],
    nu_0 = nu[0],
    nu_1 = nu[1],
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

Additionalparam_fctn <- function(data, stateVec, regimeVec, paramList) {
  # epsilon
  paramList$epsilon <- SdEpsilon_fctn(data = data, stateVec = stateVec, regimeVec = regimeVec)
  # xi
  xiVec <- SdXiSwitch_fctn(data = data, stateVec = stateVec, regimeVec = regimeVec, paramList = paramList)
  paramList$xi$xi_0 <- xiVec[1]
  paramList$xi$xi_1 <- xiVec[2]
  # nu_1
  paramList$nu_1 <- Nu_1_fctn(stateVec = stateVec, regimeVec = regimeVec, paramList = paramList)
  # omega
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- SdOmega_fctn(stateVec = stateVec, seasPos = 3)
  # p and q
  transProbs <- TransProbs_fctn(regimeVec = regimeVec, paramList = paramList)
  paramList$Probs$q <- transProbs[1]
  paramList$Probs$p <- transProbs[2]
  
  return(paramList)
}


#' @description Produces a vector of additional parameters from a list
#' @param paramList standardized list of additional parameters
#' @return vector with relevant parameters

outputParam_fctn <- function(paramList) {
  return(c(
    "xi_0" = paramList$xi$xi_0, "xi_1" = paramList$xi$xi_1, "omega" = paramList$omega$omega_0, "epsilon" = paramList$epsilon,
    "nu_1" = paramList$nu_1, "q" = paramList$Probs$q, "p" = paramList$Probs$p
  ))
}
