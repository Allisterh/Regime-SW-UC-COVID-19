#' @description Function to constrain a parameter vector to the relevant parameter space for the specification with
#' deterministic seasonality
#' @param par vector of unconstrained parameters
#' @return a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:2] <- exp(par[1:2])
  # Defines regimes by nu_1 < 0
  constrPar[3] <- -exp(par[3])
  # Constrains probabilities to >=0 & <=1
  constrPar[4:5] <- 1 / (1 + exp(par[4:5]))
  return(constrPar)
}


#' @description Function that reverses the ParConstrain_fctn to correct for the SE calculation for the specification with
#' deterministic seasonality
#' @param par vector of constrained parameters
#' @return a diagonal matrix with the correction factors

InvParConstrain_fctn <- function(par) {
  correction <- rep(1, length(par))
  # Assures positive sd of system innovations
  correction[1:3] <- diag(numDeriv::jacobian(exp, par[1:3]))
  # Defines regimes by nu_1 < 0
  correction[3] <- -correction[3]
  # Constrains probabilities to >=0 & <=1
  correction[4:5] <- diag(numDeriv::jacobian(function(x) 1 / (1 + exp(x)), par[4:5]))
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
  xi <- par[1] # Sd of innovations to the trend component
  epsilon <- par[2] # Sd of measurement equation residual
  nu_1 <- par[3] # Additional drift for down-turning regime 1
  q <- par[4] # Probability of staying in regime 0
  p <- par[5] # Probability of staying in regime 1
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- paramList$xi$xi_1 <- xi
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- 0
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


#' @description Function that constructs the model specific random parameter table
#' @param nRandom number of random parameter vectors
#' @return tibble with the random parameter vectors

ThetaRand_fctn <- function(nRandom) {
  thetaRand <- tibble(
    xi = log(runif(nRandom, 0.01, 0.25)),
    epsilon = log(runif(nRandom, 0.01, 0.3)),
    nu_1 = log(runif(nRandom, 0.001, 0.3)),
    q = runif(nRandom, -5, 5),
    p = runif(nRandom, -5, 5)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi = sigma[xi],
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
  paramList$xi$xi_0 <- paramList$xi$xi_1 <- SdXi_fctn(data = data, stateVec = stateVec, regimeVec = regimeVec, paramList = paramList)
  # nu_1
  paramList$nu_1 <- Nu_1_fctn(stateVec = stateVec, regimeVec = regimeVec, paramList = paramList)
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
    "xi" = paramList$xi$xi_0, "epsilon" = paramList$epsilon,
    "nu_1" = paramList$nu_1, "q" = paramList$Probs$q, "p" = paramList$Probs$p
  ))
}
