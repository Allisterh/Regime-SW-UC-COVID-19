#' @description Function to constrain a parameter vector to the relevant parameter space for the model specification
#' with a MA(7) cyclical component
#' @param par vector of unconstrained parameters
#' @returns a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:3] <- exp(par[1:3])
  # Defines regimes by nu_1 < 0
  constrPar[4] <- -exp(par[4])
  # Assures invertible MA(7) cycle
  
  # Constrains probabilities to >=0 & <=1
  constrPar[12:13] <- 1 / (1 + exp(par[12:13]))
  return(constrPar)
}


#' @description Function that reverses the ParConstrain_fctn to correct for the SE calculation for the model specification
#' with a MA(7) cyclical component
#' @param par vector of constrained parameters
#' @returns a diagonal matrix with the correction factors

InvParConstrain_fctn <- function(par) {
  correction <- rep(1, length(par))
  # Assures positive sd of system innovations
  correction[1:4] <- diag(numDeriv::jacobian(exp, par[1:4]))
  # Defines regimes by nu_1 < 0
  correction[4] <- -correction[4]
  # Constrains probabilities to >=0 & <=1
  correction[12:13] <- diag(numDeriv::jacobian(function(x) 1 / (1 + exp(x)), par[12:13]))
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
  omega <- par[2] # Sd of innovations to the seasonal component
  eta <- par[3] # Sd of innovations to the cyclical component
  nu_1 <- par[4] # Additional drift for down-turning regime 1
  phi_1 <- par[5] # MA covariates
  phi_2 <- par[6]
  phi_3 <- par[7]
  phi_4 <- par[8]
  phi_5 <- par[9]
  phi_6 <- par[10]
  phi_7 <- par[11]
  q <- par[12] # Probability of staying in regime 0
  p <- par[13] # Probability of staying in regime 1
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- paramList$xi$xi_1 <- xi
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- omega
  paramList$epsilon <- 0
  paramList$nu_1 <- nu_1
  paramList$Probs$q <- q
  paramList$Probs$p <- p
  paramList$eta$eta_0 <- paramList$eta$eta_1 <- eta
  paramList$phi <- c(phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7)
  
  return(paramList)
}


#' @description Function that sets up a list with system matrices for the filter
#' @param paramList list with additional parameters
#' @return a list with the system matrices

SystemMat_fctn <- function(paramList) {
  # Transition matrix
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  Tt_cycle <- rbind(rep(0, 8), cbind(diag(7), rep(0, 7)))
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  Tt <- as.matrix(bdiag(Tt_trend, Tt_cycle, Tt_seasDeterm, Tt_seasUr))
  # allows the seasonal unit root to load on gamma_t+1
  Tt[11, 17] <- 1
  Dimens <- NCOL(Tt)
  # Measurement matrix
  Z <- matrix(c(1, 0, 1, paramList$phi, 1, rep(0, Dimens - 11)), nr = 1)
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, 2), 1, rep(0, Dimens - 3)),
    c(rep(0, 16), 1, rep(0, Dimens - 17))
  )
  # Create the output list
  systemList <- list()
  systemList$Tt <- Tt
  systemList$Z <- Z
  systemList$R <- R
  systemList$R <- R
  
  return(systemList)
}


#' @description Function that constructs the model specific random parameter table for the model specification with a
#' MA(7) cyclical component
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
    phi_3 = runif(nRandom, -7, 7),
    phi_4 = runif(nRandom, -7, 7),
    phi_5 = runif(nRandom, -7, 7),
    phi_6 = runif(nRandom, -7, 7),
    phi_7 = runif(nRandom, -7, 7),
    q = runif(nRandom, -5, 5),
    p = runif(nRandom, -5, 5)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the model specification with with a
#' MA(7) cyclical component
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
    phi_3 = Phi[3],
    phi_4 = Phi[4],
    phi_5 = Phi[5],
    phi_6 = Phi[6],
    phi_7 = Phi[7],
    p = italic(p),
    q = italic(q)
  )
  return(labelMap)
}
