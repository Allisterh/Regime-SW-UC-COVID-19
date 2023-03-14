#' @description Function to constrain a parameter vector to the relevant parameter space for the original model specification
#' @param par vector of unconstrained parameters
#' @returns a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive semi definite var-cov variance of system innovations
  constrPar[1:3] <- exp(par[1:3])
  constrPar[4] <- par[4] / (1 + abs(par[4]))
  # Defines regimes by nu_1 < 0
  constrPar[5] <- -exp(par[5])
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
  xi <- par[1] # Sd of innovations to the trend component
  omega <- par[2] # Sd of innovations to seasonal component
  epsilon <- par[3] # Sd of measurement equation residual
  varrho <- par[4] # Correlation between trend and regime
  nu_1 <- par[5] # Additional drift for down-turning regime 1
  beta_0 <- par[6] # Intercept for the Probit Markov chain specification 
  beta_1 <- par[7] # Covariate for the Probit Markov chain specification
  # Create the output list
  paramList <- list()
  paramList$xi$xi_0 <- paramList$xi$xi_1 <- xi
  paramList$omega$omega_0 <- paramList$omega$omega_1 <- omega
  paramList$epsilon <- epsilon
  paramList$nu_1 <- nu_1
  paramList$beta$beta_0 <- beta_0
  paramList$beta$beta_1 <- beta_1
  paramList$varrho <- varrho
  paramList$eta$eta_0 <- paramList$eta$eta_1 <- 0
  
  return(paramList)
}


#' @description Function that sets up a list with system matrices for the filter
#' @param paramList list with the additional parameters
#' @return a list with the system matrices

SystemMat_fctn <- function(paramList = NULL) {
  # Transition matrix
  Tt_trend <- Tt_trend_fctn()
  Tt_seasDeterm <- Tt_seasDeterm_fctn()
  Tt_seasUr <- Tt_seasUr_fctn()
  Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_seasUr))
  # allows the seasonal unit root to load on gamma_t+1
  Tt[3, 9] <- 1
  Dimens <- NCOL(Tt)
  # Measurement matrix
  Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, Dimens)),
    c(rep(0, 8), 1, rep(0, Dimens - 9))
  )
  # Create the output list
  systemList <- list()
  systemList$Tt <- Tt
  systemList$Z <- Z
  systemList$R <- R
  
  return(systemList)
}


#' @description Function that constructs the model specific random parameter table for the original model specification
#' @param nRandom number of random parameter vectors
#' @return tibble with the random parameter vectors

ThetaRand_fctn <- function(nRandom) {
  thetaRand <- tibble(
    xi = log(runif(nRandom, 0.01, 0.25)),
    omega = log(runif(nRandom, 0, 0.3)),
    epsilon = log(runif(nRandom, 0.01, 0.3)),
    varrho = runif(nRandom, -5, 5),
    nu_1 = log(runif(nRandom, 0.001, 0.3)),
    beta_0 = runif(nRandom, -1, 1),
    beta_1 = runif(nRandom, 2.5, 4.5)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the original model specification
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi = sigma[xi],
    omega = sigma[omega],
    epsilon = sigma[epsilon],
    varrho = sigma[varrho],
    nu_0 = nu[0],
    nu_1 = nu[1],
    p = italic(p),
    q = italic(q)
  )
  return(labelMap)
}


#' @description Function that converts Probit AR covariates to transition probabilities
#' @param coefMat Matrix with the ML parameters
#' @return Matrix where AR covariates are switched for probabilities

ProbitProb_fctn <- function(coefMat) {
  probitCoefs <- colnames(coefMat)[str_detect(colnames(coefMat), "beta")]
  coefMatReduced <- coefMat[,!(colnames(coefMat) %in% probitCoefs)]
  probitCoefMat <- coefMat[,probitCoefs]
  probsMat <- apply(probitCoefMat, 1, function(coefs){
    beta_0 <- coefs[1]
    beta_1 <- coefs[2]
    q <- pnorm(-beta_0)
    p <- 1 - pnorm(-beta_0 - beta_1)
    return(c(q, p))
  }) %>%
    Transp()
  colnames(probsMat) <- c("q", "p")
  return(cbind(coefMatReduced, probsMat))
}
