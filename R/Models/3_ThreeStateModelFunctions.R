#' @description Function to constrain a parameter vector to the relevant parameter space for the three state model specification
#' @param par vector of unconstrained parameters
#' @returns a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:3] <- exp(par[1:3])
  # Defines regimes by nu_1 < 0
  constrPar[4] <- -exp(par[4])
  # Constrains probabilities to >=0 & <=1
  posProbs <- 1 / (1 + exp(-par[5:10] + .5))
  z <- c(posProbs[1] + posProbs[2], posProbs[3] + posProbs[4], posProbs[5] + posProbs[6])
  transZ <- 1 / (1 + exp(-z + .5))
  constrPar[5] <- transZ[1] - posProbs[2]
  constrPar[6] <- posProbs[2]
  constrPar[7] <- transZ[2] - posProbs[4]
  constrPar[8] <- posProbs[4]
  constrPar[9] <- transZ[3] - posProbs[6]
  constrPar[10] <- posProbs[6]
  return(constrPar)
}


#' @description Function that reverses the ParConstrain_fctn to correct for the SE calculation for the three state model specification
#' @param par vector of constrained parameters
#' @returns a diagonal matrix with the correction factors

InvParConstrain_fctn <- function(par) {
  correction <- rep(1, length(par))
  # Assures positive sd of system innovations
  correction[1:4] <- diag(numDeriv::jacobian(exp, par[1:4]))
  # Defines regimes by nu_1 < 0
  correction[4] <- -correction[4]
  # Constrains probabilities to >=0 & <=1
  correction[5:10] <- diag(numDeriv::jacobian(function(x) 1 / (1 + exp(x)), par[5:10]))
  correctionMat <- diag(correction)
  return(correctionMat)
}


#' @description Function that constructs the model specific random parameter table for the three state model specification
#' @param nRandom number of random parameter vectors
#' @return tibble with the random parameter vectors

ThetaRand_fctn <- function(nRandom) {
  thetaRand <- tibble(
    xi = log(runif(nRandom, 0.01, 0.25)),
    omega = log(runif(nRandom, 0, 0.3)),
    epsilon = log(runif(nRandom, 0.01, 0.3)),
    nu_1 = log(runif(nRandom, 0.001, 0.3)),
    p_22 = runif(nRandom, -7, -2),
    p_01 = runif(nRandom, -7, -2),
    p_02 = runif(nRandom, -7, -2),
    p_11 = runif(nRandom, -7, -2),
    p_12 = runif(nRandom, -7, -2),
    p_21 = runif(nRandom, -7, -2)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the three state model specification
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi = sigma[xi],
    omega = sigma[omega],
    epsilon = sigma[epsilon],
    nu_0 = nu[0],
    nu_1 = nu[1],
    p_01 = P["01"],
    p_02 = P["02"],
    p_11 = P[11],
    p_12 = P[11],
    p_21 = P[21],
    p_22 = P[22]
  )
  return(labelMap)
}
