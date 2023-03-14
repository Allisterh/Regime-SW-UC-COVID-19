#' @description Function to constrain a parameter vector to the relevant parameter space for the three state model specification
#' @param par vector of unconstrained parameters
#' @returns a vector of constrained parameters

ParConstrain_fctn <- function(par) {
  constrPar <- par <- as.numeric(par)
  # Assures positive sd of system innovations
  constrPar[1:2] <- exp(par[1:2])
  # Defines regimes by nu_1 < 0
  constrPar[3] <- -exp(par[3])
  # Constrains probabilities to >=0 & <=1
  posProbs <- 1 / (1 + exp(-par[4:9]))
  constrPar[4:9] <- posProbs
  phi_1 <- par[10] / (1 + abs(par[10]))
  phi_2 <- (1 - abs(phi_1)) * par[11] / (1 + abs(par[11])) + abs(phi_1) - phi_1^2
  constrPar[10] <- 2 * phi_1
  constrPar[11] <- -(phi_1^2 + phi_2)
  return(constrPar)
}


#' @description Function that constructs the model specific random parameter table for the three state model specification
#' @param nRandom number of random parameter vectors
#' @return tibble with the random parameter vectors

ThetaRand_fctn <- function(nRandom) {
  thetaRand <- tibble(
    xi = log(runif(nRandom, 0.01, 0.25)),
    eta = log(runif(nRandom, 0.01, 0.3)),
    nu_1 = log(runif(nRandom, 0.001, 0.3)),
    p_00 = runif(nRandom, -7, -2),
    p_01 = runif(nRandom, -7, -2),
    p_10 = runif(nRandom, -7, -2),
    p_11 = runif(nRandom, -7, -2),
    p_20 = runif(nRandom, -7, -2),
    p_21 = runif(nRandom, -7, -2),
    phi_1 = runif(nRandom, -7, 7),
    phi_2 = runif(nRandom, -7, 7)
  )
  return(thetaRand)
}


#' @description Function that constructs the model specific label map for the three state model specification
#' @return label map

LabelMap_fctn <- function() {
  labelMap <- expression(
    xi = sigma[xi],
    eta = sigma[eta],
    nu_0 = nu[0],
    nu_1 = nu[1],
    p_00 = P["00"],
    p_01 = P["01"],
    p_10 = P["10"],
    p_11 = P["11"],
    p_20 = P["20"],
    p_21 = P["21"],
    phi_1 = Phi[1],
    phi_2 = Phi[2]
  )
  return(labelMap)
}
