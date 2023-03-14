#' @description Function that constructs the transition matrix governing
#' the trend component
#' @return a square transition matrix

Tt_trend_fctn <- function() {
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  return(Tt_trend)
}


#' @description Function that constructs the transition matrix governing
#' a deterministic seasonal component
#' @return a square transition matrix

Tt_seasDeterm_fctn <- function() {
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  return(Tt_seasDeterm)
}


#' @description Function that constructs the transition matrix governing
#' an AR(2) cyclical component
#' @return a square transition matrix

Tt_AR2cycle_fctn <- function(phiVec) {
  Tt_cycle <- matrix(c(phiVec, 1, 0), byrow = T, nc = 2)
  return(Tt_cycle)
}


#' @description Function that constructs the transition matrix governing
#' a unit root seasonal process (deterministic part of the seasonal unit root
#' has to be adapted)
#' @return a square transition matrix

Tt_seasUr_fctn <- function() {
  Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  return(Tt_seasUr)
}
