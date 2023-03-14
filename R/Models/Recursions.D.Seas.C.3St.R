#' Function that computes the log likelihood produced by the Kim filter forecast error decomposition
#' for the three state model specification
#' @param par vector of unconstrained parameters
#' @param data vector of the data
#' @param Ini initial value for the trend component
#' @returns the log likelihood

Likelihood_fctn <- function(par, data, Ini) {
  #---------------------------------------------------------------------------------------#
  # Initialize filter
  #---------------------------------------------------------------------------------------#

  # Reverse parameter transformation
  par <- ParConstrain_fctn(par)

  xi <- par[1] # Sd of innovations to drift component
  # omega <- par[2] # Sd of innovations to seasonal component
  omega <- 0
  # epsilon <- par[2] # Sd of measurement equation residual
  epsilon <- 0
  # epsilonSq <- epsilon^2
  epsilonSq <- 1e-5
  eta <- par[2]
  nu_1 <- par[3] # Additional drift for down-turning regime 1
  P <- matrix(NA, 3, 3)
  phi_1 <- par[10]
  phi_2 <- par[11]
  
  # Notation here is like in Hamilton (1994) : P[2, 1] gives Pr(0 --> 1)
  P[, 1] <- c(par[4], par[5], 1 - par[4] - par[5]) # P(St =  | St-1 = 0)
  P[, 2] <- c(par[6], par[7], 1 - par[6] - par[7]) # P(St =  | St-1 = 1)
  P[, 3] <- c(par[8], par[9], 1 - par[8] - par[9]) # P(St =  | St-1 = 1)
  if(any(c(par[4] + par[5], par[6] + par[7], par[8] + par[9]) > 1)){
    return(1e6 * max(par[4] + par[5], par[6] + par[7], par[8] + par[9]))
  }

  # Transition matrix
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  Tt_cycle <- matrix(c(phi_1, phi_2, 1, 0), byrow = T, nc = 2)
  
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  # Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  # Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_seasUr))
  # Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm))
  Tt <- as.matrix(bdiag(Tt_trend, Tt_cycle, Tt_seasDeterm))
  
  Tt_noDrift <- Tt
  Tt_noDrift[1, 2] <- 0
  # allows the seasonal unit root to load on gamma_t+1
  # Tt[3, 9] <- 1
  TranspTt <- Transp(Tt)
  
  # Dimension of state vector
  Dimens <- NCOL(Tt)
  
  # Number of time periods
  periods <- NROW(data)
  
  # Measurement matrix
  # Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)
  Z <- matrix(c(1, 0, 1, 0, 1, rep(0, Dimens - 5)), 1, Dimens)
  
  TranspZ <- Transp(Z)
  
  # Initial unconditional regime probabilities (steady state unconditional probabilities)
  A <- matrix(rbind(diag(3) - P, rep(1, 3)), nc = 3)
  Pi <- Inverse(Transp(A) %*% A) %*% t(A) %*% matrix(diag(4)[, 4], ncol = 1)
  Pr_ct_0 <- Pi[1]
  Pr_ct_1 <- Pi[2]
  Pr_ct_2 <- Pi[3]
  
  # Var-Cov Matrix of innovations in transition equation
  Q <- diag(c(xi^2, eta^2, omega^2))
  
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, 2), 1, rep(0, Dimens - 3)),
    c(rep(0, Dimens))
  )
  
  transpR <- Transp(R)
  # Var-Cov Matrix of innovations in transition equation
  expandQ <- R %*% Q %*% transpR
  
  # Vector holding additional Drift in Regime 1
  lambda <- matrix(c(nu_1, rep(0, Dimens - 1)), Dimens)
  
  # Initial values for state vector
  a_t_clt_00 <- a_t_clt_10 <- a_t_clt_20 <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens)
  a_t_clt_01 <- a_t_clt_11 <- a_t_clt_21 <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens) + lambda
  a_t_clt_02 <- a_t_clt_12 <- a_t_clt_22 <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens)
  
  # Initial variance matrix for state vector (Exact diffuse initialization)
  P_t_clt_00 <- P_t_clt_10 <- P_t_clt_20 <- P_t_clt_01 <- P_t_clt_11 <- P_t_clt_21 <- P_t_clt_02 <- P_t_clt_12 <- P_t_clt_22 <- diag(Dimens) * 10000

  # Initializes log-likelihood to record values
  logLik_T <- 0

  #---------------------------------------------------------------------------------------#
  # Kim filter to set up likelihood function
  #---------------------------------------------------------------------------------------#

  for (i in 1:periods) {
    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # One step ahead prediction error with state vector prediction from
    # (t|t-1)
    v_t_00 <- as.numeric(data[i] - Z %*% a_t_clt_00)
    v_t_10 <- as.numeric(data[i] - Z %*% a_t_clt_10)
    v_t_20 <- as.numeric(data[i] - Z %*% a_t_clt_20)
    v_t_01 <- as.numeric(data[i] - Z %*% a_t_clt_01)
    v_t_11 <- as.numeric(data[i] - Z %*% a_t_clt_11)
    v_t_21 <- as.numeric(data[i] - Z %*% a_t_clt_21)
    v_t_02 <- as.numeric(data[i] - Z %*% a_t_clt_02)
    v_t_12 <- as.numeric(data[i] - Z %*% a_t_clt_12)
    v_t_22 <- as.numeric(data[i] - Z %*% a_t_clt_22)

    # Variance of prediction error
    F_t_00 <- as.numeric(Z %*% P_t_clt_00 %*% TranspZ + epsilonSq)
    F_t_10 <- as.numeric(Z %*% P_t_clt_10 %*% TranspZ + epsilonSq)
    F_t_20 <- as.numeric(Z %*% P_t_clt_20 %*% TranspZ + epsilonSq)
    F_t_01 <- as.numeric(Z %*% P_t_clt_01 %*% TranspZ + epsilonSq)
    F_t_11 <- as.numeric(Z %*% P_t_clt_11 %*% TranspZ + epsilonSq)
    F_t_21 <- as.numeric(Z %*% P_t_clt_21 %*% TranspZ + epsilonSq)
    F_t_02 <- as.numeric(Z %*% P_t_clt_02 %*% TranspZ + epsilonSq)
    F_t_12 <- as.numeric(Z %*% P_t_clt_12 %*% TranspZ + epsilonSq)
    F_t_22 <- as.numeric(Z %*% P_t_clt_22 %*% TranspZ + epsilonSq)
    
    # Updating step
    # (t|t, S_t, S_t-1)
    # State vector
    a_t_ct_00 <- a_t_clt_00 + P_t_clt_00 %*% TranspZ %*% (1 / F_t_00 * v_t_00)
    a_t_ct_10 <- a_t_clt_10 + P_t_clt_10 %*% TranspZ %*% (1 / F_t_10 * v_t_10)
    a_t_ct_20 <- a_t_clt_20 + P_t_clt_20 %*% TranspZ %*% (1 / F_t_20 * v_t_20)
    a_t_ct_01 <- a_t_clt_01 + P_t_clt_01 %*% TranspZ %*% (1 / F_t_01 * v_t_01)
    a_t_ct_11 <- a_t_clt_11 + P_t_clt_11 %*% TranspZ %*% (1 / F_t_11 * v_t_11)
    a_t_ct_21 <- a_t_clt_21 + P_t_clt_21 %*% TranspZ %*% (1 / F_t_21 * v_t_21)
    a_t_ct_02 <- a_t_clt_02 + P_t_clt_02 %*% TranspZ %*% (1 / F_t_02 * v_t_02)
    a_t_ct_12 <- a_t_clt_12 + P_t_clt_12 %*% TranspZ %*% (1 / F_t_12 * v_t_12)
    a_t_ct_22 <- a_t_clt_22 + P_t_clt_22 %*% TranspZ %*% (1 / F_t_22 * v_t_22)

    # State vector Var-Matrix
    P_t_ct_00 <- P_t_clt_00 - P_t_clt_00 %*% TranspZ %*% (1 / F_t_00) %*% Z %*% P_t_clt_00
    P_t_ct_10 <- P_t_clt_10 - P_t_clt_10 %*% TranspZ %*% (1 / F_t_10) %*% Z %*% P_t_clt_10
    P_t_ct_20 <- P_t_clt_20 - P_t_clt_20 %*% TranspZ %*% (1 / F_t_20) %*% Z %*% P_t_clt_20
    P_t_ct_01 <- P_t_clt_01 - P_t_clt_01 %*% TranspZ %*% (1 / F_t_01) %*% Z %*% P_t_clt_01
    P_t_ct_11 <- P_t_clt_11 - P_t_clt_11 %*% TranspZ %*% (1 / F_t_11) %*% Z %*% P_t_clt_11
    P_t_ct_21 <- P_t_clt_21 - P_t_clt_21 %*% TranspZ %*% (1 / F_t_21) %*% Z %*% P_t_clt_21
    P_t_ct_02 <- P_t_clt_02 - P_t_clt_02 %*% TranspZ %*% (1 / F_t_02) %*% Z %*% P_t_clt_02
    P_t_ct_12 <- P_t_clt_12 - P_t_clt_12 %*% TranspZ %*% (1 / F_t_12) %*% Z %*% P_t_clt_12
    P_t_ct_22 <- P_t_clt_22 - P_t_clt_22 %*% TranspZ %*% (1 / F_t_22) %*% Z %*% P_t_clt_22

    #-------------------#
    ### Hamilton part ###
    #-------------------#

    # Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
    # (t|t-1)
    # Again note: P[1, 2] holds P(S_t = 0 | S_t-1 = 1)
    Pr_clt_00 <- P[1, 1] * Pr_ct_0
    Pr_clt_10 <- P[1, 2] * Pr_ct_1
    Pr_clt_20 <- P[1, 3] * Pr_ct_2
    Pr_clt_01 <- P[2, 1] * Pr_ct_0
    Pr_clt_11 <- P[2, 2] * Pr_ct_1
    Pr_clt_21 <- P[2, 3] * Pr_ct_2
    Pr_clt_02 <- P[3, 1] * Pr_ct_0
    Pr_clt_12 <- P[3, 2] * Pr_ct_1
    Pr_clt_22 <- P[3, 3] * Pr_ct_2

    # Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
    y_dens_clt_00 <- dnorm(v_t_00, sd = sqrt(F_t_00))
    y_dens_clt_10 <- dnorm(v_t_10, sd = sqrt(F_t_10))
    y_dens_clt_20 <- dnorm(v_t_20, sd = sqrt(F_t_20))
    y_dens_clt_01 <- dnorm(v_t_01, sd = sqrt(F_t_01))
    y_dens_clt_11 <- dnorm(v_t_11, sd = sqrt(F_t_11))
    y_dens_clt_21 <- dnorm(v_t_21, sd = sqrt(F_t_21))
    y_dens_clt_02 <- dnorm(v_t_02, sd = sqrt(F_t_02))
    y_dens_clt_12 <- dnorm(v_t_12, sd = sqrt(F_t_12))
    y_dens_clt_22 <- dnorm(v_t_22, sd = sqrt(F_t_22))

    # Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
    # (t|t-1))
    y_dens_clt <- sum(
      Pr_clt_00 * y_dens_clt_00, Pr_clt_10 * y_dens_clt_10, Pr_clt_20 * y_dens_clt_20,
      Pr_clt_01 * y_dens_clt_01, Pr_clt_11 * y_dens_clt_11, Pr_clt_21 * y_dens_clt_21,
      Pr_clt_02 * y_dens_clt_02, Pr_clt_12 * y_dens_clt_12, Pr_clt_22 * y_dens_clt_22
    )

    # Store approximate likelihood
    loglik_t <- -log(y_dens_clt)

    # Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
    # (t|t)
    Pr_ct_00 <- Pr_clt_00 * y_dens_clt_00 / y_dens_clt
    Pr_ct_10 <- Pr_clt_10 * y_dens_clt_10 / y_dens_clt
    Pr_ct_20 <- Pr_clt_20 * y_dens_clt_20 / y_dens_clt
    Pr_ct_01 <- Pr_clt_01 * y_dens_clt_01 / y_dens_clt
    Pr_ct_11 <- Pr_clt_11 * y_dens_clt_11 / y_dens_clt
    Pr_ct_21 <- Pr_clt_21 * y_dens_clt_21 / y_dens_clt
    Pr_ct_02 <- Pr_clt_02 * y_dens_clt_02 / y_dens_clt
    Pr_ct_12 <- Pr_clt_12 * y_dens_clt_12 / y_dens_clt
    Pr_ct_22 <- Pr_clt_22 * y_dens_clt_22 / y_dens_clt

    # Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
    # Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
    Pr_ct_0 <- max(Pr_ct_00 + Pr_ct_10 + Pr_ct_20, 1e-10)
    Pr_ct_1 <- max(Pr_ct_01 + Pr_ct_11 + Pr_ct_21, 1e-10)
    Pr_ct_2 <- max(Pr_ct_02 + Pr_ct_12 + Pr_ct_22, 1e-10)

    # Sum up log-likelihood over all iterations
    logLik_T <- logLik_T + loglik_t

    # Impose constraint that p,q > 0.9
    if (i == periods) {
      logLik_T <- ifelse(round(P[1, 1], 3) <= 0.9, 10000 * (1 - P[1, 1]) + logLik_T,
        ifelse(round(P[2, 2], 3) <= 0.9, 10000 * (1 - P[2, 2]) + logLik_T,
          logLik_T
        )
      )
    }

    #------------------------#
    ### Approximation part ###
    #------------------------#

    # Approximate updated values to break exponential growth of required values
    # (t|t, S_t)
    a_t_ct_0 <- (a_t_ct_00 * Pr_ct_00 + a_t_ct_10 * Pr_ct_10 + a_t_ct_20 * Pr_ct_20) / Pr_ct_0
    a_t_ct_1 <- (a_t_ct_01 * Pr_ct_01 + a_t_ct_11 * Pr_ct_11 + a_t_ct_21 * Pr_ct_21) / Pr_ct_1
    a_t_ct_2 <- (a_t_ct_02 * Pr_ct_02 + a_t_ct_12 * Pr_ct_12 + a_t_ct_22 * Pr_ct_22) / Pr_ct_2

    P_t_ct_0 <- ((Pr_ct_00 * (P_t_ct_00 + (a_t_ct_0 - a_t_ct_00) %*% Transp(a_t_ct_0 - a_t_ct_00))) + (Pr_ct_10 * (P_t_ct_10 + (a_t_ct_0 - a_t_ct_10) %*% Transp(a_t_ct_0 - a_t_ct_10))) + (Pr_ct_20 * (P_t_ct_20 + (a_t_ct_0 - a_t_ct_20) %*% Transp(a_t_ct_0 - a_t_ct_20)))) / Pr_ct_0
    P_t_ct_1 <- ((Pr_ct_01 * (P_t_ct_01 + (a_t_ct_1 - a_t_ct_01) %*% Transp(a_t_ct_1 - a_t_ct_01))) + (Pr_ct_11 * (P_t_ct_11 + (a_t_ct_1 - a_t_ct_11) %*% Transp(a_t_ct_1 - a_t_ct_11))) + (Pr_ct_21 * (P_t_ct_21 + (a_t_ct_1 - a_t_ct_21) %*% Transp(a_t_ct_1 - a_t_ct_21)))) / Pr_ct_1
    P_t_ct_2 <- ((Pr_ct_02 * (P_t_ct_02 + (a_t_ct_2 - a_t_ct_02) %*% Transp(a_t_ct_2 - a_t_ct_02))) + (Pr_ct_12 * (P_t_ct_12 + (a_t_ct_2 - a_t_ct_12) %*% Transp(a_t_ct_2 - a_t_ct_12))) + (Pr_ct_22 * (P_t_ct_22 + (a_t_ct_2 - a_t_ct_22) %*% Transp(a_t_ct_2 - a_t_ct_22)))) / Pr_ct_2

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Prediction step with approximated updates to complete loop
    # (t+1|t)
    # Regime 0 (high drift)
    a_t_clt_00 <- Tt %*% a_t_ct_0
    a_t_clt_10 <- Tt %*% a_t_ct_1
    a_t_clt_20 <- Tt %*% a_t_ct_2
    P_t_clt_00 <- Tt %*% P_t_ct_0 %*% TranspTt + expandQ
    P_t_clt_10 <- Tt %*% P_t_ct_1 %*% TranspTt + expandQ
    P_t_clt_20 <- Tt %*% P_t_ct_2 %*% TranspTt + expandQ

    # Regime 1 (low drift)
    a_t_clt_01 <- Tt %*% a_t_ct_0 + lambda
    a_t_clt_11 <- Tt %*% a_t_ct_1 + lambda
    a_t_clt_11 <- Tt %*% a_t_ct_2 + lambda
    P_t_clt_01 <- Tt %*% P_t_ct_0 %*% TranspTt + expandQ
    P_t_clt_11 <- Tt %*% P_t_ct_1 %*% TranspTt + expandQ
    P_t_clt_21 <- Tt %*% P_t_ct_2 %*% TranspTt + expandQ

    # Regime 2 (no drift)
    a_t_clt_02 <- Tt_noDrift %*% a_t_ct_0
    a_t_clt_12 <- Tt_noDrift %*% a_t_ct_1
    a_t_clt_22 <- Tt_noDrift %*% a_t_ct_2
    P_t_clt_02 <- Tt %*% P_t_ct_0 %*% TranspTt + expandQ
    P_t_clt_12 <- Tt %*% P_t_ct_1 %*% TranspTt + expandQ
    P_t_clt_22 <- Tt %*% P_t_ct_2 %*% TranspTt + expandQ
  }

  # Final log-likelihood value
  return(logLik_T)
}


#' Function that executes the Kim filter for a specific parameter vector for the three state model specification
#' @param par vector of constrained parameters
#' @param data vector of the data
#' @param Ini initial value for the trend component
#' @returns a list of the filter output

Filter_fctn <- function(par, data, Ini) {
  #---------------------------------------------------------------------------------------#
  # Initialize filter
  #---------------------------------------------------------------------------------------#

  # Load input (provided constraining function is already reversed)
  xi <- par[1] # Sd of innovations to drift component
  # omega <- par[2] # Sd of innovations to seasonal component
  omega <- 0
  # epsilon <- par[2] # Sd of measurement equation residual
  epsilon <- 0
  epsilonSq <- epsilon^2
  eta <- par[2]
  nu_1 <- par[3] # Additional drift for down-turning regime 1
  P <- matrix(NA, 3, 3)
  P[, 1] <- c(par[4], par[5], 1 - par[4] - par[5]) # P(St =  | St-1 = 0)
  P[, 2] <- c(par[6], par[7], 1 - par[6] - par[7]) # P(St =  | St-1 = 1)
  P[, 3] <- c(par[8], par[9], 1 - par[8] - par[9]) # P(St =  | St-1 = 1)

  phi_1 <- par[10]
  phi_2 <- par[11]
  
  # Transition matrix
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  Tt_cycle <- matrix(c(phi_1, phi_2, 1, 0), byrow = T, nc = 2)
  
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  # Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  # Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_seasUr))
  Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_cycle))
  # allows the seasonal unit root to load on gamma_t+1
  # Tt[3, 9] <- 1
  Tt_noDrift <- Tt
  Tt_noDrift[1, 2] <- 0
  TranspTt <- Transp(Tt)

  # Dimension of state vector
  Dimens <- NCOL(Tt)

  # Number of time periods
  periods <- NROW(data)

  # Measurement matrix
  # Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)
  Z <- matrix(c(1, 0, 1, 0, 1, rep(0, Dimens - 5)), 1, Dimens)
  
  TranspZ <- Transp(Z)

  # Initial unconditional regime probabilities (steady state unconditional probabilities)
  A <- matrix(rbind(diag(3) - P, rep(1, 3)), nc = 3)
  Pi <- Inverse(Transp(A) %*% A) %*% t(A) %*% matrix(diag(4)[, 4], ncol = 1)
  Pr_ct_0 <- Pi[1]
  Pr_ct_1 <- Pi[2]
  Pr_ct_2 <- Pi[3]

  # Var-Cov Matrix of innovations in transition equation
  # Q <- diag(c(xi^2, omega^2))
  
  # Var-Cov Matrix of innovations in transition equation
  Q <- diag(c(xi^2, eta^2, omega^2))
  
  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, 2), 1, rep(0, Dimens - 3)),
    c(rep(0, Dimens))
  )

  # Matrix to expand Q so that it matches Var matrix of state vector
  # R <- cbind(
  #   c(1, rep(0, Dimens - 1)),
  #   c(rep(0, 8), 1, rep(0, Dimens - 9))
  # )
  # R <- cbind(
  #   c(1, rep(0, Dimens - 1)),
  #   c(rep(0, Dimens))
  # )
  transpR <- Transp(R)
  # Var-Cov Matrix of innovations in transition equation
  expandQ <- R %*% Q %*% transpR

  # Vector holding additional Drift in Regime 1
  lambda <- matrix(c(nu_1, rep(0, Dimens - 1)), Dimens)

  # Initial values for state vector
  a_t_clt_00 <- a_t_clt_10 <- a_t_clt_20 <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens)
  a_t_clt_01 <- a_t_clt_11 <- a_t_clt_21 <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens) + lambda
  a_t_clt_02 <- a_t_clt_12 <- a_t_clt_22 <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens)

  # Initial variance matrix for state vector (Exact diffuse initialization)
  P_t_clt_00 <- P_t_clt_10 <- P_t_clt_20 <- P_t_clt_01 <- P_t_clt_11 <- P_t_clt_21 <- P_t_clt_02 <- P_t_clt_12 <- P_t_clt_22 <- diag(Dimens) * 10000

  # Initializes output objects
  a_clt_array <- array(0, c(Dimens, 9, periods))
  P_clt_array <- array(0, c(Dimens, Dimens * 9, periods))
  a_ct_array <- array(0, c(Dimens, 3, periods))
  P_ct_array <- array(0, c(Dimens, Dimens * 3, periods))
  Prob_t_clt <- matrix(0, periods, 3)
  Prob_t_ct <- matrix(0, periods, 3)
  v_ct_mat <- matrix(0, periods, 9)
  F_ct_mat <- matrix(0, periods, 9)

  #---------------------------------------------------------------------------------------#
  # Kim filter to get filtered values
  #---------------------------------------------------------------------------------------#

  for (i in 1:periods) {
    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # One step ahead prediction error with state vector prediction from
    # (t|t-1)
    v_t_00 <- as.numeric(data[i] - Z %*% a_t_clt_00)
    v_t_10 <- as.numeric(data[i] - Z %*% a_t_clt_10)
    v_t_20 <- as.numeric(data[i] - Z %*% a_t_clt_20)
    v_t_01 <- as.numeric(data[i] - Z %*% a_t_clt_01)
    v_t_11 <- as.numeric(data[i] - Z %*% a_t_clt_11)
    v_t_21 <- as.numeric(data[i] - Z %*% a_t_clt_21)
    v_t_02 <- as.numeric(data[i] - Z %*% a_t_clt_02)
    v_t_12 <- as.numeric(data[i] - Z %*% a_t_clt_12)
    v_t_22 <- as.numeric(data[i] - Z %*% a_t_clt_22)

    # Variance of prediction error
    F_t_00 <- as.numeric(Z %*% P_t_clt_00 %*% TranspZ + epsilonSq)
    F_t_10 <- as.numeric(Z %*% P_t_clt_10 %*% TranspZ + epsilonSq)
    F_t_20 <- as.numeric(Z %*% P_t_clt_20 %*% TranspZ + epsilonSq)
    F_t_01 <- as.numeric(Z %*% P_t_clt_01 %*% TranspZ + epsilonSq)
    F_t_11 <- as.numeric(Z %*% P_t_clt_11 %*% TranspZ + epsilonSq)
    F_t_21 <- as.numeric(Z %*% P_t_clt_21 %*% TranspZ + epsilonSq)
    F_t_02 <- as.numeric(Z %*% P_t_clt_02 %*% TranspZ + epsilonSq)
    F_t_12 <- as.numeric(Z %*% P_t_clt_12 %*% TranspZ + epsilonSq)
    F_t_22 <- as.numeric(Z %*% P_t_clt_22 %*% TranspZ + epsilonSq)

    # Updating step
    # (t|t, S_t, S_t-1)
    # State vector
    a_t_ct_00 <- a_t_clt_00 + P_t_clt_00 %*% TranspZ %*% (1 / F_t_00 * v_t_00)
    a_t_ct_10 <- a_t_clt_10 + P_t_clt_10 %*% TranspZ %*% (1 / F_t_10 * v_t_10)
    a_t_ct_20 <- a_t_clt_20 + P_t_clt_20 %*% TranspZ %*% (1 / F_t_20 * v_t_20)
    a_t_ct_01 <- a_t_clt_01 + P_t_clt_01 %*% TranspZ %*% (1 / F_t_01 * v_t_01)
    a_t_ct_11 <- a_t_clt_11 + P_t_clt_11 %*% TranspZ %*% (1 / F_t_11 * v_t_11)
    a_t_ct_21 <- a_t_clt_21 + P_t_clt_21 %*% TranspZ %*% (1 / F_t_21 * v_t_21)
    a_t_ct_02 <- a_t_clt_02 + P_t_clt_02 %*% TranspZ %*% (1 / F_t_02 * v_t_02)
    a_t_ct_12 <- a_t_clt_12 + P_t_clt_12 %*% TranspZ %*% (1 / F_t_12 * v_t_12)
    a_t_ct_22 <- a_t_clt_22 + P_t_clt_22 %*% TranspZ %*% (1 / F_t_22 * v_t_22)

    # State vector Var-Matrix
    P_t_ct_00 <- P_t_clt_00 - P_t_clt_00 %*% TranspZ %*% (1 / F_t_00) %*% Z %*% P_t_clt_00
    P_t_ct_10 <- P_t_clt_10 - P_t_clt_10 %*% TranspZ %*% (1 / F_t_10) %*% Z %*% P_t_clt_10
    P_t_ct_20 <- P_t_clt_20 - P_t_clt_20 %*% TranspZ %*% (1 / F_t_20) %*% Z %*% P_t_clt_20
    P_t_ct_01 <- P_t_clt_01 - P_t_clt_01 %*% TranspZ %*% (1 / F_t_01) %*% Z %*% P_t_clt_01
    P_t_ct_11 <- P_t_clt_11 - P_t_clt_11 %*% TranspZ %*% (1 / F_t_11) %*% Z %*% P_t_clt_11
    P_t_ct_21 <- P_t_clt_21 - P_t_clt_21 %*% TranspZ %*% (1 / F_t_21) %*% Z %*% P_t_clt_21
    P_t_ct_02 <- P_t_clt_02 - P_t_clt_02 %*% TranspZ %*% (1 / F_t_02) %*% Z %*% P_t_clt_02
    P_t_ct_12 <- P_t_clt_12 - P_t_clt_12 %*% TranspZ %*% (1 / F_t_12) %*% Z %*% P_t_clt_12
    P_t_ct_22 <- P_t_clt_22 - P_t_clt_22 %*% TranspZ %*% (1 / F_t_22) %*% Z %*% P_t_clt_22

    # Store prediction errors and variances
    v_ct_mat[i, 1:9] <- c(v_t_00, v_t_10, v_t_20, v_t_01, v_t_11, v_t_21, v_t_02, v_t_12, v_t_22)
    F_ct_mat[i, 1:9] <- c(F_t_00, F_t_10, F_t_20, F_t_01, F_t_11, F_t_21, F_t_02, F_t_12, F_t_22)

    #-------------------#
    ### Hamilton part ###
    #-------------------#

    # Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
    # (t|t-1)
    Pr_clt_00 <- P[1, 1] * Pr_ct_0
    Pr_clt_10 <- P[1, 2] * Pr_ct_1
    Pr_clt_20 <- P[1, 3] * Pr_ct_2
    Pr_clt_01 <- P[2, 1] * Pr_ct_0
    Pr_clt_11 <- P[2, 2] * Pr_ct_1
    Pr_clt_21 <- P[2, 3] * Pr_ct_2
    Pr_clt_02 <- P[3, 1] * Pr_ct_0
    Pr_clt_12 <- P[3, 2] * Pr_ct_1
    Pr_clt_22 <- P[3, 3] * Pr_ct_2

    # Pr of S_t = j conditional on t = t-1 is given by summarizing over all regimes at t = t-1
    # (t|t-1)
    Pr_clt_0 <- Pr_clt_00 + Pr_clt_10 + Pr_clt_20
    Pr_clt_1 <- Pr_clt_01 + Pr_clt_11 + Pr_clt_21
    Pr_clt_2 <- Pr_clt_02 + Pr_clt_12 + Pr_clt_22

    # Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
    y_dens_clt_00 <- dnorm(v_t_00, sd = sqrt(F_t_00))
    y_dens_clt_10 <- dnorm(v_t_10, sd = sqrt(F_t_10))
    y_dens_clt_20 <- dnorm(v_t_20, sd = sqrt(F_t_20))
    y_dens_clt_01 <- dnorm(v_t_01, sd = sqrt(F_t_01))
    y_dens_clt_11 <- dnorm(v_t_11, sd = sqrt(F_t_11))
    y_dens_clt_21 <- dnorm(v_t_21, sd = sqrt(F_t_21))
    y_dens_clt_02 <- dnorm(v_t_02, sd = sqrt(F_t_02))
    y_dens_clt_12 <- dnorm(v_t_12, sd = sqrt(F_t_12))
    y_dens_clt_22 <- dnorm(v_t_22, sd = sqrt(F_t_22))

    # Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
    # (t|t-1))
    y_dens_clt <- sum(
      Pr_clt_00 * y_dens_clt_00, Pr_clt_10 * y_dens_clt_10, Pr_clt_20 * y_dens_clt_20,
      Pr_clt_01 * y_dens_clt_01, Pr_clt_11 * y_dens_clt_11, Pr_clt_21 * y_dens_clt_21,
      Pr_clt_02 * y_dens_clt_02, Pr_clt_12 * y_dens_clt_12, Pr_clt_22 * y_dens_clt_22
    )

    # Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
    # (t|t)
    Pr_ct_00 <- Pr_clt_00 * y_dens_clt_00 / y_dens_clt
    Pr_ct_10 <- Pr_clt_10 * y_dens_clt_10 / y_dens_clt
    Pr_ct_20 <- Pr_clt_20 * y_dens_clt_20 / y_dens_clt
    Pr_ct_01 <- Pr_clt_01 * y_dens_clt_01 / y_dens_clt
    Pr_ct_11 <- Pr_clt_11 * y_dens_clt_11 / y_dens_clt
    Pr_ct_21 <- Pr_clt_21 * y_dens_clt_21 / y_dens_clt
    Pr_ct_02 <- Pr_clt_02 * y_dens_clt_02 / y_dens_clt
    Pr_ct_12 <- Pr_clt_12 * y_dens_clt_12 / y_dens_clt
    Pr_ct_22 <- Pr_clt_22 * y_dens_clt_22 / y_dens_clt

    # Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
    # Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
    Pr_ct_0 <- max(Pr_ct_00 + Pr_ct_10 + Pr_ct_20, 1e-10)
    Pr_ct_1 <- max(Pr_ct_01 + Pr_ct_11 + Pr_ct_21, 1e-10)
    Pr_ct_2 <- max(Pr_ct_02 + Pr_ct_12 + Pr_ct_22, 1e-10)

    # Records predicted probabilities
    # (t|t-1)
    Prob_t_clt[i, 1:3] <- c(Pr_clt_0, Pr_clt_1, Pr_clt_2)

    # Records updated probability
    # (t|t)
    Prob_t_ct[i, 1:3] <- c(Pr_ct_0, Pr_ct_1, Pr_ct_2)

    #------------------------#
    ### Approximation part ###
    #------------------------#

    # Approximate updated values to break exponential growth of required values
    # (t|t, S_t)
    a_t_ct_0 <- (a_t_ct_00 * Pr_ct_00 + a_t_ct_10 * Pr_ct_10 + a_t_ct_20 * Pr_ct_20) / Pr_ct_0
    a_t_ct_1 <- (a_t_ct_01 * Pr_ct_01 + a_t_ct_11 * Pr_ct_11 + a_t_ct_21 * Pr_ct_21) / Pr_ct_1
    a_t_ct_2 <- (a_t_ct_02 * Pr_ct_02 + a_t_ct_12 * Pr_ct_12 + a_t_ct_22 * Pr_ct_22) / Pr_ct_2

    P_t_ct_0 <- ((Pr_ct_00 * (P_t_ct_00 + (a_t_ct_0 - a_t_ct_00) %*% Transp(a_t_ct_0 - a_t_ct_00))) + (Pr_ct_10 * (P_t_ct_10 + (a_t_ct_0 - a_t_ct_10) %*% Transp(a_t_ct_0 - a_t_ct_10))) + (Pr_ct_20 * (P_t_ct_20 + (a_t_ct_0 - a_t_ct_20) %*% Transp(a_t_ct_0 - a_t_ct_20)))) / Pr_ct_0
    P_t_ct_1 <- ((Pr_ct_01 * (P_t_ct_01 + (a_t_ct_1 - a_t_ct_01) %*% Transp(a_t_ct_1 - a_t_ct_01))) + (Pr_ct_11 * (P_t_ct_11 + (a_t_ct_1 - a_t_ct_11) %*% Transp(a_t_ct_1 - a_t_ct_11))) + (Pr_ct_21 * (P_t_ct_21 + (a_t_ct_1 - a_t_ct_21) %*% Transp(a_t_ct_1 - a_t_ct_21)))) / Pr_ct_1
    P_t_ct_2 <- ((Pr_ct_02 * (P_t_ct_02 + (a_t_ct_2 - a_t_ct_02) %*% Transp(a_t_ct_2 - a_t_ct_02))) + (Pr_ct_12 * (P_t_ct_12 + (a_t_ct_2 - a_t_ct_12) %*% Transp(a_t_ct_2 - a_t_ct_12))) + (Pr_ct_22 * (P_t_ct_22 + (a_t_ct_2 - a_t_ct_22) %*% Transp(a_t_ct_2 - a_t_ct_22)))) / Pr_ct_2

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Store predicted values
    # (t|t-1)
    a_clt_array[1:Dimens, 1:9, i] <- c(
      a_t_clt_00, a_t_clt_10, a_t_clt_20,
      a_t_clt_01, a_t_clt_11, a_t_clt_21,
      a_t_clt_02, a_t_clt_12, a_t_clt_22
    )

    P_clt_array[1:Dimens, 1:(Dimens * 9), i] <- c(
      P_t_clt_00, P_t_clt_10, P_t_clt_20,
      P_t_clt_01, P_t_clt_11, P_t_clt_21,
      P_t_clt_02, P_t_clt_12, P_t_clt_22
    )

    # Store updated values
    # (t|t)
    a_ct_array[1:Dimens, 1:3, i] <- c(a_t_ct_0, a_t_ct_1, a_t_ct_2)
    P_ct_array[1:Dimens, 1:(Dimens * 3), i] <- c(P_t_ct_0, P_t_ct_1, P_t_ct_2)

    # Prediction step with approximated updates to complete loop
    # (t+1|t)
    # Regime 0 (high drift)
    a_t_clt_00 <- Tt %*% a_t_ct_0
    a_t_clt_10 <- Tt %*% a_t_ct_1
    a_t_clt_20 <- Tt %*% a_t_ct_2
    P_t_clt_00 <- Tt %*% P_t_ct_0 %*% TranspTt + expandQ + 1e-4
    P_t_clt_10 <- Tt %*% P_t_ct_1 %*% TranspTt + expandQ + 1e-4
    P_t_clt_20 <- Tt %*% P_t_ct_2 %*% TranspTt + expandQ + 1e-4

    # Regime 1 (low drift)
    a_t_clt_01 <- Tt %*% a_t_ct_0 + lambda
    a_t_clt_11 <- Tt %*% a_t_ct_1 + lambda
    a_t_clt_11 <- Tt %*% a_t_ct_2 + lambda
    P_t_clt_01 <- Tt %*% P_t_ct_0 %*% TranspTt + expandQ + 1e-4
    P_t_clt_11 <- Tt %*% P_t_ct_1 %*% TranspTt + expandQ + 1e-4
    P_t_clt_21 <- Tt %*% P_t_ct_2 %*% TranspTt + expandQ + 1e-4

    # Regime 2 (no drift)
    a_t_clt_02 <- Tt_noDrift %*% a_t_ct_0
    a_t_clt_12 <- Tt_noDrift %*% a_t_ct_1
    a_t_clt_22 <- Tt_noDrift %*% a_t_ct_2
    P_t_clt_02 <- Tt %*% P_t_ct_0 %*% TranspTt + expandQ + 1e-4
    P_t_clt_12 <- Tt %*% P_t_ct_1 %*% TranspTt + expandQ + 1e-4
    P_t_clt_22 <- Tt %*% P_t_ct_2 %*% TranspTt + expandQ + 1e-4
  }

  # Output of filtered values
outputList <- list(
  a_clt = a_clt_array, # predicted state vector
  P_clt = P_clt_array, # predicted state vector var
  a_ct = a_ct_array, # updated state vector
  P_ct = P_ct_array, # updated state vector var
  Prob_clt = Prob_t_clt, # predicted regime probs
  Prob_ct = Prob_t_ct, # updated regime probs
  Pred_err = v_ct_mat, # one-step-ahead prediction error
  Pred_err_Var = F_ct_mat # Var of pred error
)
  return(outputList)
}


#' Function that executes the Kim filter for a specific parameter vector for the three state model specification
#' @param par vector of constrained parameters
#' @param filterOutput Output from the Filter_fctn
#' @param ciInterval for the smoothed trend
#' @returns a list of the smoother output

Smoother_fctn <- function(par, filterOutput, ciInterval) {
  #---------------------------------------------------------------------------------------#
  # Initialize smoother
  #---------------------------------------------------------------------------------------#

  
  # Provide system matrix
  Tt_trend <- matrix(c(1, 1, 0, 1), byrow = T, nc = 2)
  Tt_seasDeterm <- rbind(rep(-1, 6), cbind(diag(5), rep(0, 5)))
  phiVec <- par[10:11]
  Tt_cycle <- matrix(c(phiVec, 1, 0), byrow = T, nc = 2)
  
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  # Tt_seasUr <- rbind(c(rep(0, 6), 1), cbind(diag(6), rep(0, 6)))
  # Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_seasUr))
  # allows the seasonal unit root to load on gamma_t+1
  # Tt[3, 9] <- 1
  Tt <- as.matrix(bdiag(Tt_trend, Tt_seasDeterm, Tt_cycle))
  TranspTt <- Transp(Tt)

  # Dimension of state vector
  Dimens <- NCOL(Tt)

  # Number of time periods
  periods <- length(filterOutput$a_clt[1, 1, ])

  # Provide transition probabilities
  P <- matrix(NA, 3, 3)
  P[, 1] <- c(par[4], par[5], 1 - par[4] - par[5]) # P(St =  | St-1 = 0)
  P[, 2] <- c(par[6], par[7], 1 - par[6] - par[7]) # P(St =  | St-1 = 1)
  P[, 3] <- c(par[8], par[9], 1 - par[8] - par[9]) # P(St =  | St-1 = 1)

  # Realizations for t = T of Kim filter and Kim smoother are identical
  Pr_cT_0 <- filterOutput[["Prob_ct"]][periods, 1]
  Pr_cT_1 <- filterOutput[["Prob_ct"]][periods, 2]
  Pr_cT_2 <- filterOutput[["Prob_ct"]][periods, 3]

  # Initialize smoother output
  a_cT_array <- array(0, c(Dimens, 4, periods))
  P_cT_array <- array(0, c(Dimens, Dimens * 4, periods))
  Pr_cT_mat <- matrix(0, periods, 3)
  CI_a_cT_array <- array(0, c(Dimens, 2, periods))

  a_t1_cT_0 <- a_cT_array[, 1, periods] <- filterOutput[["a_ct"]][, 1, periods] %>% as.matrix()
  a_t1_cT_1 <- a_cT_array[, 2, periods] <- filterOutput[["a_ct"]][, 2, periods] %>% as.matrix()
  a_t1_cT_2 <- a_cT_array[, 3, periods] <- filterOutput[["a_ct"]][, 3, periods] %>% as.matrix()
  a_cT_array[, 3, periods] <- Pr_cT_0 * a_cT_array[, 1, periods] + Pr_cT_1 * a_cT_array[, 2, periods] + Pr_cT_2 * a_cT_array[, 3, periods]

  P_t1_cT_0 <- P_cT_array[, 1:Dimens, periods] <- filterOutput[["P_ct"]][, 1:Dimens, periods]
  P_t1_cT_1 <- P_cT_array[, (Dimens + 1):(2 * Dimens), periods] <- filterOutput[["P_ct"]][, (Dimens + 1):(2 * Dimens), periods]
  P_t1_cT_2 <- P_cT_array[, (2 * Dimens + 1):(3 * Dimens), periods] <- filterOutput[["P_ct"]][, (2 * Dimens + 1):(3 * Dimens), periods]
  P_cT_array[, (3 * Dimens + 1):(4 * Dimens), periods] <- sum(Pr_cT_0 * filterOutput[["P_ct"]][, 1:Dimens, periods],
                                                              Pr_cT_1 * filterOutput[["P_ct"]][, (Dimens + 1):(2 * Dimens), periods],
                                                              Pr_cT_2 * filterOutput[["P_ct"]][, (2 * Dimens + 1):(3 * Dimens), periods])
  
  #---------------------------------------------------------------------------------------#
  # Kim smoother (backwards) iteration
  #---------------------------------------------------------------------------------------#

  for (i in (periods - 1):5) {
    # Load necessary filtered values for each iteration
    # Updated probabilities
    # (t|t)
    Pr_ct_0 <- filterOutput[["Prob_ct"]][i, 1]
    Pr_ct_1 <- filterOutput[["Prob_ct"]][i, 2]
    Pr_ct_2 <- filterOutput[["Prob_ct"]][i, 3]

    # Filtered probabilities
    # (t+1|t)
    Pr_clt_0 <- filterOutput[["Prob_clt"]][i + 1, 1] # Note: value read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    Pr_clt_1 <- filterOutput[["Prob_clt"]][i + 1, 2]
    Pr_clt_2 <- filterOutput[["Prob_clt"]][i + 1, 3]

    # Updated (approximated) state vector
    # (t|t)
    a_t_ct_0 <- as.matrix(filterOutput[["a_ct"]][, 1, i])
    a_t_ct_1 <- as.matrix(filterOutput[["a_ct"]][, 2, i])
    a_t_ct_2 <- as.matrix(filterOutput[["a_ct"]][, 3, i])

    # Filtered state vector
    # (t+1|t)
    a_t1_ct_00 <- as.matrix(filterOutput[["a_clt"]][, 1, i + 1]) # Note: values read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    a_t1_ct_10 <- as.matrix(filterOutput[["a_clt"]][, 2, i + 1])
    a_t1_ct_20 <- as.matrix(filterOutput[["a_clt"]][, 3, i + 1])
    a_t1_ct_01 <- as.matrix(filterOutput[["a_clt"]][, 4, i + 1])
    a_t1_ct_11 <- as.matrix(filterOutput[["a_clt"]][, 5, i + 1])
    a_t1_ct_21 <- as.matrix(filterOutput[["a_clt"]][, 6, i + 1])
    a_t1_ct_02 <- as.matrix(filterOutput[["a_clt"]][, 7, i + 1])
    a_t1_ct_12 <- as.matrix(filterOutput[["a_clt"]][, 8, i + 1])
    a_t1_ct_22 <- as.matrix(filterOutput[["a_clt"]][, 9, i + 1])

    # Updated state vector Var-Matrix
    # (t|t)
    P_t_ct_0 <- filterOutput[["P_ct"]][, 1:Dimens, i]
    P_t_ct_1 <- filterOutput[["P_ct"]][, (Dimens + 1):(2 * Dimens), i]
    P_t_ct_2 <- filterOutput[["P_ct"]][, (2 * Dimens + 1):(3 * Dimens), i]

    # Filtered state vector Var-Matrix
    # (t+1|t)
    P_t1_ct_00 <- filterOutput[["P_clt"]][, 1:Dimens, i + 1] # Note: values read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    P_t1_ct_10 <- filterOutput[["P_clt"]][, (Dimens + 1):(2 * Dimens), i + 1]
    P_t1_ct_20 <- filterOutput[["P_clt"]][, (2 * Dimens + 1):(3 * Dimens), i + 1]
    P_t1_ct_01 <- filterOutput[["P_clt"]][, (3 * Dimens + 1):(4 * Dimens), i + 1]
    P_t1_ct_11 <- filterOutput[["P_clt"]][, (4 * Dimens + 1):(5 * Dimens), i + 1]
    P_t1_ct_21 <- filterOutput[["P_clt"]][, (5 * Dimens + 1):(6 * Dimens), i + 1]
    P_t1_ct_02 <- filterOutput[["P_clt"]][, (6 * Dimens + 1):(7 * Dimens), i + 1]
    P_t1_ct_12 <- filterOutput[["P_clt"]][, (7 * Dimens + 1):(8 * Dimens), i + 1]
    P_t1_ct_22 <- filterOutput[["P_clt"]][, (8 * Dimens + 1):(9 * Dimens), i + 1]

    #-------------------#
    ### Hamilton part ###
    #-------------------#

    # Compute Pr(S_t = j, S_t+1 = k) conditional on all information
    # Here Pr_cT_jk: _cT_01 implies Pr(S_t+1 = 1, S_t = 0 | T) // NOTE: probabilities conditional on t retain prev. notation --> _ct_01: conditional on t Pr(S_t = 1, S_t-1 = 0 | t)
    # (t|T)
    Pr_cT_00 <- (Pr_cT_0 * Pr_ct_0 * P[1, 1]) / Pr_clt_0
    Pr_cT_10 <- (Pr_cT_0 * Pr_ct_1 * P[1, 2]) / Pr_clt_0
    Pr_cT_20 <- (Pr_cT_0 * Pr_ct_2 * P[1, 3]) / Pr_clt_0
    Pr_cT_01 <- (Pr_cT_1 * Pr_ct_0 * P[2, 1]) / Pr_clt_1
    Pr_cT_11 <- (Pr_cT_1 * Pr_ct_1 * P[2, 2]) / Pr_clt_1
    Pr_cT_21 <- (Pr_cT_1 * Pr_ct_2 * P[2, 3]) / Pr_clt_1
    Pr_cT_02 <- (Pr_cT_2 * Pr_ct_0 * P[3, 1]) / Pr_clt_2
    Pr_cT_12 <- (Pr_cT_2 * Pr_ct_1 * P[3, 2]) / Pr_clt_2
    Pr_cT_22 <- (Pr_cT_2 * Pr_ct_2 * P[3, 3]) / Pr_clt_2

    # Compute Pr(S_t = j | all information) // NOTE: To derive at unconditional probability of Pr(S_t) on has to sum over all regimes at t = t + 1 (in difference to t = t - 1 in the filter)
    # (t|T)
    Pr_cT_0 <- max(Pr_cT_00 + Pr_cT_01 + Pr_cT_02, 1e-10)
    Pr_cT_1 <- max(Pr_cT_10 + Pr_cT_11 + Pr_cT_12, 1e-10)
    Pr_cT_2 <- max(Pr_cT_20 + Pr_cT_21 + Pr_cT_22, 1e-10)

    # Record smoothed probabilities
    Pr_cT_mat[i, 1:3] <- c(Pr_cT_0, Pr_cT_1, Pr_cT_2)

    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # Derive smoothed State vector. NOTE: again _cT_01: S_t = 0 and S_t+1 = 1 // _ct_01: S_t-1 = 0 and S_t = 1
    # (t|T)
    # P"help"_ij
    # (t|T)
    Ph_00 <- P_t_ct_0 %*% TranspTt %*% Inverse(P_t1_ct_00)
    Ph_10 <- P_t_ct_1 %*% TranspTt %*% Inverse(P_t1_ct_10)
    Ph_20 <- P_t_ct_2 %*% TranspTt %*% Inverse(P_t1_ct_20)
    Ph_01 <- P_t_ct_0 %*% TranspTt %*% Inverse(P_t1_ct_01)
    Ph_11 <- P_t_ct_1 %*% TranspTt %*% Inverse(P_t1_ct_11)
    Ph_21 <- P_t_ct_2 %*% TranspTt %*% Inverse(P_t1_ct_21)
    Ph_02 <- P_t_ct_0 %*% TranspTt %*% Inverse(P_t1_ct_02)
    Ph_12 <- P_t_ct_1 %*% TranspTt %*% Inverse(P_t1_ct_12)
    Ph_22 <- P_t_ct_2 %*% TranspTt %*% Inverse(P_t1_ct_22)

    # State vector
    a_t_cT_00 <- a_t_ct_0 + Ph_00 %*% (a_t1_cT_0 - a_t1_ct_00)
    a_t_cT_10 <- a_t_ct_1 + Ph_10 %*% (a_t1_cT_0 - a_t1_ct_10)
    a_t_cT_20 <- a_t_ct_2 + Ph_20 %*% (a_t1_cT_0 - a_t1_ct_20)
    a_t_cT_01 <- a_t_ct_0 + Ph_01 %*% (a_t1_cT_1 - a_t1_ct_01)
    a_t_cT_11 <- a_t_ct_1 + Ph_11 %*% (a_t1_cT_1 - a_t1_ct_11)
    a_t_cT_21 <- a_t_ct_2 + Ph_21 %*% (a_t1_cT_1 - a_t1_ct_21)
    a_t_cT_02 <- a_t_ct_0 + Ph_02 %*% (a_t1_cT_2 - a_t1_ct_02)
    a_t_cT_12 <- a_t_ct_1 + Ph_12 %*% (a_t1_cT_2 - a_t1_ct_12)
    a_t_cT_22 <- a_t_ct_2 + Ph_22 %*% (a_t1_cT_2 - a_t1_ct_22)

    # State vector Var-Matrix
    P_t_cT_00 <- P_t_ct_0 + Ph_00 %*% (P_t1_cT_0 - P_t1_ct_00) %*% Transp(Ph_00)
    P_t_cT_10 <- P_t_ct_1 + Ph_10 %*% (P_t1_cT_0 - P_t1_ct_10) %*% Transp(Ph_10)
    P_t_cT_20 <- P_t_ct_2 + Ph_20 %*% (P_t1_cT_0 - P_t1_ct_20) %*% Transp(Ph_20)

    P_t_cT_01 <- P_t_ct_0 + Ph_01 %*% (P_t1_cT_1 - P_t1_ct_01) %*% Transp(Ph_01)
    P_t_cT_11 <- P_t_ct_1 + Ph_11 %*% (P_t1_cT_1 - P_t1_ct_11) %*% Transp(Ph_11)
    P_t_cT_21 <- P_t_ct_2 + Ph_21 %*% (P_t1_cT_1 - P_t1_ct_21) %*% Transp(Ph_21)

    P_t_cT_02 <- P_t_ct_0 + Ph_02 %*% (P_t1_cT_2 - P_t1_ct_02) %*% Transp(Ph_02)
    P_t_cT_12 <- P_t_ct_1 + Ph_12 %*% (P_t1_cT_2 - P_t1_ct_12) %*% Transp(Ph_12)
    P_t_cT_22 <- P_t_ct_2 + Ph_22 %*% (P_t1_cT_2 - P_t1_ct_22) %*% Transp(Ph_22)


    # Confidence intervals for smoothed state vector
    # Upper
    u_a_t_cT_00 <- a_t_cT_00 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_00[1, 1])
    u_a_t_cT_10 <- a_t_cT_10 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_10[1, 1])
    u_a_t_cT_20 <- a_t_cT_20 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_20[1, 1])
    u_a_t_cT_01 <- a_t_cT_01 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_01[1, 1])
    u_a_t_cT_11 <- a_t_cT_11 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_11[1, 1])
    u_a_t_cT_21 <- a_t_cT_21 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_11[1, 1])
    u_a_t_cT_02 <- a_t_cT_02 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_02[1, 1])
    u_a_t_cT_12 <- a_t_cT_12 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_12[1, 1])
    u_a_t_cT_22 <- a_t_cT_22 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_12[1, 1])
    
    # Lower
    l_a_t_cT_00 <- a_t_cT_00 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_00[1, 1])
    l_a_t_cT_10 <- a_t_cT_10 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_10[1, 1])
    l_a_t_cT_20 <- a_t_cT_20 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_20[1, 1])
    l_a_t_cT_01 <- a_t_cT_01 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_01[1, 1])
    l_a_t_cT_11 <- a_t_cT_11 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_11[1, 1])
    l_a_t_cT_21 <- a_t_cT_21 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_21[1, 1])
    l_a_t_cT_02 <- a_t_cT_02 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_02[1, 1])
    l_a_t_cT_12 <- a_t_cT_12 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_12[1, 1])
    l_a_t_cT_22 <- a_t_cT_22 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_22[1, 1])
    
    #------------------------#
    ### Approximation part ###
    #------------------------#

    # Approximate smoothed values, so that only conditional on current regime
    # Note: In difference to the filter, probabilities are averaged over all possible regimes at t = t + 1 (not t = t - 1)
    # (t|T)
    # State vector
    a_t_cT_0 <- ((Pr_cT_00 * a_t_cT_00) + (Pr_cT_01 * a_t_cT_01) + (Pr_cT_02 * a_t_cT_02)) / Pr_cT_0
    a_t_cT_1 <- ((Pr_cT_10 * a_t_cT_10) + (Pr_cT_11 * a_t_cT_11) + (Pr_cT_12 * a_t_cT_12)) / Pr_cT_1
    a_t_cT_2 <- ((Pr_cT_20 * a_t_cT_20) + (Pr_cT_21 * a_t_cT_21) + (Pr_cT_22 * a_t_cT_22)) / Pr_cT_2

    # Also approximate confidence intervals
    u_a_t_cT_0 <- ((Pr_cT_00 * u_a_t_cT_00) + (Pr_cT_01 * u_a_t_cT_01) + (Pr_cT_02 * u_a_t_cT_02)) / Pr_cT_0
    u_a_t_cT_1 <- ((Pr_cT_10 * u_a_t_cT_10) + (Pr_cT_11 * u_a_t_cT_11) + (Pr_cT_12 * u_a_t_cT_12)) / Pr_cT_1
    u_a_t_cT_2 <- ((Pr_cT_20 * u_a_t_cT_20) + (Pr_cT_21 * u_a_t_cT_21) + (Pr_cT_22 * u_a_t_cT_22)) / Pr_cT_2

    l_a_t_cT_0 <- ((Pr_cT_00 * l_a_t_cT_00) + (Pr_cT_01 * l_a_t_cT_01) + (Pr_cT_02 * l_a_t_cT_02)) / Pr_cT_0
    l_a_t_cT_1 <- ((Pr_cT_10 * l_a_t_cT_10) + (Pr_cT_11 * l_a_t_cT_11) + (Pr_cT_12 * l_a_t_cT_12)) / Pr_cT_1
    l_a_t_cT_2 <- ((Pr_cT_20 * l_a_t_cT_20) + (Pr_cT_21 * l_a_t_cT_21) + (Pr_cT_22 * l_a_t_cT_22)) / Pr_cT_2

    # Same for State Var-Matrix
    P_t_cT_0 <- (Pr_cT_00 * (P_t_cT_00 + (a_t_cT_0 - a_t_cT_00) %*% Transp(a_t_cT_0 - a_t_cT_00)) + Pr_cT_01 * (P_t_cT_01 + (a_t_cT_0 - a_t_cT_01) %*% Transp(a_t_cT_0 - a_t_cT_01)) + Pr_cT_02 * (P_t_cT_02 + (a_t_cT_0 - a_t_cT_02) %*% Transp(a_t_cT_0 - a_t_cT_02))) / Pr_cT_0
    P_t_cT_1 <- (Pr_cT_10 * (P_t_cT_10 + (a_t_cT_1 - a_t_cT_10) %*% Transp(a_t_cT_1 - a_t_cT_10)) + Pr_cT_11 * (P_t_cT_11 + (a_t_cT_1 - a_t_cT_11) %*% Transp(a_t_cT_1 - a_t_cT_11)) + Pr_cT_12 * (P_t_cT_12 + (a_t_cT_1 - a_t_cT_12) %*% Transp(a_t_cT_1 - a_t_cT_12))) / Pr_cT_1
    P_t_cT_2 <- (Pr_cT_20 * (P_t_cT_20 + (a_t_cT_2 - a_t_cT_20) %*% Transp(a_t_cT_2 - a_t_cT_20)) + Pr_cT_21 * (P_t_cT_21 + (a_t_cT_2 - a_t_cT_21) %*% Transp(a_t_cT_2 - a_t_cT_21)) + Pr_cT_22 * (P_t_cT_22 + (a_t_cT_2 - a_t_cT_22) %*% Transp(a_t_cT_2 - a_t_cT_22))) / Pr_cT_2

    # One can again average the State Vector across all possible realizations for unconditional average
    a_t_cT <- a_t_cT_0 * Pr_cT_0 + a_t_cT_1 * Pr_cT_1 + a_t_cT_2 * Pr_cT_2
    P_t_cT <- P_t_cT_0 * Pr_cT_0 + P_t_cT_1 * Pr_cT_1 + P_t_cT_2 * Pr_cT_2
    
    u_a_t_cT <- Pr_cT_0 * u_a_t_cT_0 + Pr_cT_1 * u_a_t_cT_1 + Pr_cT_2 * u_a_t_cT_2
    l_a_t_cT <- Pr_cT_0 * l_a_t_cT_0 + Pr_cT_1 * l_a_t_cT_1 + Pr_cT_2 * l_a_t_cT_2

    # Output smoothed vectors and Var matrices
    a_cT_array[1:Dimens, 1:4, i] <- c(a_t_cT_0, a_t_cT_1, a_t_cT_2, a_t_cT)
    CI_a_cT_array[1:Dimens, 1:2, i] <- c(u_a_t_cT, l_a_t_cT)
    P_cT_array[1:Dimens, 1:(Dimens * 4), i] <- c(P_t_cT_0, P_t_cT_1, P_t_cT_2, P_t_cT)

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Update smoothed values for State Vector and Var Matrix for next iteration
    # (t + 1|T) in next iteration
    a_t1_cT_0 <- a_t_cT_0
    a_t1_cT_1 <- a_t_cT_1
    a_t1_cT_2 <- a_t_cT_2

    P_t1_cT_0 <- P_t_cT_0
    P_t1_cT_1 <- P_t_cT_1
    P_t1_cT_2 <- P_t_cT_2
  }

  return(list(
    a_cT = a_cT_array, # Smoothed state vector
    P_cT = P_cT_array, # Smoothed state vector var
    Pr_cT = Pr_cT_mat, # Smoothed regime probs
    a_CI_cT = CI_a_cT_array # Confidence intervals state vector
  ))
}
