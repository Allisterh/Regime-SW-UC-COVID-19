

Nu_1_fctn <- function(stateDraw, regimeDraw, paramList) {
  set.seed(2)
  xi <- paramList$xi$xi_1
  # Set up the priors
  zeta_0 <- 1
  meanMat_0 <- matrix(c(1, zeta_0))
  zeta_1 <- 1
  sdMat_0 <- Inverse(matrix(c(rep(0, 3), zeta_1), nc = 2))
  # Compute the SSR from the regression of of the trend on its lag and drift terms, adjusted for xi so that the error is standard normal
  yMat <- (stateDraw[1, 1, ] - lag(yMat, 1) - stateDraw[2, 1, ]) / xi
  regMat <- regimeDraw / xi
  meanMat <- Inverse(Transp(regMat) %*% regMat) %*% Transp(regMat) %*% yMat
  sdMat <- Inverse(Transp(regMat) %*% regMat)
  draw <- rmvnorm(1, mean = meanMat, sigma = xi * sdMat)
  return(draw)
}


Sdepsilon_fctn <- function(data, stateDraw, regimeDraw) {
  set.seed(2)
  # Set the priors
  zeta_0 <- 0
  delta_0 <- 1
  # Compute the SSR from the regression of the data on the drift and seasonal term
  mu <- stateDraw[1, 1, ]
  gamma <- stateDraw[3, 1, ]
  v_t <- data - mu - gamma
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  zeta_1 <- (zeta_0 + length(data))
  delta_1 <- delta_0 + SSR
  # Draw the new value
  draw <- rinvgamma(1, zeta_1 / 2, delta_1 / 2)
  return(draw)
}


SdXi_fctn <- function(data, stateDraw, regimeDraw, paramList) {
  set.seed(2)
  # Set the priors
  zeta_0 <- 0
  delta_0 <- 1
  # Compute the SSR from the regression of the trend on its lag and drift terms
  nu_1 <- paramList$nu_1
  mu_t <- stateDraw[1, 1, ]
  nu_0 <- stateDraw[2, 1, ]
  mu_lt <- lag(mu_t, n = 1)
  S_t <- regimeDraw
  v_t <- mu_t - mu_lt - nu_0 - S_t * nu_1
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  zeta_1 <- (zeta_0 + length(data))
  delta_1 <- delta_0 + SSR
  # Draw the new value
  draw <- rinvgamma(1, zeta_1 / 2, delta_1 / 2)
  return(draw)
}


SdOmega_fctn <- function(stateDraw) {
  set.seed(2)
  # Set the priors
  zeta_0 <- 0
  delta_0 <- 1
  # Compute the SSR based on the seasonal unit root (regression of the seasonal term on its determ dummies and the unit root)
  v_t <- stateDraw[3, 1, ] - stateDraw[4, 1, ] - stateDraw[5, 1, ] - stateDraw[6, 1, ] - stateDraw[7, 1, ] - stateDraw[7, 1, ] - stateDraw[9, 1, ]
  SSR <- sum(v_t^2, na.rm = T)
  # Adjust the priors
  zeta_1 <- (zeta_0 + dim(stateDraw)[3])
  delta_1 <- delta_0 + SSR
  # Draw the new value
  draw <- rinvgamma(1, zeta_1 / 2, delta_1 / 2)
  return(draw)
}


#' @description Executes the hamilton filter for a given state vector
#' @param stateVec array with draws of the state vector
#' @param paramList list with the additional parameters
#' @param endogen boolean. If TRUE, the function includes endogeneity
#' @return a matrix with the sampled regime realizations

RegimeSampling_fctn <- function(filterOutput, paramList, endogen) {
  set.seed(2)
  if (endogen == FALSE) {
    # Setting up the transition probabilities
    p <- paramList$Probs$p
    q <- paramList$Probs$q
  } else {
    beta_0 <- paramList$beta$beta_0
    beta_1 <- paramList$beta$beta_1
    varrho <- paramList$varrho
    # Transfer probit coefficients into probabilities
    p <- 1 - pnorm(-beta_0 - beta_1)
    q <- pnorm(-beta_0)
  }

  # Length of the observational period
  periods <- NROW(filterOutput)

  # Initialize output object
  S_t_mat <- matrix(NA, periods, 1)
  # Initialize the routine
  Pr_t1_cT_1 <- filterOutput[periods, 2]
  S_t_mat[periods, ] <- S_t_draw <- ifelse(runif(1) <= Pr_t1_cT_1, 1, 0)

  for (i in (periods - 1):1) {
    # Pull the filtered probability
    Pr_ct_1 <- filterOutput[i, 2]
    # Pick the relevant transition probabilities
    if (S_t_draw == 0) {
      P_vec <- c(q, 1 - q)
    } else {
      P_vec <- c(1 - p, p)
    }
    # Calculate the conditional probability for the next iteration
    Pr_t1_cT_1 <- (P_vec[2] * Pr_ct_1) / (P_vec[2] * Pr_ct_1 + P_vec[1] * (1 - Pr_ct_1))

    # Draw a regime probability
    prob <- runif(1)
    S_t_draw <- ifelse(prob <= Pr_t1_cT_1, 1, 0)
    # Record the draw
    S_t_mat[i, ] <- S_t_draw
  }
  return(S_t_mat)
}


#' @description Executes the hamilton filter for a given state vector
#' @param stateVec array with draws of the state vector
#' @param paramList list with the additional parameters
#' @param endogen boolean. If TRUE, the function includes endogeneity
#' @return the filtered output

CondHamilton_fctn <- function(stateVec, paramList, endogen) {
  periods <- dim(stateVec)[3]

  if (endogen == FALSE) {
    # Setting up the transition probabilities
    p <- paramList$Probs$p
    q <- paramList$Probs$q
  } else {
    beta_0 <- paramList$beta$beta_0
    beta_1 <- paramList$beta$beta_1
    varrho <- paramList$varrho
    # Transfer probit coefficients into probabilities
    p <- 1 - pnorm(-beta_0 - beta_1)
    q <- pnorm(-beta_0)
  }

  # Additional drift in regime 1
  nu_1 <- paramList$nu_1

  # Sd of innovations to the drift process
  xi_0 <- paramList$xi$xi_0
  xi_1 <- paramList$xi$xi_1

  # Length of the seasonal term
  seasLength <- 7

  # Initializes output object
  Prob_t_ct <- matrix(0, periods, 2)
  # and the sampler
  Prob_t_ct[1, 1] <- Pr_ct_0 <- (1 - p) / (2 - p - q)
  Prob_t_ct[1, 2] <- Pr_ct_1 <- (1 - q) / (2 - p - q)

  for (i in 2:periods) {
    # Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
    # (t|t-1)
    Pr_clt_00 <- q * Pr_ct_0
    Pr_clt_10 <- (1 - p) * Pr_ct_1
    Pr_clt_01 <- (1 - q) * Pr_ct_0
    Pr_clt_11 <- p * Pr_ct_1

    # Compute the error term of the RW trend for each regime
    mu_t <- stateVec[1, 1, i]
    mu_lt <- stateVec[1, 1, i - 1]
    nu_0 <- stateVec[2, 1, i]
    v_t_0 <- mu_t - mu_lt - nu_0
    v_t_1 <- mu_t - mu_lt - nu_0 - nu_1

    # Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
    y_dens_clt_0 <- dnorm(v_t_0, sd = xi_0)
    y_dens_clt_1 <- dnorm(v_t_1, sd = xi_1)
    # Adjust the densities in case of endogenous switching

    # Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
    # (t|t-1))
    y_dens_clt <- Pr_clt_00 * y_dens_clt_0 + Pr_clt_10 * y_dens_clt_0 + Pr_clt_01 * y_dens_clt_1 + Pr_clt_11 * y_dens_clt_1

    # Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
    # (t|t)
    Pr_ct_00 <- Pr_clt_00 * y_dens_clt_0 / y_dens_clt
    Pr_ct_10 <- Pr_clt_10 * y_dens_clt_0 / y_dens_clt
    Pr_ct_01 <- Pr_clt_01 * y_dens_clt_1 / y_dens_clt
    Pr_ct_11 <- Pr_clt_11 * y_dens_clt_1 / y_dens_clt

    # Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
    # Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
    Pr_ct_0 <- max(Pr_ct_00 + Pr_ct_10, 1e-10)
    Pr_ct_1 <- max(Pr_ct_01 + Pr_ct_11, 1e-10)

    # Records updated probability
    # (t|t)
    Prob_t_ct[i, 1:2] <- c(Pr_ct_0, Pr_ct_1)
  }

  # Output of filtered values
  return(Prob_t_ct)
}



#' @description Function to execute the forward filtering backwards sampling algorithm
#' @param filterOutput list with the updated state vector and var-cov matrix of the kalman filter
#' @param systemList list with the system matrices
#' @param paramList list with the additional parameters
#' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
#' @return the sampled state vectors

BackwardsStateSampling_fctn <- function(filterOutput, systemList, paramList, regimeVec) {
  set.seed(2)

  # unpack the filter output
  a_t_ct_array <- filterOutput$a_ct
  P_t_ct_array <- filterOutput$P_ct
  # System matrices
  Tt <- systemList$Tt
  R <- systemList$R
  # Pre transpose system matrices
  TranspTt <- Transp(Tt)
  transpR <- Transp(R)

  # Dimension of state vector and number of time periods
  Dimens <- NCOL(Tt)
  periods <- length(data)

  # Vector holding additional Drift in Regime 1
  nu_1 <- paramList$nu_1
  lambdaVec <- matrix(c(nu_1, rep(0, Dimens - 1)), nr = Dimens, nc = 1)
  lambdaNullVec <- matrix(0, nr = Dimens, nc = 1)

  # Var-Cov Matrix of innovations in transition equation
  xi_0 <- paramList$xi$xi_0
  xi_1 <- paramList$xi$xi_1
  omega_0 <- paramList$omega$omega_0
  omega_1 <- paramList$omega$omega_1
  eta_0 <- paramList$eta$eta_0
  eta_1 <- paramList$eta$eta_1
  Q_0 <- diag(c(xi_0^2, eta_0^2, omega_0^2))
  Q_1 <- diag(c(xi_1^2, eta_1^2, omega_1^2))
  expandQ_0 <- R %*% Q_0 %*% transpR
  expandQ_1 <- R %*% Q_1 %*% transpR

  # Initializes output objects
  a_t_draw_array <- array(0, c(Dimens, 1, periods))
  a_t1_cT <- a_t_ct_array[, , periods]
  P_t1_cT <- P_t_ct_array[, , periods]
  # Initialize the sampler
  a_t_draw_array[1:Dimens, 1, periods] <- a_t_draw <- Transp(rmvnorm(1, mean = a_t1_cT, sigma = P_t1_cT))

  # Execute the backwards recursions
  for (i in (periods - 1):1) {
    # Pick the updated values from the Kalman filter
    a_t_ct <- a_t_ct_array[, , i]
    P_t_ct <- P_t_ct_array[, , i]
    # Select the system matrices based on the regime
    # !!!! Check if here expansion of the Q matrix with R is correct !!!! #
    if (regimeVec[i] == 1) {
      Q <- expandQ_1
      lambda <- lambdaVec
    } else {
      Q <- expandQ_0
      lambda <- lambdaNullVec
    }
    # Compute the backwards updates
    Trans_P_t_ct_expand <- P_t_ct %*% TranspTt
    Inv_Helper <- Inverse(Tt %*% Trans_P_t_ct_expand + Q)
    a_t1_cT <- a_t_ct + Trans_P_t_ct_expand %*% Inv_Helper %*% (a_t_draw - lambda - Tt %*% a_t_ct)
    P_t1_cT <- P_t_ct - Trans_P_t_ct_expand %*% Inv_Helper %*% Tt %*% P_t_ct
    # To prevent errors
    P_t1_cT_sym <- matrix(NA, nc = Dimens, nr = Dimens)
    P_t1_cT_sym[lower.tri(P_t1_cT_sym, T)] <- P_t1_cT[lower.tri(P_t1_cT, T)]
    P_t1_cT_sym[upper.tri(P_t1_cT_sym)] <- t(P_t1_cT_sym)[upper.tri(P_t1_cT_sym)]
    # if (any(diag(P_t1_cT) < 0)) browser()
    P_t1_cT <- P_t1_cT_sym
    # Draw the state vector from a multivariate normal
    a_t_draw <- Transp(rmvnorm(1, mean = a_t1_cT, sigma = P_t1_cT))
    # Store the sampled state vectors
    a_t_draw_array[1:Dimens, 1, i] <- a_t_draw
  }
  return(a_t_draw_array)
}


#' @description Function that computes the state vector via the Kalman filter given a vector of drawn regime realizations and additional parameters
#' @param data vector of the data
#' @param Ini initial value for the trend component
#' @param systemList list with the system matrices
#' @param paramList list with the additional parameters
#' @param regimeVec vector with regime realizations. Vector takes either values 1 or 0
#' @param endogen logical. If TRUE, the function includes endogeneity
#' @return the filtered output

CondKalman_fctn <- function(data, Ini, systemList, paramList, regimeVec, endogen = FALSE) {
  # System matrices
  Tt <- systemList$Tt
  Z <- systemList$Z
  R <- systemList$R
  # Pre transpose system matrices
  TranspTt <- Transp(Tt)
  transpR <- Transp(R)
  TranspZ <- Transp(Z)

  # Dimension of state vector and number of time periods
  Dimens <- NCOL(Tt)
  periods <- length(data)

  # Residual variance in the measurement eq
  epsilon <- paramList$epsilon
  epsilonSq <- epsilon^2

  if (endogen == FALSE) {
    # Setting up the transition probabilities
    p <- paramList$Probs$p
    q <- paramList$Probs$q
  } else {
    beta_0 <- paramList$beta$beta_0
    beta_1 <- paramList$beta$beta_1
    varrho <- paramList$varrho
    # Transfer probit coefficients into probabilities
    p <- 1 - pnorm(-beta_0 - beta_1)
    q <- pnorm(-beta_0)
  }
  # Steady state probabilities
  Pr_ct_0 <- (1 - p) / (2 - p - q)
  Pr_ct_1 <- (1 - q) / (2 - p - q)

  # Vector holding additional Drift in Regime 1
  nu_1 <- paramList$nu_1
  lambda <- matrix(c(nu_1, rep(0, Dimens - 1)), Dimens)

  # Var-Cov Matrix of innovations in transition equation
  xi_0 <- paramList$xi$xi_0
  xi_1 <- paramList$xi$xi_1
  omega_0 <- paramList$omega$omega_0
  omega_1 <- paramList$omega$omega_1
  eta_0 <- paramList$eta$eta_0
  eta_1 <- paramList$eta$eta_1
  Q_0 <- diag(c(xi_0^2, eta_0^2, omega_0^2))
  Q_1 <- diag(c(xi_1^2, eta_1^2, omega_1^2))
  expandQ_0 <- R %*% Q_0 %*% transpR
  expandQ_1 <- R %*% Q_1 %*% transpR

  # Initial values for state vector
  if (regimeVec[2] == 0) {
    a_t_clt <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens)
  } else {
    a_t_clt <- matrix(c(Ini, rep(0, Dimens - 1)), Dimens) + lambda
  }

  # Initial variance matrix for state vector (Exact diffuse initialization)
  P_t_clt <- diag(Dimens) * 10000

  # Initializes output objects
  a_ct_array <- array(0, c(Dimens, 1, periods))
  P_ct_array <- array(0, c(Dimens, Dimens, periods))

  for (i in 1:periods) {
    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # One step ahead prediction error with state vector prediction from
    # (t|t-1)
    v_t <- as.numeric(data[i] - Z %*% a_t_clt)
    # Variance of prediction error
    F_t <- as.numeric(Z %*% P_t_clt %*% TranspZ + epsilonSq)

    # Updating step
    # (t|t, S_t, S_t-1)
    # State vector
    a_t_ct <- a_t_clt + P_t_clt %*% TranspZ %*% (1 / F_t * v_t)
    # State vector Var-Matrix
    P_t_ct <- P_t_clt - P_t_clt %*% TranspZ %*% (1 / F_t) %*% Z %*% P_t_clt

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Store updated values
    # (t|t)
    a_ct_array[1:Dimens, 1, i] <- a_t_ct
    P_ct_array[1:Dimens, 1:Dimens, i] <- P_t_ct

    # pull the regime realization for the next period
    if (i < periods) {
      if (regimeVec[i + 1] == 0) {
        # Prediction step with approximated updates to complete loop
        # (t+1|t)
        # Regime 0 (high drift)
        a_t_clt <- Tt %*% a_t_ct
        P_t_clt <- Tt %*% P_t_ct %*% TranspTt + expandQ_0
      } else {
        # Regime 1 (low drift)
        a_t_clt <- Tt %*% a_t_ct + lambda
        P_t_clt <- Tt %*% P_t_ct %*% TranspTt + expandQ_1
      }
    }
  }
  return(list(
    a_ct = a_ct_array, # updated state vector
    P_ct = P_ct_array # updated state vector var
  ))
}
