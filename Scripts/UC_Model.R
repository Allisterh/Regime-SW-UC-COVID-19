### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Parameter constraints ###########################################################################################################################
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Par_constrain_fctn <- function(par) {
  const_par <- par
  # Assures positive sd of system innovations
  const_par[c(1, 3, 6)] <- exp(par[c(1, 3, 6)])
  # Defines regimes by nu_1 < 0
  const_par[2] <- -exp(par[2])
  # Constrains probabilities to >=0 & <=1
  const_par[4:5] <- 1 / (1 + exp(par[4:5]))
  return(as.numeric(const_par))
}

### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Likelihood function #############################################################################################################################
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Likelihood_fctn <- function(par) {

  #---------------------------------------------------------------------------------------#
  # Initialize filter
  #---------------------------------------------------------------------------------------#

  # Reverse parameter transformation
  par <- Par_constrain_fctn(par)

  # Load input
  xi <- par[1] # Sd of innovations to drift component
  nu_1 <- par[2] # Additional drift for down-turning regime 1
  omega <- par[3] # Sd of innovations to seasonal component
  q <- par[4] # Probability of staying in regime 0
  p <- par[5] # Probability of staying in regime 1
  epsilon <- par[6] # Sd of measurement equation residual

  # Transition matrix
  Tt_temp <- matrix(c(
    1, rep(0, 7), 1, 1, rep(0, 8), -1, 1, rep(0, 6),
    -1, 0, 1, rep(0, 5), -1, 0, 0, 1, rep(0, 4),
    -1, rep(0, 3), 1, 0, 0, 0, -1, rep(0, 4), 1, 0, 0,
    -1, rep(0, 5)
  ), 8, 8)
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  Tt_seas <- rbind(
    matrix(c(rep(0, 6), 1), 1, 7),
    cbind(diag(6), rep(0, 6))
  )
  Tt <- bdiag(Tt_temp, Tt_seas)
  # allows the seasonal unit root to load on gamma_t+1
  Tt[3, 9] <- 1

  # Dimension of state vector
  Dimens <- NCOL(Tt)

  # Number of periods
  T <- Periods

  # Measurement matrix
  Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)

  # Initial unconditional regime probabilities (steady state unconditional probabilities)
  Pr_ct_0 <- (1 - p) / (2 - p - q)
  Pr_ct_1 <- (1 - q) / (2 - p - q)

  # Var-Cov Matrix of innovations in transition equation
  Q <- matrix(0, 2, 2)
  Q[1, 1] <- xi^2
  Q[2, 2] <- omega^2

  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, 8), 1, rep(0, Dimens - 9))
  )
  # Vector holding additional Drift in Regime 1
  lambda <- matrix(c(nu_1, rep(0, Dimens - 1)), Dimens)

  # Initial values for state vector
  a_t_clt_00 <- a_t_clt_10 <- matrix(c(Trend_ini, rep(0, Dimens - 1)), Dimens)
  a_t_clt_01 <- a_t_clt_11 <- matrix(c(Trend_ini, rep(0, Dimens - 1)), Dimens) + lambda

  # Initial variance matrix for state vector (Exact diffuse initialization)
  P_t_clt_00 <- P_t_clt_10 <- P_t_clt_01 <- P_t_clt_11 <- diag(Dimens) * 10000

  # Initializes log-likelihood to record values
  log_lik_T <- 0

  #---------------------------------------------------------------------------------------#
  # Kim filter to set up likelihood function
  #---------------------------------------------------------------------------------------#

  for (i in 6:T) {

    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # One step ahead prediction error with state vector prediction from
    # (t|t-1)
    v_t_00 <- Data_vec[i] - Z %*% a_t_clt_00
    v_t_10 <- Data_vec[i] - Z %*% a_t_clt_10
    v_t_01 <- Data_vec[i] - Z %*% a_t_clt_01
    v_t_11 <- Data_vec[i] - Z %*% a_t_clt_11

    # Variance of prediction error
    F_t_00 <- Z %*% P_t_clt_00 %*% t(Z) + epsilon^2
    F_t_10 <- Z %*% P_t_clt_10 %*% t(Z) + epsilon^2
    F_t_01 <- Z %*% P_t_clt_01 %*% t(Z) + epsilon^2
    F_t_11 <- Z %*% P_t_clt_11 %*% t(Z) + epsilon^2

    # Updating step
    # (t|t, S_t, S_t-1)
    # State vector
    a_t_ct_00 <- a_t_clt_00 + P_t_clt_00 %*% t(Z) %*% Inverse(F_t_00) %*% v_t_00
    a_t_ct_10 <- a_t_clt_10 + P_t_clt_10 %*% t(Z) %*% Inverse(F_t_10) %*% v_t_10
    a_t_ct_01 <- a_t_clt_01 + P_t_clt_01 %*% t(Z) %*% Inverse(F_t_01) %*% v_t_01
    a_t_ct_11 <- a_t_clt_11 + P_t_clt_11 %*% t(Z) %*% Inverse(F_t_11) %*% v_t_11

    # State vector Var-Matrix
    P_t_ct_00 <- P_t_clt_00 - P_t_clt_00 %*% t(Z) %*% Inverse(F_t_00) %*% Z %*% P_t_clt_00
    P_t_ct_10 <- P_t_clt_10 - P_t_clt_10 %*% t(Z) %*% Inverse(F_t_10) %*% Z %*% P_t_clt_10
    P_t_ct_01 <- P_t_clt_01 - P_t_clt_01 %*% t(Z) %*% Inverse(F_t_01) %*% Z %*% P_t_clt_01
    P_t_ct_11 <- P_t_clt_11 - P_t_clt_11 %*% t(Z) %*% Inverse(F_t_11) %*% Z %*% P_t_clt_11

    #-------------------#
    ### Hamilton part ###
    #-------------------#

    # Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
    # (t|t-1)
    Pr_clt_00 <- q %*% Pr_ct_0
    Pr_clt_10 <- (1 - p) %*% Pr_ct_1
    Pr_clt_01 <- (1 - q) %*% Pr_ct_0
    Pr_clt_11 <- p %*% Pr_ct_1

    # Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
    y_dens_clt_00 <- dnorm(v_t_00, mean = 0, sd = sqrt(F_t_00))
    y_dens_clt_10 <- dnorm(v_t_10, mean = 0, sd = sqrt(F_t_10))
    y_dens_clt_01 <- dnorm(v_t_01, mean = 0, sd = sqrt(F_t_01))
    y_dens_clt_11 <- dnorm(v_t_11, mean = 0, sd = sqrt(F_t_11))

    # Sum up joint densities of y_t and regimes to integrate out regime dependencies (receive density of y_t conditional on all information at
    # (t|t-1))
    y_dens_clt <- Pr_clt_00 * y_dens_clt_00 + Pr_clt_10 * y_dens_clt_10 + Pr_clt_01 * y_dens_clt_01 + Pr_clt_11 * y_dens_clt_11

    # Store approximate likelihood
    log_lik_t <- -log(y_dens_clt)

    # Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
    # (t|t)
    Pr_ct_00 <- Pr_clt_00 * y_dens_clt_00 / y_dens_clt
    Pr_ct_10 <- Pr_clt_10 * y_dens_clt_10 / y_dens_clt
    Pr_ct_01 <- Pr_clt_01 * y_dens_clt_01 / y_dens_clt
    Pr_ct_11 <- Pr_clt_11 * y_dens_clt_11 / y_dens_clt

    Pr_ct_0 <- Pr_ct_00 + Pr_ct_10
    Pr_ct_1 <- Pr_ct_01 + Pr_ct_11

    # Sum up log-likelihood over all iterations
    log_lik_T <- log_lik_T + log_lik_t

    # Impose constraint that p,q > 0.9
    if (i == T) {
      log_lik_T <- ifelse(round(p, 3) <= 0.9, 10000 * (1 - p) + log_lik_T,
        ifelse(round(q, 3) <= 0.9, 10000 * (1 - q) + log_lik_T,
          log_lik_T
        )
      )
    }

    #------------------------#
    ### Approximation part ###
    #------------------------#

    # Approximate updated values to break exponential growth of required values
    # (t|t, S_t)
    a_t_ct_0 <- (a_t_ct_00 %*% Pr_ct_00 + a_t_ct_10 %*% Pr_ct_10) %*% solve(Pr_ct_0)
    a_t_ct_1 <- (a_t_ct_01 %*% Pr_ct_01 + a_t_ct_11 %*% Pr_ct_11) %*% solve(Pr_ct_1)

    Pr_ct_00 <- as.numeric(Pr_ct_00)
    Pr_ct_10 <- as.numeric(Pr_ct_10)
    Pr_ct_01 <- as.numeric(Pr_ct_01)
    Pr_ct_11 <- as.numeric(Pr_ct_11)
    Pr_ct_0 <- as.numeric(Pr_ct_0)
    Pr_ct_1 <- as.numeric(Pr_ct_1)

    P_t_ct_0 <- ((Pr_ct_00 * (P_t_ct_00 + (a_t_ct_0 - a_t_ct_00) %*% t(a_t_ct_0 - a_t_ct_00))) + (Pr_ct_10 * (P_t_ct_10 + (a_t_ct_0 - a_t_ct_10) %*% t(a_t_ct_0 - a_t_ct_10)))) / Pr_ct_0
    P_t_ct_1 <- ((Pr_ct_01 * (P_t_ct_01 + (a_t_ct_1 - a_t_ct_01) %*% t(a_t_ct_1 - a_t_ct_01))) + (Pr_ct_11 * (P_t_ct_11 + (a_t_ct_1 - a_t_ct_11) %*% t(a_t_ct_1 - a_t_ct_11)))) / Pr_ct_1

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Prediction step with approximated updates to complete loop
    # (t+1|t)
    # Regime 0 (high drift)
    a_t_clt_00 <- Tt %*% a_t_ct_0
    a_t_clt_10 <- Tt %*% a_t_ct_1
    P_t_clt_00 <- Tt %*% P_t_ct_0 %*% t(Tt) + R %*% Q %*% t(R)
    P_t_clt_10 <- Tt %*% P_t_ct_1 %*% t(Tt) + R %*% Q %*% t(R)

    # Regime 1 (low drift)
    a_t_clt_01 <- Tt %*% a_t_ct_0 + lambda
    a_t_clt_11 <- Tt %*% a_t_ct_1 + lambda
    P_t_clt_01 <- Tt %*% P_t_ct_0 %*% t(Tt) + R %*% Q %*% t(R)
    P_t_clt_11 <- Tt %*% P_t_ct_1 %*% t(Tt) + R %*% Q %*% t(R)
  }

  # Final log-likelihood value
  return(log_lik_T)
}

### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Kim Filter ######################################################################################################################################
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Filter_fctn <- function(par, data) {

  #---------------------------------------------------------------------------------------#
  # Initialize filter
  #---------------------------------------------------------------------------------------#

  # Load input (provided constraining function is already reversed)
  xi <- par[1] # Sd of innovations to drift component
  nu_1 <- par[2] # Additional drift for down-turning regime 1
  omega <- par[3] # Sd of innovations to seasonal component
  q <- par[4] # Probability of staying in regime 0
  p <- par[5] # Probability of staying in regime 1
  epsilon <- par[6] # Sd of measurement equation residual

  # Transition matrix
  Tt_temp <- matrix(c(
    1, rep(0, 7), 1, 1, rep(0, 8), -1, 1, rep(0, 6),
    -1, 0, 1, rep(0, 5), -1, 0, 0, 1, rep(0, 4),
    -1, rep(0, 3), 1, 0, 0, 0, -1, rep(0, 4), 1, 0, 0,
    -1, rep(0, 5)
  ), 8, 8)
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  Tt_seas <- rbind(
    matrix(c(rep(0, 6), 1), 1, 7),
    cbind(diag(6), rep(0, 6))
  )
  Tt <- bdiag(Tt_temp, Tt_seas)
  # allows the seasonal unit root to load on gamma_t+1
  Tt[3, 9] <- 1

  # Dimension of state vector
  Dimens <- NCOL(Tt)

  # Number of periods
  T <- Periods

  # Measurement matrix
  Z <- matrix(c(1, 0, 1, rep(0, Dimens - 3)), 1, Dimens)

  # Initial unconditional regime probabilities (steady state unconditional probabilities)
  Pr_ct_0 <- (1 - p) / (2 - p - q)
  Pr_ct_1 <- (1 - q) / (2 - p - q)

  # Var-Cov Matrix of innovations in transition equation
  Q <- matrix(0, 2, 2)
  Q[1, 1] <- xi^2
  Q[2, 2] <- omega^2

  # Matrix to expand Q so that it matches Var matrix of state vector
  R <- cbind(
    c(1, rep(0, Dimens - 1)),
    c(rep(0, 8), 1, rep(0, Dimens - 9))
  )

  # Vector holding additional Drift in Regime 1
  lambda <- matrix(c(nu_1, rep(0, Dimens - 1)), Dimens)

  # Initial values for state vector
  a_t_clt_00 <- a_t_clt_10 <- matrix(c(Trend_ini, rep(0, Dimens - 1)), Dimens)
  a_t_clt_01 <- a_t_clt_11 <- matrix(c(Trend_ini, rep(0, Dimens - 1)), Dimens) + lambda

  # Initial variance matrix for state vector (Exact diffuse initialization)
  P_t_clt_00 <- P_t_clt_10 <- P_t_clt_01 <- P_t_clt_11 <- diag(Dimens) * 10000

  # Initializes output objects
  a_clt_array <- array(0, c(Dimens, 4, T))
  P_clt_array <- array(0, c(Dimens, Dimens * 4, T))
  a_ct_array <- array(0, c(Dimens, 2, T))
  P_ct_array <- array(0, c(Dimens, Dimens * 2, T))
  Prob_t_clt <- matrix(0, T, 2)
  Prob_t_ct <- matrix(0, T, 2)
  v_ct_mat <- matrix(0, T, 4)
  F_ct_mat <- matrix(0, T, 4)

  #---------------------------------------------------------------------------------------#
  # Kim filter to get filtered values
  #---------------------------------------------------------------------------------------#

  for (i in 6:T) {

    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # One step ahead prediction error with state vector prediction from
    # (t|t-1)
    v_t_00 <- data[i] - Z %*% a_t_clt_00
    v_t_10 <- data[i] - Z %*% a_t_clt_10
    v_t_01 <- data[i] - Z %*% a_t_clt_01
    v_t_11 <- data[i] - Z %*% a_t_clt_11

    # Variance of prediction error
    F_t_00 <- Z %*% P_t_clt_00 %*% t(Z) + epsilon^2
    F_t_10 <- Z %*% P_t_clt_10 %*% t(Z) + epsilon^2
    F_t_01 <- Z %*% P_t_clt_01 %*% t(Z) + epsilon^2
    F_t_11 <- Z %*% P_t_clt_11 %*% t(Z) + epsilon^2

    # Updating step
    # (t|t)
    a_t_ct_00 <- a_t_clt_00 + P_t_clt_00 %*% t(Z) %*% Inverse(F_t_00) %*% v_t_00
    a_t_ct_10 <- a_t_clt_10 + P_t_clt_10 %*% t(Z) %*% Inverse(F_t_10) %*% v_t_10
    a_t_ct_01 <- a_t_clt_01 + P_t_clt_01 %*% t(Z) %*% Inverse(F_t_01) %*% v_t_01
    a_t_ct_11 <- a_t_clt_11 + P_t_clt_11 %*% t(Z) %*% Inverse(F_t_11) %*% v_t_11

    P_t_ct_00 <- P_t_clt_00 - P_t_clt_00 %*% t(Z) %*% Inverse(F_t_00) %*% Z %*% P_t_clt_00
    P_t_ct_10 <- P_t_clt_10 - P_t_clt_10 %*% t(Z) %*% Inverse(F_t_10) %*% Z %*% P_t_clt_10
    P_t_ct_01 <- P_t_clt_01 - P_t_clt_01 %*% t(Z) %*% Inverse(F_t_01) %*% Z %*% P_t_clt_01
    P_t_ct_11 <- P_t_clt_11 - P_t_clt_11 %*% t(Z) %*% Inverse(F_t_11) %*% Z %*% P_t_clt_11

    # Store prediction errors and variances
    v_ct_mat[i, 1:4] <- c(v_t_00, v_t_10, v_t_01, v_t_11)
    F_ct_mat[i, 1:4] <- c(F_t_00, F_t_10, F_t_01, F_t_11)

    #-------------------#
    ### Hamilton part ###
    #-------------------#

    # Computes Pr(S_t = j, S_t-1 = i) conditional on information at time t-1
    # (t|t-1)
    Pr_clt_00 <- q %*% Pr_ct_0
    Pr_clt_10 <- (1 - p) %*% Pr_ct_1
    Pr_clt_01 <- (1 - q) %*% Pr_ct_0
    Pr_clt_11 <- p %*% Pr_ct_1

    # Pr of S_t = j conditional on t = t-1 is given by summarizing over all regimes at t = t-1
    # (t|t-1)
    Pr_clt_0 <- Pr_clt_00 + Pr_clt_10
    Pr_clt_1 <- Pr_clt_01 + Pr_clt_11

    # Density of y_t conditional on states (given by prediction error and prediction error variance from Kalman part)
    y_dens_clt_00 <- dnorm(v_t_00, mean = 0, sd = sqrt(F_t_00))
    y_dens_clt_10 <- dnorm(v_t_10, mean = 0, sd = sqrt(F_t_10))
    y_dens_clt_01 <- dnorm(v_t_01, mean = 0, sd = sqrt(F_t_01))
    y_dens_clt_11 <- dnorm(v_t_11, mean = 0, sd = sqrt(F_t_11))

    # Sum up joint densities of y_t and states to integrate out regime-condition (receive density of y_t conditional on all information at
    # (t|t-1))
    y_dens_clt <- Pr_clt_00 * y_dens_clt_00 + Pr_clt_10 * y_dens_clt_10 + Pr_clt_01 * y_dens_clt_01 + Pr_clt_11 * y_dens_clt_11

    # Update Pr(S_t = j, S_t-1 = i) now conditional on information at time t
    # (t|t)
    Pr_ct_00 <- Pr_clt_00 * y_dens_clt_00 / y_dens_clt
    Pr_ct_10 <- Pr_clt_10 * y_dens_clt_10 / y_dens_clt
    Pr_ct_01 <- Pr_clt_01 * y_dens_clt_01 / y_dens_clt
    Pr_ct_11 <- Pr_clt_11 * y_dens_clt_11 / y_dens_clt

    # Sum up updated probabilities over possible realizations at t-1 to get regime probability at t only conditional on information at time t
    # Floor at 1e-10 makes filter more robust and guarantees the scalar Pr_ct_x to be invertible
    Pr_ct_0 <- max(Pr_ct_00 + Pr_ct_10, 1e-10)
    Pr_ct_1 <- max(Pr_ct_01 + Pr_ct_11, 1e-10)

    # Records predicted probabilities
    # (t|t-1)
    Prob_t_clt[i, 1:2] <- c(Pr_clt_0, Pr_clt_1)

    # Records updated probability
    # (t|t)
    Prob_t_ct[i, 1:2] <- c(Pr_ct_0, Pr_ct_1)

    #------------------------#
    ### Approximation part ###
    #------------------------#

    # Approximate updated values to break exponential growth of required values
    # (t|t, S_t)
    a_t_ct_0 <- (a_t_ct_00 %*% Pr_ct_00 + a_t_ct_10 %*% Pr_ct_10) %*% solve(Pr_ct_0)
    a_t_ct_1 <- (a_t_ct_01 %*% Pr_ct_01 + a_t_ct_11 %*% Pr_ct_11) %*% solve(Pr_ct_1)

    Pr_ct_00 <- as.numeric(Pr_ct_00)
    Pr_ct_10 <- as.numeric(Pr_ct_10)
    Pr_ct_01 <- as.numeric(Pr_ct_01)
    Pr_ct_11 <- as.numeric(Pr_ct_11)
    Pr_ct_0 <- as.numeric(Pr_ct_0)
    Pr_ct_1 <- as.numeric(Pr_ct_1)

    P_t_ct_0 <- ((Pr_ct_00 * (P_t_ct_00 + (a_t_ct_0 - a_t_ct_00) %*% t(a_t_ct_0 - a_t_ct_00))) + (Pr_ct_10 * (P_t_ct_10 + (a_t_ct_0 - a_t_ct_10) %*% t(a_t_ct_0 - a_t_ct_10)))) / Pr_ct_0
    P_t_ct_1 <- ((Pr_ct_01 * (P_t_ct_01 + (a_t_ct_1 - a_t_ct_01) %*% t(a_t_ct_1 - a_t_ct_01))) + (Pr_ct_11 * (P_t_ct_11 + (a_t_ct_1 - a_t_ct_11) %*% t(a_t_ct_1 - a_t_ct_11)))) / Pr_ct_1

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Store predicted values
    # (t|t-1)
    a_clt_array[1:Dimens, 1:4, i] <- c(a_t_clt_00, a_t_clt_10, a_t_clt_01, a_t_clt_11)
    P_clt_array[1:Dimens, 1:(Dimens * 4), i] <- c(P_t_clt_00, P_t_clt_10, P_t_clt_01, P_t_clt_11)

    # Store updated values
    # (t|t)
    a_ct_array[1:Dimens, 1:2, i] <- c(a_t_ct_0, a_t_ct_1)
    P_ct_array[1:Dimens, 1:(Dimens * 2), i] <- c(P_t_ct_0, P_t_ct_1)

    # Prediction step with approximated updates to complete loop
    # (t+1|t)
    # Regime 0 (high drift)
    a_t_clt_00 <- Tt %*% a_t_ct_0
    a_t_clt_10 <- Tt %*% a_t_ct_1
    P_t_clt_00 <- Tt %*% P_t_ct_0 %*% t(Tt) + R %*% Q %*% t(R)
    P_t_clt_10 <- Tt %*% P_t_ct_1 %*% t(Tt) + R %*% Q %*% t(R)

    # Regime 1 (low drift)
    a_t_clt_01 <- Tt %*% a_t_ct_0 + lambda
    a_t_clt_11 <- Tt %*% a_t_ct_1 + lambda
    P_t_clt_01 <- Tt %*% P_t_ct_0 %*% t(Tt) + R %*% Q %*% t(R)
    P_t_clt_11 <- Tt %*% P_t_ct_1 %*% t(Tt) + R %*% Q %*% t(R)
  }

  # Output of filtered values
  return(list(
    a_clt = a_clt_array, # predicted state vector
    P_clt = P_clt_array, # predicted state vector var
    a_ct = a_ct_array, # updated state vector
    P_ct = P_ct_array, # updated state vector var
    Prob_clt = Prob_t_clt, # predicted regime probs
    Prob_ct = Prob_t_ct, # updated regime probs
    Pred_err = v_ct_mat, # one-step-ahead prediction error
    Pred_err_Var = F_ct_mat # Var of pred error
  ))
}

### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Kim Smoother ####################################################################################################################################
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Smoother_fctn <- function(par, Filter_output, level) {

  #---------------------------------------------------------------------------------------#
  # Initialize smoother
  #---------------------------------------------------------------------------------------#

  # Provide system matrix
  Tt_temp <- matrix(c(
    1, rep(0, 7), 1, 1, rep(0, 8), -1, 1, rep(0, 6),
    -1, 0, 1, rep(0, 5), -1, 0, 0, 1, rep(0, 4),
    -1, rep(0, 3), 1, 0, 0, 0, -1, rep(0, 4), 1, 0, 0,
    -1, rep(0, 5)
  ), 8, 8)
  # Add the transition eqn for the seasonal unit root in the seasonal errors
  Tt_seas <- rbind(
    matrix(c(rep(0, 6), 1), 1, 7),
    cbind(diag(6), rep(0, 6))
  )
  Tt <- bdiag(Tt_temp, Tt_seas)
  # allows the seasonal unit root to load on gamma_t+1
  Tt[3, 9] <- 1

  # Dimension of state vector
  Dimens <- NCOL(Tt)

  # Number of periods
  T <- Periods

  # Provide transition probabilities
  q <- par[4] # Probability of staying in regime 0
  p <- par[5] # Probability of staying in regime 1

  # Realizations for t = T of Kim filter and Kim smoother are identical
  Pr_cT_0 <- Filter_output[["Prob_ct"]][T, 1]
  Pr_cT_1 <- Filter_output[["Prob_ct"]][T, 2]

  # Initialize smoother output
  a_cT_array <- array(0, c(Dimens, 3, T))
  P_cT_array <- array(0, c(Dimens, Dimens * 2, T))
  Pr_cT_mat <- matrix(0, T, 2)
  CI_a_cT_array <- array(0, c(Dimens, 2, T))

  a_t1_cT_0 <- a_cT_array[, 1, T] <- Filter_output[["a_ct"]][, 1, T] %>% as.matrix(., 2)
  a_t1_cT_1 <- a_cT_array[, 2, T] <- Filter_output[["a_ct"]][, 2, T] %>% as.matrix(., 2)
  a_cT_array[, 3, T] <- Pr_cT_0 * a_cT_array[1, 1, T] + Pr_cT_1 * a_cT_array[1, 2, T]

  P_t1_cT_0 <- P_cT_array[, 1:Dimens, T] <- Filter_output[["P_ct"]][, 1:Dimens, T]
  P_t1_cT_1 <- P_cT_array[, 16:(2 * Dimens), T] <- Filter_output[["P_ct"]][, 16:(2 * Dimens), T]

  #---------------------------------------------------------------------------------------#
  # Kim smoother (backwards) iteration
  #---------------------------------------------------------------------------------------#

  for (i in (T - 1):5) {

    # Load necessary filtered values for each iteration
    # Updated probabilities
    # (t|t)
    Pr_ct_0 <- Filter_output[["Prob_ct"]][i, 1]
    Pr_ct_1 <- Filter_output[["Prob_ct"]][i, 2]

    # Filtered probabilities
    # (t+1|t)
    Pr_clt_0 <- Filter_output[["Prob_clt"]][i + 1, 1] # Note: value read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    Pr_clt_1 <- Filter_output[["Prob_clt"]][i + 1, 2]

    # Updated (approximated) state vector
    # (t|t)
    a_t_ct_0 <- Filter_output[["a_ct"]][, 1, i] %>% as.matrix(., 2)
    a_t_ct_1 <- Filter_output[["a_ct"]][, 2, i] %>% as.matrix(., 2)

    # Filtered state vector
    # (t+1|t)
    a_t1_ct_00 <- Filter_output[["a_clt"]][, 1, i + 1] %>% as.matrix(., 2) # Note: values read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    a_t1_ct_10 <- Filter_output[["a_clt"]][, 2, i + 1] %>% as.matrix(., 2)
    a_t1_ct_01 <- Filter_output[["a_clt"]][, 3, i + 1] %>% as.matrix(., 2)
    a_t1_ct_11 <- Filter_output[["a_clt"]][, 4, i + 1] %>% as.matrix(., 2)

    # Updated state vector Var-Matrix
    # (t|t)
    P_t_ct_0 <- Filter_output[["P_ct"]][, 1:Dimens, i]
    P_t_ct_1 <- Filter_output[["P_ct"]][, (Dimens + 1):(2 * Dimens), i]

    # Filtered state vector Var-Matrix
    # (t+1|t)
    P_t1_ct_00 <- Filter_output[["P_clt"]][, 1:Dimens, i + 1] # Note: values read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    P_t1_ct_10 <- Filter_output[["P_clt"]][, (Dimens + 1):(2 * Dimens), i + 1]
    P_t1_ct_01 <- Filter_output[["P_clt"]][, (2 * Dimens + 1):(3 * Dimens), i + 1]
    P_t1_ct_11 <- Filter_output[["P_clt"]][, (3 * Dimens + 1):(4 * Dimens), i + 1]

    #-------------------#
    ### Hamilton part ###
    #-------------------#

    # Compute Pr(S_t = j, S_t+1 = k) conditional on all information
    # Here Pr_cT_jk: _cT_01 implies Pr(S_t+1 = 1, S_t = 0 | T) // NOTE: probabilities conditional on t retain prev. notation --> _ct_01: conditional on t Pr(S_t = 1, S_t-1 = 0 | t)
    # (t|T)
    Pr_cT_00 <- (Pr_cT_0 * Pr_ct_0 * q) / Pr_clt_0
    Pr_cT_10 <- (Pr_cT_0 * Pr_ct_1 * (1 - p)) / Pr_clt_0
    Pr_cT_01 <- (Pr_cT_1 * Pr_ct_0 * (1 - q)) / Pr_clt_1
    Pr_cT_11 <- (Pr_cT_1 * Pr_ct_1 * p) / Pr_clt_1

    # Compute Pr(S_t = j | all information) // NOTE: To derive at unconditional probability of Pr(S_t) on has to sum over all regimes at t = t + 1 (in difference to t = t - 1 in the filter)
    # (t|T)
    Pr_cT_0 <- max(Pr_cT_00 + Pr_cT_01, 1e-10)
    Pr_cT_1 <- max(Pr_cT_10 + Pr_cT_11, 1e-10)

    # Record smoothed probabilities
    Pr_cT_mat[i, 1:2] <- c(Pr_cT_0, Pr_cT_1)

    #---------------------#
    ### Kalman part 1/2 ###
    #---------------------#

    # Derive smoothed State vector. NOTE: again _cT_01: S_t = 0 and S_t+1 = 1 // _ct_01: S_t-1 = 0 and S_t = 1
    # (t|T)
    # P"help"_ij
    # (t|T)
    Ph_00 <- P_t_ct_0 %*% t(Tt) %*% Inverse(P_t1_ct_00)
    Ph_10 <- P_t_ct_1 %*% t(Tt) %*% Inverse(P_t1_ct_10)
    Ph_01 <- P_t_ct_0 %*% t(Tt) %*% Inverse(P_t1_ct_01)
    Ph_11 <- P_t_ct_1 %*% t(Tt) %*% Inverse(P_t1_ct_11)

    # State vector
    a_t_cT_00 <- a_t_ct_0 + Ph_00 %*% (a_t1_cT_0 - a_t1_ct_00)
    a_t_cT_10 <- a_t_ct_1 + Ph_10 %*% (a_t1_cT_0 - a_t1_ct_10)
    a_t_cT_01 <- a_t_ct_0 + Ph_01 %*% (a_t1_cT_1 - a_t1_ct_01)
    a_t_cT_11 <- a_t_ct_1 + Ph_11 %*% (a_t1_cT_1 - a_t1_ct_11)

    # State vector Var-Matrix
    P_t_cT_00 <- P_t_ct_0 + Ph_00 %*% (P_t1_cT_0 - P_t1_ct_00) %*% t(Ph_00)
    P_t_cT_10 <- P_t_ct_1 + Ph_10 %*% (P_t1_cT_0 - P_t1_ct_10) %*% t(Ph_10)
    P_t_cT_01 <- P_t_ct_0 + Ph_01 %*% (P_t1_cT_1 - P_t1_ct_01) %*% t(Ph_01)
    P_t_cT_11 <- P_t_ct_1 + Ph_11 %*% (P_t1_cT_1 - P_t1_ct_11) %*% t(Ph_11)

    # Confidence intervals for smoothed state vector
    # Upper
    u_a_t_cT_00 <- a_t_cT_00 + qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_00[1, 1])
    u_a_t_cT_10 <- a_t_cT_10 + qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_10[1, 1])
    u_a_t_cT_01 <- a_t_cT_01 + qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_01[1, 1])
    u_a_t_cT_11 <- a_t_cT_11 + qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_11[1, 1])
    # Lower
    l_a_t_cT_00 <- a_t_cT_00 - qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_00[1, 1])
    l_a_t_cT_10 <- a_t_cT_10 - qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_10[1, 1])
    l_a_t_cT_01 <- a_t_cT_01 - qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_01[1, 1])
    l_a_t_cT_11 <- a_t_cT_11 - qnorm(1 - (1 - level) / 2) * sqrt(P_t_cT_11[1, 1])

    #------------------------#
    ### Approximation part ###
    #------------------------#

    # Approximate smoothed values, so that only conditional on current regime
    # Note: In difference to the filter, probabilities are averaged over all possible regimes at t = t + 1 (not t = t - 1)
    # (t|T)
    # State vector
    a_t_cT_0 <- ((Pr_cT_00 * a_t_cT_00) + (Pr_cT_01 * a_t_cT_01)) / Pr_cT_0
    a_t_cT_1 <- ((Pr_cT_10 * a_t_cT_10) + (Pr_cT_11 * a_t_cT_11)) / Pr_cT_1

    # Also approximate confidence intervals
    u_a_t_cT_0 <- ((Pr_cT_00 * u_a_t_cT_00) + (Pr_cT_01 * u_a_t_cT_01)) / Pr_cT_0
    u_a_t_cT_1 <- ((Pr_cT_10 * u_a_t_cT_10) + (Pr_cT_11 * u_a_t_cT_11)) / Pr_cT_1

    l_a_t_cT_0 <- ((Pr_cT_00 * l_a_t_cT_00) + (Pr_cT_01 * l_a_t_cT_01)) / Pr_cT_0
    l_a_t_cT_1 <- ((Pr_cT_10 * l_a_t_cT_10) + (Pr_cT_11 * l_a_t_cT_11)) / Pr_cT_1

    # Same for State Var-Matrix
    P_t_cT_0 <- (Pr_cT_00 * (P_t_cT_00 + (a_t_cT_0 - a_t_cT_00) %*% t(a_t_cT_0 - a_t_cT_00)) + Pr_cT_01 * (P_t_cT_01 + (a_t_cT_0 - a_t_cT_01) %*% t(a_t_cT_0 - a_t_cT_01))) / Pr_cT_0
    P_t_cT_1 <- (Pr_cT_10 * (P_t_cT_10 + (a_t_cT_1 - a_t_cT_10) %*% t(a_t_cT_1 - a_t_cT_10)) + Pr_cT_11 * (P_t_cT_11 + (a_t_cT_1 - a_t_cT_11) %*% t(a_t_cT_1 - a_t_cT_11))) / Pr_cT_1

    # One can again average the State Vector across all possible realizations for unconditional average
    a_t_cT <- a_t_cT_0 * Pr_cT_0 + a_t_cT_1 * Pr_cT_1

    u_a_t_cT <- Pr_cT_0 * u_a_t_cT_0 + Pr_cT_1 * u_a_t_cT_1
    l_a_t_cT <- Pr_cT_0 * l_a_t_cT_0 + Pr_cT_1 * l_a_t_cT_1

    # Output smoothed vectors and Var matrices
    a_cT_array[1:Dimens, 1:3, i] <- c(a_t_cT_0, a_t_cT_1, a_t_cT)
    CI_a_cT_array[1:Dimens, 1:2, i] <- c(u_a_t_cT, l_a_t_cT)
    P_cT_array[1:Dimens, 1:(Dimens * 2), i] <- c(P_t_cT_0, P_t_cT_1)

    #---------------------#
    ### Kalman part 2/2 ###
    #---------------------#

    # Update smoothed values for State Vector and Var Matrix for next iteration
    # (t + 1|T) in next iteration
    a_t1_cT_0 <- a_t_cT_0
    a_t1_cT_1 <- a_t_cT_1
    P_t1_cT_0 <- P_t_cT_0
    P_t1_cT_1 <- P_t_cT_1
  }

  return(list(
    a_cT = a_cT_array, # Smoothed state vector
    P_cT = P_cT_array, # Smoothed state vector var
    Pr_cT = Pr_cT_mat, # Smoothed regime probs
    a_CI_cT = CI_a_cT_array # Confidence intervals state vector
  ))
}

### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Supplementary functions #########################################################################################################################
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

#---------------------------------------------------------------------------------------#
# C++ function to invert matrices
#---------------------------------------------------------------------------------------#

Rcpp::cppFunction("arma::mat Inverse(const arma::mat & x) {
                  return arma::inv(x);
                  }", depends = "RcppArmadillo")

#---------------------------------------------------------------------------------------#
# Block diagonal matrix construction
#---------------------------------------------------------------------------------------#

bdiag <- function(...) {
  if (nargs() == 1) {
    x <- as.list(...)
  } else {
    x <- list(...)
  }
  n <- length(x)
  if (n == 0) {
    return(NULL)
  }
  x <- lapply(x, function(y) {
    if (length(y)) {
      as.matrix(y)
    } else {
      stop("Zero-length component in x")
    }
  })
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1, ]
  cc <- d[2, ]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1, -1] <- rcum[-n]
  ind[2, ] <- rcum
  ind[3, -1] <- ccum[-n]
  ind[4, ] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) {
    imat[
      (y[1] + 1):y[2],
      (y[3] + 1):y[4]
    ]
  }, imat = imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
}
