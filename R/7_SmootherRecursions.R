#' @description Wrapper function for the Kim smoother
#' @param param vector of constrained parameters
#' @param filterOutput Output from the Kim filter
#' @param ciInterval for the smoothed trend
#' @param endogen boolean
#' @return the smoother output

KimSmoother_fctn <- function(param, filterOutput, ciInterval, endogen) {
  if (modelSpec != "ThreeState") {
    paramList <- ParList_fctn(par = param, constrainPar = FALSE)
    systemList <- SystemMat_fctn(paramList)
    Output <- SmootherRecursions_fctn(
      filterOutput = filterOutput, systemList = systemList,
      paramList = paramList, ciInterval, endogen
    )
  } else {
    Output <- Smoother_fctn(par = param, filterOutput = filterOutput, ciInterval = ciInterval)
  }
  return(Output)
}


#' @description Function that executes the Kim filter
#' @param filterOutput Output from the Filter_fctn
#' @param systemList list with system matrices
#' @param paramList list with additional parameters
#' @param ciInterval for the smoothed trend
#' @param endogen boolean
#' @return a list of the smoother output

SmootherRecursions_fctn <- function(filterOutput, systemList, paramList, ciInterval, endogen) {
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
  periods <- length(filterOutput$a_clt[1, 1, ])

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

  # Realizations for t = T of Kim filter and Kim smoother are identical
  Pr_cT_0 <- filterOutput[["Prob_ct"]][periods, 1]
  Pr_cT_1 <- filterOutput[["Prob_ct"]][periods, 2]

  # Initialize smoother output
  a_cT_array <- array(0, c(Dimens, 3, periods))
  P_cT_array <- array(0, c(Dimens, Dimens * 3, periods))
  Pr_cT_mat <- matrix(0, periods, 2)
  CI_a_cT_array <- array(0, c(Dimens, 2, periods))

  a_t1_cT_0 <- a_cT_array[, 1, periods] <- filterOutput[["a_ct"]][, 1, periods] %>% as.matrix()
  a_t1_cT_1 <- a_cT_array[, 2, periods] <- filterOutput[["a_ct"]][, 2, periods] %>% as.matrix()
  a_cT_array[, 3, periods] <- Pr_cT_0 * a_cT_array[, 1, periods] + Pr_cT_1 * a_cT_array[, 2, periods]

  P_t1_cT_0 <- P_cT_array[, 1:Dimens, periods] <- filterOutput[["P_ct"]][, 1:Dimens, periods]
  P_t1_cT_1 <- P_cT_array[, (Dimens + 1):(2 * Dimens), periods] <- filterOutput[["P_ct"]][, (Dimens + 1):(2 * Dimens), periods]
  P_cT_array[, (2 * Dimens + 1):(3 * Dimens), periods] <- Pr_cT_0 * P_t1_cT_0 + Pr_cT_1 * P_t1_cT_1
  
  #---------------------------------------------------------------------------------------#
  # Kim smoother (backwards) iteration
  #---------------------------------------------------------------------------------------#

  for (i in (periods - 1):1) {
    # Load necessary filtered values for each iteration
    # Updated probabilities
    # (t|t)
    Pr_ct_0 <- filterOutput[["Prob_ct"]][i, 1]
    Pr_ct_1 <- filterOutput[["Prob_ct"]][i, 2]

    # Filtered probabilities
    # (t+1|t)
    Pr_clt_0 <- filterOutput[["Prob_clt"]][i + 1, 1] # Note: value read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    Pr_clt_1 <- filterOutput[["Prob_clt"]][i + 1, 2]

    # Updated (approximated) state vector
    # (t|t)
    a_t_ct_0 <- as.matrix(filterOutput[["a_ct"]][, 1, i], 2)
    a_t_ct_1 <- as.matrix(filterOutput[["a_ct"]][, 2, i], 2)

    # Filtered state vector
    # (t+1|t)
    a_t1_ct_00 <- as.matrix(filterOutput[["a_clt"]][, 1, i + 1], 2) # Note: values read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    a_t1_ct_10 <- as.matrix(filterOutput[["a_clt"]][, 2, i + 1], 2)
    a_t1_ct_01 <- as.matrix(filterOutput[["a_clt"]][, 3, i + 1], 2)
    a_t1_ct_11 <- as.matrix(filterOutput[["a_clt"]][, 4, i + 1], 2)

    # Updated state vector Var-Matrix
    # (t|t)
    P_t_ct_0 <- filterOutput[["P_ct"]][, 1:Dimens, i]
    P_t_ct_1 <- filterOutput[["P_ct"]][, (Dimens + 1):(2 * Dimens), i]

    # Filtered state vector Var-Matrix
    # (t+1|t)
    P_t1_ct_00 <- filterOutput[["P_clt"]][, 1:Dimens, i + 1] # Note: values read in as (t+1|t), equivalent to (t|t-1) with lead = 1
    P_t1_ct_10 <- filterOutput[["P_clt"]][, (Dimens + 1):(2 * Dimens), i + 1]
    P_t1_ct_01 <- filterOutput[["P_clt"]][, (2 * Dimens + 1):(3 * Dimens), i + 1]
    P_t1_ct_11 <- filterOutput[["P_clt"]][, (3 * Dimens + 1):(4 * Dimens), i + 1]

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
    Ph_00 <- P_t_ct_0 %*% TranspTt %*% Inverse(P_t1_ct_00)
    Ph_10 <- P_t_ct_1 %*% TranspTt %*% Inverse(P_t1_ct_10)
    Ph_01 <- P_t_ct_0 %*% TranspTt %*% Inverse(P_t1_ct_01)
    Ph_11 <- P_t_ct_1 %*% TranspTt %*% Inverse(P_t1_ct_11)

    # State vector
    a_t_cT_00 <- a_t_ct_0 + Ph_00 %*% (a_t1_cT_0 - a_t1_ct_00)
    a_t_cT_10 <- a_t_ct_1 + Ph_10 %*% (a_t1_cT_0 - a_t1_ct_10)
    a_t_cT_01 <- a_t_ct_0 + Ph_01 %*% (a_t1_cT_1 - a_t1_ct_01)
    a_t_cT_11 <- a_t_ct_1 + Ph_11 %*% (a_t1_cT_1 - a_t1_ct_11)

    # State vector Var-Matrix
    P_t_cT_00 <- P_t_ct_0 + Ph_00 %*% (P_t1_cT_0 - P_t1_ct_00) %*% Transp(Ph_00)
    P_t_cT_10 <- P_t_ct_1 + Ph_10 %*% (P_t1_cT_0 - P_t1_ct_10) %*% Transp(Ph_10)
    P_t_cT_01 <- P_t_ct_0 + Ph_01 %*% (P_t1_cT_1 - P_t1_ct_01) %*% Transp(Ph_01)
    P_t_cT_11 <- P_t_ct_1 + Ph_11 %*% (P_t1_cT_1 - P_t1_ct_11) %*% Transp(Ph_11)

    # Confidence intervals for smoothed state vector
    # Upper
    u_a_t_cT_00 <- a_t_cT_00 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_00[1, 1])
    u_a_t_cT_10 <- a_t_cT_10 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_10[1, 1])
    u_a_t_cT_01 <- a_t_cT_01 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_01[1, 1])
    u_a_t_cT_11 <- a_t_cT_11 + qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_11[1, 1])
    # Lower
    l_a_t_cT_00 <- a_t_cT_00 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_00[1, 1])
    l_a_t_cT_10 <- a_t_cT_10 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_10[1, 1])
    l_a_t_cT_01 <- a_t_cT_01 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_01[1, 1])
    l_a_t_cT_11 <- a_t_cT_11 - qnorm(1 - (1 - ciInterval) / 2) * sqrt(P_t_cT_11[1, 1])

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
    P_t_cT_0 <- (Pr_cT_00 * (P_t_cT_00 + (a_t_cT_0 - a_t_cT_00) %*% Transp(a_t_cT_0 - a_t_cT_00)) + Pr_cT_01 * (P_t_cT_01 + (a_t_cT_0 - a_t_cT_01) %*% Transp(a_t_cT_0 - a_t_cT_01))) / Pr_cT_0
    P_t_cT_1 <- (Pr_cT_10 * (P_t_cT_10 + (a_t_cT_1 - a_t_cT_10) %*% Transp(a_t_cT_1 - a_t_cT_10)) + Pr_cT_11 * (P_t_cT_11 + (a_t_cT_1 - a_t_cT_11) %*% Transp(a_t_cT_1 - a_t_cT_11))) / Pr_cT_1

    # One can again average the State Vector across all possible realizations for unconditional average
    a_t_cT <- a_t_cT_0 * Pr_cT_0 + a_t_cT_1 * Pr_cT_1
    P_t_cT <- P_t_cT_0 * Pr_cT_0 + P_t_cT_1 * Pr_cT_1

    u_a_t_cT <- Pr_cT_0 * u_a_t_cT_0 + Pr_cT_1 * u_a_t_cT_1
    l_a_t_cT <- Pr_cT_0 * l_a_t_cT_0 + Pr_cT_1 * l_a_t_cT_1

    # Output smoothed vectors and Var matrices
    a_cT_array[1:Dimens, 1:3, i] <- c(a_t_cT_0, a_t_cT_1, a_t_cT)
    CI_a_cT_array[1:Dimens, 1:2, i] <- c(u_a_t_cT, l_a_t_cT)
    P_cT_array[1:Dimens, 1:(Dimens * 3), i] <- c(P_t_cT_0, P_t_cT_1, P_t_cT)

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
