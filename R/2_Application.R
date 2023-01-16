#' @description Wrapper function for the Kim filter
#' @param data vector of the data
#' @param Ini initial value for the trend component
#' @param outLogLik logical. If TRUE, the function only returns the log likelihood and no filtered
#' values. If FALSE, the function only returns the log likelihood
#' @param endogen logical. If TRUE, the function includes endogeneity
#' @return either the filtered output or the log likelihood

KimFilter_fctn <- function(param, data, Ini, outLogLik, endogen) {
  if (modelSpec != "ThreeState") {
    # Set up the filter
    paramList <- ParList_fctn(par = param, constrainPar = outLogLik)
    systemList <- SystemMat_fctn(paramList)
    # Run the filter
    Output <- FilterRecursions_fctn(
      data, Ini, systemList, paramList, outLogLik,
      endogen
    )
    # Rcpp function returns a list
    if (outLogLik == TRUE) Output <- Output[[1]]
  } else {
    # Three state model is not implemented in the Rcpp function
    Output <- ifelse(outLogLik == TRUE, Likelihood_fctn(param, data, Ini),
      Filter_fctn(param, data, Ini)
    )
  }
  return(Output)
}


#' @description Function that executes the filter and smoothing recursions and collects the output
#' @param model that is being run
#' @param data vector
#' @param trendIni initial value for the trend component
#' @param dataTib tibble holding the dates
#' @returns list with the filtered and smoothed output

ModelApplication_fctn <- function(model, data, trendIni, dataTib) {
  endogenInd <- ifelse(model == "Endogen", TRUE, FALSE)
  # Collect data, model functions and parameter vector
  gridSearchPath <- paste0("Output/Output_", model, "/GridSearch/")
  thetaVec <- read_rds(file = paste0(gridSearchPath, "thetaVec.rds"))
  seList <- read_rds(file = paste0(gridSearchPath, "thetaSEList.rds"))
  mlParams <- ParConstrain_fctn(thetaVec)
  # Apply Kim Filter recursions
  filterOutput <- KimFilter_fctn(param = mlParams, data = data, Ini = trendIni, outLogLik = FALSE, endogen = endogenInd)
  # Apply Smoother recursions
  smootherOutput <- KimSmoother_fctn(param = mlParams, filterOutput = filterOutput, ciInterval = .9, endogen = endogenInd)

  #---------------------------------------------------------------------------------------#
  # Prepare output
  #---------------------------------------------------------------------------------------#

  # Compute standard errors for the ML estimates
  SEVec <- c(SE_fctn(theta = thetaVec, hessianMat = seList$Hessian), seList$nu_0_sd)
  # Collect the output used for figures 3 and 4 as well as model diagnostics
  nRegimes <- NCOL(filterOutput[["a_ct"]][, , 1])
  periods <- length(smootherOutput[["a_cT"]][1, 1, ])
  outputTib <- tibble(
    Trend = smootherOutput[["a_cT"]][1, nRegimes + 1, ],
    Smoother_Pr_S_0 = smootherOutput[["Pr_cT"]][, 1],
    Filter_Pr_S_0 = filterOutput[["Prob_clt"]][, 1],
    UTrendCI = smootherOutput[["a_CI_cT"]][1, 1, ],
    LTrendCI = smootherOutput[["a_CI_cT"]][1, 2, ],
    Date = dataTib$Date,
    Data = data
  ) %>%
    mutate(
      Diff = Data - Trend
    ) #%>%
    #slice(-c(1:5, periods))
  if (nRegimes == 3) {
    outputTib$Filter_Pr_S_2 <- filterOutput[["Prob_clt"]][, 3]
    outputTib$Smoother_Pr_S_2 <- smootherOutput[["Pr_cT"]][, 3]
  }
  nowcastingPath <- paste0("Output/Output_", model, "/Nowcasting/")
  if (dir.exists(nowcastingPath) & length(list.files(path = nowcastingPath)) != 0) {
    nowcastingOutput <- read_rds(file = paste0(nowcastingPath, "nowcastingOutput.rds"))[, c("Date", "Nowcasting_Pr_S_0")] %>%
      rename("Filter_Pr_S_0" = Nowcasting_Pr_S_0)
    outputTib <- outputTib %>%
      dplyr::select(-Filter_Pr_S_0) %>%
      left_join(., nowcastingOutput, by = "Date")
  }
  # Tibbles indicating when the infections up-turning regime is considered to be active
  smootherRegimeIndicatorTib <- RegimeIdent_fctn(outputTib = outputTib, probVar = "Smoother_Pr_S_0", threshold = 0.4)
  filterRegimeIndicatorTib <- RegimeIdent_fctn(outputTib = outputTib, probVar = "Filter_Pr_S_0", threshold = 0.4)

  return(list(
    filterOutput = filterOutput,
    outputTib = outputTib,
    SE = SEVec,
    smootherOutput = smootherOutput,
    smootherRegimeIndicatorTib = smootherRegimeIndicatorTib,
    filterRegimeIndicatorTib = filterRegimeIndicatorTib
  ))
}


#' @description Function that computes SE based on the Hessian matrix corrected for the optimization parameter
#' constraining function
#' @param theta vector of ML estimates
#' @param hessianMat Hessian of theta
#' @export Vector of the estimated standard errors

SE_fctn <- function(theta, hessianMat) {
  correctionMat <- InvParConstrain_fctn(par = theta)
  CovMat <- correctionMat %*% Inverse(hessianMat) %*% Transp(correctionMat)
  SEVec <- sqrt(diag(CovMat))
  names(SEVec) <- names(theta)
  return(SEVec)
}


#' @description Helper function to create a table of estimated up-turning regime periods
#' @param outputTib tibble with the necessary data
#' @param probVar name of the variable denoting the regime probabilities
#' @param threshold minimum probability in order to constitute an estimated up-turning probability
#' @return tibble with the start and end dates of the estimated up-turning regime

RegimeIdent_fctn <- function(outputTib, probVar, threshold) {
  # Identify the estimated up-turning regimes
  regimePeriods <- outputTib[outputTib[, probVar] > (threshold), c("Date", probVar)] %>%
    mutate(
      End = Date + 1,
      lagEnd = as.Date(lag(End)),
      leadStart = as.Date(lead(Date)),
      type = ifelse(is.na(lagEnd) | Date != lagEnd & End == leadStart, "Start",
        ifelse(is.na(leadStart) | Date == lagEnd & End != leadStart, "Last",
          ifelse(Date == lagEnd & End == leadStart, "Mid", "Lonely")
        )
      ),
      type = ifelse(is.na(leadStart) & is.na(type), "Lonely", type)
    ) %>%
    filter(type != "Mid") %>%
    mutate(End = lead(Date)) %>%
    filter(type != "Last")
  regimePeriods[regimePeriods$type == "Lonely", "End"] <- regimePeriods[regimePeriods$type == "Lonely", "Date"] + 1
  return(regimePeriods)
}
