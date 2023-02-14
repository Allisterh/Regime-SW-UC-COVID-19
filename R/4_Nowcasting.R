#' @description Function to run the nowcasting exercise
#' @param model that underlies the grid search
#' @param data vector for which the likelihood values are to be computed
#' @param trendIni inital value for the trend component
#' @param dataTib tibble with the dates
#' @param iniPeriod number of periods used to build an initial data history
#' @param estimInterval interval in between new estimations of theta
#' @param gridSearchInterval interval in between parameter vector grid searches
#' @param newRun logical. If FALSE, pick up from previous calculations
#' @export Table in the ´Model´_Output/Nowcasting folder

Nowcasting_fctn <- function(model, data, trendIni, dataTib, iniPeriod, estimInterval, gridSearchInterval,
                            newRun = FALSE) {
  endogenInd <- ifelse(model == "Endogen", TRUE, FALSE)
  # Create an output folder
  path <- paste0(getwd(), "/Output/Output_", model, "/Nowcasting/")
  if (!dir.exists(path)) {
    dir.create(path)
  }
  nPeriods <- length(data)
  # Load a parameter vector to get the number of estimated parameters
  thetaVec <- read_rds(file = paste0("Output/Output_", model, "/GridSearch/thetaVec.rds"))
  nRegimes <- ifelse(length(thetaVec) == 10, 3, 2)
  # Define the time points for parameter estimation and further grid searches
  paramEstim <- seq(iniPeriod + 1 + estimInterval, nPeriods, estimInterval)
  gridSearches <- seq(iniPeriod + 1, nPeriods, gridSearchInterval)
  eventPeriods_temp <- unique(c(paramEstim, gridSearches, nPeriods))
  eventPeriods <- eventPeriods_temp[order(eventPeriods_temp)]

  if (!("nowcastingOutput.rds" %in% list.files(path = path)) & newRun == FALSE) {
    cat("No previous results stored. Starting a new nowcasting simulation run. \n")
    newRun <- TRUE
  }
  if (newRun == FALSE) {
    # load the previously calculated time periods and pick up from there
    nowcastingOutput <- read_rds(file = paste0(path, "nowcastingOutput.rds"))
    nowcastingVec <- nowcastingOutput$Nowcasting_Pr_S_0
    paramMat <- dplyr::select(nowcastingOutput, -c(Date, Nowcasting_Pr_S_0)) %>%
      as.matrix()
    prevI <- NROW(nowcastingOutput)
    eventPeriods <- eventPeriods[eventPeriods > prevI]
    cat(paste(prevI, "periods already calculated. Proceeding with event period", eventPeriods[1]), "\n")
  } else {
    # Initialize the temporary output objects when starting a new run
    paramMat <- matrix(NA, nc = length(thetaVec) + 1, nr = iniPeriod)
    nowcastingVec <- c(rep(NA, iniPeriod))
    prevI <- iniPeriod
  }

  # Execute the routine for the remainder of the observational period
  for (i in eventPeriods) {
    # If time period i coincides with a scheduled grid search, run the "light" grid search
    if (i %in% gridSearches) {
      cat(paste("Running a robustification grid search at event period", i, "\n"))
      gridSearchOutput <- GridSearch_fctn(
        model = model, data = data[1:i], trendIni = trendIni, nRandom = 5e04, stepsGrid = 4, storeOutput = F,
        messages = F
      )
      ParamVec <- gridSearchOutput
    } else {
      # If i does not coincide with a grid search, simply run the numerical optimization with the previous
      # parameters as starting values
      optimResult <- tryCatch(optim(paramMat[NROW(paramMat), -1],
        fn = KimFilter_fctn, hessian = FALSE, method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-06), data = data[1:i], Ini = trendIni, outLogLik = TRUE, endogen = endogenInd
      ), error = function(e) NA)
      # Construct the parameter vector including the log likelihood
      ParamVec <- c(optimResult$value, optimResult$par)
    }
    # Store the parameters
    paramMat <- rbind(paramMat, matrix(rep(ParamVec, i - prevI), byrow = T, nc = NCOL(paramMat)))
    mlParams <- ParConstrain_fctn(paramMat[NROW(paramMat), -1])
    # Run the filter recursions
    filterOutput <- KimFilter_fctn(par = mlParams, data = data[1:i], Ini = trendIni, outLogLik = FALSE, endogen = endogenInd)
    nowcastingVals <- filterOutput$Prob_clt[(prevI + 1):i, 2]
    # Store the nowcasted infection numbers
    nowcastingVec <- c(nowcastingVec, nowcastingVals)
    prevI <- i
    # Construct a proper output
    rownames(paramMat) <- NULL
    colnames(paramMat) <- c("ll", names(thetaVec))
    nowcastingOutput <- dataTib[1:length(nowcastingVec), ] %>%
      dplyr::select(Date) %>%
      mutate(Nowcasting_Pr_S_0 = nowcastingVec) %>%
      bind_cols(paramMat)
    # Store the output
    write_rds(nowcastingOutput, file = paste0(path, "nowcastingOutput.rds"))
    cat(paste("\r Calculation of event period", i, "finished"))
  }
  write.xlsx(nowcastingOutput, file = paste0(path, "NowcastingOutput.xlsx"))
}
