#' @description Function that runs a Monte Carlo simulation exercise for the frequentist estimation routine
#' @param model under consideration
#' @param nPeriods number of time periods for which to simulate the DGP
#' @param nIters number of DGPs to simulate
#' @param nRandomGrid number of parameter vectors to draw in the initial grid search
#' @param stepsGrid number of grid points to search over in the second grid search
#' @param newGridSerach boolean. If TRUE, a new initial grid search is started.
#' @return tibble with moments and quantiles of the estimated parameters
#' @export .xlsx table with moments and quantiles of the estimated parameters

KimMCStudy_fctn <- function(model, nPeriods, nIters, nRandomGrid, stepsGrid, newGridSearch = FALSE) {
  set.seed(2)

  # Create an output folder
  path <- paste0(getwd(), "/Output/Output_", model, "/MonteCarlo/")
  if (!dir.exists(path)) {
    dir.create(path)
  }

  # Initialize output objects
  paramNames <- colnames(ThetaRand_fctn(nRandom = 1))
  nParams <- length(paramNames)
  paramMat <- matrix(NA, nr = nIters, nc = nParams)
  colnames(paramMat) <- paramNames

  # Initial iteration
  if (newGridSearch == FALSE) {
    if ("iniThetaVec.rds" %in% list.files(path = path)) {
      iniThetaVec <- read_rds(file = paste0(path, "iniThetaVec.rds"))
    } else {
      cat("No previous initial parameter vector stored\n")
      newGridSearch <- TRUE
    }
  }

  if (newGridSearch == TRUE) {
    cat("Performing the initial grid search\n")
    dataVec <- simData_fctn(epsilonInd = TRUE, nPeriods = nPeriods)[[1]]
    gridSearchOutput <- GridSearch_fctn(
      model = modelSpec, data = dataVec, trendIni = round(dataVec[1]),
      nRandomGrid = nRandomGrid, stepsGrid = stepsGrid, storeOutput = FALSE, messages = FALSE
    )
    iniThetaVec <- gridSearchOutput[-1]
    write_rds(iniThetaVec, file = paste0(path, "iniThetaVec.rds"))
  }
  paramMat[1, ] <- ParConstrain_fctn(iniThetaVec)

  # Construct the simulated data
  dataMat <- sapply(1:(nIters - 1), function(x, nPeriods) {
    simData_fctn(nPeriods = nPeriods)
  }, nPeriods = nPeriods)

  # Estimate the parameters
  cat("Running the monte carlo simulation study\n")
  paramMatOptim <- t(pbapply(dataMat, 2, function(x, thetaVec) {
    optimOutput <- tryCatch(
      optim(thetaVec,
        fn = KimFilter_fctn, hessian = FALSE, method = "Nelder-Mead",
        control = list(maxit = 5e3, reltol = 1e-6), data = x, Ini = round(x[1]), outLogLik = TRUE, endogen = FALSE
      ),
      error = function(e) NA
    )
    return(optimOutput$par)
  }, thetaVec = iniThetaVec))
  paramMat[-1, ] <- t(apply(paramMatOptim, 1, ParConstrain_fctn))

  # Obtain nu_0
  nu_0_vec <- apply(paramMat, 1, function(x, data) {
    filterOutput <- KimFilter_fctn(param = x, data = data, Ini = round(x[1]), outLogLik = FALSE, endogen = FALSE)
    smootherOutput <- tryCatch(
      KimSmoother_fctn(
        param = x, filterOutput = filterOutput, ciInterval = .9, endogen = FALSE
      ),
      error = function(e) list("a_cT" = array(NA, dim = c(2, 3, nPeriods)))
    )
    nu_1 <- smootherOutput$a_cT[2, 3, nPeriods]
    return(nu_1)
  }, data = dataVec)

  paramMatFinal <- cbind(paramMat, "nu_0" = nu_0_vec)

  # Summarize the results
  meanVec <- t(apply(paramMatFinal, 2, mean, na.rm = TRUE))
  sdVec <- t(apply(paramMatFinal, 2, sd, na.rm = TRUE))
  LCIVec <- meanVec - sdVec * 1.64 / sqrt(length(meanVec))
  UCIVec <- meanVec + sdVec * 1.64 / sqrt(length(meanVec))
  quantileMat <- t(apply(paramMatFinal, 2, function(x) {
    quantile(x, probs = c(0.025, .25, .5, .75, .975), na.rm = TRUE)
  }))

  # Construct and save the output
  outputTib <- tibble(
    mean = c(meanVec),
    sd = c(sdVec),
    L95CI = c(LCIVec),
    U95CI = c(UCIVec)
  ) %>%
    cbind(., quantileMat) %>%
    round(3) %>%
    rownames_to_column() %>%
    slice(match(names(LabelMap_fctn()), rowname))
  write.xlsx(outputTib, file = paste0(path, "MonteCarloOutput_", nPeriods, "_", nIters, ".xlsx"))
  write_rds(paramMatFinal, file = paste0(path, "MonteCarloParamMat_", nPeriods, "_", nIters, ".rds"))
  return(outputTib)
}




#' @description Function that runs a Monte Carlo simulation exercise for the bayesian estimation routine
#' @param model under consideration
#' @param trendIni inital value for the trend component
#' @param nPeriods number of time periods for which to simulate the DGP
#' @param nIters number of DGPs to simulate
#' @return tibble with moments and quantiles of the estimated parameters
#' @export .xlsx table with moments and quantiles of the estimated parameters

GibbsMCStudy_fctn <- function(model = modelSpec, trendIni = trendIni, nPeriods, nIters) {
  set.seed(12)

  # Create an output folder
  path <- paste0(getwd(), "/Output/Output_", model, "/GibbsMonteCarlo/")
  if (!dir.exists(path)) {
    dir.create(path)
  }

  # Initialize output objects
  paramNames <- colnames(ThetaRand_fctn(nRandom = 1))
  nParams <- length(paramNames)
  paramMat <- matrix(NA, nr = nIters, nc = nParams)
  colnames(paramMat) <- paramNames

  # Construct the simulated data
  dataMat <- sapply(1:nIters, function(x, nPeriods) {
    simData_fctn(nPeriods = nPeriods)[[1]]
  }, nPeriods = nPeriods)

  # Estimate the parameters
  cat("Running the monte carlo simulation study\n")

  outputMat <- pbapply(dataMat, 2, function(x, modelSpec, trendIni, dataTib) {
    outputList <- GibbsSampler_fctn(
      model = modelSpec, data = x, trend = trendIni, nIters = 5e3, nBurnin = 2e3,
      dataTib = dataTib, thinIntervall = 2, thinIntervallReverse = TRUE, storeOutput = FALSE,
      progressBar = FALSE
    )
    paramMat <- outputList$additionalParam
    meanVec <- apply(paramMat, 2, mean)
    sdVec <- apply(paramMat, 2, sd)
    quantileMat <- apply(paramMat, 2, function(x) quantile(x, probs = c(.025, .25, .5, .75, .975)))
    outputMat <- rbind(meanVec, sdVec, quantileMat)
    return(outputMat)
  }, modelSpec = modelSpec, trendIni = trendIni, dataTib = dataTib)

  # Bring the output into a proper format
  outputArray <- array(NA, dim = c(7, nParams, nIters))
  for (i in 0:(nIters - 1)) {
    outputArray[, , i + 1] <- outputMat[(1 + i * 7):(i * 7 + 7), ]
  }
  meanMat <- apply(outputArray, c(1, 2), mean)
  colnames(meanMat) <- paramNames
  rownames(meanMat) <- c("Mean", "Sd", ".025", ".25", ".5", ".75", ".975")

  # Construct and save the output
  write_rds(outputArray, file = paste0(path, "MonteCarloOutputArray_", nPeriods, "_", nIters, ".rds"))
  write.xlsx(as_tibble(meanMat), file = paste0(path, "MonteCarloOutput_", nPeriods, "_", nIters, ".xlsx"))
  return(as_tibble(meanMat))
}
