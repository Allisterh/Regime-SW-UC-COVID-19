#' @description Function to create distribution plots as part of the grid search output
#' @param model that underlies the grid search
#' @param data vector with the noisy measurement
#' @param trendIni inital value for the trend component
#' @param nRandomGrid number of parameter vectors to draw in the initial grid search
#' @param stepsGrid number of grid points to search over in the second grid search
#' @param storeOutput boolean If FALSE, no output is stored and only the best fitting parameter
#'        vector returned.
#' @param messages boolean If FALSE, no messages are printed
#' @export if storeOutput == TRUE Plots and tables in the ´Model´_Output/GridSearch folder
#' @return if storeOutput == FALSE, returns the best fitting parameter vector as well as the
#' associated likelihood and information criteria

GridSearch_fctn <- function(model, data, trendIni, nRandomGrid, stepsGrid, storeOutput = TRUE, messages = TRUE) {
  # Create an output folder
  path <- paste0(getwd(), "/Output/Output_", model, "/GridSearch/")
  if (!dir.exists(path) & storeOutput == TRUE) {
    dir.create(path)
  }
  notImplementedModels <- c("D.Seas.3.St.", "UR.Seas.MAC", "UR.Seas.En")
  if (modelSpec %in% notImplementedModels) {
    stop(
      "The Kim filter is not yet implemented for the model specification ",
      modelSpec, "\n"
    )
  }

  #-------------------------------------------------------------------------------------#
  # Initial random grid search                                                       ####
  #-------------------------------------------------------------------------------------#

  set.seed(2)
  # Build the random grid
  thetaRand <- ThetaRand_fctn(nRandomGrid = nRandomGrid)
  paramNames <- colnames(thetaRand)
  # Check if the model is using a Probit specification for the regime process
  probitIndicator <- ifelse(any(str_detect(paramNames, "beta")) == TRUE, TRUE, FALSE)
  # Search over parameter grid
  if (messages == TRUE) cat("Step 1/3: Begin initial random grid search \n")
  resultsRandUnfltrd <- pbapply(thetaRand, 1, function(theta, data, Ini, endogen) {
    return(c(
      tryCatch(KimFilter_fctn(param = theta, data = data, Ini = Ini, outLogLik = TRUE, endogen = endogen),
        error = function(e) NA
      ),
      theta
    ))
  }, data = dataVec, Ini = trendIni, endogen = probitIndicator) %>%
    t()

  # Collect results
  resultsRandUnfltrd <- resultsRandUnfltrd[order(resultsRandUnfltrd[, 1]), ]
  resultsRand <- resultsRandUnfltrd[!is.na(resultsRandUnfltrd[, 1]), ]
  if (storeOutput == TRUE) {
    write_rds(resultsRand, file = paste0(path, "resultsRand.rds"))
    # Analyse best thetas
    resultsRandCnstrnd <- t(apply(resultsRand[, -1], 1, ParConstrain_fctn))
    colnames(resultsRandCnstrnd) <- names(thetaRand)
    # Transform the Probit AR coefficients to probabilities if necessary
    if (probitIndicator == TRUE) resultsRandCnstrnd <- ProbitProb_fctn(resultsRandCnstrnd)
    PlotParameterDist_fctn(resultsRandCnstrnd[1:50, ],
      title = "Step1_RandGrid", path = path,
      model = model
    )
  }

  #-------------------------------------------------------------------------------------#
  # Second grid search                                                               ####
  #-------------------------------------------------------------------------------------#

  # Build grid with the most extreme realizations of the 50 best thetas
  thetaGrid <- ThetaGrid_fctn(resultsRand, stepsGrid)
  # Search over parameter grid
  if (messages == TRUE) cat("Step 2/3: Begin second grid search \n")
  resultsGridUnfltrd <- pbapply(thetaGrid, 1, function(theta, data, Ini, endogen) {
    return(c(
      tryCatch(KimFilter_fctn(param = theta, data = data, Ini = Ini, outLogLik = TRUE, endogen = endogen),
        error = function(e) NA
      ),
      theta
    ))
  }, data = dataVec, Ini = trendIni, endogen = probitIndicator) %>%
    t()
  # Collect results
  resultsGridUnfltrd <- resultsGridUnfltrd[order(resultsGridUnfltrd[, 1]), ]
  resultsGrid <- resultsGridUnfltrd[!is.na(resultsGridUnfltrd[, 1]), ]
  if (storeOutput == TRUE) {
    write_rds(resultsGrid, file = paste0(path, "resultsGrid.rds"))
    # Analyse best thetas
    resultsGridCnstrnd <- t(apply(resultsGrid[, -1], 1, ParConstrain_fctn))
    colnames(resultsGridCnstrnd) <- paramNames
    if (probitIndicator == TRUE) resultsGridCnstrnd <- ProbitProb_fctn(resultsGridCnstrnd)
    PlotParameterDist_fctn(resultsGridCnstrnd[1:50, ],
      title = "Step2_SecondGrid",
      path = path, model = model
    )
  }

  #-------------------------------------------------------------------------------------#
  # Final grid search                                                                ####
  #-------------------------------------------------------------------------------------#

  thetaFinal <- resultsGrid[1:50, -1]
  # Optimize over 50 best thetas from the second step
  if (messages == TRUE) cat("Step 3/3: Begin final optimization \n")
  resultsOptim <- pbapply(thetaFinal, 1, function(theta, data, Ini, endogen) {
    return(c(
      tryCatch(optim(theta,
        fn = KimFilter_fctn, hessian = FALSE, method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-06), data = data, Ini = Ini, outLogLik = TRUE,
        endogen = endogen
      ), error = function(e) NA),
      theta
    ))
  }, data = dataVec, Ini = trendIni, endogen = probitIndicator) %>%
    t()
  # Collect the ML parameters from the list of optim results
  resultsOptimUnfltrd <- sapply(resultsOptim, function(x, nparams) {
    if (is.list(x)) {
      outputVec <- c(x$value, x$par)
    } else {
      oututVec <- rep(NA, nparams)
    }
  }, nparams = NCOL(thetaRand) + 1) %>%
    t()
  colnames(resultsOptimUnfltrd)[1] <- "ll"
  resultsOptimUnfltrd <- resultsOptimUnfltrd[order(resultsOptimUnfltrd[, "ll"]), ]
  resultsOptimMat <- resultsOptimUnfltrd[!is.na(resultsOptimUnfltrd[, 1]), ]
  # Select the best parameter vector and save results
  thetaVec <- resultsOptimMat[1, -1]
  if (storeOutput == TRUE) {
    write_rds(resultsOptimMat, file = paste0(path, "resultsOptimMat.rds"))
    write_rds(thetaVec, file = paste0(path, "thetaVec.rds"))
  }
  # Compute nu_0 estimates based on the Kim smoother
  nu_0_mat <- matrix(ncol = 2, nrow = nrow(resultsOptimMat))
  nu_0_mat <- apply(resultsOptimMat, MARGIN = 1, function(x, data, trendIni, endogen) {
    mleParam <- ParConstrain_fctn(x[-1])
    filterOutput <- KimFilter_fctn(param = mleParam, data = data, Ini = trendIni, outLogLik = FALSE, endogen = endogen)
    smootherOutput <- tryCatch(
      KimSmoother_fctn(
        param = mleParam, filterOutput = filterOutput, ciInterval = .9, endogen = endogen
      ),
      error = function(e) NA
    )
    Dimens <- NROW(filterOutput$a_ct[, , 1])
    nRegimes <- NCOL(filterOutput$a_ct[, , 1])
    nu_0 <- tryCatch(
      c(
        smootherOutput[["a_cT"]][2, nRegimes + 1, length(data)],
        sqrt(smootherOutput[["P_cT"]][2, Dimens * nRegimes + 2, length(data)])
      ),
      error = function(e) rep(NA, 2)
    )
    return(nu_0)
  }, data = data, trendIni = trendIni, endogen = probitIndicator) %>%
    Transp()
  colnames(nu_0_mat) <- c("Estimate", "nu_0")

  # Reverse the parameter constraining function
  resultsOptimCnstrndMat <- cbind(
    -resultsOptimMat[, 1], apply(
      resultsOptimMat[, -1],
      1, ParConstrain_fctn
    ) %>% Transp()
  )
  colnames(resultsOptimCnstrndMat) <- c("ll", paramNames)
  if (probitIndicator == TRUE) resultsOptimCnstrndMat <- ProbitProb_fctn(resultsOptimCnstrndMat)
  resultsOptimCnstrndMat <- cbind(resultsOptimCnstrndMat, "nu_0" = nu_0_mat[, 1])
  # Compute the model selection criteria
  InfCritMat <- t(apply(resultsOptimCnstrndMat, MARGIN = 1, InfCritWrapper_fctn, nPeriods = length(data)))
  resultsOptimCnstrndMat <- cbind(resultsOptimCnstrndMat, InfCritMat)
  # Calculate the Hessian matrix for the best fitting parameter vector
  Hessian <- optimHess(thetaVec, KimFilter_fctn, data = data, Ini = trendIni, outLogLik = TRUE, endogen = probitIndicator)
  seList <- list(Hessian = Hessian, nu_0_sd = nu_0_mat[1, 2])

  if (storeOutput == TRUE) {
    write_rds(resultsOptimCnstrndMat, file = paste0(path, "resultsOptimCnstrndMat.rds"))
    # Save to excel
    write.xlsx(as_tibble(cbind(Index = 1:NROW(resultsOptimCnstrndMat), round(resultsOptimCnstrndMat, 3))),
      file = paste0(path, "resultsOptimCnstrndMat.xlsx")
    )
    write_rds(seList, file = paste0(path, "thetaSeList.rds"))
  } else {
    # If no output is to be stored, return the best parameter vector
    return(resultsOptimMat[1, ])
  }
}


#' @description Function that constructs a parameter grid
#' @param resultsRand matrix with the output of the initial grid search
#' @param stepsGrid number grid points in between the lower and upper limits
#' @return tibble with the parameter grid

ThetaGrid_fctn <- function(resultsRand, stepsGrid) {
  UGridLimits <- apply(resultsRand[1:50, -1], MARGIN = 2, max)
  LGridLimits <- apply(resultsRand[1:50, -1], MARGIN = 2, min)
  limitsMat <- rbind(UGridLimits, LGridLimits)
  thetaGridMat <- apply(limitsMat, MARGIN = 2, function(x, stepsGrid) {
    parSupport <- seq(x[2], x[1], length.out = stepsGrid)
    return(parSupport)
  }, stepsGrid)
  colnames(thetaGridMat) <- colnames(resultsRand)[-1]
  # Expand the grid
  thetaGridTib <- thetaGridMat %>%
    as_tibble() %>%
    expand.grid() %>%
    as_tibble()
  return(thetaGridTib)
}


#' @description Function to compute the AIC, HQ and BIC information criteria
#' @param logLik log likelihood (times -1)
#' @param nParams number of estimated parameters of the model
#' @param nDiffus number of parameters that are initialized diffusely
#' @param nPeriods number of time periods
#' @return vector of information criteria

InfCrit_fctn <- function(logLik, nParams, nDiffus, nPeriods) {
  # Define the different penalties
  Penalties <- c("AIC" = 2, "SC" = log(nPeriods), "HQ" = 2 * log(log(nPeriods)))
  # Compute the information criterions
  critValues <- sapply(Penalties, function(Penalty, logLik, nParams, nDiffus, nPeriods) {
    (-2 * logLik + Penalty * (nParams + nDiffus)) / nPeriods
  }, logLik, nParams, nDiffus, nPeriods)
  return(critValues)
}


#' @description Function that applies the InfCrit_fctn() to a row of the grid search output
#' @param outputVec row vector of the grid search output
#' @param nPeriods number of time periods
#' @return vector of information criteria

InfCritWrapper_fctn <- function(outputVec, nPeriods) {
  logLik <- outputVec[1]
  names(logLik) <- NULL
  nParams <- length(outputVec[-1])
  critValues <- InfCrit_fctn(logLik, nParams, nParams - 1, nPeriods)
  return(critValues)
}


#' @description Function to create distribution plots as part of the grid search output
#' @param paramterMat matrix holding the parameter vectors
#' @param title of the plot
#' @param path to the storage folder
#' @param model that underlies the grid search
#' @export distribution plot saved to the specified folder

PlotParameterDist_fctn <- function(paramterMat, title, path, model) {
  textSize <- 24
  # Prepare the data
  paramterMatLong <- paramterMat %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Estimate") %>%
    mutate(Parameter = factor(Parameter, levels = colnames(paramterMat)))
  # Set up the facet labels
  labelMap <- LabelMap_fctn()
  # Construct and save the faceted distribution plot
  Plot <- ggplot(paramterMatLong, aes(x = Parameter, y = Estimate)) +
    geom_violin() +
    facet_wrap(~Parameter,
      scales = "free",
      labeller = labeller(Parameter = LabelHelper_fctn(labelMap))
    ) +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_y_continuous(label = number_format(accuracy = 0.001)) +
    ggtitle(paste("Best 50 combinations |", title, "|", model)) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = textSize),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = textSize),
      axis.title = element_text(size = textSize + 3),
      strip.text = element_text(size = textSize),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks = element_line(linewidth = 1)
    )
  ggsave(Plot,
    filename = paste0(path, title, ".png"), width = 14, height = 8
  )
}
