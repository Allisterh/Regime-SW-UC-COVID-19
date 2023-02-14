#' @description Function that executes the Gibbs Sampling routine
#' @param model that underlies the grid search
#' @param data vector for which the likelihood values are to be computed
#' @param trendIni trendInital value for the trend component
#' @param nIters number of iterations for the MCMC sampler
#' @param nBurnin number of burn in periods that are not being stored
#' @param dataTib tibble holding the dates
#' @param thinIntervall interval at which one specific sample is not stored
#' @param thinIntervallReverse boolean. I TRUE, every iteration corresponding to the thinIntervall value is stored. If FALSE, each
#'        concerning iteration is discarded.
#' @param storeOutput boolean. If FALSE, no output is stored and only the best fitting parameter
#'        vector returned
#' @param outputName character string that is appended to the sampling output
#' @return a list with the sampling output
#' @export a list with the sampling output

GibbsSampler_fctn <- function(model = modelSpec, data = dataVec, trend = trendIni, nIters, nBurnin, dataTib, thinIntervall = NULL,
                              thinIntervallReverse = FALSE, storeOutput = TRUE, outputName = NULL) {
  if (modelSpec != "TwoRegime") {
    outputList <- GibbsOneRegime_fctn(
      model = model, data = data, trend = trend, nIters = nIters, nBurnin = nBurnin, dataTib = dataTib,
      thinIntervall = thinIntervall, thinIntervallReverse = thinIntervallReverse, storeOutput = storeOutput,
      outputName = outputName
    )
  } else {
    outputList <- GibbsTwoRegime_fctn(
      model = model, data = data, trend = trend, nIters = nIters, nBurnin = nBurnin, dataTib = dataTib,
      thinIntervall = thinIntervall, thinIntervallReverse = thinIntervallReverse, storeOutput = storeOutput,
      outputName = outputName
    )
  }

  if (storeOutput == TRUE) {
    write_rds(outputList, file = paste0(path, "samplingOutput", outputName, ".rds"))
    # Run diagnostics on the sampling output
    cat("\n", "Running diagnostics")
    if (modelSpec != "TwoRegime") {
      GibbsDiagnostics_fctn(model = model, samplingOutput = outputList, dataTib = dataTib, outputName = outputName)
    } else {
      GibbsTwoRegimeDiagnostics_fctn(model = model, samplingOutput = outputList, dataTib = dataTib, outputName = outputName)
    }
  }
  return(outputList)
}



#' @description Function that executes the Gibbs Sampling routine for models with only one regime process
#' @param model that underlies the grid search
#' @param data vector for which the likelihood values are to be computed
#' @param trendIni trendInital value for the trend component
#' @param nIters number of iterations for the MCMC sampler
#' @param nBurnin number of burn in periods that are not being stored
#' @param dataTib tibble holding the dates
#' @param thinIntervall interval at which one specific sample is not stored
#' @param thinIntervallReverse boolean. I TRUE, every iteration corresponding to the thinIntervall value is stored. If FALSE, each
#'        concerning iteration is discarded.
#' @param storeOutput boolean. If FALSE, no output is stored and only the best fitting parameter
#'        vector returned
#' @param outputName character string that is appended to the sampling output
#' @return a list with the sampling output
#' @export a list with the sampling output


GibbsOneRegime_fctn <- function(model, data, trend, nIters, nBurnin, dataTib, thinIntervall,
                                thinIntervallReverse, storeOutput, outputName) {
  # Create an output folder
  path <- paste0(getwd(), "/Output/Output_", model, "/GibbsSampling/")
  if (!dir.exists(path) & storeOutput == TRUE) {
    dir.create(path)
  }
  paramNames <- colnames(ThetaRand_fctn(nRandom = 1))
  nParams <- length(paramNames)
  # Get some information about the model to be run
  probitIndicator <- ifelse(any(str_detect(paramNames, "beta")) == TRUE, TRUE, FALSE)

  #-------------------------------------------------------------------------------------#
  # Initialize the sampler                                                           ####
  #-------------------------------------------------------------------------------------#

  set.seed(2)
  nPeriods <- length(data)
  if (is.null(thinIntervall)) {
    thinIntervallVec <- c()
    nThin <- 0
  } else if (thinIntervall <= 1) {
    stop("thinIntervall must be an integer greater than 1.\n")
  } else {
    thinIntervall <- round(thinIntervall)
    thinIntervallVec <- seq(nBurnin, nIters, thinIntervall)
    if (thinIntervallReverse == TRUE) thinIntervallVec1 <- (1:nIters)[!((1:nIters) %in% thinIntervallVec)]
    nThin <- length(thinIntervallVec) - 1
  }
  storageCounter <- 0
  # Set up the filter
  systemList <- SystemMat_fctn()
  dimens <- NROW(systemList$Tt)
  # Initialize the sampler
  paramList <- ParList_fctn(rep(0, nParams), constrainPar = TRUE)
  regimeDraw <- rep(NA, nPeriods)
  regimeDraw[1] <- rbinom(1, 1, .5)
  if (probitIndicator == FALSE) {
    q <- paramList$Probs$q <- p <- paramList$Probs$p <- .95
  }
  for (i in 2:nPeriods) {
    prob <- ifelse(regimeDraw[i - 1] == 1, p, 1 - q)
    regimeDraw[i] <- rbinom(1, 1, prob)
  }
  # Initialize the output objects
  if (nBurnin >= nIters) {
    nBurnin <- 0
    warning("Burn in period longer than sample iterations. Setting the burn in period to zero. \n")
  }

  stateDrawArray <- array(NA, dim = c(dimens, nPeriods, nIters - nBurnin - nThin))
  regimeDrawMat <- matrix(NA, nr = nIters - nBurnin - nThin, nc = nPeriods)
  paramMat <- matrix(NA, nr = nIters - nBurnin - nThin, nc = nParams)
  colnames(paramMat) <- paramNames

  progressBar <- progress_bar$new(
    format = "|:bar| :percent ~:eta",
    total = nIters,
    complete = "+",
    incomplete = " ",
    current = " ",
    clear = FALSE,
    width = 100
  )

  # modelOutputList <- ModelApplication_fctn(model = modelSpec, data = dataVec, trendIni = trendIni, dataTib = dataTib)
  # paramList <- ParList_fctn(readRDS("C:/Users/phaim/OneDrive/COVID Paper/COVID-19-paper/Output/Output_UCModel/GridSearch/thetaVec.rds"), TRUE)
  # regimeDraw <- apply(modelOutputList$filterOutput$Prob_ct, 1, function(x) S <- ifelse(x[1] > .5, 0, 1))
  # stateDraw <- modelOutputList$smootherOutput$a_cT[,3,]

  cat("Starting the Gibbs sampling routine \n")

  for (i in 1:nIters) {
    progressBar$tick()

    #-------------------------------------------------------------------------------------#
    # Draw the state vector                                                            ####
    #-------------------------------------------------------------------------------------#

    # Get the updated filter values
    updateStateVec <- CondKalman_fctn(data = data, Ini = trend, systemList = systemList, paramList = paramList, regimeVec = regimeDraw, endogen = probitIndicator)

    # Run the Forward filter backwards sampling routine
    stateDraw <- BackwardsStateSampling_fctn(filterOutput = updateStateVec, systemList = systemList, paramList = paramList, regimeVec = regimeDraw)

    #-------------------------------------------------------------------------------------#
    # Draw the regime vector                                                           ####
    #-------------------------------------------------------------------------------------#

    # Get the updated filter values
    updateRegimeVec <- CondHamilton_fctn(stateVec = stateDraw, paramList = paramList, endogen = probitIndicator)

    # Sample from the filtered values
    regimeDraw <- RegimeSampling_fctn(filterOutput = updateRegimeVec, paramList = paramList, endogen = probitIndicator)

    #-------------------------------------------------------------------------------------#
    # Draw the additional parameters                                                   ####
    #-------------------------------------------------------------------------------------#

    paramList <- Additionalparam_fctn(data = data, stateVec = stateDraw, regimeVec = regimeDraw, paramList = paramList)

    # # draw betas and varrho
    # betaVec <- Beta_fctn(regimeVec = regimeDraw)
    # paramList$beta$beta_0 <- betaVec[1]
    # paramList$beta$beta_1 <- betaVec[2]
    # # Transfer probit coefficients into probabilities
    # transProbs <- c(pnorm(-betaVec[1]), 1 - pnorm(-betaVec[1] - betaVec[2]))

    # Store the sampled values
    if (i > nBurnin) {
      if (!(i %in% thinIntervallVec)) {
        storageCounter <- storageCounter + 1
        stateDrawArray[1:dimens, 1:nPeriods, storageCounter] <- stateDraw
        regimeDrawMat[storageCounter, ] <- regimeDraw
        paramVec <- outputParam_fctn(paramList)
        paramMat[storageCounter, ] <- paramVec
      }
    }
  }
  outputList <- list(
    "stateDraw" = stateDrawArray,
    "regimeDraw" = regimeDrawMat,
    "additionalParam" = paramMat
  )
  if (storeOutput == TRUE) {
    write_rds(outputList, file = paste0(path, "samplingOutput", outputName, ".rds"))
    # Run diagnostics on the sampling output
    cat("\n", "Running diagnostics")
    GibbsDiagnostics_fctn(model = model, samplingOutput = outputList, dataTib = dataTib, outputName = outputName)
  }
  return(outputList)
}


#' @description Function that executes the Gibbs Sampling routine for models with only one regime process
#' @param model that underlies the grid search
#' @param data vector for which the likelihood values are to be computed
#' @param trendIni trendInital value for the trend component
#' @param nIters number of iterations for the MCMC sampler
#' @param nBurnin number of burn in periods that are not being stored
#' @param dataTib tibble holding the dates
#' @param thinIntervall interval at which one specific sample is not stored
#' @param thinIntervallReverse boolean. I TRUE, every iteration corresponding to the thinIntervall value is stored. If FALSE, each
#'        concerning iteration is discarded.
#' @param storeOutput boolean. If FALSE, no output is stored and only the best fitting parameter
#'        vector returned
#' @param outputName character string that is appended to the sampling output
#' @return a list with the sampling output
#' @export a list with the sampling output


GibbsTwoRegime_fctn <- function(model, data, trend, nIters, nBurnin, dataTib, thinIntervall,
                                thinIntervallReverse, storeOutput, outputName) {
  # Create an output folder
  path <- paste0(getwd(), "/Output/Output_", model, "/GibbsSampling/")
  if (!dir.exists(path) & storeOutput == TRUE) {
    dir.create(path)
  }
  paramNames <- colnames(ThetaRand_fctn(nRandom = 1))
  nParams <- length(paramNames)

  #-------------------------------------------------------------------------------------#
  # Initialize the sampler                                                           ####
  #-------------------------------------------------------------------------------------#

  set.seed(2)
  nPeriods <- length(data)
  if (is.null(thinIntervall)) {
    thinIntervallVec <- c()
    nThin <- 0
  } else if (thinIntervall <= 1) {
    stop("thinIntervall must be an integer greater than 1.\n")
  } else {
    thinIntervall <- round(thinIntervall)
    thinIntervallVec <- seq(nBurnin, nIters, thinIntervall)
    if (thinIntervallReverse == TRUE) thinIntervallVec1 <- (1:nIters)[!((1:nIters) %in% thinIntervallVec)]
    nThin <- length(thinIntervallVec) - 1
  }
  storageCounter <- 0
  # Set up the filter
  systemList <- SystemMat_fctn()
  dimens <- NROW(systemList$Tt)
  # Initialize the sampler
  paramList <- ParList_fctn(rep(0, nParams), constrainPar = TRUE)
  regimeDriftDraw <- regimeHeteroscDraw <- rep(NA, nPeriods)
  regimeDriftDraw[1] <- regimeHeteroscDraw[1] <- rbinom(1, 1, .5)
  if (probitIndicator == FALSE) {
    q <- paramList$Probs$q <- p <- paramList$Probs$p <- r <- paramList$Probs$r <- s <- paramList$Probs$s <- .95
  }
  for (i in 2:nPeriods) {
    probDrift <- ifelse(regimeDriftDraw[i - 1] == 1, p, 1 - q)
    regimeDriftDraw[i] <- rbinom(1, 1, probDrift)
    probHeterosc <- ifelse(regimeHeteroscDraw[i - 1] == 1, r, 1 - s)
    regimeHeteroscDraw[i] <- rbinom(1, 1, probHeterosc)
  }
  # Initialize the output objects
  if (nBurnin >= nIters) {
    nBurnin <- 0
    warning("Burn in period longer than sample iterations. Setting the burn in period to zero. \n")
  }

  stateDrawArray <- array(NA, dim = c(dimens, nPeriods, nIters - nBurnin - nThin))
  regimeDrawDriftMat <- matrix(NA, nr = nIters - nBurnin - nThin, nc = nPeriods)
  regimeDrawHeteroscMat <- matrix(NA, nr = nIters - nBurnin - nThin, nc = nPeriods)
  paramMat <- matrix(NA, nr = nIters - nBurnin - nThin, nc = nParams)
  colnames(paramMat) <- paramNames

  progressBar <- progress_bar$new(
    format = "|:bar| :percent ~:eta",
    total = nIters,
    complete = "+",
    incomplete = " ",
    current = " ",
    clear = FALSE,
    width = 100
  )

  # modelOutputList <- ModelApplication_fctn(model = modelSpec, data = dataVec, trendIni = trendIni, dataTib = dataTib)
  # paramList <- ParList_fctn(readRDS("C:/Users/phaim/OneDrive/COVID Paper/COVID-19-paper/Output/Output_UCModel/GridSearch/thetaVec.rds"), TRUE)
  # regimeDraw <- apply(modelOutputList$filterOutput$Prob_ct, 1, function(x) S <- ifelse(x[1] > .5, 0, 1))
  # stateDraw <- modelOutputList$smootherOutput$a_cT[,3,]

  cat("Starting the Gibbs sampling routine \n")

  for (i in 1:nIters) {
    progressBar$tick()

    #-------------------------------------------------------------------------------------#
    # Draw the state vector                                                            ####
    #-------------------------------------------------------------------------------------#

    # Get the updated filter values
    updateStateVec <- CondTwoRegimeKalman_fctn(
      data = data, Ini = trend, systemList = systemList, paramList = paramList,
      regimeDriftVec = regimeDriftDraw, regimeHeteroscVec = regimeHeteroscDraw
    )

    # Run the Forward filter backwards sampling routine
    stateDraw <- BackwardsTwoRegimeStateSampling_fctn(
      filterOutput = updateStateVec, systemList = systemList, paramList = paramList,
      regimeDriftVec = regimeDriftDraw, regimeHeteroscVec = regimeHeteroscDraw
    )

    #-------------------------------------------------------------------------------------#
    # Draw the regime vectors                                                          ####
    #-------------------------------------------------------------------------------------#

    # Get the updated filter values for the regime process governing the drift
    updateRegimeDriftVec <- CondDriftHamilton_fctn(stateVec = stateDraw, regimeVec = regimeHeteroscDraw, paramList = paramList)
    # Sample from the filtered values
    regimeDriftDraw <- RegimeSampling_fctn(filterOutput = updateRegimeDriftVec, paramList = paramList)

    # Get the updated filter values for the regime process governing the heteroscedasticity
    updateRegimeHeteroscVec <- CondHeteroscHamilton_fctn(stateVec = stateDraw, regimeVec = regimeDriftDraw, paramList = paramList)
    # Sample from the filtered values
    regimeHeteroscDraw <- RegimeSampling_fctn(filterOutput = updateRegimeHeteroscVec, paramList = paramList)

    #-------------------------------------------------------------------------------------#
    # Draw the additional parameters                                                   ####
    #-------------------------------------------------------------------------------------#

    paramList <- Additionalparam_fctn(
      data = data, stateVec = stateDraw, regimeDriftVec = regimeDriftDraw, regimeHeteroscVec = regimeHeteroscDraw,
      paramList = paramList
    )

    # Store the sampled values
    if (i > nBurnin) {
      if (!(i %in% thinIntervallVec)) {
        storageCounter <- storageCounter + 1
        stateDrawArray[1:dimens, 1:nPeriods, storageCounter] <- stateDraw
        regimeDriftDrawMat[storageCounter, ] <- regimeDriftDraw
        regimeHeteroscDrawMat[storageCounter, ] <- regimeHeteroscDraw
        paramVec <- outputParam_fctn(paramList)
        paramMat[storageCounter, ] <- paramVec
      }
    }
  }
  outputList <- list(
    "stateDraw" = stateDrawArray,
    "regimeDriftDraw" = regimeDriftDrawMat,
    "regimeHeteroscDraw" = regimeHeteroscDraw,
    "additionalParam" = paramMat
  )
  if (storeOutput == TRUE) {
    write_rds(outputList, file = paste0(path, "samplingOutput", outputName, ".rds"))
    # Run diagnostics on the sampling output
    cat("\n", "Running diagnostics")
    GibbsTwoRegimeDiagnostics_fctn(model = model, samplingOutput = outputList, dataTib = dataTib, outputName = outputName)
  }
  return(outputList)
}


#' @description Function that evaluates the Gibbs Sampler routine
#' @param model that underlies the grid search
#' @param samplingOutput Output from the GibbsSampler_fctn
#' @param dataTib tibble holding the dates
#' @param outputName character string that is appended to the sampling output
#' @export trace plots, tables with the sampling results and convergence diagnostics

GibbsDiagnostics_fctn <- function(model, samplingOutput = NULL, dataTib, outputName) {
  path <- paste0(getwd(), "/Output/Output_", model, "/GibbsSampling/")
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Load output if none is provided as an argument
  if (is.null(samplingOutput)) {
    if ("samplingOutput.rds" %in% list.files(path = path)) {
      samplingOutput <- readRDS(file = paste0(path, "samplingOutput.rds"))
    } else {
      stop("No Gibbs Sampling output stored and no output provided.\n")
    }
  }

  #-------------------------------------------------------------------------------------#
  # Create trace plots                                                               ####
  #-------------------------------------------------------------------------------------#

  # Define aesthetics
  orange <- "#FF6347"
  black <- "#000000"
  grey <- "#525252"
  blue <- "#1874CD"
  textSize <- 24
  themeElement <- theme(
    legend.position = "none",
    legend.text = element_text(size = textSize),
    axis.text = element_text(size = textSize),
    axis.title = element_text(size = textSize + 3, vjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(linewidth = 1)
  )
  themeElementFacet <- theme(
    legend.position = "none",
    axis.text.y = element_text(size = textSize),
    axis.text.x = element_text(size = textSize),
    axis.title = element_text(size = textSize + 3),
    strip.text = element_text(size = textSize),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(linewidth = 1)
  )

  # Transpose the sampling output into a long format
  paramMat <- samplingOutput$additionalParam
  labelMap <- LabelMap_fctn()
  paramTib <- as_tibble(paramMat) %>%
    mutate(Iteration = 1:n()) %>%
    pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = names(labelMap)))
  # Construct the figure
  tracePlot <- paramTib %>%
    ggplot() +
    theme_bw() +
    geom_line(aes(x = Iteration, y = Value), linewidth = 1) +
    facet_wrap(~Parameter,
      scales = "free", ncol = 2,
      labeller = labeller(Parameter = LabelHelper_fctn(labelMap))
    ) +
    scale_y_continuous(name = "") +
    themeElementFacet
  ggsave(tracePlot, filename = paste0(path, "/Trace_Plot", outputName, ".png"), height = 5, width = 14)

  #-------------------------------------------------------------------------------------#
  # CumSum plots                                                                     ####
  #-------------------------------------------------------------------------------------#

  # Prepare the data
  paramCumSumMat <- apply(paramMat, 2, function(x) {
    mean <- mean(x)
    cumSum <- cumsum(x - mean)
  })
  paramCumSumTib <- as_tibble(paramCumSumMat) %>%
    mutate(Iteration = 1:n()) %>%
    pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = names(labelMap)))
  # Construct the figure
  cumSumPlot <- paramCumSumTib %>%
    ggplot() +
    theme_bw() +
    geom_line(aes(x = Iteration, y = Value), linewidth = 1) +
    geom_hline(yintercept = 0, color = blue, linetype = 2, linewidth = 1) +
    facet_wrap(~Parameter,
      scales = "free", ncol = 2,
      labeller = labeller(Parameter = LabelHelper_fctn(labelMap))
    ) +
    scale_y_continuous(name = "") +
    themeElementFacet
  ggsave(cumSumPlot, filename = paste0(path, "/CumSum_Plot", outputName, ".png"), height = 5, width = 14)

  #-------------------------------------------------------------------------------------#
  # Produce descriptive statistics of the posteriors                                 ####
  #-------------------------------------------------------------------------------------#

  # Additional parameters
  summaryAdditionalParam <- apply(paramMat, 2, function(x) {
    mean <- mean(x)
    sd <- sd(x)
    median <- median(x)
    dens <- density(x)
    mod <- dens$x[which.max(dens$x)]
    quantiles <- quantile(x, probs = c(.025, .975))
    outPutVec <- c("Mean" = mean, "Std.Dev" = sd, "Modus" = mod, "Median" = median, quantiles[1], quantiles[2])
    return(outPutVec)
  }) %>%
    t()

  # Save to excel
  write.xlsx(as_tibble(cbind("Parameter" = rownames(summaryAdditionalParam), round(summaryAdditionalParam, 3))),
    file = paste0(path, "AdditionalParameter", outputName, ".xlsx")
  )

  # State vector
  stateVecArray <- samplingOutput$stateDraw
  # Extract the mean and quantiles for every period across all iterations
  stateMeanMat <- apply(stateVecArray, 2, function(x) {
    apply(x, 1, median)
  })
  stateQuantileMat <- apply(stateVecArray, 2, function(x) {
    apply(x, 1, quantile, probs = c(.25, .75))
  })

  # Regimes
  regimeMat <- samplingOutput$regimeDraw
  regimeMeanVec <- apply(regimeMat, 2, mean)

  #-------------------------------------------------------------------------------------#
  # Replicate figure 3                                                               ####
  #-------------------------------------------------------------------------------------#

  outputTib <- tibble(
    Trend = stateMeanMat[1, ],
    UTrendQ = stateQuantileMat[2, ],
    LTrendQ = stateQuantileMat[1, ],
    Pr_S_0 = 1 - regimeMeanVec,
    Date = dataTib$Date,
    Data = dataTib$logI
  )
  RegimeIndicatorTib <- RegimeIdent_fctn(outputTib = outputTib, probVar = "Pr_S_0", threshold = 0.4)
  colors <- c(
    "Data" = grey, "Trend" = orange, "Pr regime 1" = blue,
    ".025Trend" = black, ".975Trend" = black
  )
  regimePlot <- outputTib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Data, colour = "Data"), linewidth = 0.7) +
    geom_line(aes(x = Date, y = Trend, colour = "Trend"), linewidth = 1) +
    geom_line(aes(x = Date, y = Pr_S_0 * max(Data, na.rm = T), colour = "Pr regime 1"), linewidth = 1) +
    geom_rect(
      data = RegimeIndicatorTib, inherit.aes = F,
      aes(ymin = 0, ymax = max(outputTib$Data, na.rm = T), xmin = Date, xmax = End),
      fill = "black", alpha = 0.3
    ) +
    scale_y_continuous(
      name = expression(paste(log, "(", i[t], ")")),
      breaks = seq(0, 14, by = 2),
      sec.axis = sec_axis(~ . * (max(outputTib$Pr_S_0, na.rm = T) / max(outputTib$Data, na.rm = T)),
        name = expression(paste(Pr, " (", S[t] == 0, ")")),
        breaks = seq(0, 1, by = 0.2)
      ),
      limits = c(0, max(outputTib$Data, na.rm = T))
    ) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    themeElement
  ggsave(regimePlot, filename = paste0(path, "/Regime_plot", outputName, ".png"), height = 5, width = 14)

  #-------------------------------------------------------------------------------------#
  # Convergence criteria                                                             ####
  #-------------------------------------------------------------------------------------#

  nIters <- NROW(samplingOutput$additionalParam)
  mcmcOutput1 <- mcmc(samplingOutput$additionalParam[1:(nIters / 2), ])
  mcmcOutput2 <- mcmc(samplingOutput$additionalParam[(1 + nIters / 2):nIters, ])
  mcmcOutputList <- mcmc.list(list(mcmcOutput1, mcmcOutput2))
  # gelmanInd <- gelman.diag(mcmcOutputList)$mpsrf
  # if (gelmanInd >= 1.1) {
  #   cat("Gelman-Rubin criterion ", round(gelmanInd, 4), ". Convergence likely not archieved. \n")
  # } else {
  #   cat("Gelman-Rubin criterion ", round(gelmanInd, 4), ". Convergence archieved. \n")
  # }
}


#' @description Function that evaluates the Gibbs Sampler routine entailing two regime processes
#' @param model that underlies the grid search
#' @param samplingOutput Output from the GibbsSampler_fctn
#' @param dataTib tibble holding the dates
#' @param outputName character string that is appended to the sampling output
#' @export trace plots, tables with the sampling results and convergence diagnostics

GibbsTwoRegimeDiagnostics_fctn <- function(model, samplingOutput = NULL, dataTib, outputName) {
  path <- paste0(getwd(), "/Output/Output_", model, "/GibbsSampling/")
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Load output if none is provided as an argument
  if (is.null(samplingOutput)) {
    if ("samplingOutput.rds" %in% list.files(path = path)) {
      samplingOutput <- readRDS(file = paste0(path, "samplingOutput.rds"))
    } else {
      stop("No Gibbs Sampling output stored and no output provided.\n")
    }
  }

  #-------------------------------------------------------------------------------------#
  # Create trace plots                                                               ####
  #-------------------------------------------------------------------------------------#

  # Define aesthetics
  orange <- "#FF6347"
  black <- "#000000"
  grey <- "#525252"
  blue <- "#1874CD"
  green <- "#00b159"
  textSize <- 24
  themeElement <- theme(
    legend.position = "none",
    legend.text = element_text(size = textSize),
    axis.text = element_text(size = textSize),
    axis.title = element_text(size = textSize + 3, vjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(linewidth = 1)
  )
  themeElementFacet <- theme(
    legend.position = "none",
    axis.text.y = element_text(size = textSize),
    axis.text.x = element_text(size = textSize),
    axis.title = element_text(size = textSize + 3),
    strip.text = element_text(size = textSize),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(linewidth = 1)
  )

  # Transpose the sampling output into a long format
  paramMat <- samplingOutput$additionalParam
  labelMap <- LabelMap_fctn()
  paramTib <- as_tibble(paramMat) %>%
    mutate(Iteration = 1:n()) %>%
    pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = names(labelMap)))
  # Construct the figure
  tracePlot <- paramTib %>%
    ggplot() +
    theme_bw() +
    geom_line(aes(x = Iteration, y = Value), linewidth = 1) +
    facet_wrap(~Parameter,
      scales = "free", ncol = 2,
      labeller = labeller(Parameter = LabelHelper_fctn(labelMap))
    ) +
    scale_y_continuous(name = "") +
    themeElementFacet
  ggsave(tracePlot, filename = paste0(path, "/Trace_Plot", outputName, ".png"), height = 5, width = 14)

  #-------------------------------------------------------------------------------------#
  # CumSum plots                                                                     ####
  #-------------------------------------------------------------------------------------#

  # Prepare the data
  paramCumSumMat <- apply(paramMat, 2, function(x) {
    mean <- mean(x)
    cumSum <- cumsum(x - mean)
  })
  paramCumSumTib <- as_tibble(paramCumSumMat) %>%
    mutate(Iteration = 1:n()) %>%
    pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = names(labelMap)))
  # Construct the figure
  cumSumPlot <- paramCumSumTib %>%
    ggplot() +
    theme_bw() +
    geom_line(aes(x = Iteration, y = Value), linewidth = 1) +
    geom_hline(yintercept = 0, color = blue, linetype = 2, linewidth = 1) +
    facet_wrap(~Parameter,
      scales = "free", ncol = 2,
      labeller = labeller(Parameter = LabelHelper_fctn(labelMap))
    ) +
    scale_y_continuous(name = "") +
    themeElementFacet
  ggsave(cumSumPlot, filename = paste0(path, "/CumSum_Plot", outputName, ".png"), height = 5, width = 14)

  #-------------------------------------------------------------------------------------#
  # Produce descriptive statistics of the posteriors                                 ####
  #-------------------------------------------------------------------------------------#

  # Additional parameters
  summaryAdditionalParam <- apply(paramMat, 2, function(x) {
    mean <- mean(x)
    sd <- sd(x)
    median <- median(x)
    dens <- density(x)
    mod <- dens$x[which.max(dens$x)]
    quantiles <- quantile(x, probs = c(.025, .975))
    outPutVec <- c("Mean" = mean, "Std.Dev" = sd, "Modus" = mod, "Median" = median, quantiles[1], quantiles[2])
    return(outPutVec)
  }) %>%
    t()

  # Save to excel
  write.xlsx(as_tibble(cbind("Parameter" = rownames(summaryAdditionalParam), round(summaryAdditionalParam, 3))),
    file = paste0(path, "AdditionalParameter", outputName, ".xlsx")
  )

  # State vector
  stateVecArray <- samplingOutput$stateDraw
  # Extract the mean and quantiles for every period across all iterations
  stateMeanMat <- apply(stateVecArray, 2, function(x) {
    apply(x, 1, median)
  })
  stateQuantileMat <- apply(stateVecArray, 2, function(x) {
    apply(x, 1, quantile, probs = c(.25, .75))
  })

  # Regimes
  regimeDriftMat <- samplingOutput$regimeDriftDraw
  regimeDriftMeanVec <- apply(regimeDriftMat, 2, mean)
  regimeHeteroscMat <- samplingOutput$regimeHeteroscDraw
  regimeHeteroscMeanVec <- apply(regimeHeteroscMat, 2, mean)

  #-------------------------------------------------------------------------------------#
  # Replicate figure 3                                                               ####
  #-------------------------------------------------------------------------------------#

  outputTib <- tibble(
    Trend = stateMeanMat[1, ],
    UTrendQ = stateQuantileMat[2, ],
    LTrendQ = stateQuantileMat[1, ],
    Drift_Pr_S_0 = 1 - regimeDriftMeanVec,
    Heterosc_Pr_S_0 = 1 - regimeHeteroscMeanVec,
    Date = dataTib$Date,
    Data = dataTib$logI
  )
  RegimeIndicatorTib <- RegimeIdent_fctn(outputTib = outputTib, probVar = "Drift_Pr_S_0", threshold = 0.4)
  colors <- c(
    "Data" = grey, "Trend" = orange, "Drift Pr regime 1" = blue, "Heterosc Pr regime 1" = green,
    ".025Trend" = black, ".975Trend" = black
  )
  regimePlot <- outputTib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Data, colour = "Data"), linewidth = 0.7) +
    geom_line(aes(x = Date, y = Trend, colour = "Trend"), linewidth = 1) +
    geom_line(aes(x = Date, y = Drift_Pr_S_0 * max(Data, na.rm = T), colour = "Drift Pr regime 1"), linewidth = 1) +
    geom_line(aes(x = Date, y = Heterosc_Pr_S_0 * max(Data, na.rm = T), colour = "Heterosc Pr regime 1"), linewidth = 1) +
    geom_rect(
      data = RegimeIndicatorTib, inherit.aes = F,
      aes(ymin = 0, ymax = max(outputTib$Data, na.rm = T), xmin = Date, xmax = End),
      fill = "black", alpha = 0.3
    ) +
    scale_y_continuous(
      name = expression(paste(log, "(", i[t], ")")),
      breaks = seq(0, 14, by = 2),
      sec.axis = sec_axis(~ . * (max(outputTib$Drift_Pr_S_0, na.rm = T) / max(outputTib$Data, na.rm = T)),
        name = expression(paste(Pr, " (", S[t] == 1, ")")),
        breaks = seq(0, 1, by = 0.2)
      ),
      limits = c(0, max(outputTib$Data, na.rm = T))
    ) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    themeElement
  ggsave(regimePlot, filename = paste0(path, "/Regime_plot", outputName, ".png"), height = 5, width = 14)

  #-------------------------------------------------------------------------------------#
  # Convergence criteria                                                             ####
  #-------------------------------------------------------------------------------------#

  nIters <- NROW(samplingOutput$additionalParam)
  mcmcOutput1 <- mcmc(samplingOutput$additionalParam[1:(nIters / 2), ])
  mcmcOutput2 <- mcmc(samplingOutput$additionalParam[(1 + nIters / 2):nIters, ])
  mcmcOutputList <- mcmc.list(list(mcmcOutput1, mcmcOutput2))
  # gelmanInd <- gelman.diag(mcmcOutputList)$mpsrf
  # if (gelmanInd >= 1.1) {
  #   cat("Gelman-Rubin criterion ", round(gelmanInd, 4), ". Convergence likely not archieved. \n")
  # } else {
  #   cat("Gelman-Rubin criterion ", round(gelmanInd, 4), ". Convergence archieved. \n")
  # }
}
