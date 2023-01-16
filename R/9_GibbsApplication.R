GibbsSampler_fctn <- function(data = dataVec, Ini = trendIni, endogen = FALSE) {
  #----------------------------#
  ### Initialize the sampler ###
  #----------------------------#
  
  nPeriods <- length(data)
  # Set up the filter
  systemList <- SystemMat_fctn()
  # Initialize the sampler
  paramList <- ParList_fctn(thetaVec, TRUE)
  regimeDraw <- rbinom(nPeriods, 1, .5)
  
  #---------------------------#
  ### Draw the state vector ###
  #---------------------------#
  
  # Get the updated filter values
  updateStateVec <- CondKalman_fctn(data = data, Ini = Ini, systemList = systemList, paramList = paramList, regimeVec = regimeDraw, endogen = endogen)
  # Run the Forward filter backwards sampling routine
  stateDraw <- BackwardsStateSampling_fctn(filterOutput = updateStateVec, systemList = systemList, paramList = paramList, regimeVec = regimeDraw)
  
  #----------------------------#
  ### Draw the regime vector ###
  #----------------------------#
  
  # Get the updated filter values
  updateRegimeVec <- CondHamilton_fctn(stateVec = stateDraw, paramList = paramList, endogen = endogen)
  # Sample from the filtered values
  regimeDraw <- RegimeSampling_fctn(filterOutput = updateRegimeVec, paramList = paramList, endogen = endogen)
  
  #--------------------------------#
  ### Draw additional parameters ###
  #--------------------------------#
  
  # epsilon
  paramList$epsilon <- Sdepsilon_fctn(data = data, stateDraw = stateDraw, regimeDraw = regimeDraw)
  # xi
  paramList$xi$xi_1 <- paramList$xi$xi_0 <- SdXi_fctn(data = data, stateDraw = stateDraw, regimeDraw = regimeDraw, paramList = paramList)
  # omega
  paramList$omega <- SdOmega_fctn(stateDraw)
  # nu_1
  paramList$nu_1 <- Nu_1_fctn(stateDraw = stateDraw, regimeDraw = regimeDraw, paramList = paramList)
  
}
