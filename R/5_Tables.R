#' @description Function to create the tables 1 to 2
#' @param model for which the tables are created
#' @param modelOutput from which the tables are created
#' @export Tables saved in the ´model´_Output/Tables folder

Tables_fctn <- function(model, modelOutput) {
  path <- paste0(getwd(), "/Output/Output_", model, "/Tables")
  if (!dir.exists(path)) {
    dir.create(path)
  }

  #---------------------------------------------------------------------------------------#
  # Table 1                                                                            ####
  #---------------------------------------------------------------------------------------#

  # Load grid search results
  gridSearchResults <- read_rds(file = paste0("Output/Output_", model, "/gridSearch/resultsOptimCnstrndMat.rds"))
  # Build table 1
  mlEstimates <- gridSearchResults[1, -c(1, (NCOL(gridSearchResults) - 2):NCOL(gridSearchResults))] %>%
    enframe() %>%
    mutate(SE = modelOutput$SE) %>%
    rename(
      Parameter = name,
      Estimate = value
    ) %>%
    mutate_if(is.numeric, round, 3) %>%
    slice(match(names(LabelMap_fctn()), Parameter))
  # Derive various regime indicators reported below table 1
  if ("q" %in% colnames(gridSearchResults)) {
    avgDuration_S_1 <- 1 / (1 - gridSearchResults[1, "p"])
    avgDuration_S_0 <- 1 / (1 - gridSearchResults[1, "q"])
  } else {
    avgDuration_S_1 <- 1 / (1 - gridSearchResults[1, "p_11"])
    avgDuration_S_0 <- 1 / (1 - gridSearchResults[1, "p_00"])
  }
  avgChange_S_1 <- (1 - exp(gridSearchResults[1, "nu_0"] + gridSearchResults[1, "nu_1"])) * 100
  avgChange_S_0 <- (exp(gridSearchResults[1, "nu_0"]) - 1) * 100
  halving_S_1 <- log(0.5) / (gridSearchResults[1, "nu_0"] + (gridSearchResults[1, "nu_1"]))
  doubeling_S_0 <- log(2) / gridSearchResults[1, "nu_0"]
  Regime_Characteristics <- tibble(
    avgDuration_S_1 = avgDuration_S_1, avgDuration_S_0 = avgDuration_S_0,
    avgChange_S_1 = avgChange_S_1, avgChange_S_0 = avgChange_S_0,
    halving_S_1 = halving_S_1, doubeling_S_0 = doubeling_S_0
  ) %>%
    mutate(
      across(!contains("Change"), ~ round(.x, 0)),
      across(contains("Change"), ~ round(.x, 2))
    ) %>%
    pivot_longer(names_to = "Features", cols = everything())

  # Store everything in a .xlsx book
  mlBook <- createWorkbook()
  addWorksheet(mlBook, "mleEstimates")
  addWorksheet(mlBook, "Diagnostics")
  addWorksheet(mlBook, "regimeCharacteristics")
  writeData(mlBook, "mleEstimates", mlEstimates)
  writeData(mlBook, "Diagnostics", round(gridSearchResults[1, c(1, (NCOL(gridSearchResults) - 2):NCOL(gridSearchResults))], 3))
  writeData(mlBook, "regimeCharacteristics", Regime_Characteristics)
  saveWorkbook(mlBook, file = paste0(path, "/Tab_1_MLE_Parameters.xlsx"), overwrite = TRUE)

  #---------------------------------------------------------------------------------------#
  # Table 2                                                                            ####
  #---------------------------------------------------------------------------------------#

  upturnRegimesSmoother <- modelOutput$smootherRegimeIndicatorTib[, c(1, 3)] %>%
    mutate(Index = 1:n())
  upturnRegimesFilter <- modelOutput$filterRegimeIndicatorTib[, c(1, 3)] %>%
    mutate(Index = 1:n())
  upturnBook <- createWorkbook()
  addWorksheet(upturnBook, "Smoother")
  addWorksheet(upturnBook, "Filter")
  writeData(upturnBook, "Smoother", upturnRegimesSmoother)
  writeData(upturnBook, "Filter", upturnRegimesFilter)
  saveWorkbook(upturnBook,
    file = paste0(path, "/Tab_2_Infectionwaves.xlsx"),
    overwrite = TRUE
  )
}
