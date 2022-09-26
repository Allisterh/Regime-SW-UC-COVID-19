### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Reproduces tables ########################################################################################################################### ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Tables_fctn <- function(Model_output, Model) {
  Path <- paste0(getwd(), "/", Model, "_Tables")
  if (!dir.exists(Path)) {
    dir.create(Path)
  }

  #---------------------------------------------------------------------------------------#
  # Table 1                                                                            ####
  #---------------------------------------------------------------------------------------#

  source("Scripts/SE_from_Hessian.R")
  # Load grid search results (constrained estimates)
  Grid_search_Results <- read_rds(file = paste0(Model, "_Grid_search/RESULTS_optim_cnstrnd_mat.rds"))
  # Load hessian matrix of the best parameter vector
  SE_list <- read_rds(file = paste0(Model, "_Grid_search/SE_list.rds"))
  # Build table 1
  ML_Estimates <- Grid_search_Results[1, -c(1, 5)] %>%
    enframe() %>%
    mutate(SE = SEfromHessian(SE_list$Hessian_mat, hessian = F, silent = F)) %>%
    rename(
      Parameter = name,
      Estimate = value
    ) %>%
    bind_rows(., tibble(Parameter = "nu_0", Estimate = Grid_search_Results[1, 5], SE = SE_list$nu_0_var)) %>%
    mutate_if(is.numeric, round, 3) %>%
    slice(match(c("xi", "omega", "epsilon", "nu_0", "nu_1", "p", "q"), Parameter))

  # Derive various regime indicators reported below table 1
  Avg_Duration_S_1 <- 1 / (1 - Grid_search_Results[1, "p"])
  Avg_Duration_S_0 <- 1 / (1 - Grid_search_Results[1, "q"])
  Avg_Change_S_1 <- (1 - exp(Grid_search_Results[1, "nu_0"] + Grid_search_Results[1, "nu_1"])) * 100
  Avg_Change_S_0 <- (exp(Grid_search_Results[1, "nu_0"]) - 1) * 100
  Halving_S_1 <- log(0.5) / (Grid_search_Results[1, "nu_0"] + (Grid_search_Results[1, "nu_1"]))
  Doubeling_S_0 <- log(2) / Grid_search_Results[1, "nu_0"]
  Regime_Characteristics <- tibble(
    Avg_Duration_S_1 = Avg_Duration_S_1, Avg_Duration_S_0 = Avg_Duration_S_0,
    Avg_Change_S_1 = Avg_Change_S_1, Avg_Change_S_0 = Avg_Change_S_0,
    Halving_S_1 = Halving_S_1, Doubeling_S_0 = Doubeling_S_0
  ) %>%
    mutate(
      across(!contains("Change"), ~ round(.x, 0)),
      across(contains("Change"), ~ round(.x, 2))
    ) %>%
    pivot_longer(names_to = "Features", cols = everything())

  # Store everything in a .xlsx book
  ML_book <- createWorkbook()
  addWorksheet(ML_book, "MLE_Estimates")
  addWorksheet(ML_book, "Log_likelihood")
  addWorksheet(ML_book, "Regime_Characteristics")
  writeData(ML_book, "MLE_Estimates", ML_Estimates)
  writeData(ML_book, "Log_likelihood", round(Grid_search_Results[1, 1], 3))
  writeData(ML_book, "Regime_Characteristics", Regime_Characteristics)
  saveWorkbook(ML_book, file = paste0(Path, "/Tab_1_MLE_Parameters.xlsx"), overwrite = TRUE)

  #---------------------------------------------------------------------------------------#
  # Table 2                                                                            ####
  #---------------------------------------------------------------------------------------#

  Upturn_Regimes_Smoother <- Model_output_list$Smoother_Regime_indicator_tib[, c(1, 3)] %>%
    mutate(Index = 1:n())
  Upturn_Regimes_Filter <- Model_output_list$Filter_Regime_indicator_tib[, c(1, 3)] %>%
    mutate(Index = 1:n())
  Upturn_book <- createWorkbook()
  addWorksheet(Upturn_book, "Smoother")
  addWorksheet(Upturn_book, "Filter")
  writeData(Upturn_book, "Smoother", Upturn_Regimes_Smoother)
  writeData(Upturn_book, "Filter", Upturn_Regimes_Filter)
  saveWorkbook(Upturn_book,
    file = paste0(Path, "/Tab_2_Infectionwaves.xlsx"),
    overwrite = TRUE
  )
}
