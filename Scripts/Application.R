### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Apply the model ############################################################################################################################# ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

#---------------------------------------------------------------------------------------#
# Helper function to create a table of estimated up-turning periods                  ####
#---------------------------------------------------------------------------------------#

Regime_ident_fctn <- function(Output_tib, Prob_var, Threshold){
  
  Regime_periods <- Output_tib[Output_tib[,Prob_var] < (1 - Threshold), c("Date", Prob_var)] %>%
    mutate(End = Date + 1,
           lag_end = as.Date(lag(End)),
           lead_start = as.Date(lead(Date)),
           type = ifelse(is.na(lag_end) | Date != lag_end & End == lead_start, "Start",
                         ifelse(is.na(lead_start) | Date == lag_end & End != lead_start, "Last",
                                ifelse(Date == lag_end & End == lead_start, "Mid", "Lonely"))
           ),
           type = ifelse(is.na(lead_start) & is.na(type), "Lonely", type)) %>%
    filter(type != "Mid") %>%
    mutate(End = lead(Date)) %>%
    filter(type != "Last")
  Regime_periods[Regime_periods$type == "Lonely", "End"] <- Regime_periods[Regime_periods$type == "Lonely", "Date"] + 1
  
  return(Regime_periods)
}

#---------------------------------------------------------------------------------------#
# Function to run the model                                                          ####
#---------------------------------------------------------------------------------------#

Model_application_fctn <- function(Model){

# Collect data, model functions and parameter vector
source(paste0("Scripts/", Model, ".R"))
THETA_vec <- read_rds(file = paste0(Model, "_Grid_search/THETA_vec.rds"))
ML_params <- THETA_vec %>% Par_constrain_fctn()

# Apply Kim Filter recursions
Filter_output <- Filter_fctn(
  par = ML_params,
  data = Data_vec
)

# Apply Smoother recursions
Smoother_output <- Smoother_fctn(
  par = ML_params,
  Filter_output = Filter_output,
  level = 0.9
)

#---------------------------------------------------------------------------------------#
# Prepare output
#---------------------------------------------------------------------------------------#

# Collect the output used for figures 3 and 4 as well as model diagnostics
Output_tib <- tibble(
  Trend = Smoother_output[["a_cT"]][1, 3, ],
  Smoother_Pr_S_1 = Smoother_output[["Pr_cT"]][, 2],
  Filter_Pr_S_1 = Filter_output$Prob_clt[, 2],
  Trend_U_CI = Smoother_output[["a_CI_cT"]][1, 1, ],
  Trend_L_CI = Smoother_output[["a_CI_cT"]][1, 2, ],
  Date = Data$Date,
  Data = Data_vec
) %>%
  mutate(
    Diff = Data - Trend) %>%
  slice(-c(1:5, Periods))

# Tibbles indicating when the infections up-turning regime is considered to be active
Smoother_Regime_indicator_tib <- Regime_ident_fctn(Output_tib = Output_tib, Prob_var = "Smoother_Pr_S_1", Threshold = 0.4)
Filter_Regime_indicator_tib <- Regime_ident_fctn(Output_tib = Output_tib, Prob_var = "Filter_Pr_S_1", Threshold = 0.4)

return(list(Filter_output = Filter_output,
            Output_tib = Output_tib,
            Smoother_output = Smoother_output,
            Smoother_Regime_indicator_tib = Smoother_Regime_indicator_tib,
            Filter_Regime_indicator_tib = Filter_Regime_indicator_tib))

}
