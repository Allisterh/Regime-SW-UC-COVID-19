### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Load data and packages ##################################################################################################################### ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Init_fctn <- function(Start_Date, End_Date, Country = "US") {

  # Install and or load packages
  if (!require(tidyverse)) install.packages("tidyverse")
  if (!require(readr)) install.packages("readr")
  if (!require(lubridate)) install.packages("lubridate")
  if (!require(scales)) install.packages("scales")
  if (!require(tictoc)) install.packages("tictoc")
  if (!require(Rcpp)) install.packages("Rcpp")
  if (!require(RcppArmadillo)) install.packages("RcppArmadillo")
  if (!require(openxlsx)) install.packages("openxlsx")
  if (!require(coda)) install.packages("coda")
  if (!require(Matrix)) install.packages("Matrix")

  JHCSSE_data <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
  # Wrangle data
  Data_raw <<- JHCSSE_data %>%
    filter(`Country/Region` == Country) %>%
    select(-c(`Country/Region`, `Province/State`, Lat, Long)) %>%
    pivot_longer(cols = 1:ncol(.), values_to = "C", names_to = "Date") %>%
    mutate(
      Date = mdy(Date),
      I = C - lag(C),
      log_I = log(I + 1),
      index = 1:n()
    )

  Data <<- Data_raw %>%
    filter(Date >= Start_Date & Date <= End_Date)
  Data_vec <<- Data$log_I %>% as.matrix()
  # Number of periods
  Periods <<- nrow(Data_vec)
  Trend_ini <<- Data_raw[Data_raw[, "Date"] == as.character(ymd(Start_Date) - 1), "log_I"] %>%
    as.numeric()
  End_Date <<- End_Date
}
