### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Load data and packages ##################################################################################################################### ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Init_fctn <- function(Start_Date, End_Date, Country, Province = NULL) {

  # Install and or load packages
  if (!require(tidyverse)) install.packages("tidyverse", quiet = T)
  if (!require(readr)) install.packages("readr", quiet = T)
  if (!require(lubridate)) install.packages("lubridate", quiet = T)
  if (!require(scales)) install.packages("scales", quiet = T)
  if (!require(tictoc)) install.packages("tictoc", quiet = T)
  if (!require(Rcpp)) install.packages("Rcpp", quiet = T)
  if (!require(RcppArmadillo)) install.packages("RcppArmadillo", quiet = T)
  if (!require(openxlsx)) install.packages("openxlsx", quiet = T)
  if (!require(coda)) install.packages("coda", quiet = T)
  if (!require(Matrix)) install.packages("Matrix", quiet = T)

  JHCSSE_data <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv",
    col_types = cols()
  )
  # Wrangle data
  if (!(Country %in% JHCSSE_data$`Country/Region`)) stop("Unknown country name")

  if (nrow(JHCSSE_data[JHCSSE_data[, "Country/Region"] == Country, ]) > 1 & is.null(Province)) {
    Province <- deframe(JHCSSE_data[JHCSSE_data[, "Country/Region"] == Country, "Province/State"])
    cat(
      "Individual analysis of the following provinces is available:\n",
      paste(Province, collapse = "; "),
      "\nTo select one or multiple provinces, submit a character vector to the (Province = ) argument\n"
    )
  } else if (nrow(JHCSSE_data[JHCSSE_data[, "Country/Region"] == Country, ]) == 1 & !is.null(Province)) {
    warning("No division into provinces available")
    Province <- deframe(JHCSSE_data[JHCSSE_data[, "Country/Region"] == Country, "Province/State"])
  } else if (is.null(Province)) {
    Province <- deframe(JHCSSE_data[JHCSSE_data[, "Country/Region"] == Country, "Province/State"])
  }

  Data_raw <- JHCSSE_data %>%
    filter(`Country/Region` == Country & `Province/State` %in% Province) %>%
    select(-c(`Country/Region`, `Province/State`, Lat, Long)) %>%
    summarise_all(sum) %>%
    pivot_longer(cols = 1:ncol(.), values_to = "C", names_to = "Date") %>%
    mutate(
      Date = mdy(Date),
      I = C - lag(C),
      log_I = log(I + 1),
      index = 1:n()
    )

  if (max(Data_raw$Date) < End_Date) {
    warning(paste(
      "End date is latest observation:",
      max(Data_raw$Date)
    ))
  }
  if (min(Data_raw$Date) > Start_Date) {
    warning(paste(
      "Start date is earliest observation:",
      min(Data_raw$Date)
    ))
  }

  Data <<- Data_raw %>%
    filter(Date >= Start_Date & Date <= End_Date)
  Data_vec <<- Data$log_I %>% as.matrix()

  # Number of periods
  Periods <<- nrow(Data_vec)
  Trend_ini <<- Data_raw[Data_raw[, "Date"] == as.character(ymd(Start_Date) - 1), "log_I"] %>%
    as.numeric()
  End_Date <<- End_Date
}
