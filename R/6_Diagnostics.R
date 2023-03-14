#' @description Function that runs further model diagnostics
#' @param model for which the figures are created
#' @param dataTib with the COVID-19 infections
#' @param modelOutput from which the figures are created
#' @export Plots saved in the´model´_Diagnostics folder

Diagnostics_fctn <- function(model, dataTib, modelOutput) {
  # Extract necessary data
  outputTib <- modelOutput[["outputTib"]]
  filterOutput <- modelOutput[["filterOutput"]]
  smootherOutput <- modelOutput[["smootherOutput"]]
  # Set the storage folder
  path <- paste0(getwd(), "/Output/Output_", model, "/Diagnostics/")  
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Set figure aesthetics
  textSize <- 24
  blue <- "#1874CD"
  themeElementFacet <- theme(
    legend.position = "none",
    legend.text = element_text(size = textSize),
    axis.text = element_text(size = textSize),
    strip.text = element_text(size = textSize),
    axis.title = element_text(size = textSize + 3),
    plot.title = element_text(size = textSize - 3),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(linewidth = 1)
  )
  themeElement <- theme(
    legend.position = "none",
    legend.text = element_text(size = textSize),
    axis.text = element_text(size = textSize),
    axis.title = element_text(size = textSize + 3),
    plot.title = element_text(size = textSize - 3),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(linewidth = 1)
  )

  #---------------------------------------------------------------------------------------#
  # Trend divergence                                                                   ####
  #---------------------------------------------------------------------------------------#

  # Construct and save the plot displaying differences between smoothed trend and log COVID-19 infections
  diffPlot <- outputTib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Diff)) +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    scale_y_continuous(name = expression(paste(hat(mu[t]), " - ", log(i[t])))) +
    labs(title = "Difference between smoothed trend and log COVID-19 infections") +
    theme_bw() +
    themeElement
  ggsave(diffPlot, filename = paste0(path, "Diff_Trend.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Prediction errors                                                                  ####
  #---------------------------------------------------------------------------------------#

  # Build a tibble in long format with the one-step prediction errors classified by the respective regime switch
  predErrTib <- BuildPredErrTib_fctn(filterOutput, dataTib)
  # Construct and save the plot
  predErrPlot <- predErrTib %>%
    ggplot() +
    geom_point(aes(x = Date, y = Error), size = 0.9) +
    geom_line(aes(x = Date, y = (Mean + 2 * Stddev)), linetype = "dashed", color = blue, linewidth = 1) +
    geom_line(aes(x = Date, y = (Mean - 2 * Stddev)), linetype = "dashed", color = blue, linewidth = 1) +
    facet_wrap(~Switch,
      labeller = label_bquote(paste(
        S[t - 1], " = ", .(substr(Switch, 1, 1)), ", ",
        S[t], " = ", .(substr(Switch, 6, 6))
      ))
    ) +
    theme_bw() +
    labs(title = "One-step-ahead prediction errors") +
    scale_y_continuous(
      name = expression(v[t]),
      breaks = seq(-1.5, 1.5, 1.5)
    ) +
    scale_x_date(date_breaks = "9 months", date_labels = c("%m.%Y"), name = "") +
    themeElementFacet
  ggsave(predErrPlot, filename = paste0(path, "Pred_err.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # ACF of Prediction errors                                                           ####
  #---------------------------------------------------------------------------------------#

  # Define a confidence interval for the ACF-Plot
  ciInterval <- 0.95
  # Build a tibble in the long format with the ACF of the one-step-ahead prediction errors
  acfTib <- BuildAcfTib_fctn(predErrTib, ciInterval)
  # Construct and save the plot
  acfPlot <- acfTib %>%
    ggplot() +
    geom_bar(aes(x = Lags, y = AC), stat = "identity", width = 0.5) +
    geom_line(aes(x = Lags, y = CI), linetype = "dashed", color = "#1874CD", linewidth = 1) +
    geom_line(aes(x = Lags, y = -CI), linetype = "dashed", color = "#1874CD", linewidth = 1) +
    geom_hline(yintercept = 0) +
    facet_wrap(~Switch,
      labeller = label_bquote(paste(
        S[t - 1], " = ", .(substr(Switch, 1, 1)), ", ",
        S[t], " = ", .(substr(Switch, 6, 6))
      ))
    ) +
    theme_bw() +
    scale_y_continuous(breaks = seq(-0.1, 0.1, 0.1)) +
    labs(
      y = expression(paste("ACF ", v[t])),
      title = "ACF of one-step-ahead prediction error"
    ) +
    scale_x_continuous(breaks = seq(0, 40, by = 5)) +
    themeElementFacet
  ggsave(acfPlot, filename = paste0(path, "Pred_err_ACF.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Seasonality                                                                        ####
  #---------------------------------------------------------------------------------------#
  # Build a tibble with the smoothed seasonal component
  seasonTib <- tibble(
    Season = smootherOutput$a_cT[3, 3, ],
    Date = dataTib$Date
  ) %>%
    slice(-NROW(dataTib))
  # Construct and save the plot
  seasonPlot <- seasonTib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Season)) +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    scale_y_continuous(expression(paste(hat(gamma[t])))) +
    labs(title = "Smoothed seasonal component") +
    theme_bw() +
    themeElement
  ggsave(seasonPlot, filename = paste0(path, "Seasonal_comp.png"), height = 8, width = 14)
}


#' @description Function constructs a tibble used to sketch the filtered prediction errors
#' @param filterOutput of the Kim filter
#' @param dataTib with the COVID-19 infections
#' @return Tibble with the one-step-ahead prediction errors and their first two moments for every regime combination

BuildPredErrTib_fctn <- function(filterOutput, dataTib) {
  k <- 0
  for (i in (0:1)) {
    for (j in (0:1)) {
      k <- k + 1
      Pred_err_tib_temp <- tibble(
        Date = dataTib$Date,
        Error = filterOutput$Pred_err[, k],
        Switch = paste(j, "to", i),
        Stddev = sd(Error),
        Mean = mean(Error)
      ) %>%
        slice(-c(1:5))
      if (!exists("predErrTib")) {
        predErrTib <- Pred_err_tib_temp
      } else {
        predErrTib <- bind_rows(predErrTib, Pred_err_tib_temp) %>%
          arrange(Date, Switch)
      }
    }
  }
  return(predErrTib)
}


#' @description Function constructs a tibble used to sketch the ACF of the filtered prediction errors
#' @param predErrTib output from the BuildPredErrTib_fctn()
#' @param ciInterval preferred confidence interval
#' @return Tibble with the ACF and CI for every regime combination

BuildAcfTib_fctn <- function(predErrTib, ciInterval) {
  for (i in (0:1)) {
    for (j in (0:1)) {
      Switch <- paste(i, "to", j)
      ACF_tib_temp <- tibble(
        AC = acf(predErrTib[predErrTib[, "Switch"] == Switch, "Error"], plot = FALSE)$acf[-1, 1, 1],
        CI = qnorm((1 + ciInterval) / 2) / sqrt(acf(predErrTib[predErrTib[, "Switch"] == Switch, ]$Error,
          plot = FALSE
        )$n.used),
        Switch = predErrTib[predErrTib[, "Switch"] == Switch, ]$Switch[1]
      ) %>%
        mutate(Lags = 1:n())
      if (!exists("acfTib")) {
        acfTib <- ACF_tib_temp
      } else {
        acfTib <- bind_rows(acfTib, ACF_tib_temp) %>%
          arrange(Lags, Switch)
      }
    }
  }
  return(acfTib)
}
