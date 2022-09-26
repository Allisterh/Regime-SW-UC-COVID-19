### ---------------------------------------------------------------------------------------------------------------------------------------------###
### ---------------------------------------------------------------------------------------------------------------------------------------------###
### Diagnostics ###################################################################################################################################
### ---------------------------------------------------------------------------------------------------------------------------------------------###
### ---------------------------------------------------------------------------------------------------------------------------------------------###

Diagnostics_fctn <- function(Data_tib, Model_output, Model) {
  Output_tib <- Model_output[["Output_tib"]]
  Filter_output <- Model_output[["Filter_output"]]
  Smoother_output <- Model_output[["Smoother_output"]]

  Text_size <- 24
  blue <- "#1874CD"

  Path <- paste0("/", Model, "Diagnostics/")
  if (!dir.exists(Path)) {
    dir.create(Path)
  }

  Theme_element_facet <- theme(
    legend.position = "none",
    legend.text = element_text(size = Text_size),
    axis.text = element_text(size = Text_size),
    strip.text = element_text(size = Text_size),
    axis.title = element_text(size = Text_size + 3),
    plot.title = element_text(size = Text_size - 3),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(size = 1)
  )

  Theme_element <- theme(
    legend.position = "none",
    legend.text = element_text(size = Text_size),
    axis.text = element_text(size = Text_size),
    axis.title = element_text(size = Text_size + 3),
    plot.title = element_text(size = Text_size - 3),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(size = 1)
  )

  #---------------------------------------------------------------------------------------#
  # Trend divergence                                                                   ####
  #---------------------------------------------------------------------------------------#

  # Difference between smoothed trend and log COVID-19 infections
  Diff_plot <- Output_tib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Diff)) +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    scale_y_continuous(name = expression(paste(hat(mu[t]), " - ", log(i[t] + 1)))) +
    labs(title = "Difference between smoothed trend and log COVID-19 infections") +
    theme_bw() +
    Theme_element
  ggsave(Diff_plot, filename = paste0(Path, "Diff_Trend.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Prediction errors                                                                  ####
  #---------------------------------------------------------------------------------------#

  # Build a tibble in long format with the one-step prediction errors classified by the respective regime switch
  Pred_err_tib <- Build_Pred_Err_tib(Filter_output, Data_tib)

  Pred_err_plot <- Pred_err_tib %>%
    ggplot() +
    geom_point(aes(x = Date, y = Error), size = 0.9) +
    geom_line(aes(x = Date, y = (Mean + 2 * Stddev)), linetype = "dashed", color = blue, size = 1) +
    geom_line(aes(x = Date, y = (Mean - 2 * Stddev)), linetype = "dashed", color = blue, size = 1) +
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
    Theme_element_facet
  ggsave(Pred_err_plot, filename = paste0(Path, "Pred_err.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # ACF of Prediction errors                                                           ####
  #---------------------------------------------------------------------------------------#

  # Define a confidence interval for the ACF-Plot
  CI_interval <- 0.95

  # Build a tibble in the long format with the ACF of the one-step-ahead prediction errors
  ACF_tib <- Build_ACF_tib(Pred_err_tib, CI_interval)

  ACF_plot <- ACF_tib %>%
    ggplot() +
    geom_bar(aes(x = Lags, y = AC), stat = "identity", width = 0.5) +
    geom_line(aes(x = Lags, y = CI), linetype = "dashed", color = "#1874CD", size = 1) +
    geom_line(aes(x = Lags, y = -CI), linetype = "dashed", color = "#1874CD", size = 1) +
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
    Theme_element_facet
  ggsave(ACF_plot, filename = paste0(Path, "Pred_err_ACF.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Seasonality                                                                        ####
  #---------------------------------------------------------------------------------------#

  # Build a tibble with the smoothed seasonal component
  Season_tib <- tibble(
    Season = Smoother_output$a_cT[3, 3, ],
    Date = Data_tib$Date
  ) %>%
    slice(-Periods)

  Season_plot <- Season_tib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Season)) +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    scale_y_continuous(expression(paste(hat(gamma[t])))) +
    labs(title = "Smoothed seasonal component") +
    theme_bw() +
    Theme_element
  ggsave(Season_plot, filename = paste0(Path, "Seasonal_comp.png"), height = 8, width = 14)
}

#---------------------------------------------------------------------------------------#
# Supplementary functions                                                            ####
#---------------------------------------------------------------------------------------#

Build_Pred_Err_tib <- function(Filter_output, Data_tib) {
  k <- 0
  for (i in (0:1)) {
    for (j in (0:1)) {
      k <- k + 1
      Pred_err_tib_temp <- tibble(
        Date = Data_tib$Date,
        Error = Filter_output$Pred_err[, k],
        Switch = paste(j, "to", i),
        Stddev = sd(Error),
        Mean = mean(Error)
      ) %>%
        slice(-c(1:5))
      if (!exists("Pred_err_tib")) {
        Pred_err_tib <- Pred_err_tib_temp
      } else {
        Pred_err_tib <- bind_rows(Pred_err_tib, Pred_err_tib_temp) %>%
          arrange(Date, Switch)
      }
    }
  }
  return(Pred_err_tib)
}

Build_ACF_tib <- function(Pred_err_tib, CI_interval) {
  for (i in (0:1)) {
    for (j in (0:1)) {
      Switch <- paste(i, "to", j)
      ACF_tib_temp <- tibble(
        AC = acf(Pred_err_tib[Pred_err_tib[, "Switch"] == Switch, "Error"])$acf[-1, 1, 1],
        CI = qnorm((1 + CI_interval) / 2) / sqrt(acf(Pred_err_tib[Pred_err_tib[, "Switch"] == Switch, ]$Error)$n.used),
        Switch = Pred_err_tib[Pred_err_tib[, "Switch"] == Switch, ]$Switch[1]
      ) %>%
        mutate(Lags = 1:n())
      if (!exists("ACF_tib")) {
        ACF_tib <- ACF_tib_temp
      } else {
        ACF_tib <- bind_rows(ACF_tib, ACF_tib_temp) %>%
          arrange(Lags, Switch)
      }
    }
  }
  return(ACF_tib)
}
