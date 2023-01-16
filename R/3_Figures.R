#' @description Function that creates figures 1 to 4
#' @param model for which the figures are created
#' @param modelOutput from which the figures are created
#' @param dataTib with the COVID-19 infections
#' @param dataExten with COVID-19 infections including before the start of the observational period
#' @param start date (Figure 1)
#' @param end data (Figure 1)
#' @export figures saved in the ´model´_Output/Figures folder

Figures_fctn <- function(model, modelOutput, dataTib, dataExten, start = startDate, end = endDate) {
  # Extract necessary data
  outputTib <- modelOutput[["outputTib"]]
  smootherRegimeIndicatorTib <- modelOutput[["smootherRegimeIndicatorTib"]]
  filterRegimeIndicatorTib <- modelOutput[["filterRegimeIndicatorTib"]]
  # Set the storage folder
  path <- paste0(getwd(), "/Output/Output_", model, "/Figures")
  latexPath <- paste0("C:/Users/phaim/OneDrive/COVID Paper/Latex/Plots")
  if (!dir.exists(path)) {
    dir.create(path)
  }
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

  #---------------------------------------------------------------------------------------#
  # Figure 1                                                                           ####
  #---------------------------------------------------------------------------------------#

  colors_1 <- c("i" = orange, "log i" = grey)
  # Construct and save the plot
  covidInfections <- dataExten %>%
    filter(Date <= end & !is.na(logI) & !is.infinite(logI)) %>%
    ggplot() +
    geom_line(aes(x = Date, y = logI, colour = "log i"), linewidth = 0.75) +
    geom_line(aes(
      x = Date, y = I / (max(I, na.rm = T) / max(logI, na.rm = T)),
      colour = "i"
    ), linewidth = 0.75) +
    geom_vline(xintercept = as.Date(start), color = blue, linetype = 2, linewidth = 1) +
    scale_y_continuous(
      name = expression(paste(log, "(", i[t], ")")),
      breaks = seq(0, 14, by = 2),
      sec.axis = sec_axis(~ . * (max(dataTib$I, na.rm = T) / max(dataTib$logI, na.rm = T)),
        labels = scales::comma_format(big.mark = ".", decimal.mark = ","),
        name = expression(i[t]), breaks = seq(0, 14e+05, by = 2e+05)
      )
    ) +
    scale_color_manual(values = colors_1) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    themeElement
  ggsave(covidInfections, filename = paste0(path, "/Fig_1_Infections.png"), height = 8, width = 14)
  ggsave(covidInfections, filename = paste0(latexPath, "/Fig_1_Infections.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Figure 2                                                                           ####
  #---------------------------------------------------------------------------------------#

  # Load in the grid search results
  gridSearchResults <- read_rds(file = paste0("Output/Output_", model, "/GridSearch/resultsOptimCnstrndMat.rds"))
  # Create mapping for facet label names
  labelMap <- LabelMap_fctn()
  # Transposes matrix to long format
  gridSearchResultsLong <- gridSearchResults[, -c(1, (NCOL(gridSearchResults) - 2):NCOL(gridSearchResults))] %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Estimate") %>%
    mutate(Parameter = factor(Parameter, levels = names(labelMap)))

  # Construct and save the plot
  nParams <- NCOL(gridSearchResults) - 4
  gridSearch <- ggplot(gridSearchResultsLong, aes(x = Parameter, y = Estimate)) +
    geom_violin() +
    geom_point(
      data = gridSearchResultsLong[1:nParams, ], inherit.aes = F,
      aes(x = Parameter, y = Estimate),
      color = "red", size = 3
    ) +
    facet_wrap(~Parameter,
      scales = "free", nrow = 1,
      labeller = labeller(Parameter = LabelHelper_fctn(labelMap))
    ) +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_y_continuous(label = number_format(accuracy = 0.001)) +
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = textSize),
      axis.title = element_text(size = textSize + 3),
      strip.text = element_text(size = textSize),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks = element_line(linewidth = 1)
    )
  ggsave(gridSearch, filename = paste0(path, "/Fig_2_Grid_search.png"), height = 5, width = 14)
  ggsave(gridSearch, filename = paste0(latexPath, "/Fig_2_Grid_search.png"), height = 5, width = 14)

  #---------------------------------------------------------------------------------------#
  # Figure 3                                                                           ####
  #---------------------------------------------------------------------------------------#

  colors_3 <- c("Data" = grey, "Trend" = orange, "Pr regime 1" = blue, "Pr regime 2" = black)
  # Construct and save the plot
  nRegime <- 1 + NCOL(outputTib[, str_detect(colnames(outputTib), "Smoother_Pr_S")])
  regimeProbs <- outputTib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Data, colour = "Data"), linewidth = 0.7) +
    geom_line(aes(x = Date, y = Trend, colour = "Trend"), linewidth = 1) +
    geom_line(aes(x = Date, y = Smoother_Pr_S_0 * max(Data, na.rm = T), colour = "Pr regime 1"), linewidth = 1) +
    geom_rect(
      data = smootherRegimeIndicatorTib, inherit.aes = F,
      aes(ymin = 0, ymax = max(outputTib$Data, na.rm = T), xmin = Date, xmax = End),
      fill = "black", alpha = 0.3
    ) +
    scale_y_continuous(
      name = expression(paste(log, "(", i[t], ")")),
      breaks = seq(0, 14, by = 2),
      sec.axis = sec_axis(~ . * (max(outputTib$Smoother_Pr_S_0, na.rm = T) / max(outputTib$Data, na.rm = T)),
        name = expression(paste(Pr, " (", S[t] == 0, ")")),
        breaks = seq(0, 1, by = 0.2)
      ),
      limits = c(0, max(outputTib$Data, na.rm = T))
    ) +
    scale_color_manual(values = colors_3) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    themeElement
  if (nRegime == 3) {
    regimeProbs <- regimeProbs +
      geom_line(aes(x = Date, y = Smoother_Pr_S_2 * max(Data, na.rm = T), colour = "Pr regime 2"), linewidth = 1)
  }
  ggsave(regimeProbs, filename = paste0(path, "/Fig_3_Regimes.png"), height = 8, width = 14)
  ggsave(regimeProbs, filename = paste0(latexPath, "/Fig_3_Regimes.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Figure 4                                                                           ####
  #---------------------------------------------------------------------------------------#

  colors_4 <- c(Smoother = blue, Filter = orange)
  # Construct and save the plot
  nowcasting <- outputTib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Smoother_Pr_S_0, color = "Smoother"), linewidth = 1) +
    geom_line(aes(x = Date, y = Filter_Pr_S_0, color = "Filter"), linewidth = 1) +
    geom_rect(
      data = filterRegimeIndicatorTib, inherit.aes = F,
      aes(ymin = 0, ymax = 1, xmin = Date, xmax = End),
      fill = "black", alpha = 0.3
    ) +
    geom_rect(
      data = smootherRegimeIndicatorTib, inherit.aes = F,
      aes(ymin = -.03, ymax = 0, xmin = Date, xmax = End),
      fill = "black", alpha = 1
    ) +
    scale_color_manual(values = colors_4) +
    scale_y_continuous(
      name = expression(paste(Pr, " (", S[t] == 0, ")")),
      breaks = seq(0, 1, 0.25)
    ) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    themeElement
  ggsave(nowcasting, filename = paste0(path, "/Fig_4_Nowcasting.png"), height = 8, width = 14)
  ggsave(nowcasting, filename = paste0(latexPath, "/Fig_4_Nowcasting.png"), height = 8, width = 14)
}


#' @description Helper function that applies a label map to facet names
#' @param mapping of class expression

LabelHelper_fctn <- function(mapping) {
  as_labeller(
    function(x) {
      as.list(mapping[x])
    },
    default = identity
  )
}
