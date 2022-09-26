### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### Reproduce plots ############################################################################################################################# ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

Plots_fctn <- function(Data_tib, Model_output, Model, Start = Start_Date, End = End_Date) {
  Output_tib <- Model_output[["Output_tib"]]
  Smoother_Regime_indicator_tib <- Model_output[["Smoother_Regime_indicator_tib"]]
  Filter_Regime_indicator_tib <- Model_output[["Filter_Regime_indicator_tib"]]

  Path <- paste0(getwd(), "/", Model, "_Figures")
  if (!dir.exists(Path)) {
    dir.create(Path)
  }

  # Define aesthetics
  orange <- "#FF6347"
  black <- "#000000"
  grey <- "#525252"
  blue <- "#1874CD"

  Text_size <- 24

  Theme_element <- theme(
    legend.position = "none",
    legend.text = element_text(size = Text_size),
    axis.text = element_text(size = Text_size),
    axis.title = element_text(size = Text_size + 3, vjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(size = 1)
  )

  #---------------------------------------------------------------------------------------#
  # Figure 1                                                                           ####
  #---------------------------------------------------------------------------------------#

  Colors_1 <- c("i" = orange, "log i" = grey)

  Covid_infections <- Data_tib %>%
    filter(Date <= End & !is.na(log_I)) %>%
    ggplot() +
    geom_line(aes(x = Date, y = log_I, colour = "log i"), size = 0.75) +
    geom_line(aes(
      x = Date, y = I / (max(I, na.rm = T) / max(log_I, na.rm = T)),
      colour = "i"
    ), size = 0.75) +
    geom_vline(xintercept = as.Date(Start), color = blue, linetype = 2, size = 1) +
    scale_y_continuous(
      name = expression(paste(log, "(", i[t] + 1, ")")),
      breaks = seq(0, 14, by = 2),
      sec.axis = sec_axis(~ . * (max(Data_tib$I, na.rm = T) / max(Data_tib$log_I, na.rm = T)),
        labels = scales::comma_format(big.mark = ".", decimal.mark = ","),
        name = expression(i[t]), breaks = seq(0, 14e+05, by = 2e+05)
      )
    ) +
    scale_color_manual(values = Colors_1) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    Theme_element
  ggsave(Covid_infections, filename = paste0(Path, "/Fig_1_Infections.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Figure 2                                                                           ####
  #---------------------------------------------------------------------------------------#

  Grid_search_Results <- read_rds(file = paste0(Model, "_Grid_search/RESULTS_optim_cnstrnd_mat.rds"))

  # Transpond matrix to long format
  Grid_search_Result_long <- Grid_search_Results[, -1] %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Estimate") %>%
    mutate(Parameter = factor(Parameter, levels = c("xi", "omega", "epsilon", "nu_1", "nu_0", "p", "q")))

  # Create mapping for facet label names
  Label_map <- expression(
    xi = sigma[xi],
    omega = sigma[omega],
    epsilon = sigma[epsilon],
    nu_0 = nu[0],
    nu_1 = nu[1],
    p = italic(p),
    q = italic(q)
  )

  Label_helper <- function(Mapping) {
    as_labeller(function(x) {
      as.list(Mapping[x])
    },
    default = identity
    )
  }

  Grid_search <- ggplot(Grid_search_Result_long, aes(x = Parameter, y = Estimate)) +
    geom_violin() +
    geom_point(
      data = Grid_search_Result_long[1:7, ], inherit.aes = F,
      aes(x = Parameter, y = Estimate),
      color = "red", size = 3
    ) +
    facet_wrap(~Parameter,
      scales = "free", nrow = 1,
      labeller = labeller(Parameter = Label_helper(Label_map))
    ) +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_y_continuous(label = number_format(accuracy = 0.001)) +
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = Text_size),
      axis.title = element_text(size = Text_size + 3),
      strip.text = element_text(size = Text_size),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks = element_line(size = 1)
    )
  ggsave(Grid_search, filename = paste0(Path, "/Fig_2_Grid_search.png"), height = 5, width = 14)

  #---------------------------------------------------------------------------------------#
  # Figure 3                                                                           ####
  #---------------------------------------------------------------------------------------#

  Colors_3 <- c("Data" = grey, "Trend" = orange, "Pr regime 1" = blue)

  Regime_probs <- Output_tib %>%
    ggplot() +
    geom_line(aes(x = Date, y = Data, colour = "Data"), size = 0.7) +
    geom_line(aes(x = Date, y = Trend, colour = "Trend"), size = 1) +
    geom_line(aes(x = Date, y = Smoother_Pr_S_1 * max(Data, na.rm = T), colour = "Pr regime 1"), size = 1) +
    geom_rect(
      data = Smoother_Regime_indicator_tib, inherit.aes = F,
      aes(ymin = 0, ymax = max(Output_tib$Data, na.rm = T), xmin = Date, xmax = End),
      fill = "black", alpha = 0.3
    ) +
    scale_y_continuous(
      name = expression(paste(log, "(", i[t] + 1, ")")),
      breaks = seq(0, 14, by = 2),
      sec.axis = sec_axis(~ . * (max(Output_tib$Smoother_Pr_S_1, na.rm = T) / max(Output_tib$Data, na.rm = T)),
        name = expression(paste(Pr, " (", S[t] == 1, ")")),
        breaks = seq(0, 1, by = 0.2)
      ),
      limits = c(0, max(Output_tib$Data, na.rm = T))
    ) +
    scale_color_manual(values = Colors_3) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    Theme_element
  ggsave(Regime_probs, filename = paste0(Path, "/Fig_3_Regimes.png"), height = 8, width = 14)

  #---------------------------------------------------------------------------------------#
  # Figure 4                                                                           ####
  #---------------------------------------------------------------------------------------#

  Colors_4 <- c(Smoother = blue, Filter = orange)

  Nowcasting <- Output_tib %>%
    ggplot() +
    geom_line(aes(x = Date, y = 1 - Smoother_Pr_S_1, color = "Smoother"), size = 1) +
    geom_line(aes(x = Date, y = 1 - Filter_Pr_S_1, color = "Filter"), size = 1) +
    geom_rect(
      data = Filter_Regime_indicator_tib, inherit.aes = F,
      aes(ymin = 0, ymax = 1, xmin = Date, xmax = End),
      fill = "black", alpha = 0.3
    ) +
    scale_color_manual(values = Colors_4) +
    scale_y_continuous(name = expression(paste(Pr, " (", S[t] == 0, ")"))) +
    theme_bw() +
    labs(colour = "") +
    scale_x_date(date_breaks = "5 months", date_labels = c("%m.%Y"), name = "") +
    Theme_element
  ggsave(Nowcasting, filename = paste0(Path, "/Fig_4_Nowcasting.png"), height = 8, width = 14)
}
