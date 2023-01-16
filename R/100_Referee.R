#' @description Function that simulates one seasonal unit root
#' @param omegaSd Standard deviation of the innovations to the seasonal unit root
#' @param periods to be simulated
#' @export plot to the specified folder

SeasSim_fctn <- function(omegaSd = 1, periods = 1000) {
  path <- paste0(getwd(), "/Referee_Output/")
  if (!dir.exists(path)) {
    dir.create(path)
  }

  # Simulate the seasonal unit root process
  set.seed(1)
  xVec <- rep(0, periods + 7)
  omega <- rnorm(periods + 7, sd = omegaSd)
  for (i in 8:(periods + 7)) xVec[i] <- xVec[i - 7] + omega[i]
  xVec <- tail(xVec, -7)

  # Simulate the deterministic seasonality (not really necessary)
  # gammaDeterm <- c(.0005, 0.0877, -0.1027, -0.1182, 0.0257, 0.1448, -0.0378)
  # gammaDetermVec <- rep_len(gammaDeterm, periods)

  simTib <- tibble(
    Determ = gammaDetermVec,
    Stochas = xVec,
    Gamma = gammaDetermVec + xVec,
    Periods = 1:periods
  )

  # Define aesthetics
  grey <- "#525252"
  textSize <- 24
  # Construct and save the figure
  seasSimulation <- simTib %>%
    ggplot() +
    geom_line(aes(x = Periods, y = Gamma), linewidth = .5, color = grey) +
    theme_bw() +
    theme(
      plot.title = element_text(size = textSize),
      axis.text = element_text(size = textSize),
      axis.title = element_text(size = textSize + 3, vjust = 1),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks = element_line(linewidth = 1)
    ) +
    scale_y_continuous(name = expression(x[t])) +
    scale_x_continuous(name = "") +
    ggtitle(expression(paste(
      "Simulation of the seasonal unit root ", x[t], " with Std.N innovations ",
      omega[t], " (eq. 2)"
    )))
  ggsave(seasSimulation, filename = paste0(path, "SeasSim.png"), height = 8, width = 14)
}
