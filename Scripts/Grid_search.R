### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### 3-step grid search ########################################################################################################################## ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------------------------------------------------------------- ###

#---------------------------------------------------------------------------------------#
# Function to create distribution plots                                              ####
#---------------------------------------------------------------------------------------#

Plot_parameter_dist <- function(Paramter_mat, Title, Path, Model) {
  Text_size <- 24

  Paramter_mat_long <- Paramter_mat %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Estimate") %>%
    mutate(Parameter = factor(Parameter, levels = c("xi", "omega", "epsilon", "nu_1", "p", "q")))

  Label_map <- expression(
    epsilon = sigma[epsilon],
    nu_1 = nu[1],
    omega = sigma[omega],
    p = italic(p),
    q = italic(q),
    xi = sigma[xi]
  )
  
  Label_helper <- function(Mapping) {
    as_labeller(function(x) {
      as.list(Mapping[x])
    },
    default = identity
    )
  }

  Plot <- ggplot(Paramter_mat_long, aes(x = Parameter, y = Estimate)) +
    geom_violin() +
    facet_wrap(~Parameter,
      scales = "free",
      labeller = labeller(Parameter = Label_helper(Label_map))
    ) +
    theme_bw() +
    scale_x_discrete(name = "") +
    scale_y_continuous(label = number_format(accuracy = 0.001)) +
    labs(title = paste("Best 50 combinations |", Title, "|", Model)) +
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
  ggsave(Plot,
    filename = paste0(Path, "/", Title, ".png"), device = "png",
    width = 14, height = 8
  )
}

#---------------------------------------------------------------------------------------#
# Grid search function                                                               ####
#---------------------------------------------------------------------------------------#

Grid_search_fctn <- function(Model, n_random, steps_grid) {
  source(paste0("Scripts/", Model, ".R"))
  
  Path <- paste0(getwd(), "/" , Model, "_Grid_search")
    if (!dir.exists(Path)) {
    dir.create(Path)
  }

  #-------------------------------------------------------------------------------------#
  # Initial random grid search                                                       ####
  #-------------------------------------------------------------------------------------#

  set.seed(2)
  # Build random grid
  THETA_rand <- tibble(
    xi = log(runif(n_random, 0.01, 0.25)),
    nu_1 = log(runif(n_random, 0.001, 0.3)),
    omega = log(runif(n_random, 0, 0.3)),
    q = runif(n_random, -5, 5),
    p = runif(n_random, -5, 5),
    epsilon = log(runif(n_random, 0.01, 0.3))
  )

  cat("Step 1/3: Begin initial random grid search \n")
  tic("          Initial random grid search")
  RESULTS_rand_unfltrd <- apply(THETA_rand, 1, function(x) {
    return(c(
      tryCatch(Likelihood_fctn(x), error = function(e) NA),
      x
    ))
  }) %>%
    t()
  toc()

  # Collect results
  RESULTS_rand_unfltrd <- RESULTS_rand_unfltrd[order(RESULTS_rand_unfltrd[, 1]), ]
  RESULTS_rand <- RESULTS_rand_unfltrd[!is.na(RESULTS_rand_unfltrd[, 1]), ]
  write_rds(RESULTS_rand, file = paste0(Path, "/RESULTS_rand.rds"))
  
  # Analyse best thetas
  RESULTS_rand_cnstrnd <- t(apply(RESULTS_rand[, -1], 1, Par_constrain_fctn))
  colnames(RESULTS_rand_cnstrnd) <- names(THETA_rand)
  Plot_parameter_dist(RESULTS_rand_cnstrnd[1:50,], Title = "Step1_RandGrid", Path = Path, 
                      Model = Model)

  #-------------------------------------------------------------------------------------#
  # Second grid search                                                               ####
  #-------------------------------------------------------------------------------------#

  # Build grid with the most extreme realizations of the 50 best thetas
  Grid_U_limits <- apply(RESULTS_rand[1:50, -1], MARGIN = 2, max)
  Grid_L_limits <- apply(RESULTS_rand[1:50, -1], MARGIN = 2, min)

  THETA_grid <- tibble(
    xi = seq(Grid_L_limits[1], Grid_U_limits[1], length.out = steps_grid),
    nu_1 = seq(Grid_L_limits[2], Grid_U_limits[2], length.out = steps_grid),
    omega = seq(Grid_L_limits[3], Grid_U_limits[3], length.out = steps_grid),
    q = seq(Grid_L_limits[4], Grid_U_limits[4], length.out = steps_grid),
    p = seq(Grid_L_limits[5], Grid_U_limits[5], length.out = steps_grid),
    epsilon = seq(Grid_L_limits[6], Grid_U_limits[6], length.out = steps_grid)
  ) %>%
    expand.grid()

  # Search over parameter grid
  cat("Step 2/3: Begin second grid search \n")
  tic("          Second grid search")
  RESULTS_grid_unfltrd <- apply(THETA_grid, 1, function(x) {
    return(c(
      tryCatch(Likelihood_fctn(x), error = function(e) NA),
      x
    ))
  }) %>%
    t()
  toc()

  # Collect results
  RESULTS_grid_unfltrd <- RESULTS_grid_unfltrd[order(RESULTS_grid_unfltrd[, 1]), ]
  RESULTS_grid <- RESULTS_grid_unfltrd[!is.na(RESULTS_grid_unfltrd[, 1]), ]
  write_rds(RESULTS_grid, file = paste0(Path, "/RESULTS_grid.rds"))
  
  # Analyse best thetas
  RESULTS_grid_cnstrnd <- t(apply(RESULTS_grid[, -1], 1, Par_constrain_fctn))
  colnames(RESULTS_grid_cnstrnd) <- names(THETA_rand)
  Plot_parameter_dist(RESULTS_grid_cnstrnd[1:50,], Title = "Step2_SecondGrid", 
                      Path = Path, Model = Model)

  #-------------------------------------------------------------------------------------#
  # Final grid search                                                                 ####
  #-------------------------------------------------------------------------------------#

  Theta_final <- RESULTS_grid[1:50, -1]

  # Optimize over 50 best thetas from the second step
  cat("Step 3/3: Begin final optimization \n")
  tic("          Final grid search")
  RESULTS_optim <- apply(Theta_final, 1, function(x) {
    return(c(
      tryCatch(optim(x,
        fn = Likelihood_fctn, hessian = TRUE, method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-06)
      ), error = function(e) NA),
      x
    ))
  }) %>%
    t()
  toc()

  # Collect parameters from the list of optim results
  RESULTS_optim_unfltrd <- matrix(NA, nrow = length(RESULTS_optim), ncol = 7)
  for (i in (1:length(RESULTS_optim))) {
    if (is.list(RESULTS_optim[[i]])) {
      RESULTS_optim_unfltrd[i, 1] <- RESULTS_optim[[i]]$value
      RESULTS_optim_unfltrd[i, 2:7] <- RESULTS_optim[[i]]$par
    } else {
      RESULTS_optim_unfltrd[i, ] <- NA
    }
  }
  colnames(RESULTS_optim_unfltrd) <- c("ll", names(THETA_rand))
  RESULTS_optim_unfltrd <- RESULTS_optim_unfltrd[order(RESULTS_optim_unfltrd[, "ll"]), ]
  RESULTS_optim_mat <- RESULTS_optim_unfltrd[!is.na(RESULTS_optim_unfltrd[, 1]), ]
  write_rds(RESULTS_optim_mat, file = paste0(Path, "/RESULTS_optim_mat.rds"))

  # Best parameters
  THETA_vec <- RESULTS_optim_mat[1, -1]
  write_rds(THETA_vec, file = paste0(Path, "/THETA_vec.rds"))

  # Compute nu_0 estimates
  nu_0_mat <- matrix(ncol = 2, nrow = nrow(RESULTS_optim_mat))
  for (i in 1:nrow(RESULTS_optim_mat)) {
    MLE_param <- Par_constrain_fctn(RESULTS_optim_mat[i, -1])
    Filter_output <- Filter_fctn(par = MLE_param, data = Data_vec)
    Smoother_output <- Smoother_fctn(par = MLE_param, Filter_output, level = 0.9)
    nu_0_mat[i, ] <- c(Smoother_output[["a_cT"]][2, 1, Periods], sqrt(Smoother_output[["P_cT"]][2, 2, Periods]))
  }

  # Apply constraining function to the matrix of parameter estimates
  RESULTS_optim_cnstrnd_mat <- cbind(
    -RESULTS_optim_mat[, 1], apply(
      RESULTS_optim_mat[, -1],
      1, Par_constrain_fctn
    ) %>% t(),
    nu_0_mat[, 1]
  )
  colnames(RESULTS_optim_cnstrnd_mat) <- c("ll", names(THETA_rand), "nu_0")
  write_rds(RESULTS_optim_cnstrnd_mat, file = paste0(Path, "/RESULTS_optim_cnstrnd_mat.rds"))
  # Save to excel
  write.xlsx(as_tibble(RESULTS_optim_cnstrnd_mat), file = paste0(Path, "/RESULTS_optim_cnstrnd_mat.xlsx"))
  
  # Save Hessian matrix of best theta to compute standard errors in table 1 (append estimate for nu_0 derived from the state vector)
  for (index in 1:length(RESULTS_optim)) {
    if (RESULTS_optim[[index]]$value == RESULTS_optim_mat[1, 1]) {
      break
    }
  }
  SE_list <- list(Hessian_mat = RESULTS_optim[[index]]$hessian, nu_0_var = nu_0_mat[1, 2])
  write_rds(SE_list, file = paste0(Path, "/SE_list.rds"))
}
