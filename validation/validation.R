# Clear workspace and set working directory
rm(list = ls()) 
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
source(file = "bTBwl_func.R") 

# Simulation parameters
years <- 3
times <- seq(from = 0, to = years * 12, by = 12 / 365) 
seedQuarter <- 1 
prop_superSpreader <- 0.05 
pct <- 0.02 
verbose <- 0 
reps <- 50 
sizes <- c(10, 50, 100, 250, 500) 
infType <- "seeded" 
runtype <- "test_" 
name_out <- paste0(runtype, infType, "_", pct, '-', prop_superSpreader, '-') 

# Paths and file saving
pth <- file.path(getwd()) 
save_runs <- TRUE 
save_plots <- TRUE 

# Testing flags
test_mode <- TRUE
test_birth <- TRUE
test_death_n <- TRUE
test_death_h <- TRUE
test_disease <- FALSE

# Parallelization setup
n.cores <- floor(detectCores() * (3 / 4)) 
cl <- makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { 
  library(deSolve) 
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(rstudioapi)
  library(foreach)
  library(doParallel)
  library(RColorBrewer)
  library(scales)
})

# Main simulation loop
foreach(size = iter(sizes)) %dopar% {
  
  # Initialize parameters
  size <- as.integer(size)
  pars <- parameter_set_wl(
    k = size, 
    scenario = infType, 
    initial_exposed = pct * size, 
    SS_prop = prop_superSpreader,
    start_quarter = seedQuarter, 
    test = test_mode,
    birth = test_birth,
    death_h = test_death_h,
    death_n = test_death_n,
    disease = test_disease,
    verbose = verbose
  )
  
  # Deterministic model
  X0_full <- c(
    S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
    sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"])
  )
  tic <- Sys.time()
  ode(
    func = SEI_model_full,
    y = X0_full,
    times = times,
    parms = pars,
    method = "rk4"
  ) %>%
    as.data.frame() -> out
  toc <- Sys.time()
  print(toc - tic)
  
  out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
  data <- out %>% gather(variable, value, -time)
  data$value[is.nan(data$value)] <- 0 
  data$value[data$value < 0] <- 0
  
  # Stochastic model
  initial_state <- data.frame(
    S_0 = pars["S_0"], E1_0 = pars["E1_0"], I_0 = pars["I_0"], 
    SuperS_0 = pars["SuperS_0"], SuperE1_0 = pars["SuperE1_0"], SuperI_0 = pars["SuperI_0"]
  )
  population_parameters <- data.frame(
    K = pars["K"], eta_hunt = pars["eta_hunt"], eta_nat = pars["eta_nat"], 
    theta = pars["theta"], gamma = pars["gamma"], alpha_max = pars["alpha_max"], 
    ksi = pars["ksi"], omega = pars["omega"], s = pars["s"], alpha = pars["alpha"]
  )
  disease_parameters <- data.frame(
    beta = pars["beta"], area = pars["area"], p1 = pars["p1"], 
    p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], 
    p2_q4 = pars["p2_q4"], phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], 
    sigma1_rate = pars["sigma1_rate"]
  )
  param_vals <- merge(population_parameters, disease_parameters)
  
  tic <- Sys.time()
  sto_out <- wildlife_model(
    n_reps = reps, parameters = param_vals, initial_state = initial_state, 
    nyears = years, seed_quarter = as.integer(pars["start_q"]), 
    verbose = verbose, batch_name = "test", type = 'c'
  )
  toc <- Sys.time()
  print(toc - tic)
  
  # Format data
  data_sto <- sto_out[, c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>%
    pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  models_res <- list(
    as.data.frame(data), 
    as.data.frame(sto_out), 
    as.data.frame(data_sto)
  )
  
  # Plot results
  sto_by_class <- split(as.data.frame(models_res[[3]]), models_res[[3]]$variable)
  det_by_class <- split(as.data.frame(models_res[[1]]), models_res[[1]]$variable)
  
  SEIcols <- RColorBrewer::brewer.pal(11, "Spectral")[c(1, 2, 4, 5, 8, 9, 10)]
  
  # Define and save plots
  Nplot <- ggplot() + 
    geom_line(data = sto_by_class[["N"]], aes(x = time, y = value, color = variable, group = rep), size = 1, alpha = .9) + 
    geom_line(data = det_by_class[["N"]], aes(x = time, y = value), size = 1) +
    scale_color_manual(values = SEIcols[7]) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(breaks = seq(0, years * 12, 12), limits = c(0, years * 12))
  
  Splot <- ggplot() + 
    geom_line(data = sto_by_class[["S"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
    geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[6]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  Eplot <- ggplot() + 
    geom_line(data = sto_by_class[["E1"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
    geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[4]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  Iplot <- ggplot() + 
    geom_line(data = sto_by_class[["I"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
    geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[2]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  SSSplot <- ggplot() + 
    geom_line(data = sto_by_class[["SS_S"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
    geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[5]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  ESSplot <- ggplot() + 
    geom_line(data = sto_by_class[["SS_E1"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
    geom_line(data = det_by_class[["sE"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[3]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  ISSplot <- ggplot() + 
    geom_line(data = sto_by_class[["SS_I"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
    geom_line(data = det_by_class[["sI"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[1]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
  print(Nplot)
  
  if (save_plots) {
    jpeg(filename = file.path(pth, "validation/figures", paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg")))
    print(Nplot)
    dev.off()
    
    jpeg(filename = file.path(pth, "validation/figures", paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg")))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
  }
  
  # Save outputs
  if (save_runs) {
    save(models_res, file = file.path(pth, "test_runs", paste0(name_out, size, '-', Sys.Date(), '.RData')))
  }
}

stopCluster(cl)