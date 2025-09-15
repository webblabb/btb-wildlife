system("g++ -L/usr/lib/x86_64-linux-gnu src/bTB_wildlifeModel_DTMC.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu src/bTB_wildlifeModel_CTMC.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")

####################
## Initialization ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "R/bTBwl_func.R")

years = 10
times = seq(from=0, to=years*12, by=12/365)
seedQuarter = 1
prop_superSpreader = 0.05
pct = 0.02
verbose = 0
reps = 500
sizes = c(10, 50, 100, 250, 500)
infType = "seeded"
runtype = "stdPres_"
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-')
pth = getwd()
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
save_runs = T
save_plots = T

for (size in sizes) {
  models_res <- vector(mode = "list", length = 3)
  
  pars <- parameter_set_wl(k = size, 
                           scenario = infType, 
                           initial_exposed = pct*size, 
                           SS_prop = prop_superSpreader,
                           start_quarter = seedQuarter, 
                           test = test_mode,
                           birth = test_birth,
                           death_h = test_death_h,
                           death_n = test_death_n,
                           disease = test_disease,
                           verbose = verbose)
  
  # deterministic model
  X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
              sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
  
  tic = Sys.time()
  out <- ode(
    func = SEI_model_full,
    y = X0_full,
    times = times,
    parms = pars,
    method = "rk4"
  ) %>% as.data.frame()
  toc = Sys.time()
  print(toc - tic)
  
  out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
  data <- out %>% tidyr::gather(variable, value, -time)
  remove(out)
  data$value[is.nan(data$value)] <- 0 
  data$value[which(data$value < 0)] <- 0
  
  # stochastic model
  initial_state <- data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters <- data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"])
  disease_parameters <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  param_vals <- data.frame(merge(population_parameters, disease_parameters))
  
  tic = Sys.time()
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = param_vals,
    initial_state = initial_state,
    nyears = years,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = verbose,
    batch_name = "test",
    type = 'c'
  )
  toc = Sys.time()
  print(toc - tic)
  remove(initial_state, population_parameters, disease_parameters, param_vals)
  
  # Only keep the columns you need for long-format and plotting
  sto_long <- sto_out |>
    dplyr::select(rep, tstep, time, N, S, SS_S, E1, SS_E1, I, SS_I) |>
    tidyr::pivot_longer(cols = c(N, S, SS_S, E1, SS_E1, I, SS_I),
                        names_to = "variable", values_to = "value")
  
  # Remove the big stochastic output as soon as possible
  rm(sto_out)
  gc()
  
  # Split for plotting, but only keep the variables you need
  sto_by_class <- split(sto_long, sto_long$variable)
  det_by_class <- split(data, data$variable)
  
  SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
  
  # You can define a quick plotting function to avoid repetition
  plot_class <- function(var, sto_color, det_var=NULL, det_data=NULL, ylab=NULL) {
    p <- ggplot() +
      geom_line(data = sto_by_class[[var]], aes(x = time, y = value, color = variable, group = rep),
                size = 1, alpha = .7) 
    if (!is.null(det_var) && !is.null(det_data)) {
      p <- p + geom_line(data = det_data[[det_var]], aes(x = time, y = value), size = 1)
    }
    p + scale_color_manual(values = sto_color) +
      scale_y_continuous(limits = c(0, NA), name = ylab) +
      scale_x_continuous(breaks = seq(0, years * 12, 12), limits = c(0, years * 12))
  }
  
  # Example: If you have det_by_class, pass it to the plotting function as det_data
  Nplot   <- plot_class("N", SEIcols[7], "N", det_by_class, "Total N")
  Splot   <- plot_class("S", SEIcols[6], "S", det_by_class, "Susceptible")
  Eplot   <- plot_class("E1", SEIcols[4], "E", det_by_class, "Exposed")
  Iplot   <- plot_class("I", SEIcols[2], "I", det_by_class, "Infectious")
  SSSplot <- plot_class("SS_S", SEIcols[5], "sS", det_by_class, "Super S")
  ESSplot <- plot_class("SS_E1", SEIcols[3], "sE", det_by_class, "Super Exposed")
  ISSplot <- plot_class("SS_I", SEIcols[1], "sI", det_by_class, "Super Infectious")
  
  # If no deterministic model, just plot stochastic
  Nplot   <- plot_class("N", SEIcols[7], ylab = "Total N")
  Splot   <- plot_class("S", SEIcols[6], ylab = "Susceptible")
  Eplot   <- plot_class("E1", SEIcols[4], ylab = "Exposed")
  Iplot   <- plot_class("I", SEIcols[2], ylab = "Infectious")
  SSSplot <- plot_class("SS_S", SEIcols[5], ylab = "Super S")
  ESSplot <- plot_class("SS_E1", SEIcols[3], ylab = "Super Exposed")
  ISSplot <- plot_class("SS_I", SEIcols[1], ylab = "Super Infectious")
  
  if (exists("save_plots") && save_plots) {
    jpeg(filename = paste0("figures/cont_Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggpubr::ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol = 2, nrow = 3))
    dev.off()
    
    jpeg(filename = paste0("figures/cont_Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(Nplot)
    dev.off()
  }
  
  # Save only what's needed, not the whole list
  if (exists("save_runs") && save_runs) {
    save(sto_long, file = paste0('data/raw/cont_', name_out, size, '-', Sys.Date(), '.RData'))
  }
  
  # Remove everything except what you want to keep for later
  rm(list = setdiff(ls(), c("name_out", "size", "years", "SEIcols", "save_plots", "save_runs")))
  gc()
  
  data_sto <- sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>%
    tidyr::pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  models_res[[1]] <- as.data.frame(data)
  remove(data)
  models_res[[2]] <- as.data.frame(sto_out)
  remove(sto_out)
  models_res[[3]] <- as.data.frame(data_sto)
  remove(data_sto)
  
  gc()
  
  ###########
  ## Plots ##
  ###########
  
  #class comparison plots
  sto_by_class <- split(as.data.frame(models_res[[3]]), as.data.frame(models_res[[3]])$variable)
  det_by_class <- split(as.data.frame(models_res[[1]]), as.data.frame(models_res[[1]])$variable)
  
  SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
  
  Nplot <- ggplot() + 
    geom_line(data = sto_by_class[["N"]], aes(x = time, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
    geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[7]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
  
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
  
  # print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
  # print(Nplot)
  
  if(save_plots){
    jpeg(filename = paste0("figures/cont_Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("figures/cont_Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(Nplot)
    dev.off()
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    save( models_res, file = paste0('data/raw/cont_', name_out, size, '-', Sys.Date(), '.RData'))
  }
  remove(models_res) #clear from workspace
  ##################
  
  gc()
}

##########################################################################################
##########################################################################################
############### Pre-simulations for determination of monthly Lambda values ###############
##########################################################################################
##########################################################################################

####################
## Initialization ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "R/bTBwl_func.R")
years = 3
seedQuarter = 1
prop_superSpreader = 0.0
pct = 0.02
verbose = 0
reps = 200
sizes = c(10, 50, 100, 250, 500, 750, 1000)
infType = "spillover"
runtype = "stdLambda_"
name_out = paste0(runtype,infType,"_",pct,'-',prop_superSpreader,'-')
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
save_runs = TRUE
save_plots = TRUE
test_mode = FALSE; test_birth = FALSE; test_death_n = FALSE; test_death_h = FALSE; test_disease = FALSE;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum

for (size in sizes) {
  ###########
  ## Runs ##
  ##########
  lambda_res <- vector(mode = "list", length = 2) # storage for single herd size results - deterministic, true stochastic, and formatted stochastic data
  lambda_out <- vector(mode = "list", length = 1)
  
  pars <- parameter_set_wl(k = size, 
                           scenario = infType, 
                           initial_exposed = pct*size, 
                           SS_prop = prop_superSpreader,
                           start_quarter = seedQuarter, 
                           test = test_mode,
                           birth = test_birth,
                           death_h = test_death_h,
                           death_n = test_death_n,
                           disease = test_disease,
                           verbose = verbose)
  
  # stochastic model
  initial_state <- data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters <- data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  param_vals <- data.frame(merge(population_parameters, disease_parameters))
  
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'c')
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state, population_parameters, disease_parameters, param_vals)
  
  # assign data frames to list, then clear data from workspace
  lambda_res[[2]] <- as.data.frame(sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% tidyr::pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value"))
  sto_out <- sto_out[,c("time", "Month", "Lambda", "rep")]
  lambda_res[[1]] <- as.data.frame(sto_out)
  lambda_out[[1]] <- data.frame(min = getLambda_vec(data = lambda_res[[1]], type = 'min'), 
                                max = getLambda_vec(data = lambda_res[[1]], type = 'max'),
                                mean = getLambda_vec(data = lambda_res[[1]], type = 'mean'),
                                median = getLambda_vec(data = lambda_res[[1]], type = 'median'),
                                month = c(1:12))
  remove(sto_out)
  
  ###########
  ## Plots ##
  ###########
  lambdaVec_plot <- ggplot(data = lambda_out[[1]]) +
    geom_ribbon(aes(x=month, ymax=max,ymin=min ), color = "darkgray", fill = "aquamarine", alpha=.5) +
    geom_line(aes(x=month, y=median), color = "blue") + 
    geom_line(aes(x=month, y=mean), color = "red")
  
  lambdaVal_plot <- ggplot(data = lambda_res[[1]]) + 
    geom_line(aes(x=time, y=Lambda, group = rep))
  lambdaVal_plot2 <- ggplot(data = lambda_res[[1]]) + 
    geom_line(aes(x=Month, y=Lambda, group = rep))
  
  if(save_plots){
    jpeg(filename = paste0("figures/lambda_preSim_Rplot_LambdaVec_", name_out, size, ".jpeg"))
    print(lambdaVec_plot)
    dev.off()
    
    jpeg(filename = paste0("figures/lambda_preSim_Rplot_LambdaVal_", name_out, size, ".jpeg"))
    print(ggarrange(lambdaVal_plot, lambdaVal_plot2, ncol=1, nrow=2))
    dev.off()
  }
  remove(lambdaVec_plot, lambdaVal_plot, lambdaVal_plot2)
  
  ##################
  ## save outputs ##
  ##################
  if(save_runs){
    save(lambda_res, file = paste0("data/raw/lambda_preSim_",name_out, size, '-', Sys.Date(), '.RData'))
    save(lambda_out, file = paste0('data/raw/vals_', name_out, size, '.RData'))
  }
  remove(lambda_res)
  gc()
}

##########################################################################################
##########################################################################################
#######################     Multiple Herd Sizes -- Discrete time    ######################
##########################################################################################
##########################################################################################

rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "R/bTBwl_func.R")
years = 7
scaled = TRUE
times = seq(from=0, to=years*12, by=.5)
seedQuarter = 1
prop_superSpreader = 0.1
pct = 0.02
verbose = 0
reps = 500
sizes = c(10, 50, 100, 250, 500, 750)
infType = "seeded"
runtype = "stdD_pub_"
type_of_integral = 3
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
lambda_name_out = paste0('stdLambda_', infType,"_", pct,'-', prop_superSpreader, '-')
lambda_factor = 1.2
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/"
test_mode = FALSE; test_birth = FALSE; test_death_n = FALSE; test_death_h = FALSE; test_disease = FALSE;
save_runs = TRUE
save_plots = TRUE

for (size in sizes) {
  models_res <- vector(mode = "list", length = 3)
  
  pars <- parameter_set_wl(k = size, 
                           scenario = infType, 
                           initial_exposed = pct*size, 
                           SS_prop = prop_superSpreader,
                           start_quarter = seedQuarter, 
                           test = test_mode,
                           birth = test_birth,
                           death_h = test_death_h,
                           death_n = test_death_n,
                           disease = test_disease,
                           verbose = verbose)
  
  # deterministic model
  X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
              sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
  tic = Sys.time()
  out <- ode(
    func = SEI_model_full,
    y = X0_full,
    times = times,
    parms = pars,
    method = "rk4"
  ) %>% as.data.frame()
  toc = Sys.time()
  print(toc - tic)
  
  out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
  data <- out %>% tidyr::gather(variable,value,-time)
  remove(out)
  data$value[is.nan(data$value)] <- 0 
  data$value[which(data$value < 0)] <- 0
  
  # stochastic model
  initial_state <- data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters <- data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  param_vals <- data.frame(merge(population_parameters, disease_parameters))
  
  # get lambda
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 5, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
  toc = Sys.time()
  print(toc-tic)
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  remove(sto_out)
  
  # run
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state, population_parameters, disease_parameters, param_vals)
  
  data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% tidyr::pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  models_res[[1]] <- as.data.frame(data)
  remove(data)
  models_res[[2]] <- as.data.frame(sto_out)
  remove(sto_out)
  models_res[[3]] <- as.data.frame(data_sto)
  remove(data_sto)
  
  ###########
  ## Plots ##
  ###########
  sto_by_class <- split(as.data.frame(models_res[[3]]), as.data.frame(models_res[[3]])$variable)
  det_by_class <- split(as.data.frame(models_res[[1]]), as.data.frame(models_res[[1]])$variable)
  SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
  
  # susceptible
  sto_by_class[['S']]$value <- sto_by_class[["S"]]$value
  det_by_class[['S']]$value <- det_by_class[["S"]]$value
  # susceptible SS
  sto_by_class[['SS_S']]$value <- sto_by_class[["SS_S"]]$value
  det_by_class[['sS']]$value <- det_by_class[['sS']]$value
  # exposed -- merged SS
  sto_by_class[['E1']]$value <- sto_by_class[["E1"]]$value + sto_by_class[["SS_E1"]]$value
  det_by_class[['E']]$value <- det_by_class[["E"]]$value + det_by_class[["sE"]]$value
  # infected -- merged SS
  sto_by_class[['I']]$value <- sto_by_class[["I"]]$value + sto_by_class[["SS_I"]]$value
  det_by_class[['I']]$value <- det_by_class[["I"]]$value + det_by_class[["sI"]]$value
  
  Nplot <- ggplot() + 
    geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
    geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
    scale_color_manual(values=SEIcols[7]) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_x_continuous(breaks=seq(0,years*12,1), limits = c(0, years*12))
  
  Splot <- ggplot() + 
    geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
    geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
    scale_color_manual(values=SEIcols[6]) +
    theme(text = element_text(size = 16)) +
    scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
    scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  Eplot <- ggplot() + 
    geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
    geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
    scale_color_manual(values=SEIcols[4]) +
    theme(text = element_text(size = 16)) +
    scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
    scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  Iplot <- ggplot() + 
    geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
    geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
    scale_color_manual(values=SEIcols[2]) +
    theme(text = element_text(size = 16)) +
    scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
    scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  SSSplot <- ggplot() + 
    geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) + 
    geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
    scale_color_manual(values=SEIcols[5]) +
    theme(text = element_text(size = 16)) +
    scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
    scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
  
  print(ggpubr::ggarrange(Splot, Eplot,  Iplot, SSSplot, ncol=3, nrow=2))
  print(Nplot)
  
  if(scaled){
    # Proportion susceptible
    sto_by_class[['S']]$value <- sto_by_class[["S"]]$value/sto_by_class[["N"]]$value
    det_by_class[['S']]$value <- det_by_class[["S"]]$value/det_by_class[["N"]]$value
    # Proportion susceptible SS
    sto_by_class[['SS_S']]$value <- sto_by_class[["SS_S"]]$value/sto_by_class[["N"]]$value
    det_by_class[['sS']]$value <- det_by_class[['sS']]$value/det_by_class[["N"]]$value
    # Proportion exposed -- merged SS
    sto_by_class[['E1']]$value <- sto_by_class[["E1"]]$value/sto_by_class[["N"]]$value 
    det_by_class[['E']]$value <- det_by_class[["E"]]$value/det_by_class[["N"]]$value 
    # Proportion infected -- merged SS
    sto_by_class[['I']]$value <- sto_by_class[["I"]]$value/sto_by_class[["N"]]$value 
    det_by_class[['I']]$value <- det_by_class[["I"]]$value/det_by_class[["N"]]$value 
    
    Splot_scl <- ggplot() + 
      geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
      geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=SEIcols[6]) +
      theme(text = element_text(size = 16)) +
      scale_y_continuous(name = 'Proportion S' ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
    Eplot_scl <- ggplot() + 
      geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
      geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=SEIcols[4]) +
      theme(text = element_text(size = 16)) +
      scale_y_continuous(name = 'Proportion E1' ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
    Iplot_scl <- ggplot() + 
      geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
      geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=SEIcols[2]) +
      theme(text = element_text(size = 16)) +
      scale_y_continuous(name = 'Proportion I' ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
    SSSplot_scl <- ggplot() + 
      geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) + 
      geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=SEIcols[5]) +
      theme(text = element_text(size = 16)) +
      scale_y_continuous(name = expression('Proportion S'['SS']) ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12))
  }
  
  if(save_plots){
    plot <- ggpubr::ggarrange(Splot, Eplot,  Iplot, SSSplot, ncol=3, nrow=2)
    plot_scl <- ggpubr::ggarrange(Splot_scl, Eplot_scl,  Iplot_scl, SSSplot_scl, ncol=3, nrow=2)
    plot <- ggpubr::annotate_figure(plot, top = ggpubr::text_grob("Discrete time stochastic model outbreak trajectories", color = "black", face = "bold", size = 18))
    plot_scl <- ggpubr::annotate_figure(plot_scl, top = ggpubr::text_grob("Discrete time stochastic model outbreak trajectories", color = "black", face = "bold", size = 18))
    jpeg(filename = paste0("figures/disc_Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(plot)
    dev.off()
    jpeg(filename = paste0("figures/disc_Rplot_class_scl_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(plot_scl)
    dev.off()
    remove(plot, plot_scl)
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  
  ##################
  ## save outputs ##
  ##################
  if(save_runs){
    save(models_res, file = paste0('data/raw/disc_', name_out, size, '-', Sys.Date(), '.RData'))
  }
  remove(models_res)
  gc()
}


# ####################
# ##    Fadeout     ##
# ####################
# # takes about 8 min to run 20x21
# rm(list = ls())
# setwd(dirname(getActiveDocumentContext()$path))
# source(file = "R/bTBwl_func.R")
# 
# years = 20
# times = seq(from=0, to=years*12, by=12/365) # daily time steps
# seedQuarter = 1
# prop_superSpreader = 0
# verbose = 0
# reps = 500
# lambda_factor = 1.2
# sizes = seq(from = 0, to = 750, length.out = 31)[-1]
# pct = seq(from = 0, to = .5, length.out = 21) # pct inf - fadeout
# infType = "seeded"
# runtype = "fade_"
# type_of_integral = 3
# name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
# pth = getwd()
# save_runs = FALSE
# save_plots = FALSE
# save_plots_final = TRUE
# fade_data = matrix(data = NA, nrow = length(sizes), ncol = length(pct))
# fadeTime_data = matrix(data = NA, nrow = length(sizes), ncol = length(pct))
# 
# for (i in seq_along(sizes)) {
#   data_over_ranges = vector(mode = 'list', length = length(pct))
#   fade_data_out = matrix(data = NA, nrow = 2, ncol = length(pct))
#   size = as.integer(sizes[i])
#   
#   for (j in seq_along(pct)) {
#     pars <- parameter_set_wl(k = size, 
#                              scenario = infType, 
#                              initial_exposed = pct[j] * size, 
#                              SS_prop = 0,
#                              start_quarter = seedQuarter, 
#                              verbose = verbose)
#     
#     # stochastic model
#     initial_state <- data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
#     population_parameters <- data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
#     disease_parameters <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
#     param_vals <- data.frame(merge(population_parameters, disease_parameters))
#     
#     # get lambda
#     tic = Sys.time()
#     sto_out <- wildlife_model(
#       n_reps = 200, 
#       parameters = param_vals, 
#       initial_state = initial_state, 
#       nyears = 5, 
#       seed_quarter = as.integer(pars["start_q"]), 
#       verbose = -3, 
#       batch_name = "test", 
#       type = 'c'
#     )
#     toc = Sys.time()
#     print(toc - tic)
#     lambda_out <- getLambda_vec(data = sto_out, type = 'max')
#     remove(sto_out)
#     
#     # run with new lambda
#     tic = Sys.time()
#     sto_out <- wildlife_model(
#       n_reps = reps, 
#       parameters = param_vals, 
#       initial_state = initial_state, 
#       nyears = years, 
#       seed_quarter = as.integer(pars["start_q"]), 
#       verbose = verbose, 
#       batch_name = "test", 
#       type = 'd', 
#       lambda = lambda_out * lambda_factor, 
#       integrate_type = type_of_integral
#     )
#     toc = Sys.time()
#     print(toc - tic)
#     remove(initial_state, population_parameters, disease_parameters, param_vals)
#     
#     # assign data frames to list, then clear data from workspace
#     data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>%
#       pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
#     sto_out <- as.data.frame(sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", "fadeout", "fadeout time")])
#     data_over_ranges[[j]] <- aggregate(data = subset(sto_out, select = -c(rep)), .~tstep, FUN = mean)
#     A <- aggregate(data = subset(sto_out, select = c(rep, `fadeout time`)), .~rep, FUN = max) 
#     FOT <- mean(A$`fadeout time`[A$`fadeout time` != 0])
#     if (is.nan(FOT)) { FOT = 1e10 }
#     remove(sto_out, A)
#     
#     # Plots (optional, unchanged)
#     if (save_plots) {
#       X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]),
#                   sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
#       tic = Sys.time()
#       out <- ode(
#         func = SEI_model_full,
#         y = X0_full,
#         times = times,
#         parms = pars,
#         method = "rk4"
#       ) %>% as.data.frame()
#       toc = Sys.time()
#       print(toc - tic)
#       
#       out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
#       data <- out %>% gather(variable, value, -time)
#       remove(out)
#       #remove NaN or negative values
#       data$value[is.nan(data$value)] <- 0 
#       data$value[which(data$value < 0)] <- 0
#       
#       # class comparison plots (unchanged)
#       sto_by_class <- split(as.data.frame(data_sto), as.data.frame(data_sto)$variable)
#       det_by_class <- split(as.data.frame(data), as.data.frame(data)$variable)
#       SEIcols <- RColorBrewer::brewer.pal(11, "Spectral")[c(1,2,4,5,8,9,10)]
#       
#       Nplot <- ggplot() + 
#         geom_line(data = sto_by_class[["N"]], aes(x = tstep, y = value, color = variable, group = rep), size = 1, alpha = .9) + 
#         geom_line(data = det_by_class[["N"]], aes(x = time, y = value), size = 1) +
#         scale_color_manual(values = SEIcols[7]) +
#         scale_y_continuous(limits = c(0, NA)) +
#         scale_x_continuous(breaks = seq(0, years * 12, 12), limits = c(0, years * 12))
#       
#       Splot <- ggplot() + 
#         geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
#         geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1) +
#         scale_color_manual(values=SEIcols[6]) +
#         scale_y_continuous(limits = c(0,NA)) +
#         scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
#       
#       Eplot <- ggplot() + 
#         geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
#         geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1) +
#         scale_color_manual(values=SEIcols[4]) +
#         scale_y_continuous(limits = c(0,NA)) +
#         scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
#       
#       Iplot <- ggplot() + 
#         geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
#         geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1) +
#         scale_color_manual(values=SEIcols[2]) +
#         scale_y_continuous(limits = c(0,NA)) +
#         scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
#       
#       SSSplot <- ggplot() + 
#         geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
#         geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1) +
#         scale_color_manual(values=SEIcols[5]) +
#         scale_y_continuous(limits = c(0,NA)) +
#         scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
#       
#       ESSplot <- ggplot() + 
#         geom_line(data = sto_by_class[["SS_E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
#         geom_line(data = det_by_class[["sE"]], aes(x = time, y=value), size = 1) +
#         scale_color_manual(values=SEIcols[3]) +
#         scale_y_continuous(limits = c(0,NA)) +
#         scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
#       
#       ISSplot <- ggplot() + 
#         geom_line(data = sto_by_class[["SS_I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
#         geom_line(data = det_by_class[["sI"]], aes(x = time, y=value), size = 1) +
#         scale_color_manual(values=SEIcols[1]) +
#         scale_y_continuous(limits = c(0,NA)) +
#         scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
#       
#       # setwd(paste0(pth, "Disc_runs/fadeout_runs/validation_plots/"))
#       
#       jpeg(filename = paste0("figures/Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
#       print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
#       dev.off()
#       
#       jpeg(filename = paste0("figures/Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
#       print(Nplot)
#       dev.off()
#       
#       # setwd(pth)
#       remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot, det_by_class, sto_by_class, data, data_sto)
#       
#     }
#     ###########
#     
#     # save fadeout results to data matrix
#     p_fade <- data_over_ranges[[j]]$fadeout[length(data_over_ranges[[j]]$fadeout)] # extract mean fade value from last averaged row
#     fade_data_out[1,j] <- p_fade 
#     fade_data_out[2,j] <- FOT
#   }
#   
#   ##################
#   ## save outputs ##
#   ##################
#   
#   if(save_runs){
#     name = c(paste0('data_over_ranges',i))
#     assign(x=name, value = data_over_ranges)
#     save( list = name[1], file = paste0('data/raw/fadeout_dat',  size, '-', Sys.Date(), '.RData'))
#   }
#   
#   ##################
#   fade_data_out
# }

####################
##    Fadeout     ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "R/bTBwl_func.R")

years = 20
times = seq(from=0, to=years*12, by=12/365) # daily time steps
seedQuarter = 1
prop_superSpreader = 0
verbose = 0
reps = 500
lambda_factor = 1.2
sizes = seq(from = 0, to = 750, length.out = 31)[-1]
pct = seq(from = 0, to = .5, length.out = 21) # pct inf - fadeout
infType = "seeded"
runtype = "fade_"
type_of_integral = 3
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
pth = getwd()
save_runs = FALSE
save_plots = FALSE
save_plots_final = TRUE

fade_data     = matrix(NA, nrow=length(sizes), ncol=length(pct))
fadeTime_data = matrix(NA, nrow=length(sizes), ncol=length(pct))

for (i in seq_along(sizes)) {
  size = as.integer(sizes[i])
  for (j in seq_along(pct)) {
    pars <- parameter_set_wl(
      k = size,
      scenario = infType,
      initial_exposed = pct[j] * size,
      SS_prop = 0,
      start_quarter = seedQuarter,
      verbose = verbose
    )
    # stochastic model
    initial_state <- data.frame(
      S_0 = pars["S_0"], E1_0 = pars["E1_0"], I_0 = pars["I_0"],
      SuperS_0 = pars["SuperS_0"], SuperE1_0 = pars["SuperE1_0"], SuperI_0 = pars["SuperI_0"]
    )
    population_parameters <- data.frame(
      K = pars["K"], eta_hunt = pars["eta_hunt"], eta_nat = pars["eta_nat"], theta = pars["theta"],
      gamma = pars["gamma"], alpha_max = pars["alpha_max"], ksi = pars["ksi"], omega = pars["omega"],
      s = pars["s"], alpha = pars["alpha"]
    )
    disease_parameters <- data.frame(
      beta = pars["beta"], area = pars["area"], p1 = pars["p1"],
      p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], p2_q4 = pars["p2_q4"],
      phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], sigma1_rate = pars["sigma1_rate"]
    )
    param_vals <- data.frame(merge(population_parameters, disease_parameters))
    
    # get lambda
    sto_out <- wildlife_model(
      n_reps = 200,
      parameters = param_vals,
      initial_state = initial_state,
      nyears = 5,
      seed_quarter = as.integer(pars["start_q"]),
      verbose = -3,
      batch_name = "test",
      type = 'c'
    )
    lambda_out <- getLambda_vec(data = sto_out, type = 'max')
    rm(sto_out)
    
    # run with new lambda
    sto_out <- wildlife_model(
      n_reps = reps,
      parameters = param_vals,
      initial_state = initial_state,
      nyears = years,
      seed_quarter = as.integer(pars["start_q"]),
      verbose = verbose,
      batch_name = "test",
      type = 'd',
      lambda = lambda_out * lambda_factor,
      integrate_type = type_of_integral
    )
    rm(initial_state, population_parameters, disease_parameters, param_vals)
    
    # extract fadeout/fadeout time
    # Use tail() to get the last tstep for fadeout
    fadeout_last <- tail(sto_out$fadeout, 1)
    fade_data[i, j] <- fadeout_last
    
    # Get mean fadeout time for reps that did not fadeout instantly (nonzero)
    FOT_vals <- tapply(sto_out$`fadeout time`, sto_out$rep, max)
    FOT <- mean(FOT_vals[FOT_vals != 0], na.rm=TRUE)
    if (is.nan(FOT) || length(FOT) == 0) FOT = 1e10
    fadeTime_data[i, j] <- FOT
    
    rm(sto_out, FOT_vals)
    gc()
  }
}

# fade_data and fadeTime_data are now matrices with [length(sizes), length(pct)]
# Optionally, save results:
if (save_runs) {
  save(fade_data, fadeTime_data, file = paste0("data/raw/fadeout_matrix_", Sys.Date(), ".RData"))
}
