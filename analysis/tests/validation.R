##########################################################################################
##########################################################################################
#######################    Multiple Herd Sizes -- test mode output   #####################
##########################################################################################
##########################################################################################

####################
## Initialization ##
####################
rm(list = ls()) #clear workspace for memory
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))
source(file = "bTBwl_func.R") #load model and supplementary functions/libraries to the workspace
years = 3 #number of years for the simulation to run
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1 #staring quarter for the simulation
prop_superSpreader = 0.05 #percent of SS individuals in the population 
pct = 0.02 #initial percent of the population exposed
verbose = 0 #model diagnostic outputs - note that any message will likely cause the data to import incorrectly
reps = 50 #stochastic simulation replicates
sizes = c(10, 50, 100, 250, 500) #vector of herd sizes
infType = "seeded" # infection either 'seeded' or 'spillover' - changes forces of infection
runtype = "test_" # type of run - std for standard, SS for super susceptible, H for hunting, 
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-') #output data and plots name
pth = getwd() # "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
#must also have output directories 'runs/validation_plots' and 'lambda_preSim/lambda_plots'
save_runs = T #T if runs should be saved to the appropriate folder
save_plots = T #T if plots should be saved to the appropriate folder.
test_mode = T; test_birth = T; test_death_n = T; test_death_h = T; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { 
  library(deSolve); 
  library(tidyverse);
  library(ggplot2);
  library(ggpubr);
  library(grid);
  library(rstudioapi);
  library(foreach);
  library(doParallel);
  library(RColorBrewer);
  library(scales); })
#####################

foreach(i = 1:length(sizes)) %dopar% {
  ###########
  ## Runs ##
  ##########
  
  setwd(pth)
  models_res <- vector(mode = "list", length = 3) # storage for single herd size results - deterministic, true stochastic, and formatted stochastic data
  
  size=as.integer(sizes[i])
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
                           verbose = verbose) #set parameters dependent 
  
  #deterministic model
  X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
              sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
  
  tic = Sys.time()
  ode(
    func = SEI_model_full,
    y = X0_full,
    times = times,
    parms = pars,
    method = "rk4"
  ) %>%
    as.data.frame() -> out
  toc = Sys.time()
  print(toc-tic)
  
  out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
  #out <- out[, colSums(out != 0) > 0] # remove classes if they are all 0
  data <- out %>% gather(variable,value,-time)
  remove(out)
  #remove NaN or negative values
  data$value[is.nan(data$value)] <- 0 
  data$value[which(data$value < 0)] <- 0
  
  
  #stochastic model
  initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  
  param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'c')
  #n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears = years; seed_quarter = as.integer(pars["start_q"]); verbose = verbose; batch_name = "test"; type = 'c';
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state,population_parameters,disease_parameters,param_vals)
  
  
  data_sto <- sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  #assign data frames to list
  #then clear data from workspace
  
  models_res[[1]] <- as.data.frame(data)
  remove(data)
  models_res[[2]] <- as.data.frame(sto_out)
  remove(sto_out)
  models_res[[3]] <- as.data.frame(data_sto)
  remove(data_sto)
  
  ##########
  
  
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
  
  print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
  print(Nplot)
  
  if(save_plots){
    #setwd(paste0(pth,"test_runs/validation_plots/"))
    
    jpeg(filename = paste0("tests/results/Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("tests/results/Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(Nplot)
    dev.off()
    
    #setwd(pth)
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    
    save( models_res, file = paste0('tests/results/', name_out, size, '-', Sys.Date(), '.RData'))
    
  }
  remove(models_res) #clear from workspace
  ##################
  
}

stopCluster(cl)

##########################################################################################
##########################################################################################
#######################    Multiple Herd Sizes -- Continuous Time   ######################
##########################################################################################
##########################################################################################

####################
## Initialization ##
####################
rm(list = ls()) #clear workspace for memory
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R") #load model and supplementary functions/libraries to the workspace
years = 10 #number of years for the simulation to run
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1 #staring quarter for the simulation
prop_superSpreader = 0.05 #percent of SS individuals in the population 
pct = 0.02 #initial percent of the population exposed
verbose = 0 #model diagnostic outputs - note that any message will likely cause the data to import incorrectly
reps = 500 #stochastic simulation replicates
sizes = c(10, 50, 100, 250, 500) #vector of herd sizes
infType = "seeded" # infection either 'seeded' or 'spillover' - changes forces of infection
runtype = "stdPres_" # type of run - std for standard, SS for super susceptible, H for hunting, 
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-') #output data and plots name
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
#must also have output directories 'runs/validation_plots' and 'lambda_preSim/lambda_plots'
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
save_runs = T #T if runs should be saved to the appropriate folder
save_plots = T #T if plots should be saved to the appropriate folder.
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { 
  library(deSolve); 
  library(tidyverse);
  library(ggplot2);
  library(ggpubr);
  library(grid);
  library(rstudioapi);
  library(foreach);
  library(doParallel);
  library(RColorBrewer);
  library(scales); })
#####################

foreach(i = 1:length(sizes)) %dopar% {
  ###########
  ## Runs ##
  ##########
  
  setwd(pth)
  models_res <- vector(mode = "list", length = 3) # storage for single herd size results - deterministic, true stochastic, and formatted stochastic data
  
  size=as.integer(sizes[i])
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
                           verbose = verbose) #set parameters dependent 
  
  #deterministic model
  X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
              sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
  
  tic = Sys.time()
  ode(
    func = SEI_model_full,
    y = X0_full,
    times = times,
    parms = pars,
    method = "rk4"
  ) %>%
    as.data.frame() -> out
  toc = Sys.time()
  print(toc-tic)
  
  out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
  #out <- out[, colSums(out != 0) > 0] # remove classes if they are all 0
  data <- out %>% gather(variable,value,-time)
  remove(out)
  #remove NaN or negative values
  data$value[is.nan(data$value)] <- 0 
  data$value[which(data$value < 0)] <- 0
  
  
  #stochastic model
  initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  
  param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'c')
  #n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears = years; seed_quarter = as.integer(pars["start_q"]); verbose = verbose; batch_name = "test"; type = 'c';
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state,population_parameters,disease_parameters,param_vals)
  
  
  data_sto <- sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  #assign data frames to list
  #then clear data from workspace
  
  models_res[[1]] <- as.data.frame(data)
  remove(data)
  models_res[[2]] <- as.data.frame(sto_out)
  remove(sto_out)
  models_res[[3]] <- as.data.frame(data_sto)
  remove(data_sto)
  
  ##########
  
  
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
  
  print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
  print(Nplot)
  
  if(save_plots){
    #setwd(paste0(pth,"Cont_runs/validation_plots/"))
    
    jpeg(filename = paste0("/results/cont_Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("/results/cont_Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(Nplot)
    dev.off()
    
    setwd(pth)
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    
    save( models_res, file = paste0('Cont_runs/', name_out, size, '-', Sys.Date(), '.RData'))
    
  }
  remove(models_res) #clear from workspace
  ##################
  
}

stopCluster(cl)

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
source(file = "bTBwl_func.R")
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
#must also have output directories 'runs/validation_plots' and 'lambda_preSim/lambda_plots'
save_runs = T
save_plots = T
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { 
  library(deSolve); 
  library(tidyverse);
  library(ggplot2);
  library(ggpubr);
  library(grid);
  library(rstudioapi);
  library(foreach);
  library(doParallel);
  library(RColorBrewer);
  library(scales); })
#####################
foreach(i = 1:length(sizes)) %dopar% {
  
  ###########
  ## Runs ##
  ##########
  
  setwd(pth)
  lambda_res <- vector(mode = "list", length = 2) # storage for single herd size results - deterministic, true stochastic, and formatted stochastic data
  lambda_out <- vector(mode = "list", length = 1)
  
  size=as.integer(sizes[i])
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
                           verbose = verbose) #set parameters dependent 
  
  #stochastic model
  initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  
  param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'c')
  #n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears = years; seed_quarter = as.integer(pars["start_q"]); verbose = verbose; batch_name = "test"; type = 'c';
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state,population_parameters,disease_parameters,param_vals)
  
  
  
  
  #assign data frames to list
  #then clear data from workspace
  lambda_res[[2]] <- as.data.frame(sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value"))
  sto_out <- sto_out[,c("time", "Month", "Lambda", "rep")]
  lambda_res[[1]] <- as.data.frame(sto_out) 
  lambda_out[[1]] <- data.frame(min = getLambda_vec(data = lambda_res[[1]], type = 'min'), 
                                max = getLambda_vec(data = lambda_res[[1]], type = 'max'),
                                mean = getLambda_vec(data = lambda_res[[1]], type = 'mean'),
                                median = getLambda_vec(data = lambda_res[[1]], type = 'median'),
                                month = c(1:12))
  remove(sto_out)
  ##########
  
  
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
    #setwd(paste0(pth, "lambda_preSim/lambda_plots/"))
    
    jpeg(filename = paste0("lambda_preSim_Rplot_LambdaVec_", name_out, size, ".jpeg"))
    print(lambdaVec_plot)
    dev.off()
    
    jpeg(filename = paste0("lambda_preSim_Rplot_LambdaVal_", name_out, size, ".jpeg"))
    print(ggarrange(lambdaVal_plot, lambdaVal_plot2, ncol=1, nrow=2))
    dev.off()
    
    setwd(pth)
  }
  remove(lambdaVec_plot, lambdaVal_plot, lambdaVal_plot2)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    setwd(paste0(pth, "lambda_preSim/"))
    
    save( lambda_res, file = paste0( name_out, size, '-', Sys.Date(), '.RData'))
    save( lambda_out, file = paste0('vals_', name_out, size, '.RData'))
    
    setwd(pth)
  }
  remove(lambda_res) #clear from workspace
  ##################
  
}

stopCluster(cl)

##########################################################################################
##########################################################################################
#######################     Multiple Herd Sizes -- Discrete time    ######################
##########################################################################################
##########################################################################################


####################
## Initialization ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
years = 7
scaled = T
#times = seq(from=0, to=years*12, by=12/365) #daily time steps
times = seq(from=0, to=years*12, by=.5) #daily time steps
seedQuarter = 1
prop_superSpreader = 0.1
pct = 0.02
verbose = 0
reps = 500
sizes = c(10, 50, 100, 250, 500, 750)
infType = "seeded"
runtype = "stdD_pub_"
type_of_integral = 3;
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
lambda_name_out = paste0('stdLambda_', infType,"_", pct,'-', prop_superSpreader, '-')
lambda_factor = 1.2
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
#must also have output directories 'runs/validation_plots' and 'lambda_preSim/lambda_plots'
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
save_runs = T
save_plots = T
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { 
  library(deSolve); 
  library(tidyverse);
  library(ggplot2);
  library(ggpubr);
  library(grid);
  library(rstudioapi);
  library(foreach);
  library(doParallel);
  library(RColorBrewer);
  library(scales); })
#####################

foreach(i = 1:length(sizes)) %dopar% {
  ###########
  ## Runs ##
  ##########
  
  setwd(pth)
  models_res <- vector(mode = "list", length = 3) # storage for single herd size results - deterministic, true stochastic, and formatted stochastic data
  
  size=as.integer(sizes[i])
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
                           verbose = verbose) #set parameters dependent 
  
  #deterministic model
  X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
              sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
  
  tic = Sys.time()
  ode(
    func = SEI_model_full,
    y = X0_full,
    times = times,
    parms = pars,
    method = "rk4"
  ) %>%
    as.data.frame() -> out
  toc = Sys.time()
  print(toc-tic)
  
  out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
  #out <- out[, colSums(out != 0) > 0] # remove classes if they are all 0
  data <- out %>% gather(variable,value,-time)
  remove(out)
  #remove NaN or negative values
  data$value[is.nan(data$value)] <- 0 
  data$value[which(data$value < 0)] <- 0
  
  
  #stochastic model
  initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  
  param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
  #get lambda
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 5, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
  toc = Sys.time()
  print(toc-tic)
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  remove(sto_out)
  
  #run
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
  #n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears = years; seed_quarter = as.integer(pars["start_q"]); verbose = verbose; batch_name = "test"; type = 'c';
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state,population_parameters,disease_parameters,param_vals)
  
  
  data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  #assign data frames to list
  #then clear data from workspace
  
  models_res[[1]] <- as.data.frame(data)
  remove(data)
  models_res[[2]] <- as.data.frame(sto_out)
  remove(sto_out)
  models_res[[3]] <- as.data.frame(data_sto)
  remove(data_sto)
  
  ##########
  
  
  ###########
  ## Plots ##
  ###########
  
  #class comparison plots
  sto_by_class <- split(as.data.frame(models_res[[3]]), as.data.frame(models_res[[3]])$variable)
  det_by_class <- split(as.data.frame(models_res[[1]]), as.data.frame(models_res[[1]])$variable)
  
  SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
  
  #susceptible
  sto_by_class[['S']]$value <- sto_by_class[["S"]]$value
  det_by_class[['S']]$value <- det_by_class[["S"]]$value
  
  #susceptible SS
  sto_by_class[['SS_S']]$value <- sto_by_class[["SS_S"]]$value
  det_by_class[['sS']]$valu <- det_by_class[['sS']]$value
  
  #exposed -- merged SS
  sto_by_class[['E1']]$value <- sto_by_class[["E1"]]$value + sto_by_class[["SS_E1"]]$value
  det_by_class[['E']]$value <- det_by_class[["E"]]$value + det_by_class[["sE"]]$value
  
  #infected -- merged SS
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
  
  
  print(ggarrange(Splot, Eplot,  Iplot, SSSplot, ncol=3, nrow=2))
  print(Nplot)
  
  if(scaled){
    #Proportion susceptible
    sto_by_class[['S']]$value <- sto_by_class[["S"]]$value/sto_by_class[["N"]]$value
    det_by_class[['S']]$value <- det_by_class[["S"]]$value/det_by_class[["N"]]$value
    
    #Proportion susceptible SS
    sto_by_class[['SS_S']]$value <- sto_by_class[["SS_S"]]$value/sto_by_class[["N"]]$value
    det_by_class[['sS']]$value <- det_by_class[['sS']]$value/det_by_class[["N"]]$value
    
    
    #Proportion exposed -- merged SS
    sto_by_class[['E1']]$value <- sto_by_class[["E1"]]$value/sto_by_class[["N"]]$value 
    det_by_class[['E']]$value <- det_by_class[["E"]]$value/det_by_class[["N"]]$value 
    
    #Proportion infected -- merged SS
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
    # setwd(paste0(pth, "Disc_runs/validation_plots/"))
    
    plot <- ggarrange(Splot, Eplot,  Iplot, SSSplot, ncol=3, nrow=2)
    plot_scl <- ggarrange(Splot_scl, Eplot_scl,  Iplot_scl, SSSplot_scl, ncol=3, nrow=2)
    
    plot <- annotate_figure(plot, top = text_grob("Discrete time stochastic model outbreak trajectories", 
                                                  color = "black", face = "bold", size = 18))
    plot_scl <- annotate_figure(plot_scl, top = text_grob("Discrete time stochastic model outbreak trajectories", 
                                                          color = "black", face = "bold", size = 18))
    
    jpeg(filename = paste0("results/disc_Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(plot)
    dev.off()
    
    jpeg(filename = paste0("results/disc_Rplot_class_scl_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(plot_scl)
    dev.off()
    
    #jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    #print(Nplot)
    #dev.off()
    
    remove(plot, plot_scl)
    setwd(pth)
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    save( models_res, file = paste0('Disc_runs/', name_out, size, '-', Sys.Date(), '.RData'))
  }
  remove(models_res) #clear from workspace
  ##################
  
}

stopCluster(cl)