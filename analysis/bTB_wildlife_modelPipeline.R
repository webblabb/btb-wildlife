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
source(file = "R/bTBwl_func.R") #load model and supplementary functions/libraries to the workspace
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
pth = getwd() #"/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
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
    setwd(paste0(pth,"/results/"))
    
    jpeg(filename = paste0("cont_Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("cont_Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
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
    setwd(paste0(pth,"/data/"))
    
    save( models_res, file = paste0('cont_', name_out, size, '-', Sys.Date(), '.RData'))
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

##########################################################################################
##########################################################################################
#######################               Other Plot Types              ######################
##########################################################################################
##########################################################################################

####################
##    Fadeout     ##
####################
#takes about 8 min to run 20x21
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
years = 20
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
prop_superSpreader = 0
verbose = 0
reps = 500
lambda_factor = 1.2
sizes = seq(from = 0, to = 750, length.out = 31)[-1]
#range_vec = c(2,4,6,8,10) #num inf - fadeout
pct = seq(from = 0, to = .5, length.out = 21) #pct inf - fadeout
infType = "seeded"
runtype = "fade_"
type_of_integral = 3;
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
pth = getwd() # "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
save_runs = F
save_plots = F
save_plots_final = T
fade_data = matrix(data = NA, nrow = length(sizes), ncol = length(pct))
fadeTime_data = matrix(data = NA, nrow = length(sizes), ncol = length(pct))
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

R <- foreach(i = 1:length(sizes), .combine='cbind', .inorder=TRUE) %dopar% {

  data_over_ranges = vector( mode = 'list', length = length(pct))
  fade_data_out = matrix(data = NA, nrow = 2, ncol = length(pct))
  
  for(j in 1:length(pct)){
    
    ###########
    ## Runs ##
    ##########
    
    setwd(pth)
    size=as.integer(sizes[i])
    pars <- parameter_set_wl(k = size, 
                             scenario = infType, 
                             initial_exposed = pct[j]*size, 
                             SS_prop = 0,
                             start_quarter = seedQuarter, 
                             verbose = verbose) #set parameters dependent 
    
    
    
    #stochastic model
    initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
    population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
    disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #get lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = 200, parameters = param_vals, initial_state = initial_state, nyears = 5, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
    toc = Sys.time()
    print(toc-tic)
    lambda_out <- getLambda_vec(data = sto_out, type = 'max')
    remove(sto_out)
    
    #run with new lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
    toc = Sys.time()
    print(toc-tic)
    remove(initial_state,population_parameters,disease_parameters,param_vals)
    
    #assign data frames to list
    #then clear data from workspace
    data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
    sto_out <- as.data.frame(sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", 'fadeout', 'fadeout time') ])
    data_over_ranges[[j]] <- aggregate(data=subset(sto_out, select=-c(rep)), .~tstep, FUN = mean) #determine average trajectory
    A <- aggregate(data=subset(sto_out, select=c(rep, `fadeout time`)), .~rep, FUN = max) 
    FOT <- mean( A$`fadeout time`[A$`fadeout time`!=0] )
    if(is.nan(FOT)){FOT=1e10}
    remove(sto_out, A)
    ##########
    
    ###########
    ## Plots ##
    ###########
    if(save_plots){
      
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
      
      #class comparison plots
      sto_by_class <- split(as.data.frame(data_sto), as.data.frame(data_sto)$variable)
      det_by_class <- split(as.data.frame(data), as.data.frame(data)$variable)
      SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
      
      Nplot <- ggplot() + 
        geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
        geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[7]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Splot <- ggplot() + 
        geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[6]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Eplot <- ggplot() + 
        geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[4]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Iplot <- ggplot() + 
        geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[2]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      SSSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[5]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ESSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sE"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[3]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ISSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sI"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[1]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      setwd(paste0(pth, "Disc_runs/fadeout_runs/validation_plots/"))
      
      jpeg(filename = paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
      dev.off()
      
      jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(Nplot)
      dev.off()
      
      setwd(pth)
      remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot, det_by_class, sto_by_class, data, data_sto)
      
    }
    ###########
    
    # save fadeout results to data matrix
    p_fade <- data_over_ranges[[j]]$fadeout[length(data_over_ranges[[j]]$fadeout)] # extract mean fade value from last averaged row
    fade_data_out[1,j] <- p_fade 
    fade_data_out[2,j] <- FOT
  }
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    name = c(paste0('data_over_ranges',i))
    assign(x=name, value = data_over_ranges)
    save( list = name[1], file = paste0('Disc_runs/fadeout_runs/FO_dat',  size, '-', Sys.Date(), '.RData'))
  }
  
  ##################
  fade_data_out
}

stopCluster(cl)

####################
## generate plots ##
####################
setwd(paste0(pth, "Disc_runs/fadeout_runs/validation_plots/"))

fade_data = matrix(data = R[1,], nrow = length(sizes), ncol = length(pct), byrow = T)
fadeTime_data = matrix(data = R[2,], nrow = length(sizes), ncol = length(pct), byrow = T)
rownames(fade_data) <- sizes
colnames(fade_data) <- pct
rownames(fadeTime_data) <- sizes
colnames(fadeTime_data) <- pct

fade_data_m <- melt(fade_data)
fadeTime_data_m <- melt(fadeTime_data)
fadeTime_data_m$value[fadeTime_data_m$value>1e6] <- NA # remove entries where no fadeout occured

names(fade_data_m) <- c('herd size','proportion initially exposed','fadeout probability')
names(fadeTime_data_m) <- c('herd size','proportion initially exposed','fadeout time')

#file names for saving/loading data
fade_data_name = "data/fade_PlotDat.RData"
fadeTime_data_name = "data/fadeTime_PlotDat.RData"

#save data to plot directory if it does not exist
if(!exists(fade_data_name) && !exists(fadeTime_data_name)){
  save(x = fade_data_m, file = fade_data_name)
  save(x = fadeTime_data_m, file = fadeTime_data_name)
}

#load data from plot directory if not in workspace and also in directory
if(!("fade_data_m" %in% ls()) && !("fadeTime_data_m" %in% ls()) ){
  if(!exists(fade_data_name) && !exists(fadeTime_data_name)){
    load(file = fade_data_name)
    load(file = fadeTime_data_name)
  }
}

fade_data_m$`fadeout probability`[fade_data_m$`fadeout probability` == 0] <- NA


fade_prob_map <- ggplot(data = fade_data_m, aes(x=`herd size`, y=`proportion initially exposed`, fill=`fadeout probability`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fade_data_m[3]), 1)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_prob_map <- annotate_figure(fade_prob_map, top = text_grob("Fadeout probability with fixed proportion initially exposed", 
                                      color = "black", face = "bold", size = 18))
fade_prob_map

fade_time_map <- ggplot(data = fadeTime_data_m, aes(x=`herd size`, y=`proportion initially exposed`, fill=`fadeout time`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fadeTime_data_m[3]), 1)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  labs(fill='fadeout\ntime (months)') +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_time_map <- annotate_figure(fade_time_map, top = text_grob("Fadeout time with fixed proportion initially exposed", 
                                      color = "black", face = "bold", size = 18))
fade_time_map


fade_prob <- ggplot(data = fade_data_m, aes(x=`proportion initially exposed`, y=`fadeout probability`, group=`herd size`, color=`herd size`)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral")

fade_time <- ggplot(data = fadeTime_data_m, aes(x=`proportion initially exposed`, y=`fadeout time`, group=`herd size`, color=`herd size`)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral")

if(save_plots_final){
  
  jpeg(filename = paste0("results/Fadeout_maps", '-', Sys.Date(), ".jpeg"))
  print(ggarrange(fade_prob_map, fade_time_map,ncol=2, nrow=1))
  dev.off()
  
  jpeg(filename = paste0("results/FadeoutProb_maps", '-', Sys.Date(), ".jpeg"))
  print(fade_prob_map)
  dev.off()
  
  jpeg(filename = paste0("results/FadeoutTime_maps", '-', Sys.Date(), ".jpeg"))
  print(fade_time_map)
  dev.off()
  
  jpeg(filename = paste0("results/Fadeout_plots", '-', Sys.Date(), ".jpeg"))
  print(ggarrange(fade_prob, fade_time, ncol=2, nrow=1))
  dev.off()
  
}
setwd(pth)
####################

################
## Regression ##
################

#hunt fadeout: No interaction regression 
setwd(paste0(pth,'/data/'))
single.model.fadePct = as.formula( `fadeout probability` ~ `proportion initially exposed` + `herd size`)
lm.fadePct = lm(single.model.fadePct, data = as.data.frame(scale(fade_data_m)))
write.csv((lm.fadePct), file=paste0("lmOutput_FadePct_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.fadePct)
tab_model(lm.fadePct, file = paste0("reg_FadePct_", infType, Sys.Date(), ".doc"))
setwd(pth)

################

##########################################################################################
##########################################################################################

####################
## Fadeout by Num ##
####################
#takes about 8 min to run 20x21
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
years = 20
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
prop_superSpreader = 0
verbose = 0
reps = 500
lambda_factor = 1.2
sizes = seq(from = 0, to = 750, length.out = 31)[-1]
#range_vec = c(2,4,6,8,10) #num inf - fadeout
pct = seq(from = 0, to = 25, by = 1) #num inf - fadeout
infType = "seeded"
runtype = "fadeNum_"
type_of_integral = 3;
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
save_runs = T
save_plots = F
save_plots_final = T
fade_data = matrix(data = NA, nrow = length(sizes), ncol = length(pct))
fadeTime_data = matrix(data = NA, nrow = length(sizes), ncol = length(pct))
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

R <- foreach(i = 1:length(sizes), .combine='cbind', .inorder=TRUE) %dopar% {
  
  data_over_ranges = vector( mode = 'list', length = length(pct))
  fade_data_out = matrix(data = NA, nrow = 2, ncol = length(pct))
  
  for(j in 1:length(pct)){
    
    ###########
    ## Runs ##
    ##########
    
    setwd(pth)
    size=as.integer(sizes[i])
    pars <- parameter_set_wl(k = size, 
                             scenario = infType, 
                             initial_exposed = pct[j], 
                             SS_prop = 0,
                             start_quarter = seedQuarter, 
                             verbose = verbose) #set parameters dependent 
    
    
    
    #stochastic model
    initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
    population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
    disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #get lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = 200, parameters = param_vals, initial_state = initial_state, nyears = 5, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
    toc = Sys.time()
    print(toc-tic)
    lambda_out <- getLambda_vec(data = sto_out, type = 'max')
    remove(sto_out)
    
    #run with new lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
    toc = Sys.time()
    print(toc-tic)
    remove(initial_state,population_parameters,disease_parameters,param_vals)
    
    #assign data frames to list
    #then clear data from workspace
    data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
    sto_out <- as.data.frame(sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", 'fadeout', 'fadeout time') ])
    data_over_ranges[[j]] <- aggregate(data=subset(sto_out, select=-c(rep)), .~tstep, FUN = mean) #determine average trajectory
    A <- aggregate(data=subset(sto_out, select=c(rep, `fadeout time`)), .~rep, FUN = max) 
    FOT <- mean( A$`fadeout time`[A$`fadeout time`!=0] )
    if(is.nan(FOT)){FOT=1e10}
    remove(sto_out, A)
    ##########
    
    ###########
    ## Plots ##
    ###########
    if(save_plots){
      
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
      
      #class comparison plots
      sto_by_class <- split(as.data.frame(data_sto), as.data.frame(data_sto)$variable)
      det_by_class <- split(as.data.frame(data), as.data.frame(data)$variable)
      SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
      
      Nplot <- ggplot() + 
        geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
        geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[7]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Splot <- ggplot() + 
        geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[6]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Eplot <- ggplot() + 
        geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[4]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Iplot <- ggplot() + 
        geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[2]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      SSSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[5]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ESSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sE"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[3]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ISSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sI"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[1]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      setwd(paste0(pth, "Disc_runs/fadeout_runs/validation_plots/"))
      
      jpeg(filename = paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
      dev.off()
      
      jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(Nplot)
      dev.off()
      
      setwd(pth)
      remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot, det_by_class, sto_by_class, data, data_sto)
      
    }
    ###########
    
    # save fadeout results to data matrix
    p_fade <- data_over_ranges[[j]]$fadeout[length(data_over_ranges[[j]]$fadeout)] # extract mean fade value from last averaged row
    fade_data_out[1,j] <- p_fade 
    fade_data_out[2,j] <- FOT
  }
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    name = c(paste0('data_over_ranges',i))
    assign(x=name, value = data_over_ranges)
    save( list = name[1], file = paste0('Disc_runs/fadeout_runs/FOnum_dat_',  size, '-', Sys.Date(), '.RData'))
  }
  
  ##################
  fade_data_out
}

stopCluster(cl)

####################
## generate plots ##
####################
setwd(paste0(pth, "Disc_runs/fadeout_runs/validation_plots/"))

fade_data = matrix(data = R[1,], nrow = length(sizes), ncol = length(pct), byrow = T)
fadeTime_data = matrix(data = R[2,], nrow = length(sizes), ncol = length(pct), byrow = T)
rownames(fade_data) <- sizes
colnames(fade_data) <- pct
rownames(fadeTime_data) <- sizes
colnames(fadeTime_data) <- pct

fade_Ndata_m <- melt(fade_data)
fadeTime_Ndata_m <- melt(fadeTime_data)
fadeTime_Ndata_m$value[fadeTime_Ndata_m$value>1e6] <- NA # remove entries where no fadeout occured

names(fade_Ndata_m) <- c('herd size','number initially exposed','fadeout probability')
names(fadeTime_Ndata_m) <- c('herd size','number initially exposed','fadeout time')

#names for saving/loading data
fade_data_name = "fade_numPlotDat.RData"
fadeTime_data_name = "fadeTime_numPlotDat.RData"

#save data to plot directory if it does not exist
if(!exists(fade_data_name) && !exists(fadeTime_data_name)){
  save(fade_Ndata_m, file = fade_data_name)
  save(fadeTime_Ndata_m, file = fadeTime_data_name)
}

#load data from plot directory if not in workspace and also in directory
if(!("fade_Ndata_m" %in% ls()) && !("fadeTime_Ndata_m" %in% ls()) ){
  if(!exists(fade_data_name) && !exists(fadeTime_data_name)){
    load(file = fade_data_name)
    load(file = fadeTime_data_name)
  }
}

fade_Ndata_m$`fadeout probability`[fade_Ndata_m$`fadeout probability` == 0] <- NA

fade_prob_map <- ggplot(data = fade_Ndata_m, aes(x=`herd size`, y=`number initially exposed`, fill=`fadeout probability`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fade_Ndata_m[3]), 1)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_prob_map <- annotate_figure(fade_prob_map, top = text_grob("Fadeout probability with fixed number initially exposed", 
                                                                color = "black", face = "bold", size = 18))
fade_prob_map

fade_time_map <- ggplot(data = fadeTime_Ndata_m, aes(x=`herd size`, y=`number initially exposed`, fill=`fadeout time`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fadeTime_Ndata_m[3]), 1)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  labs(fill='fadeout\ntime (months)') +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_time_map <- annotate_figure(fade_time_map, top = text_grob("Fadeout time with fixed number initially exposed", 
                                                                color = "black", face = "bold", size = 18))
fade_time_map

fade_prob <- ggplot(data = fade_Ndata_m, aes(x=`number initially exposed`, y=`fadeout probability`, group=`herd size`, color=`herd size`)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral")

fade_time <- ggplot(data = fadeTime_Ndata_m, aes(x=`number initially exposed`, y=`fadeout time`, group=`herd size`, color=`herd size`)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral")


if(save_plots_final){
  jpeg(filename = paste0("Fadeout_maps", '_num-', Sys.Date(), ".jpeg"))
  print(ggarrange(fade_prob_map, fade_time_map, ncol=2, nrow=1))
  dev.off()
  
  jpeg(filename = paste0("FadeoutProb_maps", '_num-', Sys.Date(), ".jpeg"))
  print(fade_prob_map)
  dev.off()
  
  jpeg(filename = paste0("FadeoutTime_maps", '_num-', Sys.Date(), ".jpeg"))
  print(fade_time_map)
  dev.off()
  
  jpeg(filename = paste0("Fadeout_plots", '_num-', Sys.Date(), ".jpeg"))
  print(ggarrange(fade_prob, fade_time, ncol=2, nrow=1))
  dev.off()
}
setwd(pth)
####################

################
## Regression ##
################

setwd(paste0(pth,'Disc_runs/fadeout_runs/regression/'))
single.model.fadeNum = as.formula( `fadeout probability` ~ `number initially infected` + `herd size`)
lm.fadeNum = lm(single.model.fadeNum, data = as.data.frame(scale(fade_Ndata_m)))
write.csv(lm.fadeNum, file=paste0("lmOutput_FadeNum_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.fadeNum)
tab_model(lm.fadeNum, file = paste0("reg_FadeNum_", infType, Sys.Date(), ".doc"))
setwd(pth)

################

##########################################################################################
##########################################################################################

####################
## super susceptible ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_DTMC.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_CTMC.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")
years = 1
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
# prop_superSpreader = 0.05
pct = 0
verbose = 0
lambda_factor = 1.2
reps = 500
sizes = seq(from = 0, to = 750, length.out = 31)[-1]
range_vec = seq(from = 0, to = 1, length.out = 21) # prop SS - epi. extent
infType = "spillover"
runtype = "SS_pub_"
type_of_integral = 3;
#name_out = paste0(runtype, infType,"_", pct,'-', range_vec[j], '-it', type_of_integral, '-')
#lambda_name_out = paste0('stdLambda_', infType,"_", pct,'-', prop_superSpreader, '-')
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files

save_runs = T
save_plots = F
save_plots_final = T
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

R <- foreach(i = 1:length(sizes), .combine='cbind', .inorder=TRUE) %dopar% {
  size=sizes[i]
  data_over_ranges = vector( mode = 'list', length = length(range_vec))
  ss_stats_out = matrix(data = NA, nrow = 4, ncol = length(range_vec))
  rownames(ss_stats_out) <- c('min', 'median', 'mean', 'max')
  nm <- c("tstep", "pct_inf", "num_inf", "pct_SS")
  ss_data_out = data.frame( matrix(data = NA, nrow=0, ncol=length(nm)) )
  names(ss_data_out) <- nm
  
  for(j in 1:length(range_vec)){
    name_out = paste0(runtype, infType,"_", pct,'-', range_vec[j], '-it', type_of_integral, '-')
    ###########
    ## Runs ##
    ##########
    
    setwd(pth)
    size=as.integer(sizes[i])
    pars <- parameter_set_wl(k = size, 
                             scenario = infType, 
                             initial_exposed = pct*size, 
                             SS_prop = range_vec[j],
                             start_quarter = seedQuarter, 
                             verbose = verbose) #set parameters dependent 

    #stochastic model
    initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
    population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
    disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #get lambda vals
    #get lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 3, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
    toc = Sys.time()
    print(toc-tic)
    lambda_out <- getLambda_vec(data = sto_out, type = 'max')
    remove(sto_out)
    
    #run model with lambdas
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
    toc = Sys.time()
    print(toc-tic)
    remove(initial_state,population_parameters,disease_parameters,param_vals)
    
    #assign data frames to list
    #then clear data from workspace
    data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
    sto_out <- as.data.frame(sto_out[,c( "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", 'fadeout') ])
    sto_out <- aggregate(data=sto_out, .~tstep, FUN = mean)
    data_over_ranges[[j]] <- sto_out
    pct_prev <- sto_out[,"Total Infected"]/sto_out[,'N']
    ss_data_out <- rbind(ss_data_out , cbind(sto_out[,"tstep"], sto_out[,"Total Infected"], pct_prev, range_vec[j]))
    ss_stats_out[,j] <- c(min(pct_prev), median(pct_prev), mean(pct_prev), max(pct_prev))
    remove(sto_out)
    ##########
    
    ###########
    ## Plots ##
    ###########
    if(save_plots){
      
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
      
      #class comparison plots
      sto_by_class <- split(as.data.frame(data_sto), as.data.frame(data_sto)$variable)
      det_by_class <- split(as.data.frame(data), as.data.frame(data)$variable)
      SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
      
      Nplot <- ggplot() + 
        geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
        geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[7]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,1), limits = c(0, years*12))
      
      Splot <- ggplot() + 
        geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[6]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Eplot <- ggplot() + 
        geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[4]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Iplot <- ggplot() + 
        geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[2]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      SSSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[5]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ESSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sE"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[3]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ISSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sI"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[1]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
    
      setwd(paste0(pth, "Disc_runs/ss_runs/validation_plots/"))
      
      jpeg(filename = paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
      dev.off()
      
      jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(Nplot)
      dev.off()
      
      setwd(pth)
      remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot, det_by_class, sto_by_class, data, data_sto)
      
    }
    ###########
  }
    
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    #name = c(paste0('data_over_ranges',i))
    #assign(x=name, value = data_over_ranges)
    #save( list = name[1], file = paste0('Disc_runs/ss_runs/ss_dat',  size, '-', Sys.Date(), '.RData'))
    
    save( data_over_ranges, file = paste0('Disc_runs/ss_runs/ss_dat',  size, '-', Sys.Date(), '.RData'))
  }
  names(ss_data_out) <- nm
  return(list(ss_data_out, ss_stats_out))
  
  ##################
}

stopCluster(cl)

#####################
## Post Processing ##
#####################
# reassemble data because foreach is annoying,lol
st_nm <- c('min', 'median', 'mean', 'max')
ss_stat_full <- data.frame( matrix(data = NA, nrow=0, ncol=length(st_nm)) )
names(ss_stat_full) <- st_nm #dataframe for stat metrics

nm <- c("tstep", "number infected", "percent infected", "percent super susceptible", "size")
ss_data_full <- data.frame( matrix(data = NA, nrow=0, ncol=length(nm)) )
names(ss_data_full) <- nm #dataframe for epidemic extent data

setwd(paste0(pth, "Disc_runs/ss_runs/validation_plots/"))

for(i in 1:(length(R)/2)){
  ss_data_full <- rbind( ss_data_full, cbind(as.data.frame(R[[2*i-1]]), sizes[i]) ) #epi extent data
  ss_stat_full <- rbind(ss_stat_full,as.data.frame(t(R[[2*i]]))) #stat metrics - transpose b/c some idiot made the output as row vectors - someone should really fix that :)
}
names(ss_data_full) <- nm

#names of output data files for saving/reloading data
ss_data_name <- 'ss_plotData.RData'
ss_stat_name <- 'ss_plotStat.RData'
  
#save data to plot directory if it does not exist
if(!exists(ss_data_name) && !exists(ss_stat_name)){
  save(ss_data_full, file = ss_data_name)
  save(ss_stat_full, file = ss_stat_name)
}

#load data from plot directory if not in workspace and also in directory
if(!("ss_data_full" %in% ls()) && !("ss_stat_full" %in% ls()) ){
  if(!exists(ss_data_name) && !exists(ss_stat_name)){
    load(file = ss_data_name)
    load(file = ss_stat_name)
  }
}
#####################

################
# metric plots #
################
min_mat <- matrix(data = ss_stat_full$min, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(min_mat) <- sizes
colnames(min_mat) <- range_vec

median_mat <- matrix(data = ss_stat_full$median, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(median_mat) <- sizes
colnames(median_mat) <- range_vec

mean_mat <- matrix(data = ss_stat_full$mean, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(mean_mat) <- sizes
colnames(mean_mat) <- range_vec

max_mat <- matrix(data = ss_stat_full$max, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(max_mat) <- sizes
colnames(max_mat) <- range_vec

melty_min_mat <- melt(min_mat)
melty_median_mat <- melt(median_mat)
melty_mean_mat <- melt(mean_mat)
melty_max_mat <- melt(max_mat)

names(melty_min_mat) <- c('herd size','proportion of super susceptibles','min. infection prevalence')
names(melty_median_mat) <- c('herd size','proportion of super susceptibles','median infection prevalence')
names(melty_mean_mat) <- c('herd size','proportion of super susceptibles','mean infection prevalence')
names(melty_max_mat) <- c('herd size','proportion of super susceptibles','max. infection prevalence')

min_map <- ggplot(data = melty_min_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`min. infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_median_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral")

median_map <- ggplot(data = melty_median_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`median infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_median_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral")
median_map <- annotate_figure(median_map, top = text_grob("Median infection prevalence\nwith constant force of infection", 
                                                                color = "black", face = "bold", size = 18))
median_map

mean_map <- ggplot(data = melty_mean_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`mean infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_mean_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral")
mean_map <- annotate_figure(mean_map, top = text_grob("Mean infection prevalence with constant force of infection", 
                                                                color = "black", face = "bold", size = 18))
mean_map

max_map <- ggplot(data = melty_max_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`max. infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_median_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral")
################

################
## Regression ##
################

#ss proportion: No interaction regression 
setwd(paste0(pth,'Disc_runs/ss_runs/regression/'))
single.model.ss_prev = as.formula( `mean inf. prevalence` ~ `proportion of super susceptibles` + `herd size`)
lm.ss_prev = lm(single.model.ss_prev, data = as.data.frame(scale(melty_mean_mat)))
write.csv(lm.ss_prev, file=paste0("lmOutput_mean_ssPrev_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.ss_prev)
tab_model(lm.ss_prev, file = paste0("reg_ssPrev_", infType, Sys.Date(), ".doc"))
setwd(pth)
################

######################
# epi. extent plots #
######################
#select some sizes to view in plots
setwd(paste0(pth, "Disc_runs/ss_runs/validation_plots/"))

sizes_select <- sizes[c(2, 5, 8, 12, 15, 18)]
pct_select <- range_vec[c(1,5,9,14,13,17,21)]
ss_data_full <- ss_data_full[ss_data_full$size %in% sizes_select,]
ss_data_full <- ss_data_full[ss_data_full$`percent super susceptible` %in% pct_select,]

p1 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[1],], 
             aes(x=tstep, y=`percent infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
p2 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[2],], 
             aes(x=tstep, y=`percent infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
p3 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[3],], 
             aes(x=tstep, y=`percent infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
p4 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[4],], 
             aes(x=tstep, y=`percent infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
p5 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[5],], 
             aes(x=tstep, y=`percent infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
p6 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[6],], 
             aes(x=tstep, y=`percent infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()


n1 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[1],], 
             aes(x=tstep, y=`number infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
n2 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[2],], 
             aes(x=tstep, y=`number infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
n3 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[3],], 
             aes(x=tstep, y=`number infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
n4 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[4],], 
             aes(x=tstep, y=`number infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
n5 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[5],], 
             aes(x=tstep, y=`number infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()
n6 <- ggplot(data = ss_data_full[ss_data_full$size == sizes_select[6],], 
             aes(x=tstep, y=`number infected`, group = `percent super susceptible`, color = `percent super susceptible`)) +
  geom_line() +
  geom_point()



######################
if(save_plots_final){
  jpeg(filename = paste0("SS_statMaps", '_rng-', Sys.Date(), ".jpeg"))
  print(ggarrange(min_map, max_map, ncol=1, nrow=2))
  dev.off()
  
  jpeg(filename = paste0("SS_statMaps", '-', Sys.Date(), ".jpeg"))
  print(ggarrange(median_map, mean_map, ncol=1, nrow=2))
  dev.off()
  
  jpeg(filename = paste0("SS_meanMaps", '-', Sys.Date(), ".jpeg"))
  print(mean_map)
  dev.off()
  
  jpeg(filename = paste0("SS_medianMaps", '_mid-', Sys.Date(), ".jpeg"))
  print(median_map)
  dev.off()
  
  jpeg(filename = paste0("Fadeout_maps", '_pct-', Sys.Date(), ".jpeg"))
  print(ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3))
  dev.off()
  
  jpeg(filename = paste0("SSepi_plots", '_num-', Sys.Date(), ".jpeg"))
  print(ggarrange(n1, n2, n3, n4, n5, n6, ncol=2, nrow=3))
  dev.off()
}
setwd(pth)

####################
## Hunting  Rates ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
years = 20
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
prop_superSpreader = 0.00
pct = .025 # explore low medium and high prevalence c(.025, .1, .25)
case = "low"
verbose = 0
reps = 300
sizes = seq(from = 0, to = 750, length.out = 31)[-1]
range_vec = seq(from = 0, to = .5, length.out = 21) # prop SS - epi. extent
#range_vec = c(2,4,6,8,10) #num inf - fadeout
#range_vec = c(0, .02, .05, .1, .25, .5) #pct inf - fadeout
#range_vec = c(0, .1, .25, .5, .75, .9, 1) #harvest prop -- prev. estimate, fadeout
infType = "seeded"
runtype = "hunt_lowPrev_"
type_of_integral = 3;
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
lambda_factor = 1.2
save_runs = T
save_plots = F
save_plots_final = T
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

R <- foreach(i = 1:length(sizes), .combine='cbind', .inorder=TRUE) %dopar% {
  
  nm <- c("tstep", "mean", "median")
  nm_stat <- c( "mean", "median")
  nm_fade <- c("FO")
  
  fade_data_out = matrix(data = NA, nrow = 1, ncol = length(range_vec)) # fadeout data vector
  hunt_data_out = data.frame( matrix(data = NA, nrow=0, ncol=length(nm)) ) # aggregate simulation data with select cols
  prevEst_out = matrix(data = NA, nrow = 2, ncol = length(range_vec)) # prev estimate (metrics) vectors
  rownames(prevEst_out) <- nm_stat
  names(fade_data_out) <- nm_fade
  names(hunt_data_out) <- nm
  
  data_over_ranges = vector( mode = 'list', length = length(range_vec))
  
  for(j in 1:length(range_vec)){
    
    ###########
    ## Runs ##
    ##########
    
    setwd(pth)
    size=as.integer(sizes[i])
    pars <- parameter_set_wl(k = size, 
                             scenario = infType, 
                             initial_exposed = pct*size, 
                             SS_prop = 0,
                             start_quarter = seedQuarter, 
                             verbose = verbose) #set parameters dependent 
    
    
    
    #stochastic model
    initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
    population_parameters<-data.frame(K=pars["K"], eta_hunt=range_vec[j], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
    disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #get lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 3, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
    toc = Sys.time()
    print(toc-tic)
    lambda_out <- getLambda_vec(data = sto_out, type = 'max')
    remove(sto_out)
    
    #run with new lambda
    tic = Sys.time()
    sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
    toc = Sys.time()
    print(toc-tic)
    remove(initial_state,population_parameters,disease_parameters,param_vals)
    
    #assign data frames to list
    #then clear data from workspace
    data_sto <- sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
    sto_out <- as.data.frame(sto_out[,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", 'Total Hunt', "Total Infected", 'quarter', 'Infected Hunted', 'fadeout', 'fadeout time', 'Hunt Prevalence') ])
    data_over_ranges[[j]] <- aggregate(data=subset(sto_out, select=-c(rep)), .~tstep, FUN = mean) #determine average trajectory
    p_fade <- data_over_ranges[[j]]$fadeout[length(data_over_ranges[[j]]$fadeout)] # extract mean fade value from last averaged row
    
    sto_out <- as.data.frame(sto_out[,c("rep", "tstep", "N", 'Total Hunt', "Total Infected", 'quarter', 'Infected Hunted', 'fadeout', 'fadeout time', 'Hunt Prevalence') ])
    sto_out <- sto_out[sto_out$quarter==4,] # extract only data with hunting
    prevEst_mean <- mean(sto_out$`Hunt Prevalence`)
    prevEst_median <- median(sto_out$`Hunt Prevalence`)
    
    # calculate aggregate metrics
    med_Estim <- aggregate(data=subset(sto_out, select=c("tstep", 'Hunt Prevalence')), .~tstep, FUN = median)
    mean_Estim <- aggregate(data=subset(sto_out, select=c("tstep", 'Hunt Prevalence')), .~tstep, FUN = mean)
    
    remove(sto_out)
    ##########
    
    ###########
    ## Plots ##
    ###########
    if(save_plots){
      
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
      
      #class comparison plots
      sto_by_class <- split(as.data.frame(data_sto), as.data.frame(data_sto)$variable)
      det_by_class <- split(as.data.frame(data), as.data.frame(data)$variable)
      SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
      
      Nplot <- ggplot() + 
        geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
        geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[7]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Splot <- ggplot() + 
        geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[6]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Eplot <- ggplot() + 
        geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[4]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      Iplot <- ggplot() + 
        geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7) + 
        geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[2]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      SSSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[5]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ESSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sE"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[3]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      ISSplot <- ggplot() + 
        geom_line(data = sto_by_class[["SS_I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6) + 
        geom_line(data = det_by_class[["sI"]], aes(x = time, y=value), size = 1) +
        scale_color_manual(values=SEIcols[1]) +
        scale_y_continuous(limits = c(0,NA)) +
        scale_x_continuous(breaks=seq(0,years*12,12), limits = c(0, years*12))
      
      setwd(paste0(pth, "Disc_runs/fadeout_runs/validation_plots/"))
      
      jpeg(filename = paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
      dev.off()
      
      jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
      print(Nplot)
      dev.off()
      
      setwd(pth)
      remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot, det_by_class, sto_by_class, data, data_sto)
      
    }
    ###########
    
    # save fadeout results to data matrix
    
    fade_data_out[1,j] <- p_fade 
    hunt_data_out <- rbind( hunt_data_out, data.frame(tstep=med_Estim$tstep, `harvest rate` = range_vec[j], mean=mean_Estim$`Hunt Prevalence`, median=med_Estim$`Hunt Prevalence`) )
    prevEst_out[1,j] <- prevEst_mean
    prevEst_out[2,j] <- prevEst_median
  }
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    save( data_over_ranges, file = paste0('Disc_runs/hunt_runs/hunt_dat_',  size, '_', case, '-', Sys.Date(), '.RData'))
  }
  
  ##################
  return(list(fade_data_out, hunt_data_out, prevEst_out))
}

stopCluster(cl)

#####################
## Post Processing ##
#####################
# reassemble data because foreach is annoying,lol
nm <- c("tstep", "percent harvested", "mean", "median", "size")
nm_stat <- c( "mean", "median")
nm_fade <- c("FO")

fade_data_full = data.frame( matrix(data = NA, nrow=0, ncol=1) ) # fadeout data vector
hunt_data_full = data.frame( matrix(data = NA, nrow=0, ncol=length(nm)) ) # aggregate simulation data with select cols
prevEst_full = data.frame( matrix(data = NA, nrow = 0, ncol = 2) ) # prev estimate (metrics) vectors


setwd(paste0(pth, "Disc_runs/hunt_runs/validation_plots/"))

for(i in 1:(length(R)/3)){
  
  fade_data_full <- rbind( fade_data_full, as.data.frame( R[[3*i-2]] ) ) #
  hunt_data_full <- rbind( hunt_data_full, cbind( as.data.frame(R[[3*i-1]]), sizes[i] ) ) #
  prevEst_full <- rbind( prevEst_full, t( as.data.frame(R[[3*i]]) ) ) # 
  
}

names(prevEst_full) <- nm_stat
names(hunt_data_full) <- nm
fade_data_full <- as.matrix(fade_data_full)
rownames(fade_data_full) <- sizes
colnames(fade_data_full) <- range_vec


hunt_data_name <- paste0('hunt_plotData_', case, '.RData')
hunt_stat_name <- paste0('hunt_plotStat_', case, '.RData')
hunt_fade_name <- paste0('hunt_plotFade_', case, '.RData')

#save data to plot directory if it does not exist
if(!exists(hunt_data_name) && !exists(hunt_stat_name) && !exists(hunt_fade_name)){
  save(fade_data_full, file = hunt_fade_name)
  save(hunt_data_full, file = hunt_data_name)
  save(prevEst_full, file = hunt_stat_name)
}

#load data from plot directory if not in workspace and also in directory
if(!("fade_data_full" %in% ls()) && !("hunt_data_full" %in% ls()) && !("prevEst_full" %in% ls())){
  if(!exists(hunt_data_name) && !exists(hunt_stat_name) && !exists(hunt_fade_name)){
    load(file = hunt_fade_name)
    load(file = hunt_data_name)
    load(file = hunt_stat_name)
  }
}

mean_mat <- matrix(data = prevEst_full$mean, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(mean_mat) <- sizes
colnames(mean_mat) <- range_vec

median_mat <- matrix(data = prevEst_full$median, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(median_mat) <- sizes
colnames(median_mat) <- range_vec
#####################

######################
## Fadeout heatmaps ##
######################
fade_data_melt <- melt(fade_data_full)
names(fade_data_melt) <- c('herd size', 'percent harvested', 'fadeout probability')
fade_data_melt$`fadeout probability`[fade_data_melt$`fadeout probability` == 0] <- NA

fade_prob_map <- ggplot(data = fade_data_melt, aes(x=`herd size`, y=`percent harvested`, fill=`fadeout probability`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fade_data_melt[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_prob_map <- annotate_figure(fade_prob_map, top = text_grob("Fadeout probability across different hunter harvest rates", 
                                                          color = "black", face = "bold", size = 18))
fade_prob_map

######################

################
## Regression ##
################

#hunt fadeout: No interaction regression 
setwd(paste0(pth,'Disc_runs/hunt_runs/regression/'))
single.model.hunt_fade = as.formula( `fadeout probability` ~ `percent harvested` + `herd size`)
lm.hunt_fade = lm(single.model.hunt_fade, data = as.data.frame(scale(fade_data_melt)))
write.csv(lm.hunt_fade, file=paste0("lmOutput_huntFade_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.hunt_fade)
tab_model(lm.hunt_fade, file = paste0("reg_huntFade_", infType, Sys.Date(), ".doc"))
setwd(pth)
################

#####################
## Prev. Est. maps ##
#####################
median_mat_melt <- melt(median_mat)
names(median_mat_melt) <- c('herd size', 'percent harvested', 'median estimated prevalence ratio')
mean_mat_melt <- melt(mean_mat)
names(mean_mat_melt) <- c('herd size', 'percent harvested', 'mean estimated prevalence ratio')

median_map <- ggplot(data = median_mat_melt, aes(x=`herd size`, y=`percent harvested`, fill=`median estimated prevalence ratio`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(median_mat_melt[3]), 18)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral")
median_map <- annotate_figure(median_map, top = text_grob("Median estimated prevalence ratio", 
                                                      color = "black", face = "bold", size = 18))
median_map

mean_map <- ggplot(data = mean_mat_melt, aes(x=`herd size`, y=`percent harvested`, fill=`mean estimated prevalence ratio`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(mean_mat_melt[3]), 18)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral")
mean_map <- annotate_figure(mean_map, top = text_grob("Mean estimated prevalence ratio", 
                                                      color = "black", face = "bold", size = 18))
mean_map

#####################

################
## Regression ##
################
#mean hunt ratio: No interaction regression 
setwd(paste0(pth,'Disc_runs/hunt_runs/regression/'))
single.model.mean_est = as.formula( `mean estimated prevalence` ~ `percent harvested` + `herd size`)
lm.mean_est = lm(single.model.mean_est, data = as.data.frame(scale(mean_mat_melt)))
write.csv(lm.mean_est, file=paste0("lmOutputMeanEst_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.mean_est)
tab_model(lm.mean_est, file = paste0("reg_huntMean_", infType, Sys.Date(), ".doc"))
setwd(pth)
################

##########################
## Prev. Est. over time ##
##########################
names(hunt_data_full) <- c("tstep", "percent harvested", "mean observed prevalence", "median observed prevalence", "size")
sizes_select <- sizes[c(2, 5, 8, 12, 15, 18)]
harvest_select <- range_vec[c(1,5,9,14,13,17,21)]

hunt_data_full <- hunt_data_full[hunt_data_full$size %in% sizes_select,]
hunt_data_full <- hunt_data_full[hunt_data_full$`percent harvested` %in% pct_select,]

p1 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[1],], 
             aes(x=tstep, y=`mean observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
p2 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[2],], 
             aes(x=tstep, y=`mean observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
p3 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[3],], 
             aes(x=tstep, y=`mean observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
p4 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[4],], 
             aes(x=tstep, y=`mean observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
p5 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[5],], 
             aes(x=tstep, y=`mean observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
p6 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[6],], 
             aes(x=tstep, y=`mean observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()


n1 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[1],], 
             aes(x=tstep, y=`median observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
n2 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[2],], 
             aes(x=tstep, y=`median observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
n3 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[3],], 
             aes(x=tstep, y=`median observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
n4 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[4],], 
             aes(x=tstep, y=`median observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
n5 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[5],], 
             aes(x=tstep, y=`median observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()
n6 <- ggplot(data = hunt_data_full[hunt_data_full$size == sizes_select[6],], 
             aes(x=tstep, y=`median observed prevalence`, group = `percent harvested`, color = `percent harvested`)) +
  geom_line() +
  geom_point()



##########################
setwd(paste0(pth, "Disc_runs/hunt_runs/validation_plots/"))
if(save_plots_final){
  jpeg(filename = paste0("hunt_fadeMaps", '_', case,'-', Sys.Date(), ".jpeg"))
  print(fade_prob_map)
  dev.off()
  
  jpeg(filename = paste0("hunt_statMaps", '_', case,'-', Sys.Date(), ".jpeg"))
  print(ggarrange(median_map, mean_map, ncol=1, nrow=2))
  dev.off()
  
  jpeg(filename = paste0("hunt_medianMaps", '_', case,'-', Sys.Date(), ".jpeg"))
  print(median_map)
  dev.off()
  
  jpeg(filename = paste0("hunt_meanMaps", '_', case,'-', Sys.Date(), ".jpeg"))
  print(mean_map)
  dev.off()
  
  jpeg(filename = paste0("hunt_prevMean", '_', case,'-', Sys.Date(), ".jpeg"))
  print(ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3))
  dev.off()
  
  jpeg(filename = paste0("hunt_prevMedian", '_', case,'-', Sys.Date(), ".jpeg"))
  print(ggarrange(n1, n2, n3, n4, n5, n6, ncol=2, nrow=3))
  dev.off()
}
setwd(pth)
##########################


##########################################################################################
##########################################################################################
#######################            Sensitivity Analysis             ######################
##########################################################################################
##########################################################################################

######################
## Setup Parameters ##
######################
rm(list = ls())
library(rstudioapi);
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
N_LHS_sets = 1000
size_range = c(50,1000)
years = 4
prop_superSpreader = 0.1
pct = .02
verbose = 0
reps = 500
def_seed = 43
infection = "seeded"
infType = "seeded_q1"
fix_q = 1 # set to NA or desired quarter
runtype = "LHS_"
type_of_integral = 3;
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
par_file = paste0('bTBwl_LHSpars_',infType)
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
lambda_factor = 1.2
save_pars = T
save_runs = T
save_plots = F
save_plots_final = T
recalc = F
averaged = T
##parameters excluded from analysis: 
# omega - birth pulse timing, verbose, integral type
######################

###############################
## Generate or load LHS pars ##
###############################
# see if lh is in workspace or sensitivity_analysis folder; otherwise, generate new LHS set
if(recalc){
  print('calculating parameter sets')
  lh <- LHS_ParSets(k_lim = size_range, 
                    scenario = infection,
                    infected = pct,
                    pct = T,
                    SS_prop = prop_superSpreader,
                    verbose = 0,
                    Num_LHS_sets = N_LHS_sets,
                    file.path = '/home/webblab/Documents/Brandon/bTB_wildlife_code/sensitivity_analysis/parameters/',
                    file.out = save_pars,
                    file.name = par_file,
                    seed_q = fix_q,
                    seed = def_seed)
}
if( ("lh" %in% ls()) ){
  print('found in workspace')
}else if( paste0('bTBwl_LHSpars_', infType, '_summary.csv') %in% list.files("./sensitivity_analysis/parameters/") ){
  print('loading from file')
  lh <- read.csv(file = paste0('./sensitivity_analysis/parameters/bTBwl_LHSpars_', infType, '.csv'))
}else{
  print('recalculating parameter sets')
  lh <- LHS_ParSets(k_lim = size_range, 
                  scenario = infection,
                  infected = pct,
                  pct = T,
                  SS_prop = prop_superSpreader,
                  verbose = 0,
                  Num_LHS_sets = N_LHS_sets,
                  file.path = '/home/webblab/Documents/Brandon/bTB_wildlife_code/sensitivity_analysis/parameters/',
                  file.out = save_pars,
                  file.name = par_file,
                  seed_q = fix_q,
                  seed = def_seed)
}
###############################

#########################################
## run model and generate summary file ##
#########################################
# run with LHS parameters and generate summary file

if(def_seed == 43 & N_LHS_sets == 1000 & dim(lh)[1] == N_LHS_sets){
  lh <- lh[-c(363,367),] #remove these 2 entries due to excessive runtime issues -- getting massive vals for lambda O(10^15)
}

TIC <- Sys.time()

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


R <- foreach(i = 1:dim(lh)[1], .combine='rbind', .inorder=TRUE) %dopar% {
  pars <- lh[i,]
  
  initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=0) 
  disease_parameters<-data.frame(beta=pars["beta"], area=1, p1=1, p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  
  param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
  #get lambda
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
  toc = Sys.time()
  print(toc-tic)
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  remove(sto_out)
  
  #run with new lambda
  # verbose of -19 used to take minimal output necessary for sensitivity runs
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = -19, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
  toc = Sys.time()
  print(toc-tic)
  remove(initial_state,population_parameters,disease_parameters,param_vals)
  if(save_runs){
    write.csv(sto_out, file = paste0('./sensitivity_analysis/sens_runs/sens_run_',infType, '_', i, '.csv'))
  }
  
  # mean hunt prev
  prev_est = aggregate(data=subset(sto_out[sto_out$quarter==4,], select=c("rep", 'Hunt Prevalence')), .~rep, FUN = mean)
  prev_est <- prev_est$`Hunt Prevalence`
  #truncate to last row only
  sto_out <- sto_out[sto_out$tstep == max(sto_out$tstep),]
  sto_out$`Hunt Prevalence` <- prev_est
  sto_out <- sto_out[,c('N','Total Infected','fadeout','fadeout time','Hunt Prevalence')]
  remove(prev_est)
  #merge with pars
  #workaround to save memory
  pars <- pars[rep(seq_len(nrow(pars)), reps),]
  #pars<-replicate(reps, i)
  sto_out <- cbind(sto_out, pars)
  remove(pars)
  return(sto_out)
  remove(sto_out)
}

stopCluster(cl)

TOC <- Sys.time()
TOC-TIC

c('K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 'alpha_max', 'ksi', 
  'omega', 's', 'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi', 
  'sigma1_mean', 'sigma1_rate', 'start_q')

#R <- cbind(R, lh[rep(seq_len(nrow(lh)), each = reps), ])
# scale parameter values
bTB_wl_scaled <- as.data.frame(R %>% mutate_at(vars('K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                                        'alpha_max', 'ksi', 'omega', 's', 
                                        'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 
                                        'phi', 'sigma1_mean', 'sigma1_rate', 'start_q'),
                                   list(~scale(as.vector(.)))))

# save summary file and rescaled summary
write.csv(R, file = paste0("./sensitivity_analysis/summary_files/LHS_summary_", infType, ".csv" ), row.names = F)
remove(R)
write.csv(bTB_wl_scaled, file = paste0("./sensitivity_analysis/summary_files/LHS_scaled_summary_", infType, ".csv" ), row.names = F)
#########################################

#####################
## Monotonic Plots ##
#####################
# load in data if not in workspace
# define infType variable and averaged 
# infType = 'seeded_q1'; remove(bTB_wl_scaled)
averaged = T
setwd(pth)

if( ("bTB_wl_scaled" %in% ls()) ){
  print('found in workspace')
}else if( paste0('LHS_scaled_summary_', infType, '.csv') %in% list.files("./sensitivity_analysis/summary_files/") ){
  print('loading from file')
  bTB_wl_scaled <- read.csv(file = paste0('./sensitivity_analysis/summary_files/LHS_scaled_summary_', infType, '.csv'))
}
mono_plots = F
bTB_wl_scaled <- bTB_wl_scaled[,!(names(bTB_wl_scaled) %in% c('N', "pars", "S_0", "E1_0", "SuperS_0", "SuperE1_0", "I_0", "SuperI_0"))]
names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max', 'ksi', 'omega', 's', 
                          'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi', 'sigma1_mean', 'sigma1_rate', 
                          'start_q')

if(averaged){
  bTB_wl_agg <- aggregate.data.frame(x = bTB_wl_scaled, by = list(rep(1:dim(lh)[1], times=1, each=reps)), FUN = mean)
  bTB_wl_agg <- bTB_wl_agg[,names(bTB_wl_agg) %in% names(bTB_wl_scaled)]
  bTB_wl_scaled <- bTB_wl_agg
  infType <- paste0(infType, '_avg')
  write.csv(bTB_wl_scaled, file = paste0("./sensitivity_analysis/summary_files/LHS_scaled_summary_", infType, ".csv" ), row.names = F)
}
# dont expand - just run it...
#####################

setwd('./sensitivity_analysis/monotonic_plots/')
if(mono_plots){
jpeg(filename = paste0("monotonic_K_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)

p1=ggplot(bTB_wl_scaled, aes(x = K, y = `Total Infected`)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3))

p2=ggplot(bTB_wl_scaled, aes(x = K, y = fadeout)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3))

p3=ggplot(bTB_wl_scaled, aes(x = K, y = `fadeout time`)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3))

p4=ggplot(bTB_wl_scaled, aes(x = K, y = `Hunt Prevalence`)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3))

print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))

dev.off()
} #K
if(mono_plots){
  jpeg(filename = paste0("monotonic_eta_hunt_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = eta_hunt, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = eta_hunt, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = eta_hunt, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = eta_hunt, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #eta_hunt
if(mono_plots){
  jpeg(filename = paste0("monotonic_eta_nat_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = eta_nat, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = eta_nat, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = eta_nat, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = eta_nat, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #eta_nat
if(mono_plots){
  jpeg(filename = paste0("monotonic_theta_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = theta, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = theta, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = theta, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = theta, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #theta
if(mono_plots){
  jpeg(filename = paste0("monotonic_gamma_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = gamma, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = gamma, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = gamma, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = gamma, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #gamma
if(mono_plots){
  jpeg(filename = paste0("monotonic_alpha_max_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = alpha_max, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = alpha_max, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = alpha_max, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = alpha_max, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #alpha_max
if(mono_plots){
  jpeg(filename = paste0("monotonic_ksi_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = ksi, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = ksi, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = ksi, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = ksi, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #ksi
if(mono_plots){
  jpeg(filename = paste0("monotonic_omega_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = omega, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = omega, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = omega, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = omega, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #omega
if(mono_plots){
  jpeg(filename = paste0("monotonic_s_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = s, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = s, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = s, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = s, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #s
if(mono_plots){
  jpeg(filename = paste0("monotonic_beta_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = beta, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = beta, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = beta, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = beta, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #beta
if(mono_plots){
  jpeg(filename = paste0("monotonic_p2_q1_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = p2_q1, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = p2_q1, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = p2_q1, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = p2_q1, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #p2_q1
if(mono_plots){
  jpeg(filename = paste0("monotonic_p2_q2_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = p2_q2, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = p2_q2, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = p2_q2, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = p2_q2, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #p2_q2
if(mono_plots){
  jpeg(filename = paste0("monotonic_p2_q3_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = p2_q3, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = p2_q3, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = p2_q3, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = p2_q3, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #p2_q3
if(mono_plots){
  jpeg(filename = paste0("monotonic_p2_q4_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = p2_q4, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = p2_q4, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = p2_q4, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = p2_q4, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #p2_q4
if(mono_plots){
  jpeg(filename = paste0("monotonic_phi_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = phi, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = phi, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = phi, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = phi, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #phi
if(mono_plots){
  jpeg(filename = paste0("monotonic_sigma1_mean_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = sigma1_mean, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = sigma1_mean, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = sigma1_mean, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = sigma1_mean, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #sigma1_mean
if(mono_plots){
  jpeg(filename = paste0("monotonic_sigma1_rate_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = sigma1_rate, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = sigma1_rate, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = sigma1_rate, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = sigma1_rate, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #sigma1_rate
if(mono_plots){
  jpeg(filename = paste0("monotonic_start_q_", infType, ".jpeg"), width =840, height = 840, units = 'px', res = 100)
  
  p1=ggplot(bTB_wl_scaled, aes(x = start_q, y = `Total Infected`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p2=ggplot(bTB_wl_scaled, aes(x = start_q, y = fadeout)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p3=ggplot(bTB_wl_scaled, aes(x = start_q, y = `fadeout time`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  p4=ggplot(bTB_wl_scaled, aes(x = start_q, y = `Hunt Prevalence`)) + 
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3))
  
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  
  dev.off()
} #start_q

setwd(pth)
#####################

###########
##  PRCC ##
###########
library(epiR)
response <- c('Total Infected','fadeout','fadeout time','Hunt Prevalence')
parameters <- names(bTB_wl_scaled[,!is.na(bTB_wl_scaled[1,])])
parameters <- parameters[!(parameters %in% response)]
parmeters_sig <- parameters
p_vals <- data.frame(parameters = parameter_names)

setwd('./sensitivity_analysis/PRCC_results/')
#True prevalence PRCC model
# 1. create model; 2. run model, 3. save model output
TruPrev <- bTB_wl_scaled[ c(parameters, response[1]) ]
prcc.TruPrev <- epi.prcc(TruPrev, sided.test = 2)
write.csv( prcc.TruPrev, file = paste0("transSensitivityPRCC_truPrev_", infType,".csv"), row.names = F )
p_vals$TruPrev <- prcc.TruPrev$p.value
#p_vals$min <- prcc.TruPrev$p.value


#fadeout PRCC model
Fadeout <- bTB_wl_scaled[c(parameters, response[2])]
prcc.Fadeout<-epi.prcc(Fadeout, sided.test = 2)
write.csv( prcc.Fadeout, file = paste0("transSensitivityPRCC_fadeout_", infType,".csv"), row.names = F )
p_vals$Fadeout <- prcc.Fadeout$p.value
#p_vals$min <- pmin(p_vals$min, prcc.Fadeout$p.value)
#p_vals$max <- pmax(p_vals$min, prcc.Fadeout$p.value)


#True prevalence PRCC model
FadeoutTime <- bTB_wl_scaled[ c(parameters, response[3]) ]
prcc.FadeoutTime <- epi.prcc(FadeoutTime, sided.test = 2)
write.csv( prcc.FadeoutTime, file = paste0("transSensitivityPRCC_truPrev_", infType,".csv"), row.names = F )
p_vals$FadeoutTime <- prcc.FadeoutTime$p.value
#p_vals$min <- pmin(p_vals$min, prcc.FadeoutTime$p.value)
#p_vals$max <- pmax(p_vals$max, prcc.FadeoutTime$p.value)

#fadeout PRCC model
HuntPrev <- bTB_wl_scaled[c(parameters, response[4])]
prcc.HuntPrev <- epi.prcc(HuntPrev, sided.test = 2)
write.csv( prcc.HuntPrev, file = paste0("transSensitivityPRCC_fadeout_", infType,".csv"), row.names = F )
p_vals$HuntPrev <- prcc.HuntPrev$p.value
#p_vals$min <- pmin(p_vals$min, prcc.HuntPrev$p.value)
#p_vals$max <- pmax(p_vals$max, prcc.HuntPrev$p.value)

#Plots
parameter_names <- c('carrying.capacity','hunt.mortality','base.mortality','dens.dep.asymetry','dens.dep.mortality','max.birth.rate','proportion.SS','birth.timing','birth.duration','transmission.rate','SS.contact.factor','latency.mean','latency.rate')

# parameter_names <- as.expression(c('carrying.capacity (K)', paste('hunt.mortality (', expression(eta),')'), paste('base.mortality (', expression(eta),')'), 
#                                    paste('dens.dep.asymetry (', expression(theta),')'), paste('dens.dep.mortality (', expression(gamma),')'),
#                                    'max.birth.rate (A)', paste('prop.SS (', expression(xi),')'), paste('birth.timing (', expression(omega),')'), 'birth.duration (s)',
#                                    paste('trans.rate (', expression(beta),')'), paste('SS.contact.factor (', expression(phi),')'), 
#                                    paste('latency.mean (', expression(sigma),')'), paste('latency.rate (', expression(sigma),')'))) %>% as.vector()

#Setting the order of the parameter
#prcc.TruPrev$parameters <- parameters
#prcc.Fadeout$parameters <- parameters
#prcc.FadeoutTime$parameters <- parameters
#prcc.HuntPrev$parameters <- parameters

#Setting the order of the parameter
prcc.TruPrev$parameters <- parameter_names
prcc.Fadeout$parameters <- parameter_names
prcc.FadeoutTime$parameters <- parameter_names
prcc.HuntPrev$parameters <- parameter_names


#ordering the number infected results by gamma value, and setting all other models to match the order of the number infected model
prcc.TruPrev <- prcc.TruPrev[order(prcc.TruPrev$est),]
prcc.Fadeout <- prcc.Fadeout[order(match(prcc.Fadeout$parameters,prcc.TruPrev$parameters)),]
prcc.FadeoutTime <- prcc.FadeoutTime[order(match(prcc.FadeoutTime$parameters,prcc.TruPrev$parameters)),]
prcc.HuntPrev <- prcc.HuntPrev[order(match(prcc.HuntPrev$parameters,prcc.TruPrev$parameters)),]
p_vals <- p_vals[order(match(p_vals$parameters,prcc.TruPrev$parameters)),]

prcc.trans <- cbind(TruPrev = prcc.TruPrev, Fadeout = prcc.Fadeout$est, FadeoutTime = prcc.FadeoutTime$est, HuntPrev = prcc.HuntPrev$est)
prcc.trans <- prcc.trans[order(prcc.trans$TruPrev.est),]



prcc.trans <- prcc.trans[,c("TruPrev.est", "TruPrev.parameters", "Fadeout", "FadeoutTime", "HuntPrev")]
prcc.trans.long <- melt(prcc.trans, id = c("TruPrev.parameters"))

p_vals.long <- melt(p_vals, id = c("parameters"))

colnames(prcc.trans.long)[1] <- "parameters"
prcc.trans.long$parameters2 <- factor(prcc.trans.long$parameters, levels = prcc.trans$TruPrev.parameters)
prcc.trans.long$time_ordered = factor(prcc.trans.long$variable, levels=c('TruPrev.est','Fadeout', 'FadeoutTime', 'HuntPrev'))


#Renaming facets
variable_names <- list(
  "TruPrev.est" = "True Prevalence",
  "Fadeout" = "Fadeout",
  "FadeoutTime" = "Fadeout Time" ,
  "HuntPrev" = "Observed Prevalence"
)

variable_labeller <- function(variable,value){
  return(variable_names[value])
}


significant=T
if(significant){
  prcc.trans.long$value[p_vals.long$value >= .05] <- NA
}

blue<-c('#2171b5')
#plot PRCC results

sens <- ggplot() + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14)) +
  geom_col(data = prcc.trans.long, aes(y = value, x = parameters2, fill=blue)) +
  #scale_x_discrete(position = "top") +
  #scale_y_discrete(limits=c(0, 0.25, 0.5), expand=c(0.02,0.02)) +
  scale_y_discrete(limits=c(-1, -.5, 0, .5, 1)) +
  labs(y = NULL, x = NULL) + 
  coord_flip() + 
  facet_grid(~time_ordered, labeller= variable_labeller) +
  scale_fill_manual(values = blue) + 
  theme(panel.spacing.y = unit(1.5, "lines")) + 
  theme(strip.text = element_text(size = 12)) 
sens <- annotate_figure(sens, top = text_grob("PRCC Analysis", 
                                                      color = "black", face = "bold", size = 18))
sens

jpeg(filename = paste0("PRCC_Trans_", infType,".jpeg"), width = 1440, height = 840, units = 'px', res = 100)
  print(sens)
dev.off()

####################
## Logistic Model ##
####################
library(sjPlot)
library(sjmisc)
library(sjlabelled)

single.model = as.formula( paste( 'RESPONSE~', paste(parameters, collapse = '+') ) )
p_vals_lm <- data.frame(parameters = parameter_names)
#TruPrev: No interaction regression 
single.model.TruPrev = update.formula(single.model, `Total Infected` ~ .)
lm.TruPrev.PRCC = lm(single.model.TruPrev, data = bTB_wl_scaled)
write.csv(lm.TruPrev.PRCC, file=paste0("lmOutputTruPrev_", infType, ".csv"), row.names = F)
summary(lm.TruPrev.PRCC)
tab_model(lm.TruPrev.PRCC)

values.in.PRCC <- lm.TruPrev.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC)) #regression results scaled by the largest abs regression value, premsize here

lmTruPrevScale<-data.frame(values.in.PRCC.scaled, row.names = parameters)

p_vals_lm$TruPrev <- summary(lm.TruPrev.PRCC)$coefficients[-1,4]


#Fadeout: No interaction regression 
single.model.Fadeout = update.formula(single.model, fadeout ~ .)
lm.Fadeout.PRCC = lm(single.model.Fadeout, data = bTB_wl_scaled)
write.csv(lm.Fadeout.PRCC, file=paste0("lmOutputFadeout_", infType, ".csv"), row.names = F)
summary(lm.Fadeout.PRCC)
tab_model(lm.Fadeout.PRCC)

values.in.PRCC <- lm.Fadeout.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC))  #regression results scaled by the largest abs regression value

lmFadeoutScale<-data.frame(values.in.PRCC.scaled,row.names = parameters)

p_vals_lm$Fadeout <- summary(lm.Fadeout.PRCC)$coefficients[-1,4]

#Fadeout Time: No interaction regression 
single.model.FadeoutTime = update.formula(single.model, `fadeout time` ~ .)
lm.FadeoutTime.PRCC = lm(single.model.FadeoutTime, data = bTB_wl_scaled)
write.csv(lm.FadeoutTime.PRCC, file=paste0("lmOutputFadeoutTime_", infType, ".csv"), row.names = F)
summary(lm.FadeoutTime.PRCC)
tab_model(lm.FadeoutTime.PRCC)

values.in.PRCC <- lm.FadeoutTime.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC))  #regression results scaled by the largest abs regression value

lmFadeoutTimeScale<-data.frame(values.in.PRCC.scaled,row.names = parameters)

p_vals_lm$FadeoutTime <- summary(lm.FadeoutTime.PRCC)$coefficients[-1,4]

#Hunt Prevalence: No interaction regression 
single.model.HuntPrev = update.formula(single.model, `Hunt Prevalence` ~ .)
lm.HuntPrev.PRCC = lm(single.model.HuntPrev, data = bTB_wl_scaled)
write.csv(lm.HuntPrev.PRCC, file=paste0("lmOutputHuntPrev_", infType, ".csv"), row.names = F)
summary(lm.HuntPrev.PRCC)
tab_model(lm.HuntPrev.PRCC)

values.in.PRCC <- lm.HuntPrev.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC))  #regression results scaled by the largest abs regression value

lmHuntPrevScale<-data.frame(values.in.PRCC.scaled,row.names = parameters)
p_vals_lm$HuntPrev <- summary(lm.HuntPrev.PRCC)$coefficients[-1,4]

#Plot 
lmTruPrevScale$parameters <- parameter_names
lmFadeoutTimeScale$parameters <- parameter_names
lmHuntPrevScale$parameters <- parameter_names
lmFadeoutScale$parameters <- parameter_names

#ordering values for plot
lmTruPrevScale <-  lmTruPrevScale[order(match(lmTruPrevScale$parameters,       prcc.TruPrev$parameters)),]
lmFadeoutTimeScale <-lmFadeoutTimeScale[order(match(lmFadeoutTimeScale$parameters,   prcc.TruPrev$parameters)),]
lmHuntPrevScale <-  lmHuntPrevScale[order(match(lmHuntPrevScale$parameters,       prcc.TruPrev$parameters)),]
lmFadeoutScale <- lmFadeoutScale[order(match(lmFadeoutScale$parameters,     prcc.TruPrev$parameters)),]

p_vals_lm <- p_vals_lm[order(match(p_vals_lm$parameters,     prcc.TruPrev$parameters)),]


lmScale <- cbind(TruPrev=lmTruPrevScale, FadeoutTime=lmFadeoutTimeScale$values.in.PRCC.scaled, HuntPrev=lmHuntPrevScale$values.in.PRCC.scaled, Fadeout=lmFadeoutScale$values.in.PRCC.scaled)
 # uncomment the line below if you want lm plots in ascending order -- i.e. probably not identical to PRCC ordering
#lmScale <- lmScale[order(lmScale$TruPrev.values.in.PRCC.scaled),]
lmScale <- lmScale[,c("TruPrev.values.in.PRCC.scaled", "TruPrev.parameters", "Fadeout", "FadeoutTime", "HuntPrev")]
lmScale.long <- melt(lmScale,id=c("TruPrev.parameters"))

p_vals_lm.long <- melt(p_vals_lm,id=c("parameters"))

colnames(lmScale.long)[1] <- "parameters"
lmScale.long$parameters2 <- factor(lmScale.long$parameters,levels=lmScale$TruPrev.parameters)
lmScale.long$time_ordered <- factor(lmScale.long$variable,levels=c('TruPrev.values.in.PRCC.scaled','Fadeout','FadeoutTime','HuntPrev'))

#Renaming facets
variable_names2 <- list(
  "TruPrev.values.in.PRCC.scaled" = "True Prevalence",
  "Fadeout" = "Fadeout",
  "FadeoutTime" = "Fadeout Time" ,
  "HuntPrev" = "Observed Prevalence"
)

variable_labeller <- function(variable,value){
  return(variable_names2[value])
}


if(significant){
  lmScale.long$value[p_vals_lm.long$value >= .05] <- NA
}

#plot single model results
reg <- ggplot() + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14)) +
  geom_col(data = lmScale.long, aes(y = value, x = parameters2, fill = blue)) + 
  scale_y_discrete(limits=c(-1, -.5, 0, .5, 1)) +
  coord_flip() + 
  facet_grid(~time_ordered, labeller= variable_labeller) +
  scale_fill_manual (values = blue) + 
  theme(panel.spacing.y = unit(1.5, "lines")) + 
  theme(strip.text = element_text(size = 12))
reg <- annotate_figure(reg, top = text_grob("Regression Analysis", 
                                              color = "black", face = "bold", size = 18))
reg

jpeg(filename = paste0("SingleModel_Trans_", infType, ".jpeg"), width = 1440, height = 840, units = 'px', res = 100)
reg
dev.off()


####################

