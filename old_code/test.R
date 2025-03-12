##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
years = 10
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
prop_superSpreader = 0.05
pct = 0.02
verbose = 0
reps = 50
size = 50
infType = "seeded"
runtype = "CprelimTest_"
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum

##############################
sto_res <- vector(mode = "list", length = length(size))
sto_full_res <- vector(mode = "list", length = length(size))
det_res <- vector(mode = "list", length = length(size))

setwd(dirname(getActiveDocumentContext()$path))


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

data_sto <- sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")

#assign data frames to list
sto_res[[1]] <- as.data.frame(data_sto)
sto_full_res[[1]] <- as.data.frame(sto_out)
det_res[[1]] <- as.data.frame(data)

remove(data_sto)
remove(sto_out)
remove(data)

##############################

###########
## Plots ##
###########

#class comparison plots
sto_by_class <- split(as.data.frame(sto_res[[1]]), as.data.frame(sto_res[[1]])$variable)
det_by_class <- split(as.data.frame(det_res[[1]]), as.data.frame(det_res[[1]])$variable)

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

ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3)
Nplot



#save outputs
save( sto_res, file = paste0('data_sto_', runtype, Sys.Date(), '.RData'))
save( sto_full_res, file = paste0('data_full_res_', runtype, Sys.Date(), '.RData'))
save( det_res, file = paste0('data_det_', runtype, Sys.Date(), '.RData'))


setwd(paste0(dirname(getActiveDocumentContext()$path), "/validation_plots/"))
jpeg(filename = paste0("Rplot_class_", runtype, size, ".jpeg"))
print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
dev.off()
jpeg(filename = paste0("Rplot_N_", runtype, size, ".jpeg"))
print(Nplot)
dev.off()





##########################################################################################
##########################################################################################
#######################    Multiple Herd Sizes -- Continuous Time   ######################
##########################################################################################
##########################################################################################

####################
## Initialization ##
####################
rm(list = ls())
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")
years = 10
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
prop_superSpreader = 0.05
pct = 0.02
verbose = 0
reps = 50
sizes = c(10, 50, 100, 250, 500)
infType = "seeded"
runtype = "std_"
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-')
save_runs = T
save_plots = T
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
#####################

foreach(i = 1:length(sizes)) %do% {
  ###########
  ## Runs ##
  ##########
  
  setwd(dirname(getActiveDocumentContext()$path))
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
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/runs/validation_plots/"))
    
    jpeg(filename = paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(Nplot)
    dev.off()
    
    setwd(dirname(getActiveDocumentContext()$path))
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################

  if(save_runs){
    save( models_res, file = paste0('runs/', name_out, size, '-', Sys.Date(), '.RData'))
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
years = 15
seedQuarter = 1
prop_superSpreader = 0.05
pct = 0.02
verbose = 0
reps = 250
sizes = c(10, 50, 100, 250, 500)
infType = "seeded"
runtype = "stdLambda_"
name_out = paste0(runtype,infType,"_",pct,'-',prop_superSpreader,'-')
save_runs = T
save_plots = T
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)

#####################
foreach(i = 1:length(sizes)) %do% {
  
  ###########
  ## Runs ##
  ##########
  
  setwd(dirname(getActiveDocumentContext()$path))
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
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/lambda_preSim/lambda_plots/"))
    
    jpeg(filename = paste0("Rplot_LambdaVec_", name_out, size, ".jpeg"))
    print(lambdaVec_plot)
    dev.off()
    
    jpeg(filename = paste0("Rplot_LambdaVal_", name_out, size, ".jpeg"))
    print(ggarrange(lambdaVal_plot, lambdaVal_plot2, ncol=1, nrow=2))
    dev.off()
    
    setwd(dirname(getActiveDocumentContext()$path))
  }
  remove(lambdaVec_plot, lambdaVal_plot, lambdaVal_plot2)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/lambda_preSim/"))
    
    save( lambda_res, file = paste0( name_out, size, '-', Sys.Date(), '.RData'))
    save( lambda_out, file = paste0('vals_', name_out, size, '.RData'))
    
    setwd(dirname(getActiveDocumentContext()$path))
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
years = 10
times = seq(from=0, to=years*12, by=12/365) #daily time steps
seedQuarter = 1
prop_superSpreader = 0.05
pct = 0.02
verbose = 0
reps = 50
sizes = c(10, 50, 100, 250, 500)
infType = "seeded"
runtype = "stdD_"
type_of_integral = 3;
name_out = paste0(runtype, infType,"_", pct,'-', prop_superSpreader, '-it', type_of_integral, '-')
lambda_name_out = paste0('stdLambda_', infType,"_", pct,'-', prop_superSpreader, '-')
save_runs = T
save_plots = T
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
#####################

foreach(i = 1:length(sizes)) %do% {
  ###########
  ## Runs ##
  ##########
  
  setwd(dirname(getActiveDocumentContext()$path))
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
  
  load(paste0('lambda_preSim/vals_', lambda_name_out, size, '.RData'))
  
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out[[1]]$max, integrate_type = type_of_integral)
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
  
  print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
  print(Nplot)
  
  if(save_plots){
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/runs/validation_plots/"))
    
    jpeg(filename = paste0("Rplot_class_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("Rplot_N_", name_out, size, '-', Sys.Date(), ".jpeg"))
    print(Nplot)
    dev.off()
    
    setwd(dirname(getActiveDocumentContext()$path))
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    save( models_res, file = paste0('runs/', name_out, size, '-', Sys.Date(), '.RData'))
  }
  remove(models_res) #clear from workspace
  ##################
  
}

stopCluster(cl)




