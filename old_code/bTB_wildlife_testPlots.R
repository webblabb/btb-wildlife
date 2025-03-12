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
name_out = paste0(runtype,infType,"_",pct,'-',prop_superSpreader,'-')
save_runs = T
save_plots = T
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
sto_res <- vector(mode = "list", length = length(sizes))
sto_full_res <- vector(mode = "list", length = length(sizes))
det_res <- vector(mode = "list", length = length(sizes))
#####################

for(i in 1:length(sizes)){
  
  ###########
  ## Runs ##
  ##########
  
  setwd(dirname(getActiveDocumentContext()$path))
  
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
  sto_res[[i]] <- as.data.frame(data_sto)
  sto_full_res[[i]] <- as.data.frame(sto_out)
  det_res[[i]] <- as.data.frame(data)
  
  remove(data_sto)
  remove(sto_out)
  remove(data)
  
  ##########
  
  
  ###########
  ## Plots ##
  ###########
  
  #class comparison plots
  sto_by_class <- split(as.data.frame(sto_res[[i]]), as.data.frame(sto_res[[i]])$variable)
  det_by_class <- split(as.data.frame(det_res[[i]]), as.data.frame(det_res[[i]])$variable)
  
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
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/validation_plots/"))
    
    jpeg(filename = paste0("Rplot_class_", name_out, size, ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("Rplot_N_", name_out, size, ".jpeg"))
    print(Nplot)
    dev.off()
    
    setwd(dirname(getActiveDocumentContext()$path))
  }
  
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  if(i == length(sizes)){
    
    if(save_runs){
      save( sto_res, file = paste0('data_sto_', name_out, size, Sys.Date(), '.RData'))
      save( sto_full_res, file = paste0('data_full_res_', name_out, size, Sys.Date(), '.RData'))
      save( det_res, file = paste0('data_det_', name_out, size, Sys.Date(), '.RData'))
    }
    
  }  
  ##################
  
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
source(file = "bTBwl_func.R")
years = 30
seedQuarter = 1
prop_superSpreader = 0.05
pct = 0.02
verbose = 0
reps = 300
sizes = c(10, 50, 100, 250, 500)
infType = "seeded"
runtype = "stdLambda_"
name_out = paste0(runtype,infType,"_",pct,'-',prop_superSpreader,'-')
save_runs = T
save_plots = T
verbose = 0 #verbose of -1 forces a .25 month timestep maximum
sto_res <- vector(mode = "list", length = length(sizes))
lambda_res <- vector(mode = "list", length = length(sizes))
lambda_vecs <- vector(mode = "list", length = length(sizes))

#####################

for(i in 1:length(sizes)){
  
  ###########
  ## Runs ##
  ##########
  
  setwd(dirname(getActiveDocumentContext()$path))
  
  size=as.integer(sizes[i])
  pars <- parameter_set_wl(k = size, 
                           scenario = infType, 
                           initial_exposed = pct*size, 
                           SS_prop = prop_superSpreader,
                           start_quarter = seedQuarter, 
                           test = F,
                           birth = F,
                           death_h = F,
                           death_n = F,
                           disease = F,
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
  
  data_sto <- sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  #assign data frames to list
  sto_res[[i]] <- as.data.frame(data_sto)
  lambda_res[[i]] <- as.data.frame(sto_out[,c("time", "Month", "Lambda", "rep")])
  lambda_vecs[[i]] <- getLambda_vec(data = sto_out[,c("Month", "Lambda")], type = '')
  
  sto_out <- sto_out[,c("time", "Month", "Lambda", "rep")]
  remove(data_sto)
  
  
  ##########

  
  ###########
  ## Plots ##
  ###########
  
  dat_metrics <- data.frame(min = getLambda_vec(data = sto_out, type = 'min'), 
                         max = getLambda_vec(data = sto_out, type = 'max'),
                         mean = getLambda_vec(data = sto_out, type = 'mean'),
                         median = getLambda_vec(data = sto_out, type = 'median'),
                         month = c(1:12)
                         )
  
  lambdaVec_plot <- ggplot(data = dat_metrics) +
    geom_ribbon(aes(x=month, ymax=max,ymin=min ), color = "darkgray", fill = "aquamarine", alpha=.5) +
    geom_line(aes(x=month, y=median), color = "blue") + 
    geom_line(aes(x=month, y=mean), color = "red")
  
  lambdaVal_plot <- ggplot(data = sto_out) + 
    geom_line(aes(x=time, y=Lambda, group = rep))
  lambdaVal_plot2 <- ggplot(data = sto_out) + 
    geom_line(aes(x=Month, y=Lambda, group = rep))
  
  
  if(save_plots){
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/lambda_preSim/lambda_plots/"))
    
    jpeg(filename = paste0("Rplot_LambdaVec_", name_out, size, ".jpeg"))
    print(lambdaVec_plot)
    dev.off()
    
    jpeg(filename = paste0("Rplot_LambdaVal_", name_out, size, ".jpeg"))
    print(ggarrange(lambdaVal_plot, lambdaVal_plot2, ncol=1, nrow=1))
    dev.off()
    
    setwd(dirname(getActiveDocumentContext()$path))
  }
  
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  if(i == length(sizes)){
    
    if(save_runs){
      save( sto_res, file = paste0('lambda_preSim/data_sto_', name_out, size, Sys.Date(), '.RData'))
      save( lambda_res, file = paste0('lambda_preSim/data_lambda_res_', name_out, size, Sys.Date(), '.RData'))
      save( lambda_vecs, file = paste0('lambda_preSim/data_lambda_vecs_', name_out, size, '.RData'))
    }
    
  }  
  ##################
  
}



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
runtype = "std_"
name_out = paste0(runtype,infType,"_",pct,'-',prop_superSpreader,'-')
save_runs = T
save_plots = T
test_mode = F; test_birth = F; test_death_n = F; test_death_h = F; test_disease = F;
verbose = 0 
sto_res <- vector(mode = "list", length = length(sizes))
sto_full_res <- vector(mode = "list", length = length(sizes))
det_res <- vector(mode = "list", length = length(sizes))
#####################

for(i in 1:length(sizes)){
  
  ###########
  ## Runs ##
  ##########
  
  setwd(dirname(getActiveDocumentContext()$path))
  
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
  #remove NaN or negative values
  data$value[is.nan(data$value)] <- 0 
  data$value[which(data$value < 0)] <- 0
  
  
  #stochastic model
  initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
  population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
  disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
  
  param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
  
  list.files(path = './lambda_preSim/', pattern = paste0('data_lambda_vecs_', name_out, size))
  lambda_data <- load()
  tic = Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd')
  #n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears = years; seed_quarter = as.integer(pars["start_q"]); verbose = verbose; batch_name = "test"; type = 'c';
  toc = Sys.time()
  print(toc-tic)
  
  data_sto <- sto_out[,c("rep", "tstep", "time", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>% pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  #assign data frames to list
  sto_res[[i]] <- as.data.frame(data_sto)
  sto_full_res[[i]] <- as.data.frame(sto_out)
  det_res[[i]] <- as.data.frame(data)
  
  remove(data_sto)
  remove(sto_out)
  remove(data)
  
  ##########
  
  
  ###########
  ## Plots ##
  ###########
  
  #class comparison plots
  sto_by_class <- split(as.data.frame(sto_res[[i]]), as.data.frame(sto_res[[i]])$variable)
  det_by_class <- split(as.data.frame(det_res[[i]]), as.data.frame(det_res[[i]])$variable)
  
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
    setwd(paste0(dirname(getActiveDocumentContext()$path), "/validation_plots/"))
    
    jpeg(filename = paste0("Rplot_class_", name_out, size, ".jpeg"))
    print(ggarrange(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, ncol=2, nrow=3))
    dev.off()
    
    jpeg(filename = paste0("Rplot_N_", name_out, size, ".jpeg"))
    print(Nplot)
    dev.off()
    
    setwd(dirname(getActiveDocumentContext()$path))
  }
  
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  if(i == length(sizes)){
    
    if(save_runs){
      save( sto_res, file = paste0('data_sto_', name_out, size, Sys.Date(), '.RData'))
      save( sto_full_res, file = paste0('data_full_res_', name_out, size, Sys.Date(), '.RData'))
      save( det_res, file = paste0('data_det_', name_out, size, Sys.Date(), '.RData'))
    }
    
  }  
  ##################
  
}

