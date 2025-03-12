#parameter setting method#
#k takes positive integer values, and scenario is 1 or 2#
#scenario 1 has a constant force of infection whereas 2 has ?10? unreactive exposed individuals
parameter_set_wl_test <-function(k, SS_prop = 0, scenario, SS_inf = F, inf_0 = 10){
  
  parameter_names=c('S_0', 'E1_0', 'I_0', 'SuperS_0', 'SuperE1_0', 'SuperI_0', 'K', 'eta_hunt', 'eta_q1', 'eta_q2', 'eta_q3', 'eta_q4', 'alpha_q1', 'alpha_q2', 'alpha_q3', 'alpha_q4', 'ksi', 'beta', 'p1', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi', 'sigma1_mean', 'sigma1_rate')
  #initialize parameters to be set
  if(SS_inf && SS_prop <= 0 ){print('Warning: incompatible arguments passed. SS_inf should only be enabled if SS_prop is nonzero')}
  if(k<10){print('insufficient starting population: need k >= 10. No output produced.'); print(k); break}
  if(k<inf_0){print('number of initial infections exceeds start population size'); print(paste0('k = ',k)); print(paste0('inf_0 = ', inf_0)); break}
  K = k
  S_0 = k
  E1_0 = 0
  SuperE1_0 = 0
  eta_q1 = .1
  eta_q2 = .1
  eta_q3 = .1
  eta_q4 = .1
  eta_hunt = 1.5
  alpha_q1 = 0
  alpha_q2 = .05
  alpha_q3 = 0
  alpha_q4  =0
  ksi =  0
  beta = 2
  p1 = .1
  p2_q1 = .05
  p2_q2 = .1
  p2_q3 = 0
  p2_q4 = .01
  sigma1_mean = 3.5
  sigma1_rate = .2
  
  
  
  if(scenario == 1 ){
    E1_0 = inf_0
    S_0 = S_0 - inf_0
    p2_q1 = 0
    p2_q2 = 0
    p2_q3 = 0
    p2_q4 = 0
  }
  #set SS proportion and contact factor
  
  phi = 5 #increased farm contact factor
  ksi = SS_prop
  
  #calculate number of super spreaders and update S_0
  nSS=floor(S_0*SS_prop)
  SuperS_0 = nSS
  S_0 = S_0 - nSS 
  
  #make intitally exposed individuals SS instead if desired
  if(E1_0 != 0 && SS_inf == T){
    E1_0 = 0
    SuperE1_0 = inf_0
  }
  
  
  
  #parameters not influenced by size, type, or infection introduction are hard coded in. All others set in loop.
  paramvec=c(S_0=S_0, E1_0=E1_0, I_0=0, SuperS_0 = nSS, SuperE1_0 = SuperE1_0, SuperI_0 = 0,
             K=k, eta_hunt=eta_hunt, eta_q1=eta_q1, eta_q2=eta_q2, eta_q3=eta_q3, eta_q4=eta_q4, alpha_q1=alpha_q1, alpha_q2=alpha_q2, alpha_q3=alpha_q3, alpha_q4=alpha_q4, ksi=ksi,
             beta=beta, p1=p1, p2_q1=p2_q1, p2_q2=p2_q2, p2_q3=p2_q3, p2_q4=p2_q4, phi=phi, sigma1_mean=sigma1_mean, sigma1_rate=sigma1_rate)
  if( min(paramvec) < 0 ){print('negative parameter value detected. No output produced'); print(paramvec); break}
  return(paramvec)  #output vector of parameters for a run of given farm size type and infection scenario
}


###########################################################################################
#wrapper method code######################################################################
#########################################################################################

setwd("~/Documents/Brandon/bTB_wildlife_code/") #set working ddirectory to model code location
wildlife_model<-function(n_reps, parameters, initial_state, nyears, seed_quarter, batch_name, verbose = 0, version = 2){
  
  parvec = c( n_reps, #Number of replicates
              initial_state$S_0, #0, start Susc.
              initial_state$E1_0, #1, start E1
              initial_state$I_0, #2, start I
              initial_state$SuperS_0, #3
              initial_state$SuperE1_0, #4
              initial_state$SuperI_0, #5
              parameters$K, #6, Herd size (county level)
              parameters$eta_hunt, #7, Hunting induced mortality rate
              parameters$eta_q1, #8, Natural mortality rate
              parameters$eta_q2, #9, Natural mortality rate
              parameters$eta_q3, #10, Natural mortality rate
              parameters$eta_q4, #11, Natural mortality rate
              parameters$alpha_q1, #12 Seasonal birth rate
              parameters$alpha_q2, #13 Seasonal birth rate
              parameters$alpha_q3, #14 Seasonal birth rate
              parameters$alpha_q4, #15 Seasonal birth rate
              parameters$ksi, #16, proportion of newborns that are superSpreaders
              parameters$beta, #17 Density Dependent Transmission rate
              parameters$p1, #18 deer-deer contact rate
              parameters$p2_q1, #19 C-to-wildl contact rate quarter 1
              parameters$p2_q2, #20 C-to-wildl contact rate quarter 2
              parameters$p2_q3, #21 C-to-wildl contact rate quarter 3
              parameters$p2_q4, #22 C-to-wildl contact rate quarter 4
              parameters$phi, #23, increased SS farm contact rate factor
              parameters$sigma1_mean, #24, mean exposed period length
              parameters$sigma1_rate, #25, exposed rate
              nyears, #26
              seed_quarter, #27
              batch_name, #28
              verbose, #29
              version #30
  )
  strvec = format(parvec, digits = 5)
  
  #The path to the bTB cpp binary file must be set correctly in the sys call below:
  #r <- system2(command = "/home/webblab/Documents/Brandon/bTB_wildlife_code/a.out", args=strvec, stdout=TRUE) #call for v1.0 code (i.e. no super spreaders, hunting calc)
  r <- system2(command = "./wl_model_3.0.exe", args=strvec, stdout=TRUE) #call for most updated model code
  res <- read.table(text=r, header=TRUE, sep=';', check.names = FALSE)
  return (res)
}

#setting run state variables#

reps=10
years=20
popsize=500
seedQuarter=1
infType=2

root="~/Documents/Brandon/bTB_wildlife_code/" #output file location

par<-parameter_set_wl_test(k=popsize, SS_prop=0, scenario=infType, SS_inf = F, inf_0 = 10)

initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
disease_parameters<-data.frame(beta=par[18], p1=0, p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])

param_vals<-data.frame(merge(population_parameters,disease_parameters))

#run#
vsn=2
run_time = replicate(100, NA)
reps=10
years=20
popsize=500
seedQuarter=1
infType=2

tic=Sys.time()
run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 5, version = vsn)
toc=Sys.time()
toc-tic

for(i in 1:length(run_time)){
tic=Sys.time()
run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
fname='wildlife_SS_0_test.csv'
write.csv(run, paste0(root,fname))
toc=Sys.time()
run_time[i] = toc - tic
}
mean(run_time)
###testing###
#############
#n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears=years; seed_quarter=1; batch_name="test";
#cbind(reps, initial_state, param_vals, years, "test")
#X=cbind(reps, initial_state, param_vals, years, "test")
#length(X)
#############


############################################################################
############################################################################
############################################################################


#model runs pipeline
#set working directory to folder with bTB Wildlife model executable file
setwd("~/Documents/Brandon/bTB_wildlife_code")

## Make a directory for runs

if (!dir.exists("bTB_wildlife_runs")){
  dir.create("bTB_wildlife_runs")
}

#simulation duration parameters
reps = 1000 # 1000 replicates should be plenty. scaled down for testing
years = 40 #arbitrarily chosen

#initialize relevant categorical range vectors
N = c(10, 50, 100, 200, 500, 750, 1000) #population sizes
N_inf0 = c(1:10) #number of initially infected for seeded infections
Prop_SS = c(0, .05, .1, .25, .5, .75, 1) #proportion of individuals that are superSpreaders
hunt_rates = c(0, .5, 1,  5, 10, 25, 50) 

alpha_key = letters #alphabet list applied to names. this precents mis-ordering in data reads due to inconsistencies in numbers in alpha-order.
infType = c(1,2) #input value for infection source input
infectionSource = c('seeded', 'spillover') #name for infection source
start_SS_inf = F #seeded infections will have all infections start in SS individuals if true 
                 #only enable this option if super spreaders are enabled 

#single core runs
#############################################################################################################################
#Standard Runs - SS_prop, Ninf_0 fixed
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Standard_Runs")){
  dir.create("Standard_Runs")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

tic = Sys.time()
for(i in 1:length(infType)){
 
  for(n in 1:length(N)){
    
    #set parameters for this run
    par<-parameter_set_wl_test(k=N[n], SS_prop=Prop_SS[3], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = 10) #
    
    #reformat parameters for model input
    initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
    population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
    disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #run#
    run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = 2)
    fname = paste0('bTBwl_', letters[n], '_', N[n], '_', Prop_SS[3], '_', infectionSource[i],'.csv') #name formatting: bTBwl_popsize_SSprop_infectionSource.csv
    if(start_SS_inf){
      fname = paste0('bTBwl_', letters[n], '_', N[n], '_', Prop_SS[3], '_', infectionSource[i],'_startSSinf.csv') 
    }
    write.csv(run, paste0('./bTB_wildlife_runs/Standard_Runs/',fname))
    
    
  }

}
toc = Sys.time()
toc-tic


#varied initial infected runs - SS_prop fixed
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_initialInf_Runs")){
  dir.create("Var_initialInf_Runs")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

tic = Sys.time()

for( k in 1:length(infType) ){
  
  for( n in 1:length(N) ){
    
    for( i in 1:length(N_inf0) ){
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[n], SS_prop=Prop_SS[3], scenario=infType[k], SS_inf = start_SS_inf, inf_0 = N_inf0[i]) #
    
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
    
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = 2)
      fname = paste0('bTBwl_', letters[n], letters[i], '_N', N[n], '_', Prop_SS[3], '_', N_inf0[i], '#inf_', infectionSource[k],'.csv') #name formatting: bTBwl_popsize_SSprop_E1_0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_', letters[n], letters[i], '_N', N[n], '_', Prop_SS[3], '_', N_inf0[i], '_', infectionSource[k],'_startSSinf.csv')
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_initialInf_Runs/', fname))
    }
    
  }
  
}
toc = Sys.time()
toc-tic



#varied SSprop runs - Start population and Ninf_0 fixed
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_SSprop")){
  dir.create("Var_SSprop")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

for( i in 1:length(infType) ){
  
  for( n in 1:length(Prop_SS) ){
    
    fixPop = 5
    #set parameters for this run
    par<-parameter_set_wl_test(k=N[fixPop], SS_prop=Prop_SS[n], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = 10) #
    
    #reformat parameters for model input
    #### P1 = 0 #####
    initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
    population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
    disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #run#
    run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = 2)
    fname = paste0('bTBwl_', letters[n], '_', N[fixPop], '_', as.integer(Prop_SS[n]*100), '%SS_',  infectionSource[i], '.csv') #name formatting: bTBwl_popsize_percentSS_infectionSource.csv
    if(start_SS_inf){
      fname = paste0('bTBwl_', letters[n], '_', N[fixPop], '_', as.integer(Prop_SS[n]*100), '%SS_',  infectionSource[i], '_startSSinf.csv')
    }
    write.csv(run, paste0('./bTB_wildlife_runs/Var_SSprop/',fname))
    
    
  }
  
}


#varied hunt rates runs - Start population and Ninf_0 fixed
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_huntRates")){
  dir.create("Var_huntRates")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

for(i in 1:length(infType)){
  for(k in 1:length(N)){
    for(n in 1:length(hunt_rates)){
    
      fixPop = 5
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[k], SS_prop=Prop_SS[3], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = N_inf0[10]) #
    
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=hunt_rates[n], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
    
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = 2)
      fname = paste0('bTBwl_', letters[n], '_', N[k], '_huntRate', hunt_rates[n], '_', N_inf0[10], '_', infectionSource[i], '.csv') #name formatting: bTBwl_popsize_huntRate_nInf0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_', letters[n], '_',  N[k], '_huntRate', hunt_rates[n], '_', N_inf0[10], '_', infectionSource[i], '_startSSinf.csv') #adjust name if SS start infected
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_huntRates/',fname))
    
    }
  }
  
}


#varied hunt rates runs - start with 2% infections
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_huntRates")){
  dir.create("Var_huntRates")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

for(i in 1:length(infType)){
  for(k in 1:length(N)){
    for(n in 1:length(hunt_rates)){
      
      
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[k], SS_prop=Prop_SS[3], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = ceiling(N[k]*0.02)) #designate number of infected as 2 percent of totla(rounded up)
      
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=hunt_rates[n], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
      
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
      
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = 2)
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[k], '_huntRate', hunt_rates[n], '_2%inf_', infectionSource[i], '_2%inf.csv') #name formatting: bTBwl_popsize_huntRate_nInf0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_',  N[k], '_huntRate', hunt_rates[n], '_2%inf_', infectionSource[i], '_startSSinf_2%inf.csv') #adjust name if SS start infected
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_huntRates/',fname))
      
    }
  }
  
}
#############################################################################################################################

#parallel runs -- standard runs
#############################################################################################################################
library(foreach)
library(doParallel)
detectCores()
#n.cores <- 6
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
#simulation duration parameters
reps = 1000 # 1000 replicates should be plenty. scaled down for testing
years = 40 #arbitrarily chosen
vsn = 3 #model version to be used -- currently 1 or 2\

#initialize relevant categorical range vectors
N = c(10, 50, 100, 200, 500, 750, 1000) #population sizes
N_inf0 = c(1:10) #number of initially infected for seeded infections
Prop_SS = c(0, .05, .1, .25, .5, .75, 1) #proportion of individuals that are superSpreaders
hunt_rates = c(0, .5, 1,  5, 10, 25, 50) 

alpha_key = letters #alphabet list applied to names. this precents mis-ordering in data reads due to inconsistencies in numbers in alpha-order.
infType = c(1,2) #input value for infection source input
infectionSource = c('seeded', 'spillover') #name for infection source
start_SS_inf = F #seeded infections will have all infections start in SS individuals if true 
#only enable this option if super spreaders are enabled 

tic = Sys.time()
setwd("~/Documents/Brandon/bTB_wildlife_code")
out_file_names <- foreach(n = 1:length(N), .combine = "cbind") %dopar% {
  
  for(i in 1:length(infType)){
    #set parameters for this run
    par<-parameter_set_wl_test(k=N[n], SS_prop=Prop_SS[3], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = 10) #
  
    #reformat parameters for model input
    initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
    population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
    disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
  
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
  
    #run#
    run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
    fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[n], '_', Prop_SS[3], '_', infectionSource[i],'.csv') #name formatting: bTBwl_popsize_SSprop_infectionSource.csv
    if(start_SS_inf){
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[n], '_', Prop_SS[3], '_', infectionSource[i],'_startSSinf.csv') 
    }
    write.csv(run, paste0('./bTB_wildlife_runs/Standard_Runs/',fname))
  }
  
  fname
  
}
toc = Sys.time()
toc-tic
stopCluster(cl)
out_file_names


#parallel runs -- fadeout runs by #inf
#############################################################################################################################


detectCores()
#n.cores <- 6
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)

setwd("~/Documents/Brandon/bTB_wildlife_code")
tic = Sys.time()
out_file_name <- foreach( i = 1:length(N_inf0) )  %dopar% {
  
  for( n in 1:length(N) ){
    
    for( k in 1:length(infType) ){
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[n], SS_prop=Prop_SS[3], scenario=infType[k], SS_inf = start_SS_inf, inf_0 = N_inf0[i]) #
      
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
      
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
      
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], letters[i], '_N', N[n], '_', Prop_SS[3], '_', N_inf0[i], '#inf_', infectionSource[k],'.csv') #name formatting: bTBwl_popsize_SSprop_E1_0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_v.',vsn,'_', letters[n], letters[i], '_N', N[n], '_', Prop_SS[3], '_', N_inf0[i], '#inf_', infectionSource[k],'_startSSinf.csv')
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_initialInf_Runs/', fname))
    }
    
  }
  
}

toc = Sys.time()
toc-tic
stopCluster(cl)
out_file_names



#parallel runs #fadeout runs by %inf
#############################################################################################################################


detectCores()
#n.cores <- 6
pctInf = c(.01,.02,.03,.04,.05,.06,.07,.08,.09,.1)
n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)

setwd("~/Documents/Brandon/bTB_wildlife_code")
tic = Sys.time()
out_file_name <- foreach( i = 1:length(pctInf) )  %dopar% {
  
  for( n in 1:length(N) ){
    
    for( k in 1:length(infType) ){
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[n], SS_prop=Prop_SS[3], scenario=infType[k], SS_inf = start_SS_inf, inf_0 = floor(N[n]*pctInf[i])) #
      
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
      
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
      
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], letters[i], '_N', N[n], '_', Prop_SS[3], '_', pctInf[i]*100, '%inf_', infectionSource[k],'.csv') #name formatting: bTBwl_popsize_SSprop_E1_0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_v.',vsn,'_', letters[n], letters[i], '_N', N[n], '_', Prop_SS[3], '_', pctInf[i]*100, '%inf_', infectionSource[k],'_startSSinf.csv')
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_initialInf_Runs/', fname))
    }
    
  }
  
}

toc = Sys.time()
toc-tic
stopCluster(cl)
out_file_names
#############################################################################################################################

#varied SSprop runs - Start population and Ninf_0 fixed
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_SSprop")){
  dir.create("Var_SSprop")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
tic = Sys.time()

out_file_name <- foreach( n = 1:length(Prop_SS) ) %dopar% {
  
  for( i in 1:length(infType) ){
    
    fixPop = 5
    #set parameters for this run
    par<-parameter_set_wl_test(k=N[fixPop], SS_prop=Prop_SS[n], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = 10) #
    
    #reformat parameters for model input
  
    initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
    population_parameters<-data.frame(K=par[7], eta_hunt=par[8], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
    disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
    
    param_vals<-data.frame(merge(population_parameters,disease_parameters))
    
    #run#
    run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
    fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[fixPop], '_', as.integer(Prop_SS[n]*100), '%SS_',  infectionSource[i], '.csv') #name formatting: bTBwl_popsize_percentSS_infectionSource.csv
    if(start_SS_inf){
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[fixPop], '_', as.integer(Prop_SS[n]*100), '%SS_',  infectionSource[i], '_startSSinf.csv')
    }
    write.csv(run, paste0('./bTB_wildlife_runs/Var_SSprop/',fname))
  }
}
toc = Sys.time()
toc-tic
stopCluster(cl)
out_file_names
#############################################################################################################################

#varied hunt rates runs - Start population and Ninf_0 fixed
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_huntRates")){
  dir.create("Var_huntRates")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
tic = Sys.time()

out_file_name <- foreach(k = 1:length(N)) %dopar% {
  for(i in 1:length(infType)){
    for(n in 1:length(hunt_rates)){
      
      fixPop = 5
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[k], SS_prop=Prop_SS[3], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = N_inf0[10]) #
      
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=hunt_rates[n], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
      
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
      
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[k], '_huntRate', hunt_rates[n], '_', N_inf0[10], '_', infectionSource[i], '.csv') #name formatting: bTBwl_popsize_huntRate_nInf0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_',  N[k], '_huntRate', hunt_rates[n], '_', N_inf0[10], '_', infectionSource[i], '_startSSinf.csv') #adjust name if SS start infected
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_huntRates/',fname))
      
    }
  }
}

toc = Sys.time()
toc-tic
stopCluster(cl)
out_file_names

#############################################################################################################################

#varied hunt rates runs - start with 2% infections
setwd("~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs")
if (!dir.exists("Var_huntRates")){
  dir.create("Var_huntRates")
}
setwd("~/Documents/Brandon/bTB_wildlife_code")

n.cores <- floor( detectCores()*(3/4) ) #flexible core determination
cl<-makeCluster(n.cores)
registerDoParallel(cl)
tic = Sys.time()

out_file_name <- foreach(k = 1:length(N)) %dopar% {
  for(i in 1:length(infType)){
    for(n in 1:length(hunt_rates)){
      
      
      #set parameters for this run
      par<-parameter_set_wl_test(k=N[k], SS_prop=Prop_SS[3], scenario=infType[i], SS_inf = start_SS_inf, inf_0 = ceiling(N[k]*0.02)) #designate number of infected as 2 percent of totla(rounded up)
      
      #reformat parameters for model input
      initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3], SuperS_0=par[4], SuperE1_0=par[5], SuperI_0=par[6])
      population_parameters<-data.frame(K=par[7], eta_hunt=hunt_rates[n], eta_q1=par[9], eta_q2=par[10], eta_q3=par[11], eta_q4=par[12], alpha_q1=par[13], alpha_q2=par[14], alpha_q3=par[15], alpha_q4=par[16], ksi=par[17]) 
      disease_parameters<-data.frame(beta=par[18], p1=par[19], p2_q1=par[20], p2_q2=par[21], p2_q3=par[22], p2_q4=par[23], phi=par[24], sigma1_mean=par[25], sigma1_rate=par[26])
      
      param_vals<-data.frame(merge(population_parameters,disease_parameters))
      
      #run#
      run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test", verbose = 0, version = vsn)
      fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_', N[k], '_huntRate', hunt_rates[n], '_2%inf_', infectionSource[i], '_2%inf.csv') #name formatting: bTBwl_popsize_huntRate_nInf0_infectionSource.csv
      if(start_SS_inf){
        fname = paste0('bTBwl_v.',vsn,'_', letters[n], '_',  N[k], '_huntRate', hunt_rates[n], '_2%inf_', infectionSource[i], '_startSSinf_2%inf.csv') #adjust name if SS start infected
      }
      write.csv(run, paste0('./bTB_wildlife_runs/Var_huntRates/',fname))
      
    }
  }
}

toc = Sys.time()
toc-tic
stopCluster(cl)
out_file_names

#########################
#Create plots for output#
#########################

writePlot(data_vec = c(1), totalPop_vec = NULL, herdSize = "500", infType = "seeded", dataType = "", category = "", n_reps = 0, n_years = 0, minMax = F, rep_Plots = 0, Scaled = F, Regular = T, Paired = F, as.jpeg = F, to.console = T, color = 'blue')
  
tic = Sys.time()
Fade_data = plotsData_reformat(version = vsn, type = c('fadeout','seeded.csv',''), runPath = "~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_initialInf_Runs/", scaled = F, scaleBy = NULL, herdSizes = N)
save( Fade_data, file = 'bTB_wildlife_runs/Fade_data.RData')
toc = Sys.time() #takes about 4 min
toc-tic

tic = Sys.time()
Fade_data_pct = plotsData_reformat(version = vsn, type = c('fadeout%','seeded.csv',''), runPath = "~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_initialInf_Runs/", scaled = F, scaleBy = NULL, herdSizes = N)
save( Fade_data_pct, file = 'bTB_wildlife_runs/Fade_data_pct.RData')
toc = Sys.time() #takes about 4 min
toc-tic

tic = Sys.time()
Fade_data_spill = plotsData_reformat(version = vsn, type = c('fadeout','spillover.csv',''), runPath = "~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_initialInf_Runs/", scaled = F, scaleBy = NULL, herdSizes = N)
save( Fade_data_spill, file = 'bTB_wildlife_runs/Fade_data_spillover.RData')
toc = Sys.time() #takes about 4 min
toc-tic

tic = Sys.time()
Fade_data_pct_spill = plotsData_reformat(version = vsn, type = c('fadeout%','spillover.csv',''), runPath = "~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_initialInf_Runs/", scaled = F, scaleBy = NULL, herdSizes = N)
save( Fade_data_pct_spill, file = 'bTB_wildlife_runs/Fade_data_pct_spillover.RData')
toc = Sys.time() #takes about 4 min
toc-tic

#tic = Sys.time()
#Hunt_data = plotsData_reformat(type = c('hunt','seeded.csv',''), runPath = '~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_huntRates/', scaled = F, scaleBy = NULL, herdSizes = N)
#save( Hunt_data, file = 'bTB_wildlife_runs/Hunt_data.RData')
#toc = Sys.time() #takes about 10.5 min
#toc-tic

tic = Sys.time()
Inf_data = plotsData_reformat(version = vsn, type = c('multi','seeded.csv','Total_inf'), runPath = '~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Standard_Runs/', scaled = F, scaleBy = NULL, herdSizes = N)
save( Inf_data, file = 'bTB_wildlife_runs/Inf_data.RData')
toc = Sys.time() #takes about 1 min
toc-tic



tic = Sys.time()
SS_inf_data_spillover = plotsData_reformat(version = vsn, type = c('multi','spillover.csv','Total_inf'), runPath = '~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_SSprop/', scaled = F, scaleBy = NULL, herdSizes = N)
save( SS_inf_data_spillover, file = 'bTB_wildlife_runs/SS_inf_data_spillover.RData')
toc = Sys.time() #takes about 1 min
toc-tic

tic = Sys.time()
SS_inf_data_seeded = plotsData_reformat(version = vsn, type = c('multi','seeded.csv','Total_inf'), runPath = '~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_SSprop/', scaled = F, scaleBy = NULL, herdSizes = N)
save( SS_inf_data_seeded, file = 'bTB_wildlife_runs/SS_inf_data_seeded.RData')
toc = Sys.time() #takes about 1 min
toc-tic

load(file = 'bTB_wildlife_runs/Fade_data.RData')
load(file = 'bTB_wildlife_runs/Hunt_data.RData')
load(file = 'bTB_wildlife_runs/Inf_data.RData')
load(file = 'bTB_wildlife_runs/SS_inf_data_spillover.RData')
load(file = 'bTB_wildlife_runs/SS_inf_data_seeded.RData')
#multiPlot for infections over herd size
multiPlots(data = Inf_data, data_scaled = data.frame(), pairwiseComp = 'pop. size', pairwiseVec = N, category = 'total infections', infType = 'seeded', n_reps = reps,  n_years = years, scaled = F, regular = T, paired = F, col_vec = NA )
 
#multi plots over various SS proportions
multiPlots(data = SS_inf_data_spillover, data_scaled = data.frame(), pairwiseComp = 'SS proportion', pairwiseVec = Prop_SS, category = 'total infections', infType = 'spillover', n_reps = reps,  n_years = years, scaled = F, regular = T, paired = F, col_vec = NA )
#multi plots over various SS proportions
multiPlots(data = SS_inf_data_seeded, data_scaled = data.frame(), pairwiseComp = 'SS proportion', pairwiseVec = Prop_SS, category = 'total infections', infType = 'seeded', n_reps = reps,  n_years = years, scaled = F, regular = T, paired = F, col_vec = NA )


#multi plots over various SS proportions
multiPlots(data = Hunt_mortality_data, data_scaled = data.frame(), pairwiseComp = 'hunting rate', pairwiseVec = hunt_rates, category = 'total hunted animals', infType = 'seeded', n_reps = reps,  n_years = years, scaled = F, regular = T, paired = F, col_vec = NA )

  
fadeoutPlots(data = Fade_data, regular = T, fadeTime = T, paired = F, infType = 'seeded', n_reps = reps,  n_years = years, n_inf = N_inf0, herdSizes = herdSizes, col_vec = NA)
    
fadeoutPlots(data = Fade_data_pct, regular = T, fadeTime = T, paired = F, infType = 'seeded', n_reps = reps,  n_years = years, n_inf = N_inf0, herdSizes = herdSizes, col_vec = NA)


fadeoutPlots(data = Fade_data_spill, regular = T, fadeTime = T, paired = F, infType = 'spillover', n_reps = reps,  n_years = years, n_inf = N_inf0, herdSizes = herdSizes, col_vec = NA)

fadeoutPlots(data = Fade_data_pct_spill, regular = T, fadeTime = T, paired = F, infType = 'spillover', n_reps = reps,  n_years = years, n_inf = N_inf0, herdSizes = herdSizes, col_vec = NA)

#hunt_prop_plots(data = Hunt_data, pairedData = list(NA), regular = T, paired = F, infType = 'seeded', n_reps = reps,  n_years = years, par_vals = hunt_rates, herdSizes = herdSizes, col_vec = NA)
    
  




writePlot(data_vec = run$Total_inf, totalPop_vec = NULL, herdSize = "", infType = "", dataType = "", category = "", n_reps = 1000, n_years = 40, minMax = T, rep_Plots = 10, Scaled = F, Regular = T, Paired = F, as.jpeg = F, to.console = T, color = 'blue')
  
  