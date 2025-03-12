parameter_set_wl_LHS <-function(k, SS_prop = 0, scenario, SS_inf = F, inf_0 = 10){
  
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



LHS_ParSets <- function(FarmSizeClass, FarmTypeClass, InfTypeClass, Num_LHS_sets, Model_Version=2){
  
  setwd('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis2.0/')
  #### Latin Hypercube ####
  library(lhs)
  set.seed(8) #not sure if this is relevant
  h <- as.data.frame(optimumLHS(Num_LHS_sets, 34))
  lh <- as.data.frame(h)
  
  sig1=(.21/9.75)+1
  sig2=(1.85/1.6)+1
  
  ParRangeVec=bTB_LHS_Parameters(FarmSizeClass, FarmTypeClass, InfTypeClass) #set parameters dependent 
  #on farm type, size, and infection introduction type.
  
  #parameter ranges are determined by taking 10 * the default value for the parameter. 
  
  #farm factors
  lh[,1] <- ceiling(qunif(h[,1], ParRangeVec[1], ParRangeVec[2])) # K
  lh[,2] <- qunif(h[,2], 0, 10*ParRangeVec[4]) # mu *use 10 as max for beef
  lh[,3] <- qunif(h[,3], 0, 10*ParRangeVec[5]) # q4_mu
  lh[,4] <- qunif(h[,4], 0, 10*ParRangeVec[8]) #r *use 0 as max for beef
  lh[,5] <- qunif(h[,5], 0, 10*ParRangeVec[6]) # nu *use .2 as max for beef
  lh[,6] <- qunif(h[,6], 0, 10*ParRangeVec[7]) # alpha *use 12 as max for beef
  
  #transmission factors
  #pi and beta depend on farm type
  lh[,7] <- qunif(h[,7], 0, 10*ParRangeVec[9]) # beta
  lh[,8] <- qunif(h[,8], 0, 2) # p1
  lh[,9] <- qunif(h[,9], 0, 10*ParRangeVec[11]) # p2_q1
  lh[,10] <- qunif(h[,10], 0, 10*ParRangeVec[12]) # p2_q2
  lh[,11] <- qunif(h[,11], 0, 10*ParRangeVec[13]) # p2_q3
  lh[,12] <- qunif(h[,12], 0, 10*ParRangeVec[14]) # p2_q4
  lh[,13] <- qunif(h[,13], 0, 10*ParRangeVec[10]) # p3
  lh[,14] <- qunif(h[,14], 0, 10*ParRangeVec[15]) # pi
  lh[,15] <- qunif(h[,15], 0, 1) # psi1
  lh[,16] <- qunif(h[,16], 0, 1) # psi2
  lh[,17] <- qunif(h[,17], 0, 1) # psi3
  lh[,18] <- qunif(h[,18], 0, 1) # psi4
  if(Model_Version==1){
    lh[,18]=replicate(Num_LHS_sets,0)
  }
  
  #latency factors
  #always constant
  lh[,19] <- qunif(h[,19], 0, 2) # phi_sigma_1
  lh[,20] <- qunif(h[,20], 0, 2) # phi_sigma_2
  lh[,21] <- qunif(h[,21], 0, 2) # phi_delta1
  lh[,22] <- qunif(h[,22], 0, 2) # phi_delta2
  lh[,23] <- qunif(h[,24], 0, 10*sig1 ) # sigma1
  lh[,24] <- qunif(h[,24], 0, 2.1) # sigma1_rate
  lh[,25] <- qunif(h[,25], 0, 10*sig2) # sigma2
  lh[,26] <- qunif(h[,26], 0, 18.5) # sigma2_rate
  
  #testing factors
  #always constant
  lh[,27] <- qunif(h[,27], 0, 1) # pa_0
  if(Model_Version==1){
    lh[,27]=replicate(Num_LHS_sets,0)
  }
  lh[,28] <- qunif(h[,28], 0, 1) # pa_1
  lh[,29] <- qunif(h[,29], 0, 1) # pa_2
  lh[,30] <- qunif(h[,30], 0, 1) # pa_3
  lh[,31] <- qunif(h[,31], 0, 1) # t1
  lh[,32] <- qunif(h[,32], 0, 1) # t2
  lh[,33] <- qunif(h[,33], 0, 1) # pa_t2
  
  #Seed#
  #always constant
  lh[,34] <- ceiling(qunif(h[,34], 0, 4)) #seed_q
  
  lh <- round(lh,6)
  
  
  farm_factors<-c('K','mu','q4_mu','r','nu','alpha')
  transmission_factors<-c('beta','p1','p2_q1','p2_q2','p2_q3','p2_q4','p3','pi','psi1','psi2','psi3','psi4')
  latency_factors<-c('phi_sigma_1','phi_sigma_2','phi_delta1','phi_delta2','sigma1','sigma1_rate','sigma2','sigma2_rate')
  testing_factors<-c('pa_0','pa_1','pa_2','pa_3','t1','t2','pa_t2')
  other_factors<-c('seed_q')
  
  LHpars=as.data.frame(lh) #data frame without column names in case there is an issue with reading in params
  colnames(lh) <- c(farm_factors,transmission_factors,latency_factors,testing_factors,other_factors)
  
  head(lh)
  LHS_par_table=paste0('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis',Model_Version,'.0/',FarmSizeClass, '_', FarmTypeClass, '_', InfTypeClass, '_LHS_Control_Parameters.csv');
  LHS_par_sumary=paste0('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis',Model_Version,'.0/',FarmSizeClass, '_', FarmTypeClass, '_', InfTypeClass, '_LHS_Control_Param_Summary.csv');
  #write.csv(x = lh, file = LHS_par_table, row.names = FALSE) #write out parameter data file
  lh.summary <- as.data.frame(round(apply(lh, 2, function(x) summary(x)),4)) #parameter summary file
  #write.csv(x = lh.summary, file = LHS_par_sumary) #write out above
  return(lh)
}

####################################################
#Function to run model with each LHS parameter set#
###################################################

bTB_LHS_Runs <- function(FarmSizeClass, FarmTypeClass, InfTypeClass, nreps){
  
  #compile ith replicate into appropriate format and run model with each set of LHS parameters
  setwd('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis2.0/')
  LHS_Parameters <- read.csv(file=paste0(FarmSizeClass, '_', FarmTypeClass, '_', InfTypeClass, '_LHS_Control_Parameters.csv'))
  num_LHS_parsets=nrow(LHS_Parameters)
  repCount=replicate(num_LHS_parsets,0)
  fname=''
  for(n in 1:num_LHS_parsets){
    E1U=0
    if(InfTypeClass==2){
      E1U=10
    }
    
    initial_state<-data.frame(S0=LHS_Parameters[n,1],E1U0=E1U,E1R0=0,E2R0=0,E2U0=0,IU0=0,IR0=0) 
    farm_parameters<-data.frame(K=LHS_Parameters[n,1],mu=LHS_Parameters[n,2],q4_mu=LHS_Parameters[n,3],r=LHS_Parameters[n,4],nu=LHS_Parameters[n,5],alpha=LHS_Parameters[n,6]) 
    transmission_parameters<-data.frame(beta=LHS_Parameters[n,7],p1=LHS_Parameters[n,8],p2_q1=LHS_Parameters[n,9],p2_q2=LHS_Parameters[n,10],p2_q3=LHS_Parameters[n,11],p2_q4=LHS_Parameters[n,12],p3=LHS_Parameters[n,13],pi=LHS_Parameters[n,14],psi1=LHS_Parameters[n,15],psi2=LHS_Parameters[n,16],psi3=LHS_Parameters[n,17],psi4=LHS_Parameters[n,18])
    latency_parameters<-data.frame(phi_sigma_1=LHS_Parameters[n,19],phi_sigma_2=LHS_Parameters[n,20],phi_delta1=LHS_Parameters[n,21],phi_delta2=LHS_Parameters[n,22],sigma1=LHS_Parameters[n,23],sigma1_rate=LHS_Parameters[n,24],sigma2=LHS_Parameters[n,25],sigma2_rate=LHS_Parameters[n,26])
    testing_parameters<-data.frame(pa_0=LHS_Parameters[n,27],pa_1=LHS_Parameters[n,28],pa_2=LHS_Parameters[n,29],pa_3=LHS_Parameters[n,30],t1=LHS_Parameters[n,31],t2=LHS_Parameters[n,32],pa_t2=LHS_Parameters[n,33])
    
    parameters<-data.frame(merge(merge(farm_parameters,transmission_parameters),merge(latency_parameters,testing_parameters)))
    
    LHrun=withinherd_model(nreps, parameters, initial_state, nyears=3, farm_type=FarmTypeClass, seed_quarter=LHS_Parameters[n,34], batch_name="test")
    
    root=paste0('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis2.0/')
    fname=paste0('LHS_',n,'_',FarmSizeClass, '_',FarmTypeClass, '_',InfTypeClass,'.csv')
    source=paste0(root,fname)
    write.csv(LHrun, source)
    repCount[n]=fname
  }
  
}

###########################################
#Generate LHS parameters for all scenarios#
###########################################
inf=c(1,2) #i
farms=c('beef', 'dairy') #j
sizes=c('small', 'medium','large') #k
model_version=1
for(i in 1:2){
  for(j in 1:2){
    for(k in 1:3){
      setwd(paste0('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis',model_version,'.0/'))
      LHS_ParSets(FarmSizeClass=sizes[k], FarmTypeClass=farms[j], InfTypeClass=inf[i], Num_LHS_sets=200,Model_Version = model_version)
    }
  }
}

#######################################
#Run model for all LHS parameter sets#
######################################
inf=c(1,2) #i
farms=c('beef', 'dairy') #j
sizes=c('small', 'medium','large') #k
for(i in 1:2){
  for(j in 1:2){
    for(k in 1:3){
      bTB_LHS_Runs(FarmSizeClass=sizes[k], FarmTypeClass=farms[j], InfTypeClass=inf[i], nreps=200)
    }
  }
}



for(i in 1:2){ #mode of infection
  for(j in 1:2){ #farm type
    #Plots to check monotonic assumptions
    #independent variables are: "K", "mu", "q4_mu", "nu", "alpha", "r", 
    #"beta", "p1", "p3", "p2_q1", "p2_q2", "p2_q3", "p2_q4", "pi", "psi1", "psi2", "psi3", 
    #"sigma1", "sigma1_rate", "sigma2", "sigma2_rate", "phi", "phi_delta1", "phi_delta2", 
    #"pa_1", "pa_2", "pa_3", "pa_t2", "t1", "t2", 
    #"seed_q"
    #dependent outbreak variables are: "True_Prevalence" "Observed_Prevalence" "Difference" "fade_out"
    
    library(gridExtra)
    setwd(paste0('/home/webblab/Documents/Brandon/bTB_WH/Runs/Reruns/Sensitivity_Analysis2.0/LHS_Monotonic_Plots/',farms[j],'_',inf[i],'/')) 
    names(bTB_WH_scaled)[names(bTB_WH_scaled) == "EndTinf"] <- "TruPrev"
    names(bTB_WH_scaled)[names(bTB_WH_scaled) == "EndObs"] <- "ObsPrev"
    names(bTB_WH_scaled)[names(bTB_WH_scaled) == "EndDiff"] <- "Diff"
    names(bTB_WH_scaled)[names(bTB_WH_scaled) == "FOT"] <- "FOtime"
    names(bTB_WH_scaled)[names(bTB_WH_scaled) == "Fadeout"] <- "fadeout"
    
    ########################
    #K PLOTS WITH EACH OUTBREAK VARIABLE
    jpeg(filename = "monotonic_K.jpeg", width =840, height = 840, units = 'px', res = 100)
    
    p1=ggplot(bTB_WH_scaled, aes(x=K, y=TruPrev)) + 
      geom_point()
    
    p2=ggplot(bTB_WH_scaled, aes(x=K, y=ObsPrev)) + 
      geom_point()
    
    p3=ggplot(bTB_WH_scaled, aes(x=K, y=Diff)) + 
      geom_point()
    
    p4=ggplot(bTB_WH_scaled, aes(x=K, y=FOtime)) + 
      geom_point()
    
    p5=ggplot(bTB_WH_scaled, aes(x=K, y=fadeout)) + 
      geom_point()
    grid.arrange(p1,p2,p3,p4,p5, nrow=2)
    
    dev.off()
    
    
    
    #mu PLOTS WITH EACH OUTBREAK VARIABLE
    jpeg(filename = "monotonic_mu.jpeg", width =840, height = 840, units = 'px', res = 100)
    
    p1=ggplot(bTB_WH_scaled, aes(x=mu, y=TruPrev)) + 
      
      geom_point()
    
    p2=ggplot(bTB_WH_scaled, aes(x=mu, y=ObsPrev)) + 
      geom_point()
    
    p3=ggplot(bTB_WH_scaled, aes(x=mu, y=Diff)) + 
      geom_point()
    
    p4=ggplot(bTB_WH_scaled, aes(x=mu, y=FOtime)) + 
      geom_point()
    
    p5=ggplot(bTB_WH_scaled, aes(x=mu, y=fadeout)) + 
      geom_point()
    grid.arrange(p1,p2,p3,p4,p5, nrow=2)
    
    dev.off()
    
    
    
    #q4_mu PLOTS WITH EACH OUTBREAK VARIABLE
    jpeg(filename = "monotonic_q4_mu.jpeg", width =840, height = 840, units = 'px', res = 100)
    
    p1=ggplot(bTB_WH_scaled, aes(x=q4_mu, y=TruPrev)) + 
      geom_point()
    
    
    p2=ggplot(bTB_WH_scaled, aes(x=q4_mu, y=ObsPrev)) + 
      geom_point()
    
    
    p3=ggplot(bTB_WH_scaled, aes(x=q4_mu, y=Diff)) + 
      geom_point()
    
    p4=ggplot(bTB_WH_scaled, aes(x=q4_mu, y=FOtime)) + 
      geom_point()
    
    p5=ggplot(bTB_WH_scaled, aes(x=q4_mu, y=fadeout)) + 
      geom_point()
    grid.arrange(p1,p2,p3,p4,p5, nrow=2)
    
    dev.off()
    
  }
}