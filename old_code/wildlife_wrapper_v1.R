#parameter setting method#
#k takes positive integer values, and scenario is 1 or 2#
#scenario 1 has a constant force of infection whereas 2 has ?10? unreactive exposed individuals
parameter_set_wl <-function(k,scenario){
  
  parameter_names=c('S_0','E1_0','I_0','K','eta_hunt','eta_q1','eta_q2','eta_q3','eta_q4','alpha_q1','alpha_q2','alpha_q3','alpha_q4','beta','p1','p2_q1','p2_q2','p2_q3','p2_q4','sigma1_mean','sigma1_rate')
  #initialize parameters to be set
  K = k
  S_0 = k
  E1_0 = 0
  eta_q1 = .1
  eta_q2 = .1
  eta_q3 = .1
  eta_q4 = .1
  eta_hunt = .1
  alpha_q1=0
  alpha_q2=.05
  alpha_q3=0
  alpha_q4=0
  beta = 1
  p1 = .2
  p2_q1 = .05
  p2_q2 = .1
  p2_q3 = .01
  p2_q4 = .1
  sigma1_mean=3.5
  sigma1_rate=2.6
  
  if(scenario == 1){
    E1_0 = 10
    S_0 = S_0 - 10
    p2_q1 = 0
    p2_q2 = 0
    p2_q3 = 0
    p2_q4 = 0
  }
  #parameters not influenced by size, type, or infection introduction are hard coded in. All others set in loop.
  paramvec=c(S_0=S_0,E1_0=E1_0,I_0=0, 
             K=k, eta_hunt=eta_hunt, eta_q1=eta_q1, eta_q2=eta_q2, eta_q3=eta_q3, eta_q4=eta_q4, alpha_q1=alpha_q1, alpha_q2=alpha_q2, alpha_q3=alpha_q3, alpha_q4=alpha_q4, 
             beta=beta, p1=p1, p2_q1=p2_q1, p2_q2=p2_q2, p2_q3=p2_q3, p2_q4=p2_q4, sigma1_mean=sigma1_mean, sigma1_rate=sigma1_rate)
  
  return(paramvec)  #output vector of parameters for a run of given farm size type and infection scenario
}


###########################################################################################
#wrapper method code######################################################################
#########################################################################################

setwd("~/Documents/Brandon/bTB_wildlife_code/") #set working ddirectory to model code location
wildlife_model<-function(n_reps, parameters, initial_state, nyears, seed_quarter, batch_name){
  
  parvec = c( n_reps, #Number of replicates
              initial_state$S_0, #0, start Susc.
              initial_state$E1_0, #1, start E1
              initial_state$I_0, #2, start I
              parameters$K, #3, Herd size (county level)
              parameters$eta_hunt, #4, Hunting induced mortality rate
              parameters$eta_q1, #5, Natural mortality rate
              parameters$eta_q2, #6, Natural mortality rate
              parameters$eta_q3, #7, Natural mortality rate
              parameters$eta_q4, #8, Natural mortality rate
              parameters$alpha_q1, #9 Seasonal birth rate
              parameters$alpha_q2, #10 Seasonal birth rate
              parameters$alpha_q3, #11 Seasonal birth rate
              parameters$alpha_q4, #12 Seasonal birth rate
              parameters$beta, #13 Density Dependent Transmission rate
              parameters$p1, #14 deer-deer contact rate
              parameters$p2_q1, #15 C-to-wildl contact rate quarter 1
              parameters$p2_q2, #16 C-to-wildl contact rate quarter 2
              parameters$p2_q3, #17 C-to-wildl contact rate quarter 3
              parameters$p2_q4, #18 C-to-wildl contact rate quarter 4
              parameters$sigma1_mean, #19
              parameters$sigma1_rate, #20
              nyears, #21
              seed_quarter, #22
              batch_name #23
  )
  strvec = format(parvec, digits = 5)
  
  #The path to the bTB cpp binary file must be set correctly in the sys call below:
  #r <- system2(command = "/home/webblab/Documents/Brandon/bTB_wildlife_code/a.out", args=strvec, stdout=TRUE)
  r <- system2(command = "./wh_model_v1.exe", args=strvec, stdout=TRUE)
  res <- read.table(text=r, header=TRUE, sep=';', check.names = FALSE)
  return (res)
}

#setting run state variables#

reps=50
years=10
popsize=500
seedQuarter=1
infType=1

root="~/Documents/Brandon/bTB_wildlife_code/" #output file location

par<-parameter_set_wl(k=popsize,scenario=infType)

initial_state<-data.frame(S_0=par[1], E1_0=par[2], I_0=par[3])
population_parameters<-data.frame(K=par[4], eta_hunt=par[5], eta_q1=par[6], eta_q2=par[7], eta_q3=par[8], eta_q4=par[9], alpha_q1=par[10], alpha_q2=par[11], alpha_q3=par[12], alpha_q4=par[13]) 
disease_parameters<-data.frame(beta=par[14], p1=par[15], p2_q1=par[16], p2_q2=par[17], p2_q3=par[18], p2_q4=par[19], sigma1_mean=par[20], sigma1_rate=par[21])

param_vals<-data.frame(merge(population_parameters,disease_parameters))

#run#
run<-wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears=years, seed_quarter=1, batch_name="test")
fname='wildlife_test.csv'
write.csv(run, paste0(root,fname))

###testing###
#############
#n_reps = reps; parameters = param_vals; initial_state = initial_state; nyears=years; seed_quarter=1; batch_name="test";
#cbind(reps, initial_state, param_vals, years, 1, "test")
#############
