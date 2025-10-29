library(deSolve);
library(tidyverse);
library(ggplot2);
library(ggpubr);
library(grid);
library(rstudioapi);
library(foreach);
library(doParallel);
library(RColorBrewer);
library(scales);
library(lhs)
library(reshape2)
library(stringr)
##################################################
#### time & population varied model functions ####
##################################################
# seasonal birth pulse function
alpha <- function(Amp = 1, synchrony = 0, phase = 0, start_q = 1, t){ 
  return( Amp*exp( -synchrony*( cos((pi/12)*( t + 3*(start_q-1) - phase) ) )^2 ) ) 
}

# explicit mod 12 function for sapply usage
mod12 <- function(a){ return(a%%12) }

# function for natural and density dependent mortality
eta <- function(baseline_mortality = 0, theta = 1, y_scale = 1, N_frac){
  return( baseline_mortality + y_scale*(N_frac^theta) ) # sum of mortality rates, baseline scaled upward if population exceeds carrying capacity
}

# seasonal function for hunter harvest
eta_hunt <- function(hunt_mortality = 0, start_q = 1, t){
  t <- t + ((start_q-1)*3) 
  t <- t %>% sapply(floor) %>% sapply(FUN = mod12) # shift function input t to correspond with start quarter, then convert t into month
  Ind = as.numeric(t %in% c(9, 10, 11)) # indicator for quarter = 4 - note indexing starts at 0
  return( Ind*hunt_mortality ) # sum of mortality rates, baseline scaled upward if population exceeds carrying capacity
}

# seasonal farm contact rates
p2 <- function(probs = c(0,0,0,0), start_q = 1, t){
  t <- t + ((start_q-1)*3) 
  t <- t %>% sapply(floor) %>% sapply(FUN = mod12) # shift function input t to correspond with start quarter, then convert t into month
  t <- floor(t/3) + 1 # translate month to quarter, add 1 for R indexing
  return(probs[t])
}

# parameter set method - should be consistent between ODE and C++ model input requirements
parameter_set_wl <-function(k = 10, 
                            SS_prop = 0, 
                            scenario = "spillover", 
                            initial_exposed = 0, 
                            start_quarter = 1, 
                            test = F,
                            birth = F,
                            death_h = F,
                            death_n = F,
                            disease = F,
                            verbose = 0){
  
  #initialize parameters to be set
  start_q = 1
  
  # mortality parameters
  K = k # carrying capacity
  eta_nat = 0.05/12 # natural mortality rate
  eta_hunt = - (1/3) * log(1 - 0.2)
  # eta_hunt = (.2*k)/3 # seasonal hunter harvest -- set to depend on county population size and start pop. size
  theta = 2.0 # density dependent mortality asymmetry
  gamma = 2.0*eta_nat
  
  # birth pulse parameters
  alpha_max = 1.0 # amplitude 
  ksi = SS_prop # proportion of newborn super spreaders
  omega = 10.3 # phase shift
  s = 68 # synchrony factor
  alpha_func = function(x){alpha_max*exp(-s*(cos(pi*(x-omega)/12))^2)}
  alpha_annual <- integrate(alpha_func, lower = 0, upper = 12)$value
  
  # transmission / transition parameters
  beta = 0.025 # = 0.05*0.5 #### 0.025 # 0.08333333333 # .5 # transmission rate
  A = 1 # occupiable habitat for density dependent transmission
  p1 = 1 # 0.05 probability of contact sufficient for transmission
  # Need to divide these by 0.05 so we can get rid of deer-deer contact when multiply by tranmission rate in the model
  p2_q1 = .01/0.05 # probability of infected farm contact
  p2_q2 = .015/0.05 # probability of infected farm contact
  p2_q3 = 0/0.05 # probability of infected farm contact
  p2_q4 = .0075/0.05 # probability of infected farm contact 
  # p2_q1 = .005 # probability of infected farm contact
  # p2_q2 = .01 # probability of infected farm contact
  # p2_q3 = .001 # probability of infected farm contact
  # p2_q4 = .01 # probability of infected farm contact 
  phi = 5 # farm contact scaling factor
  sigma1_mean = 12*5 # mean latency period duration
  sigma1_rate = 0.6 # rate parameter for gamma distn.
  
  # parameter flex by scenario
  # assume spillover infection is the default mode
  E0=ceiling(initial_exposed)
  if(scenario %in% c('seeded', 'Seeded')){
    if(initial_exposed>k){E0=k}
    E1_0 = E0
    S_0 = k - E0
    p2_q1 = 0
    p2_q2 = 0
    p2_q3 = 0
    p2_q4 = 0
  }
  
  else if(scenario %in% c('spillover','Spillover')){
    if(initial_exposed>k){E0=k}
    E1_0 = E0
    S_0 = k - E0
    # no other parameters to reset since this is the default
  }
  else{print(paste0('Warning: Undefined scenario selection: ', scenario))
    print('Options include: seeded, spillover')
    break
  }
  
  # initial population parameters
  nSS=floor(S_0*SS_prop)
  S_0 = S_0 - nSS
  SuperS_0 = nSS
  
  if(start_quarter %in% c(1,2,3,4)){start_q = start_quarter}
  
  #adjust appropriate parameters if test mode is enabled
  if(test == T){
    if(birth == F){
      alpha_max = 0 # amplitude 
      alpha_func = function(x){alpha_max*exp(-s*(cos(pi*(x-omega)/12))^2)}
      alpha_annual <- integrate(alpha_func, lower = 0, upper = 12)$value
    }
    if(death_h == F){
      eta_hunt = 0 # seasonal hunter harvest 
    }
    if(death_n == F){
      eta_nat = 0 # natural mortality rate
      gamma = 0
    }
    if(disease == F){
      beta = 0 # transmission rate
      sigma1_mean = 0 # mean latency period duration
      sigma1_rate = 0 # rate parameter for gamma distn.
    }
  }
  
  #parameters not influenced by size, type, or infection introduction are hard coded in. All others set in loop.
  paramvec=c(S_0 = S_0, E1_0 = E1_0, I_0 = 0, SuperS_0 = nSS, SuperE1_0 = 0, SuperI_0 = 0,
             K = K, eta_hunt = eta_hunt, eta_nat = eta_nat, theta = theta, gamma = gamma,
             alpha_max = alpha_max, ksi = ksi, omega = omega, s = s, alpha = alpha_annual,
             beta = beta,  area = A, p1 = p1, p2_q1 = p2_q1, p2_q2 = p2_q2, p2_q3 = p2_q3, p2_q4 = p2_q4, phi = phi, sigma1_mean = sigma1_mean, sigma1_rate = sigma1_rate, 
             start_q = start_q, verbose = verbose)
  
  return(paramvec)  #output vector of parameters for a run of given farm size type and infection scenario
}
# latin hypercube parameter generator
LHS_ParSets <- function(k_lim = c(10,100), 
                        scenario = "spillover", 
                        infected = 10,
                        pct = F,
                        SS_prop = 0.1,
                        verbose = 1,
                        Num_LHS_sets,
                        file.path = NA,
                        file.out = F,
                        file.name = "bTB_wl_LHS_pars",
                        seed_q = NA,
                        seed = 43){
  
  pars=parameter_set_wl(k = k_lim[1], scenario = scenario, SS_prop = SS_prop, verbose = verbose) #set parameters dependent 
  # drop <- c("S_0", "E1_0", "I_0", "SuperS_0", "SuperE1_0", "SuperI_0", 'alpha', 'p1', 'verbose', 'area')
  drop <- c("S_0", "E1_0", "I_0", "SuperS_0", "SuperE1_0", "SuperI_0", 'alpha', 'verbose', 'area')
  pars <- as.data.frame(t(pars)) 
  pars <- pars[,!(names(pars) %in% drop)]
  
  #### Latin Hypercube ####
  library(lhs)
  set.seed(seed) 
  h <- as.data.frame(improvedLHS(Num_LHS_sets, length(pars))) 
  names(h) <- names(pars) # apply parameter names to data columns
  # h <- select(h, -c("S_0", "E1_0", "I_0", "SuperS_0", "SuperE1_0", "SuperI_0"))
  lh <- as.data.frame(h)
  
  #parameter ranges are determined by taking 10 * the default value for the parameter. 
  
  #mortality parameters
  lh[,"K"] <- floor(qunif(h[,"K"], k_lim[1], k_lim[2]))
  lh[,"eta_hunt"] <- qunif(h[,"eta_hunt"], 0, 1)
  lh[,"eta_nat"] <- qunif(h[,"eta_nat"], 0, 10*as.numeric(pars["eta_nat"]))
  lh[,"theta"] <- qunif(h[,"theta"], 0, 10*as.numeric(pars["theta"]))
  lh[,"gamma"] <- qunif(h[,"gamma"], 0, 10*as.numeric(pars["gamma"]))
  
  #birth parameters
  lh[,"alpha_max"] <- qunif(h[,"alpha_max"], .9*as.numeric(pars["alpha_max"]), 1.1*as.numeric(pars["alpha_max"]))
  lh[,"ksi"] <- qunif(h[,"ksi"], 0, 1)
  lh[,"omega"] <- qunif(h[,"omega"], .9*as.numeric(pars["omega"]), 1.1*as.numeric(pars["omega"]))
  lh[,"s"] <- qunif(h[,"s"], .9*as.numeric(pars["s"]), 1.1*as.numeric(pars["s"]))

  
  #disease transition factors
  lh[,"beta"] <- qunif(h[,"beta"], 0, 10*as.numeric(pars["beta"]))

  lh[,"p2_q1"] <- qunif(h[,"p2_q1"], 0, 1)
  lh[,"p2_q2"] <- qunif(h[,"p2_q2"], 0, 1)
  lh[,"p2_q3"] <- qunif(h[,"p2_q3"], 0, 1)
  lh[,"p2_q4"] <- qunif(h[,"p2_q4"], 0, 1)
  
  lh[,"p1"] <- 1 # qunif(h[,"p1"], 0, 1)
  
  #set these to 0 instead if seeded
  if(scenario == 'seeded'){
    lh[,"p2_q1"] <- pars["p2_q1"]
    lh[,"p2_q2"] <- pars["p2_q2"]
    lh[,"p2_q3"] <- pars["p2_q3"]
    lh[,"p2_q4"] <- pars["p2_q4"]
  }
  lh[,"phi"] <- qunif(h[,"phi"], 0, 10*as.numeric(pars["phi"]))
  lh[,"sigma1_mean"] <- qunif(h[,"sigma1_mean"], 0, 10*as.numeric(pars["sigma1_mean"]))
  lh[,"sigma1_rate"] <- qunif(h[,"sigma1_rate"], 0, 10*as.numeric(pars["sigma1_rate"]))
  
  #Seed#
  lh[,"start_q"] <- ceiling(qunif(h[,"start_q"], 0, 4)) #seed_q
  if(!is.na(seed_q)){ lh[,"start_q"] <- seed_q}
  ## determine number of animals, starting at carrying capacity
  
  lh[,"S_0"] <- lh[,"K"]
  #determine infected individuals
  if(pct == T){
    lh[,"E1_0"] <- ceiling(lh[,"S_0"]*infected)
    lh[,"S_0"] <- lh[,"S_0"] - lh[,"E1_0"]
  }else{
    lh[,"E1_0"] <- as.integer(infected)
    lh[,"S_0"] <- lh[,"S_0"] - lh[,"E1_0"]
  }
  
  #determine super spreaders
  lh[,"SuperS_0"] <- floor(lh[,"ksi"]*lh[,"S_0"])
  lh[,"S_0"] <- lh[,"S_0"] - lh[,"SuperS_0"]
  
  #assign other classes to 0
  lh[,"SuperE1_0"] <- 0
  lh[,"I_0"] <- 0
  lh[,"SuperI_0"] <- 0
  #####
  if(sum(lh[,"E1_0"]<0) > 0 | 
     sum(lh[,"SuperS_0"]<0) > 0 | 
     sum(lh[,"S_0"]<0) > 0){
    print("negative value detected - assigning a value of 0")
    lh[which(lh[,"E1_0"]<0),"E1_0"] <- 0
    lh[which(lh[,"SuperS_0"]<0),"SuperS_0"] <- 0
    lh[which(lh[,"S_0"]<0),"S_0"] <- 0
  }
  lh <- round(lh,6)
  
  LHpars=as.data.frame(lh) #data frame without column names in case there is an issue with reading in params
  
  # output parameters and summary file
  if(file.out == T){
    # save file to path if supplied, otherwise write out to the current working directory
    if(is.na(file.path)){
      write.csv(x = lh, file = paste0(file.name, '.csv'), row.names=FALSE)
      lh.summary <- as.data.frame(round(apply(lh, 2, function(x) summary(x)),4)) #parameter summary file
      write.csv(x = lh.summary, file = paste0(file.name, "_summary.csv"), row.names=FALSE)
    }else{
      write.csv(x = lh, file = paste0(file.path, file.name, '.csv'), row.names=FALSE)
      lh.summary <- as.data.frame(round(apply(lh, 2, function(x) summary(x)),4)) #parameter summary file
      write.csv(x = lh.summary, file = paste0(file.path, file.name, "_summary.csv"), row.names=FALSE)
    }
  }
  
  return(lh)
}

SEI_model_full <- function (t, x, parms) {
  
  S  <- as.numeric(x["S"])
  E  <- as.numeric(x["E"])
  I  <- as.numeric(x["I"])
  sS <- as.numeric(x["sS"])
  sE <- as.numeric(x["sE"])
  sI <- as.numeric(x["sI"])
  
  # births
  Alpha_max   <- as.numeric(parms["alpha_max"])     # amplitude
  ksi         <- as.numeric(parms["ksi"])           # fraction of newborn SS
  omega       <- as.numeric(parms["omega"])         # phase
  s           <- as.numeric(parms["s"])             # synchrony
  alpha_annual<- as.numeric(parms["alpha"])         # annual integral (unused directly here)
  
  # mortality
  K     <- as.numeric(parms["K"])
  eta_n <- as.numeric(parms["eta_nat"])
  eta_h <- as.numeric(parms["eta_hunt"])
  theta <- as.numeric(parms["theta"])
  gamma <- as.numeric(parms["gamma"])
  
  # transmission / transition
  beta_wild <- as.numeric(parms["beta"])            # interpret 'beta' as deer–deer rate (month^-1)
  A         <- as.numeric(parms["area"])            # density-dependent scaling area
  p1        <- as.numeric(parms["p1"])              # kept for interface compatibility (unused)
  p2_vec    <- c(as.numeric(parms["p2_q1"]),        # unitless seasonal multipliers
                 as.numeric(parms["p2_q2"]),
                 as.numeric(parms["p2_q3"]),
                 as.numeric(parms["p2_q4"]))
  phi       <- as.numeric(parms["phi"])             # SS farm-contact multiplier
  sigma     <- 1 / as.numeric(parms["sigma1_mean"])
  if (as.numeric(parms["sigma1_mean"]) == 0) sigma <- 0
  
  Q <- as.numeric(parms["start_q"])                 # start quarter (1..4)
  
  N <- S + E + I + sS + sE + sI
  
  # convert beta_wild to beta_farm by multiiplying by beta_wild*p_q
  beta_farm_t <- function(tt) p2(probs = p2_vec, start_q = Q, t = tt) * beta_wild
  
  dSdt  <- (1 - ksi) * alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t) * N
  dSdt  <- dSdt - beta_wild * S  * (I + sI) / A               # deer–deer
  dSdt  <- dSdt - beta_farm_t(t) * S                          # farm
  dSdt  <- dSdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * S
  dSdt  <- dSdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * S # (S / N)
  
  dEdt  <-  beta_wild * S  * (I + sI) / A                     # deer–deer (S)
  dEdt  <- dEdt + beta_farm_t(t) * S                          # farm (S)
  dEdt  <- dEdt + beta_wild * sS * (I + sI) / A               # deer–deer (SS)
  dEdt  <- dEdt + beta_farm_t(t) * phi * sS                   # farm (SS boosted)
  dEdt  <- dEdt - sigma * E
  dEdt  <- dEdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * E
  dEdt  <- dEdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * E #  (E / N)
  
  dIdt  <-  sigma * E
  dIdt  <- dIdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * I
  dIdt  <- dIdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * I # (I / N)
  
  dsSdt <-  ksi * alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t) * N
  dsSdt <- dsSdt - beta_wild * sS * (I + sI) / A              # deer–deer
  dsSdt <- dsSdt - beta_farm_t(t) * phi * sS                  # farm (SS boosted)
  dsSdt <- dsSdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * sS
  dsSdt <- dsSdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * sS # (sS / N)
  
  dsEdt <-  beta_wild * sS * (I + sI) / A                     # deer–deer (SS)
  dsEdt <- dsEdt + beta_farm_t(t) * phi * sS                  # farm (SS boosted)
  dsEdt <- dsEdt - sigma * sE
  dsEdt <- dsEdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * sE
  dsEdt <- dsEdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * sE # (sE / N)
  
  dsIdt <-  sigma * sE
  dsIdt <- dsIdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * sI
  dsIdt <- dsIdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * sI # (sI / N)
  
  list(c(dSdt, dEdt, dIdt, dsSdt, dsEdt, dsIdt))
}


# # full deterministic model structures - with super spreaders
# SEI_model_full <- function (t, x, parms) {
#   
#   ## first extract the state variables
#   S <- as.numeric(x["S"])
#   E <- as.numeric(x["E"])
#   I <- as.numeric(x["I"])
#   sS <- as.numeric(x["sS"])
#   sE <- as.numeric(x["sE"])
#   sI <- as.numeric(x["sI"])
#   
#   ## now extract the parameters
#   
#   # birth pulse parameters
#   Alpha_max <- as.numeric(parms["alpha_max"]) # birth rate amplitude
#   ksi <- as.numeric(parms["ksi"]) # proportion of newborn super spreaders
#   omega <- as.numeric(parms["omega"]) # birth rate phase shift
#   s <- as.numeric(parms["s"]) # birth rate synchrony
#   alpha_annual = as.numeric(parms["alpha"]) # annual birth rate value
#   
#   # mortality rate parameters
#   K <- as.numeric(parms["K"]) # environment carrying capacity / initial pop. size
#   eta_n <- as.numeric(parms["eta_nat"]) # baseline natural mortality rate
#   eta_h <- as.numeric(parms["eta_hunt"]) # hunt mortality rate
#   theta <- as.numeric(parms["theta"]) # scale disproportionality in density dependent mortality
#   gamma <- as.numeric(parms["gamma"])
#   #gamma <- alpha_annual/12 - eta_n # density dependent mortality scalar
#   
#   
#   # disease state parameters
#   beta <- as.numeric(parms["beta"]) #transmission rate
#   A <- as.numeric(parms["area"])
#   p1 <- as.numeric(parms["p1"]) # infectious contact prob
#   p2_vec <- c(as.numeric(parms["p2_q1"]), as.numeric(parms["p2_q2"]), as.numeric(parms["p2_q3"]), as.numeric(parms["p2_q4"])) # infected farm contact prob.
#   phi = as.numeric(parms["phi"]) # farm contact scaling factor for super spreaders
#   sigma <- 1 / as.numeric(parms["sigma1_mean"]) # mean rate of transition to infectious
#   if(as.numeric(parms["sigma1_mean"]) == 0){sigma = 0} # prevent infinite transition rate
#   
#   Q <- as.numeric(parms["start_q"]) # start quarter
#   
#   # define total population
#   N <- S + E + I + sS + sE + sI
#   
#   dSdt <- (1-ksi)*alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t)*N - p1*beta*S*(I+sI)/A - p2(probs = p2_vec, start_q = Q, t)*beta*S - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K)*S - eta_hunt(hunt_mortality = eta_h, start_q = Q, t)*(S/N)  
#   
#   dEdt <- p1*beta*S*(I+sI)/A + p2(probs = p2_vec, start_q = Q, t)*beta*S - sigma*E - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K)*E - eta_hunt(hunt_mortality = eta_h, start_q = Q, t)*(E/N)  
#   
#   dIdt <- sigma*E - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K)*I - eta_hunt(hunt_mortality = eta_h, start_q = Q, t)*(I/N)  
#   
#   dsSdt <- (ksi)*alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t)*N - p1*beta*sS*(I+sI)/A - p2(probs = p2_vec, start_q = Q, t)*phi*beta*sS - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K)*sS - eta_hunt(hunt_mortality = eta_h, start_q = Q, t)*(sS/N)  
#   
#   dsEdt <- p1*beta*sS*(I+sI)/A + p2(probs = p2_vec, start_q = Q, t)*phi*beta*sS - sigma*sE - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K)*sE - eta_hunt(hunt_mortality = eta_h, start_q = Q, t)*(sE/N)  
#   
#   dsIdt <- sigma*sE - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K)*sI - eta_hunt(hunt_mortality = eta_h, start_q = Q, t)*(sI/N)  
#   
#   
#   ## combine results into a single vector
#   dxdt <- c(dSdt, dEdt, dIdt, dsSdt, dsEdt, dsIdt)
#   ## return result as a list!
#   list(dxdt)
# }

# # Frequency-dependent SEI model with super spreaders
# SEI_model_full_freq <- function(t, x, parms) {
#   
#   # State variables
#   S  <- as.numeric(x["S"])
#   E  <- as.numeric(x["E"])
#   I  <- as.numeric(x["I"])
#   sS <- as.numeric(x["sS"])
#   sE <- as.numeric(x["sE"])
#   sI <- as.numeric(x["sI"])
#   
#   # Parameters
#   Alpha_max <- as.numeric(parms["alpha_max"])
#   ksi       <- as.numeric(parms["ksi"])
#   omega     <- as.numeric(parms["omega"])
#   s         <- as.numeric(parms["s"])
#   alpha_annual <- as.numeric(parms["alpha"])
#   
#   K      <- as.numeric(parms["K"])
#   eta_n  <- as.numeric(parms["eta_nat"])
#   eta_h  <- as.numeric(parms["eta_hunt"])
#   theta  <- as.numeric(parms["theta"])
#   gamma  <- as.numeric(parms["gamma"])
#   
#   beta   <- as.numeric(parms["beta"])
#   A      <- as.numeric(parms["area"])
#   p1     <- as.numeric(parms["p1"])
#   p2_vec <- c(as.numeric(parms["p2_q1"]), as.numeric(parms["p2_q2"]),
#               as.numeric(parms["p2_q3"]), as.numeric(parms["p2_q4"]))
#   phi    <- as.numeric(parms["phi"])
#   sigma  <- 1 / as.numeric(parms["sigma1_mean"])
#   if(as.numeric(parms["sigma1_mean"]) == 0) sigma <- 0
#   
#   Q <- as.numeric(parms["start_q"])
#   
#   # Total population
#   N <- S + E + I + sS + sE + sI
#   
#   # Differential equations (frequency dependent transmission)
#   dSdt  <- (1 - ksi) * alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t) * N -
#     p1 * beta * S * ((I + sI) / N) -
#     p2(probs = p2_vec, start_q = Q, t) * beta * S * (1 / N) -
#     eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K) * S -
#     eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (S/N)
#   
#   dEdt  <- p1 * beta * S * ((I + sI) / N) +
#     p2(probs = p2_vec, start_q = Q, t) * beta * S * (1 / N) -
#     sigma * E -
#     eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K) * E -
#     eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (E/N)
#   
#   dIdt  <- sigma * E -
#     eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K) * I -
#     eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (I/N)
#   
#   dsSdt <- ksi * alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t) * N -
#     p1 * beta * sS * ((I + sI) / N) -
#     p2(probs = p2_vec, start_q = Q, t) * phi * beta * sS * (1 / N) -
#     eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K) * sS -
#     eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (sS/N)
#   
#   dsEdt <- p1 * beta * sS * ((I + sI) / N) +
#     p2(probs = p2_vec, start_q = Q, t) * phi * beta * sS * (1 / N) -
#     sigma * sE -
#     eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K) * sE -
#     eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (sE/N)
#   
#   dsIdt <- sigma * sE -
#     eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N/K) * sI -
#     eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (sI/N)
#   
#   list(c(dSdt, dEdt, dIdt, dsSdt, dsEdt, dsIdt))
# }
# Frequency-dependent SEI model with super spreaders
# Form B: explicit beta_wild (deer–deer) and beta_farm,q = p_q * beta_wild, both / N
SEI_model_full_freq <- function(t, x, parms) {
  
  S  <- as.numeric(x["S"])
  E  <- as.numeric(x["E"])
  I  <- as.numeric(x["I"])
  sS <- as.numeric(x["sS"])
  sE <- as.numeric(x["sE"])
  sI <- as.numeric(x["sI"])
  
  # births
  Alpha_max    <- as.numeric(parms["alpha_max"])
  ksi          <- as.numeric(parms["ksi"])
  omega        <- as.numeric(parms["omega"])
  s            <- as.numeric(parms["s"])
  alpha_annual <- as.numeric(parms["alpha"])
  
  K     <- as.numeric(parms["K"])
  eta_n <- as.numeric(parms["eta_nat"])
  eta_h <- as.numeric(parms["eta_hunt"])
  theta <- as.numeric(parms["theta"])
  gamma <- as.numeric(parms["gamma"])
  
  beta_wild <- as.numeric(parms["beta"])    # deer–deer transmission coefficient (month^-1)
  A         <- as.numeric(parms["area"])    # kept for parity; unused in freq form
  p1        <- as.numeric(parms["p1"])      # kept for interface compatibility (unused)
  p2_vec    <- c(as.numeric(parms["p2_q1"]),
                 as.numeric(parms["p2_q2"]),
                 as.numeric(parms["p2_q3"]),
                 as.numeric(parms["p2_q4"]))  # unitless seasonal multipliers
  phi       <- as.numeric(parms["phi"])     # SS farm-FOI multiplier
  sigma     <- 1 / as.numeric(parms["sigma1_mean"])
  if (as.numeric(parms["sigma1_mean"]) == 0) sigma <- 0
  
  Q <- as.numeric(parms["start_q"])
  
  N <- S + E + I + sS + sE + sI
  if (N <= 0) N <- 1  # guard against division-by-zero
  
  beta_farm_t <- function(tt) p2(probs = p2_vec, start_q = Q, t = tt) * beta_wild
  
  dSdt  <- (1 - ksi) * alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t) * N
  dSdt  <- dSdt - beta_wild   * S  * (I + sI) / N          # deer–deer
  dSdt  <- dSdt - beta_farm_t(t) * S        / N            # farm
  dSdt  <- dSdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * S
  dSdt  <- dSdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (S / N)
  
  dEdt  <-  beta_wild   * S  * (I + sI) / N                # deer–deer (S)
  dEdt  <- dEdt + beta_farm_t(t) * S        / N            # farm (S)
  dEdt  <- dEdt + beta_wild   * sS * (I + sI) / N          # deer–deer (SS)
  dEdt  <- dEdt + beta_farm_t(t) * phi * sS / N            # farm (SS boosted)
  dEdt  <- dEdt - sigma * E
  dEdt  <- dEdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * E
  dEdt  <- dEdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (E / N)
  
  dIdt  <-  sigma * E
  dIdt  <- dIdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * I
  dIdt  <- dIdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (I / N)
  
  dsSdt <-  ksi * alpha(Amp = Alpha_max, synchrony = s, phase = omega, start_q = Q, t) * N
  dsSdt <- dsSdt - beta_wild   * sS * (I + sI) / N         # deer–deer
  dsSdt <- dsSdt - beta_farm_t(t) * phi * sS / N           # farm (SS boosted)
  dsSdt <- dsSdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * sS
  dsSdt <- dsSdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (sS / N)
  
  dsEdt <-  beta_wild   * sS * (I + sI) / N                # deer–deer (SS)
  dsEdt <- dsEdt + beta_farm_t(t) * phi * sS / N           # farm (SS boosted)
  dsEdt <- dsEdt - sigma * sE
  dsEdt <- dsEdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * sE
  dsEdt <- dsEdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (sE / N)
  
  dsIdt <-  sigma * sE
  dsIdt <- dsIdt - eta(baseline_mortality = eta_n, theta = theta, y_scale = gamma, N_frac = N / K) * sI
  dsIdt <- dsIdt - eta_hunt(hunt_mortality = eta_h, start_q = Q, t) * (sI / N)
  
  list(c(dSdt, dEdt, dIdt, dsSdt, dsEdt, dsIdt))
}

# C++ stochastic model functions
wildlife_model <- function(n_reps, parameters, initial_state, nyears, seed_quarter, verbose, batch_name, type = 'discrete', lambda = c(), integrate_type = 0){
  #The path to the bTB cpp binary file must be set correctly in the sys call below:
  #r <- system2(command = "/home/webblab/Documents/Brandon/bTB_wildlife_code/filename.exe", args=strvec, stdout=TRUE)
  if(type %in% c('discrete','Discrete', 'D', 'd')){
    parvec = c( n_reps, #Number of replicates
                initial_state$S_0, #0, start Susc.
                initial_state$E1_0, #1, start E1
                initial_state$I_0, #2, start I
                initial_state$SuperS_0, #3
                initial_state$SuperE1_0, #4
                initial_state$SuperI_0, #5
                
                parameters$K, #6, carrying capacity (county level)
                parameters$eta_hunt, #7, Hunting induced mortality
                parameters$eta_nat, #8, Natural mortality rate
                parameters$theta, #9, density dependent asymmetry  
                parameters$gamma, #10, density dependence scalar
                
                parameters$alpha_max, #11 birth rate amplitude
                parameters$ksi, #12 proportion of newborn super spreaders
                parameters$omega, #13 birth rate phase shift
                parameters$s, #14 birth rate synchrony
                parameters$alpha, #15, annual birth rate
                
                parameters$beta, #16 Density Dependent Transmission rate
                parameters$area, #17, available habitat area
                parameters$p1, #18 deer-deer contact rate
                parameters$p2_q1, #19 C-to-wildl contact rate quarter 1
                parameters$p2_q2, #20 C-to-wildl contact rate quarter 2
                parameters$p2_q3, #21 C-to-wildl contact rate quarter 3
                parameters$p2_q4, #22 C-to-wildl contact rate quarter 4
                parameters$phi, #23, increased SS farm contact rate factor
                parameters$sigma1_mean, #24, mean exposed period length
                parameters$sigma1_rate, #25, exposed rate
                
                lambda[1], #26, L1
                lambda[2], #27, L2
                lambda[3], #28, L3
                lambda[4], #29, L4
                lambda[5], #30, L5
                lambda[6], #31, L6
                lambda[7], #32, L7
                lambda[8], #33, L8
                lambda[9], #34, L9
                lambda[10], #35, L10
                lambda[11], #36, L11
                lambda[12], #37, L12
                integrate_type, #38, birth rate integral numeric method (midpoint=1, left=2, or trapezoid=3)
                
                nyears, #39, simulation length
                seed_quarter, #40, starting quarter
                verbose, #41, error checkpoints for debugging
                batch_name #42
    )
    strvec = format(parvec, digits = 5)
    r <- system2(command = "./wl_model_DTMC.exe", args=strvec, stdout=TRUE)
  }
  else if(type %in% c('continuous','Comtinuous', 'C', 'c')){
    parvec = c( n_reps, #Number of replicates
                initial_state$S_0, #0, start Susc.
                initial_state$E1_0, #1, start E1
                initial_state$I_0, #2, start I
                initial_state$SuperS_0, #3
                initial_state$SuperE1_0, #4
                initial_state$SuperI_0, #5
                
                parameters$K, #6, carrying capacity (county level)
                parameters$eta_hunt, #7, Hunting induced mortality
                parameters$eta_nat, #8, Natural mortality rate
                parameters$theta, #9, density dependent asymmetry  
                parameters$gamma, #10, density dependence scalar
                
                parameters$alpha_max, #11 birth rate amplitude
                parameters$ksi, #12 proportion of newborn super spreaders
                parameters$omega, #13 birth rate phase shift
                parameters$s, #14 birth rate synchrony
                parameters$alpha, #15, annual birth rate
                
                parameters$beta, #16 Density Dependent Transmission rate
                parameters$area, #17, available habitat area
                parameters$p1, #18 deer-deer contact rate
                parameters$p2_q1, #19 C-to-wildl contact rate quarter 1
                parameters$p2_q2, #20 C-to-wildl contact rate quarter 2
                parameters$p2_q3, #21 C-to-wildl contact rate quarter 3
                parameters$p2_q4, #22 C-to-wildl contact rate quarter 4
                parameters$phi, #23, increased SS farm contact rate factor
                parameters$sigma1_mean, #24, mean exposed period length
                parameters$sigma1_rate, #25, exposed rate
                
                nyears, #26, simulation length
                seed_quarter, #27, starting quarter
                verbose, #28, error checkpoints for debugging
                batch_name #29
    )
    strvec = format(parvec, digits = 5)
    r <- system2(command = "./wl_model_CTMC.exe", args=strvec, stdout=TRUE)
  }
  else{print('invalid type selection - use `c` for CTMC or `d` for DTMC'); break;}
  res <- read.table(text=r, header=TRUE, sep=';', check.names = FALSE)
  return (res)
}


#extract simulation monthly lambda vector
getLambda_vec <- function(data = NA, type = ''){
  data <- as.data.frame(data) 
  lambdas = c()
  
  if(type=='min'){lambdas <- aggregate(data[, "Lambda"], list(data$Month), min)$x}
  else if(type=='mean'){lambdas <- aggregate(data[, "Lambda"], list(data$Month), mean)$x}
  else if(type=='median'){lambdas <- aggregate(data[, "Lambda"], list(data$Month), median)$x}
  else{lambdas <- aggregate(data[, "Lambda"], list(data$Month), max)$x}
  
  return(as.vector(lambdas))
}


library(deSolve); library(tidyverse); library(ggplot2); library(ggpubr);
library(grid); library(rstudioapi); library(foreach); library(doParallel);
library(RColorBrewer); library(scales);

initialize_simulation <- function(run_type, inf_type, years, herd_sizes, pct, ss_prop) {
  list(
    years = years,
    times = seq(0, years * 12, by = 12 / 365),
    seedQuarter = 1,
    prop_superSpreader = ss_prop,
    pct = pct,
    reps = 50,
    sizes = herd_sizes,
    infType = inf_type,
    runtype = run_type,
    name_out = paste0(run_type, inf_type, "_", pct, '-', ss_prop, '-'),
    pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/",
    save_runs = TRUE,
    save_plots = TRUE
  )
}

setup_parallel <- function() {
  n.cores <- floor(detectCores() * 3 / 4)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(deSolve); library(tidyverse); library(ggplot2); library(ggpubr);
    library(grid); library(rstudioapi); library(foreach); library(doParallel);
    library(RColorBrewer); library(scales) })
  return(cl)
}

run_model <- function(params, size) {
  pars <- parameter_set_wl(k = size, scenario = params$infType, initial_exposed = params$pct * size,
                           SS_prop = params$prop_superSpreader, start_quarter = params$seedQuarter)
  
  X0_full <- c(S = pars["S_0"], E = pars["E1_0"], I = pars["I_0"],
               sS = pars["SuperS_0"], sE = pars["SuperE1_0"], sI = pars["SuperI_0"])
  
  out <- as.data.frame(ode(func = SEI_model_full, y = X0_full, times = params$times, parms = pars, method = "rk4"))
  out$N <- rowSums(out[, c("S", "E", "I", "sS", "sE", "sI")])
  out %>% pivot_longer(-time, names_to = "variable", values_to = "value")
}

plot_results <- function(data, params, size) {
  SEIcols <- RColorBrewer::brewer.pal(11, "Spectral")[c(1,2,4,5,8,9,10)]
  
  ggplot(data) +
    geom_line(aes(x = time, y = value, color = variable, group = rep), size = 1, alpha = .7) +
    scale_color_manual(values = SEIcols) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(breaks = seq(0, params$years * 12, 12), limits = c(0, params$years * 12))
}

main <- function() {
  params <- initialize_simulation("test_", "seeded", 3, c(10, 50, 100, 250, 500), 0.02, 0.05)
  cl <- setup_parallel()
  
  foreach(size = params$sizes, .packages = c("deSolve", "tidyverse", "ggplot2")) %dopar% {
    data <- run_model(params, size)
    plot_results(data, params, size)
    if (params$save_runs) save(data, file = paste0('test_runs/', params$name_out, size, '-', Sys.Date(), '.RData'))
  }
  
  stopCluster(cl)
}

##############################
## Sensitivity Functions
##############################
library(tidyverse)
library(epiR)
library(broom)

#-----------------------------
# 1. PRCC Function
#-----------------------------
run_prcc <- function(data, response, parameters) {
  # select relevant columns
  dat <- data[c(parameters, response)]
  
  # run PRCC
  res <- epi.prcc(dat, sided.test = 2)
  
  res %>%
    select(var, est, p.value) %>%
    mutate(response = response)
  
}

#-----------------------------
# 2. Linear Model Function
#-----------------------------
run_lm <- function(data, response, parameters) {
  
  form <- as.formula(paste(response, "~", paste(parameters, collapse = "+")))
  fit  <- lm(form, data = data)
  
  # extract coefficients (scaled) and tidy
  coefs <- tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      scaled   = estimate / max(abs(estimate), na.rm = TRUE),
      response = response
    )
  
  # return tidy tibble
  coefs %>%
    select(response, parameter = term, estimate = scaled, p_value = p.value)
}

run_lm_int <- function(data, response, parameters) {
  
  pars_bt <- sprintf("`%s`", parameters)
  pairs_star <- combn(pars_bt, 2, FUN = function(x) paste(x, collapse = ":"))
  rhs_star   <- paste(pairs_star, collapse = " + ")
  form  <- as.formula(paste(response, "~", rhs_star))
  fit  <- lm(form, data = data)
  
  # extract coefficients (scaled) and tidy
  coefs <- tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      scaled   = estimate / max(abs(estimate), na.rm = TRUE),
      response = response
    )
  
  # return tidy tibble
  coefs %>%
    select(response, parameter = term, estimate = scaled, p_value = p.value)
}


param_labels <- c("K"            = "carrying capacity ($K$)",
                  "eta_hunt"     = "hunt mortality ($\\eta_{h,q=4}$)",
                  "eta_nat"      = "base mortality ($\\eta$)",
                  "theta"        = "density-dep asymmetry ($\\theta$)",
                  "gamma"        = "density-dep mortality ($\\gamma$)",
                  "alpha_max"    = "max. birth rate ($A$)",
                  "ksi"          = "proportion SS ($\\xi$)",
                  "omega"        = "birth timing ($\\omega$)",
                  "s"            = "birth duration ($s$)",
                  "beta"         = "transmission rate ($\\beta$)",
                  "p1"         = "deer-deer contact rate ($p_1$)",
                  "phi"          = "SS contact factor ($\\phi$)",
                  "sigma1_mean"  = "latency mean ($\\kappa$)",
                  "sigma1_rate"  = "latency rate ($\\gamma$)")

# helper to make readable interaction labels
make_label <- function(term) {
  parts <- unlist(strsplit(term, ":", fixed = TRUE))
  nice_parts <- param_labels[parts]
  if (any(is.na(nice_parts))) nice_parts[is.na(nice_parts)] <- parts[is.na(nice_parts)]
  if (length(parts) == 1) {
    TeX(nice_parts)
  } else {
    TeX(paste0(nice_parts[1], " × ", nice_parts[2]))
  }
}

theme_pub <- theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    panel.grid = element_line(color = "gray85", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 15, face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    axis.ticks = element_line(color = "black", size = 0.3),
    legend.position = "none"
  )

make_panel <- function(sto, det, var, col, ylab, years) {
  ggplot() +
    geom_line(data = sto[[var]], aes(x = tstep, y = value, group = rep),
              color = col, alpha = 0.25, linewidth = 0.7) +
    geom_line(data = det[[var]], aes(x = time, y = value),
              color = "black", linewidth = 1.3) +
    scale_y_continuous(name = ylab, limits = c(0, NA)) +
    scale_x_continuous(name = "time (months)",
                       breaks = seq(0, years * 12, 12),
                       limits = c(0, years * 12)) +
    theme_pub
}

library(ggplot2)
library(ggpubr)

# run_wl_sim <- function(
#     size,
#     infType = c("seeded", "spillover"),
#     years = 7,
#     times = seq(0, years * 12, by = .5),
#     seedQuarter = 1,
#     prop_superSpreader = 0.1,
#     pct = 0.02,
#     reps = 500,
#     lambda_factor = 1.2,
#     integrate_type = 3,
#     verbose = 0,
#     scaled = TRUE,
#     supp = TRUE,
#     save_runs = TRUE,
#     save_plots = TRUE,
#     overrides = NULL
# ) {
# 
#   pars <- parameter_set_wl(
#     k               = as.integer(size),
#     scenario        = infType,
#     initial_exposed = pct * as.integer(size),
#     SS_prop         = prop_superSpreader,
#     start_quarter   = seedQuarter,
#     test = FALSE, birth = FALSE, death_h = FALSE, death_n = FALSE,
#     disease = FALSE, verbose = verbose
#   )
#   if (!is.null(overrides) && length(overrides) > 0) {
#     if (!is.list(overrides)) {
#       if (is.null(names(overrides)) || any(names(overrides) == "")) {
#         stop("`overrides` must be a named list or named vector.")
#       }
#       overrides <- as.list(overrides)
#     }
#     for (nm in names(overrides)) if (nm %in% names(pars)) pars[[nm]] <- overrides[[nm]]
#   }
#   
#   #deterministic model
#   X0_full = c(S = as.integer(pars["S_0"]), E = as.integer(pars["E1_0"]), I = as.integer(pars["I_0"]), 
#               sS = as.integer(pars["SuperS_0"]), sE = as.integer(pars["SuperE1_0"]), sI = as.integer(pars["SuperI_0"]))
#   
#   tic = Sys.time()
#   ode(
#     func = SEI_model_full,
#     y = X0_full,
#     times = times,
#     parms = pars,
#     method = "rk4"
#   ) %>%
#     as.data.frame() -> out
#   toc = Sys.time()
#   print(toc-tic)
#   
#   out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
#   data <- out %>% gather(variable,value,-time)
#   remove(out)
#   data$value[is.nan(data$value)] <- 0 
#   data$value[which(data$value < 0)] <- 0
#   
#   
#   #stochastic model
#   initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
#   population_parameters<-data.frame(K=pars["K"], eta_hunt=pars["eta_hunt"], eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]) 
#   disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
#   
#   param_vals<-data.frame(merge(population_parameters,disease_parameters))
#   
#   #get lambda
#   tic = Sys.time()
#   sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 5, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
#   toc = Sys.time()
#   print(toc-tic)
#   lambda_out <- getLambda_vec(data = sto_out, type = 'max')
#   remove(sto_out)
#   
#   #run
#   tic = Sys.time()
#   sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = integrate_type)
#   toc = Sys.time()
#   print(toc-tic)
#   remove(initial_state,population_parameters,disease_parameters,param_vals)
#   
#   data_sto <- sto_out |> 
#     select("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I") |> 
#     pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
#   
#   models_res <- vector(mode = "list", length = 3)
#   
#   models_res[[1]] <- as.data.frame(data)
#   remove(data)
#   models_res[[2]] <- as.data.frame(sto_out)
#   remove(sto_out)
#   models_res[[3]] <- as.data.frame(data_sto)
#   remove(data_sto)
#   
#   ##########
#   
#   ###########
#   ## Plots ##
#   ###########
#   
#   #class comparison plots
#   sto_by_class <- split(as.data.frame(models_res[[3]]), as.data.frame(models_res[[3]])$variable)
#   det_by_class <- split(as.data.frame(models_res[[1]]), as.data.frame(models_res[[1]])$variable)
#   
#   SEIcols <- RColorBrewer::brewer.pal(11,"Spectral")[c(1,2,4,5,8,9,10)]
#   
#   #susceptible
#   sto_by_class[['S']]$value <- sto_by_class[["S"]]$value
#   det_by_class[['S']]$value <- det_by_class[["S"]]$value
#   
#   #susceptible SS
#   sto_by_class[['SS_S']]$value <- sto_by_class[["SS_S"]]$value
#   det_by_class[['sS']]$valu <- det_by_class[['sS']]$value
#   
#   #exposed -- merged SS
#   sto_by_class[['E1']]$value <- sto_by_class[["E1"]]$value + sto_by_class[["SS_E1"]]$value
#   det_by_class[['E']]$value <- det_by_class[["E"]]$value + det_by_class[["sE"]]$value
#   
#   #infected -- merged SS
#   sto_by_class[['I']]$value <- sto_by_class[["I"]]$value + sto_by_class[["SS_I"]]$value
#   det_by_class[['I']]$value <- det_by_class[["I"]]$value + det_by_class[["sI"]]$value
#   
#   
#   Nplot <- ggplot() +
#     geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) +
#     geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
#     scale_color_manual(values=SEIcols[7]) +
#     scale_y_continuous(limits = c(0,NA)) +
#     scale_x_continuous(breaks=seq(0,years*12,1), limits = c(0, years*12))
#   
#   Splot <- ggplot() + 
#     geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
#     geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#     scale_color_manual(values=scales::alpha(SEIcols[6], 0.1)) +
#     theme_bw() +
#     theme(text = element_text(size = 16),
#           panel.grid = element_blank()) +
#     scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
#     scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#   
#   Eplot <- ggplot() + 
#     geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
#     geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#     scale_color_manual(values=scales::alpha(SEIcols[4], 0.1)) +
#     theme_bw() +
#     theme(text = element_text(size = 16),
#           panel.grid = element_blank()) +
#     scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
#     scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#   
#   Iplot <- ggplot() + 
#     geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
#     geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#     scale_color_manual(values=scales::alpha(SEIcols[2], 0.1)) +
#     theme_bw() +
#     theme(text = element_text(size = 16),
#           panel.grid = element_blank()) +
#     scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
#     scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#   
#   SSSplot <- ggplot() + 
#     geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) + 
#     geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#     scale_color_manual(values=scales::alpha(SEIcols[5], 0.1)) +
#     theme_bw() +
#     theme(text = element_text(size = 16),
#           panel.grid = element_blank()) +
#     scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
#     scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#   
#   
#   if(scaled){
#     #Proportion susceptible
#     sto_by_class[['S']]$value <- sto_by_class[["S"]]$value/sto_by_class[["N"]]$value
#     det_by_class[['S']]$value <- det_by_class[["S"]]$value/det_by_class[["N"]]$value
#     
#     #Proportion susceptible SS
#     sto_by_class[['SS_S']]$value <- sto_by_class[["SS_S"]]$value/sto_by_class[["N"]]$value
#     det_by_class[['sS']]$value <- det_by_class[['sS']]$value/det_by_class[["N"]]$value
#     
#     
#     #Proportion exposed -- merged SS
#     sto_by_class[['E1']]$value <- sto_by_class[["E1"]]$value/sto_by_class[["N"]]$value 
#     det_by_class[['E']]$value <- det_by_class[["E"]]$value/det_by_class[["N"]]$value 
#     
#     #Proportion infected -- merged SS
#     sto_by_class[['I']]$value <- sto_by_class[["I"]]$value/sto_by_class[["N"]]$value 
#     det_by_class[['I']]$value <- det_by_class[["I"]]$value/det_by_class[["N"]]$value 
#     
#     Splot_scl <- ggplot() + 
#       geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
#       geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#       scale_color_manual(values=SEIcols[6]) +
#       theme_bw() +
#       theme(text = element_text(size = 16),
#             panel.grid = element_blank()) +
#       scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
#       scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#     
#     Eplot_scl <- ggplot() + 
#       geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
#       geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#       scale_color_manual(values=SEIcols[4]) +
#       theme_bw() +
#       theme(text = element_text(size = 16),
#             panel.grid = element_blank()) +
#       scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
#       scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#     
#     Iplot_scl <- ggplot() + 
#       geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
#       geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#       scale_color_manual(values=SEIcols[2]) +
#       theme_bw() +
#       theme(text = element_text(size = 16),
#             panel.grid = element_blank()) +
#       scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
#       scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#     
#     SSSplot_scl <- ggplot() + 
#       geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) + 
#       geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
#       scale_color_manual(values=SEIcols[5]) +
#       theme_bw() +
#       theme(text = element_text(size = 16),
#             panel.grid = element_blank()) +
#       scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
#       scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
#   }
#   
#   if(save_plots){
#     
#     plot <- ggarrange(Splot, Eplot, Iplot,SSSplot, ncol=3, nrow=2)
#     plot_scl <- ggarrange(Splot_scl, Eplot_scl,  Iplot_scl, SSSplot_scl, ncol=3, nrow=2)
#     
#     plot
#     plot_scl
#     
#     ggsave(paste0("figures/", infType,"_classes_",size,".pdf"), plot, width = 10, height = 7, dpi = 300)
#     ggsave(paste0("figures/", infType,"_classes_",size,"_scaled.pdf"), plot_scl, width = 10, height = 7, dpi = 300)
#     
#     if (supp) {
#       # compute mean trajectory for each stochastic class
#       sto_means <- lapply(sto_by_class, function(df) {
#         df %>%
#           group_by(tstep) %>%
#           summarise(mean_value = mean(value, na.rm = TRUE))
#       })
#       
#       # example for S
#       Splot <- ggplot() +
#         geom_line(data = sto_by_class[["S"]],
#                   aes(x = tstep, y = value, group = rep),
#                   color = scales::alpha(SEIcols[6], 0.4), size = 0.7, alpha = 0.6) +
#         geom_line(data = sto_means[["S"]],
#                   aes(x = tstep, y = mean_value),
#                   color = "grey50", linewidth = 0.75) +
#         geom_line(data = det_by_class[["S"]],
#                   aes(x = time, y = value),
#                   color = "black", linewidth = 1.5) +
#         theme_bw() +
#         theme(text = element_text(size = 16),
#               panel.grid = element_blank()) +
#         scale_y_continuous(name = 'Susceptible', limits = c(0, NA)) +
#         scale_x_continuous(name = 'time (months)',
#                            breaks = seq(0, years * 12, 12),
#                            limits = c(0, years * 12),
#                            expand = c(0, 0))
#       
#       # repeat pattern for E, I, and SS_S
#       Eplot <- ggplot() +
#         geom_line(data = sto_by_class[["E1"]],
#                   aes(x = tstep, y = value, group = rep),
#                   color = scales::alpha(SEIcols[4], 0.4), size = 0.7, alpha = 0.6) +
#         geom_line(data = sto_means[["E1"]],
#                   aes(x = tstep, y = mean_value),
#                   color = "grey50", linewidth = 0.75) +
#         geom_line(data = det_by_class[["E"]],
#                   aes(x = time, y = value),
#                   color = "black", linewidth = 1.5) +
#         theme_bw() +
#         theme(text = element_text(size = 16),
#               panel.grid = element_blank()) +
#         scale_y_continuous(name = 'Exposed', limits = c(0, NA)) +
#         scale_x_continuous(name = 'time (months)',
#                            breaks = seq(0, years * 12, 12),
#                            limits = c(0, years * 12),
#                            expand = c(0, 0))
#       
#       Iplot <- ggplot() +
#         geom_line(data = sto_by_class[["I"]],
#                   aes(x = tstep, y = value, group = rep),
#                   color = scales::alpha(SEIcols[2], 0.4), size = 0.7, alpha = 0.6) +
#         geom_line(data = sto_means[["I"]],
#                   aes(x = tstep, y = mean_value),
#                   color = "grey50", linewidth = 0.75) +
#         geom_line(data = det_by_class[["I"]],
#                   aes(x = time, y = value),
#                   color = "black", linewidth = 1.5) +
#         theme_bw() +
#         theme(text = element_text(size = 16),
#               panel.grid = element_blank()) +
#         scale_y_continuous(name = 'Infectious', limits = c(0, NA)) +
#         scale_x_continuous(name = 'time (months)',
#                            breaks = seq(0, years * 12, 12),
#                            limits = c(0, years * 12),
#                            expand = c(0, 0))
#       
#       SSSplot <- ggplot() +
#         geom_line(data = sto_by_class[["SS_S"]],
#                   aes(x = tstep, y = value, group = rep),
#                   color = scales::alpha(SEIcols[5], 0.4), size = 0.7, alpha = 0.6) +
#         geom_line(data = sto_means[["SS_S"]],
#                   aes(x = tstep, y = mean_value),
#                   color = "grey50", linewidth = 0.75) +
#         geom_line(data = det_by_class[["sS"]],
#                   aes(x = time, y = value),
#                   color = "black", linewidth = 1.2) +
#         theme_bw() +
#         theme(text = element_text(size = 16),
#               panel.grid = element_blank()) +
#         scale_y_continuous(name = expression('Susceptible'['SS']), limits = c(0, NA)) +
#         scale_x_continuous(name = 'time (months)',
#                            breaks = seq(0, years * 12, 12),
#                            limits = c(0, years * 12),
#                            expand = c(0, 0))
#       
#       
#       plot <- ggarrange(Splot, SSSplot, Eplot, Iplot, nrow=1)
#       
#       plot
#       
#       ggsave(paste0("figures/supp/", infType,"_classes_",size,"_supp.pdf"), plot, width = 15, height = 4, dpi = 300)
#     }
#     
#     remove(plot, plot_scl)
#   }
#   
#   remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
#   
#   if(save_runs){
#     save(models_res, file = paste0('data/raw/', name_out, size, '-', Sys.Date(), '.RData'))
#   }
#   
#   remove(models_res)
# }

## ===========================================================
##  Wrapper to run simulations easily
## ===========================================================
run_wl_sim <- function(
    size,
    infType = c("seeded", "spillover"),
    years = 7,
    times = seq(0, years * 12, by = .5),
    seedQuarter = 1,
    prop_superSpreader = 0.1,
    pct = 0.02,
    reps = 500,
    lambda_factor = 1.2,
    integrate_type = 3,
    verbose = 0,
    scaled = TRUE,           
    supp = TRUE,        
    save_runs = TRUE,
    save_plots = FALSE, 
    overrides = NULL,
    freq = FALSE,
    ode = TRUE
) {
  
  ## Parameters from your helper
  pars <- parameter_set_wl(
    k               = as.integer(size),
    scenario        = infType,
    initial_exposed = pct,
    SS_prop         = prop_superSpreader,
    start_quarter   = seedQuarter,
    test = FALSE, birth = FALSE, death_h = FALSE, death_n = FALSE,
    disease = FALSE, verbose = verbose
  )
  
  ## Inline apply_overrides (no separate function)
  if (!is.null(overrides) && length(overrides) > 0) {
    if (!is.list(overrides)) {
      if (is.null(names(overrides)) || any(names(overrides) == "")) {
        stop("`overrides` must be a named list or named vector.")
      }
      overrides <- as.list(overrides)
    }
    for (nm in names(overrides)) {
      if (nm %in% names(pars)) pars[[nm]] <- overrides[[nm]]
    }
  }
  
  ## Deterministic ODE (compact; no plotting)
  X0_full <- c(
    S  = as.integer(pars["S_0"]),
    E  = as.integer(pars["E1_0"]),
    I  = as.integer(pars["I_0"]),
    sS = as.integer(pars["SuperS_0"]),
    sE = as.integer(pars["SuperE1_0"]),
    sI = as.integer(pars["SuperI_0"])
  )
  
  if (freq) {
    pars["beta"] <- pars["beta"]*pars["K"]
  }
  
  if (ode){
    tic <- Sys.time()
    if (freq) {
      out <- ode(func = SEI_model_full_freq, y = X0_full, times = times, parms = pars, method = "rk4" ) |> as.data.frame()
      infType <- paste0("freq_",infType)
    } else {
      out <- ode(func = SEI_model_full, y = X0_full, times = times, parms = pars, method = "rk4" ) |> as.data.frame()
    }
    toc <- Sys.time(); print(toc - tic)
    
    out$N <- out$S + out$E + out$I + out$sS + out$sE + out$sI
    data <- tidyr::gather(out, variable, value, -time)
    rm(out)
    data$value[is.nan(data$value)] <- 0
    data$value[data$value < 0] <- 0
  }
  
  ## Stochastic — estimate lambda, then main run
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
    p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], p2_q4 = pars["p2_q4"],
    phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], sigma1_rate = pars["sigma1_rate"]
  )
  param_vals <- merge(population_parameters, disease_parameters)
  
  tic <- Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state,
                            nyears = 5, seed_quarter = as.integer(pars["start_q"]),
                            verbose = -3, batch_name = "prelambda", type = 'c')
  toc <- Sys.time(); print(toc - tic)
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  rm(sto_out)
  
  tic <- Sys.time()
  sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state,
                            nyears = years, seed_quarter = as.integer(pars["start_q"]),
                            verbose = verbose, batch_name = "main", type = 'd',
                            lambda = lambda_out * lambda_factor, integrate_type = integrate_type)
  toc <- Sys.time(); print(toc - tic)
  rm(initial_state, population_parameters, disease_parameters, param_vals)
  
  data_sto <- sto_out |>
    dplyr::select(rep, tstep, N, S, SS_S, E1, SS_E1, I, SS_I) |>
    tidyr::pivot_longer(c(N, S, SS_S, E1, SS_E1, I, SS_I), names_to = "variable", values_to = "value")
  
  models_res <- vector("list", 3)
  models_res[[1]] <- as.data.frame(data);
  models_res[[2]] <- as.data.frame(sto_out);
  models_res[[3]] <- as.data.frame(data_sto);
  
  if(save_plots){
    
    ###########
    ## Plots ##
    ###########
    
    #class comparison plots
    sto_by_class <- split(as.data.frame(models_res[[3]]), as.data.frame(models_res[[3]])$variable)
    if (ode) {
      det_by_class <- split(as.data.frame(models_res[[1]]), as.data.frame(models_res[[1]])$variable)
    } else {
      # dummy empty data for plotting when ode = FALSE
      det_by_class <- list(
        S = data.frame(time = numeric(0), value = numeric(0)),
        sS = data.frame(time = numeric(0), value = numeric(0)),
        E = data.frame(time = numeric(0), value = numeric(0)),
        sE = data.frame(time = numeric(0), value = numeric(0)),
        I = data.frame(time = numeric(0), value = numeric(0)),
        sI = data.frame(time = numeric(0), value = numeric(0)),
        N = data.frame(time = numeric(0), value = numeric(0))
      )
    }
    # det_by_class <- split(as.data.frame(models_res[[1]]), as.data.frame(models_res[[1]])$variable)
    
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
      scale_color_manual(values=scales::alpha(SEIcols[6], 0.1)) +
      theme_bw() +
      theme(text = element_text(size = 16),
            panel.grid = element_blank()) +
      scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
    
    Eplot <- ggplot() + 
      geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
      geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=scales::alpha(SEIcols[4], 0.1)) +
      theme_bw() +
      theme(text = element_text(size = 16),
            panel.grid = element_blank()) +
      scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
    
    Iplot <- ggplot() + 
      geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
      geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=scales::alpha(SEIcols[2], 0.1)) +
      theme_bw() +
      theme(text = element_text(size = 16),
            panel.grid = element_blank()) +
      scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
    
    SSSplot <- ggplot() + 
      geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) + 
      geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
      scale_color_manual(values=scales::alpha(SEIcols[5], 0.1)) +
      theme_bw() +
      theme(text = element_text(size = 16),
            panel.grid = element_blank()) +
      scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
      scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
    
    
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
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
        scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
      
      Eplot_scl <- ggplot() + 
        geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
        geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
        scale_color_manual(values=SEIcols[4]) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
        scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
      
      Iplot_scl <- ggplot() + 
        geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
        geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
        scale_color_manual(values=SEIcols[2]) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
        scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
      
      SSSplot_scl <- ggplot() + 
        geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) + 
        geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
        scale_color_manual(values=SEIcols[5]) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
        scale_x_continuous(name = 'time (months)', breaks=seq(0,years*12,12), limits = c(0, years*12), expand = c(0, 0))
    }
    
    plot <- ggarrange(Splot, Eplot, Iplot,SSSplot, ncol=3, nrow=2)
    plot_scl <- ggarrange(Splot_scl, Eplot_scl,  Iplot_scl, SSSplot_scl, ncol=3, nrow=2)
    
    plot
    plot_scl
    
    ggsave(paste0("figures/", infType,"_classes_",size,".pdf"), plot, width = 10, height = 7, dpi = 300)
    ggsave(paste0("figures/", infType,"_classes_",size,"_scaled.pdf"), plot_scl, width = 10, height = 7, dpi = 300)
    
    if (supp) {
      # compute mean trajectory for each stochastic class
      sto_means <- lapply(sto_by_class, function(df) {
        df %>%
          group_by(tstep) %>%
          summarise(mean_value = mean(value, na.rm = TRUE))
      })
      
      # example for S
      Splot <- ggplot() +
        geom_line(data = sto_by_class[["S"]],
                  aes(x = tstep, y = value, group = rep),
                  color = scales::alpha(SEIcols[6], 0.4), size = 0.7, alpha = 0.6) +
        geom_line(data = sto_means[["S"]],
                  aes(x = tstep, y = mean_value),
                  color = "grey50", linewidth = 0.75) +
        geom_line(data = det_by_class[["S"]],
                  aes(x = time, y = value),
                  color = "black", linewidth = 1.5) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = 'Susceptible', limits = c(0, NA)) +
        scale_x_continuous(name = 'time (months)',
                           breaks = seq(0, years * 12, 12),
                           limits = c(0, years * 12),
                           expand = c(0, 0))
      
      # repeat pattern for E, I, and SS_S
      Eplot <- ggplot() +
        geom_line(data = sto_by_class[["E1"]],
                  aes(x = tstep, y = value, group = rep),
                  color = scales::alpha(SEIcols[4], 0.4), size = 0.7, alpha = 0.6) +
        geom_line(data = sto_means[["E1"]],
                  aes(x = tstep, y = mean_value),
                  color = "grey50", linewidth = 0.75) +
        geom_line(data = det_by_class[["E"]],
                  aes(x = time, y = value),
                  color = "black", linewidth = 1.5) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = 'Exposed', limits = c(0, NA)) +
        scale_x_continuous(name = 'time (months)',
                           breaks = seq(0, years * 12, 12),
                           limits = c(0, years * 12),
                           expand = c(0, 0))
      
      Iplot <- ggplot() +
        geom_line(data = sto_by_class[["I"]],
                  aes(x = tstep, y = value, group = rep),
                  color = scales::alpha(SEIcols[2], 0.4), size = 0.7, alpha = 0.6) +
        geom_line(data = sto_means[["I"]],
                  aes(x = tstep, y = mean_value),
                  color = "grey50", linewidth = 0.75) +
        geom_line(data = det_by_class[["I"]],
                  aes(x = time, y = value),
                  color = "black", linewidth = 1.5) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = 'Infectious', limits = c(0, NA)) +
        scale_x_continuous(name = 'time (months)',
                           breaks = seq(0, years * 12, 12),
                           limits = c(0, years * 12),
                           expand = c(0, 0))
      
      SSSplot <- ggplot() +
        geom_line(data = sto_by_class[["SS_S"]],
                  aes(x = tstep, y = value, group = rep),
                  color = scales::alpha(SEIcols[5], 0.4), size = 0.7, alpha = 0.6) +
        geom_line(data = sto_means[["SS_S"]],
                  aes(x = tstep, y = mean_value),
                  color = "grey50", linewidth = 0.75) +
        geom_line(data = det_by_class[["sS"]],
                  aes(x = time, y = value),
                  color = "black", linewidth = 1.2) +
        theme_bw() +
        theme(text = element_text(size = 16),
              panel.grid = element_blank()) +
        scale_y_continuous(name = expression('Susceptible'['SS']), limits = c(0, NA)) +
        scale_x_continuous(name = 'time (months)',
                           breaks = seq(0, years * 12, 12),
                           limits = c(0, years * 12),
                           expand = c(0, 0))
      
      
      plot <- ggarrange(Splot, Eplot, Iplot, SSSplot, ncol = 3, nrow=2)
      
      plot
      
      ggsave(paste0("figures/supp/", infType,"_classes_",size,"_supp.pdf"), plot, width = 10, height = 7, dpi = 300)
    }
    
    remove(plot, plot_scl)
    remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  }
  
  if (save_runs) {
    name_out <- paste0("std_", infType, "_", pct, "-", prop_superSpreader, "-it", integrate_type, "-")
    save(models_res, file = file.path("data/raw", paste0(name_out, size, "-", Sys.Date(), ".RData")))
  }
  
  rm(data, data_sto)
  return(sto_out)
}
