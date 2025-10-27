###############################################################
## Section selector for SLURM array jobs
###############################################################
args <- commandArgs(trailingOnly = TRUE)
section_id <- ifelse(length(args) > 0, as.numeric(args[1]), 1)

cat("\n====================================\n")
cat(" Running section:", section_id, "\n")
cat(" Hostname:", Sys.info()[["nodename"]], "\n")
cat(" Timestamp:", Sys.time(), "\n")
cat("====================================\n\n")

###############################################################
## Package setup and environment
###############################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(ggplot2)
  library(ggpubr)
  library(foreach)
  library(doParallel)
  library(RColorBrewer)
  library(scales)
})

source("bTBwl_func.R")

# Configure parallel backend
# n.cores <- max(1, floor(parallel::detectCores() * 3/4))
n.cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

###############################################################
## SECTION 1: SS proportion × herd size
###############################################################
if (section_id == 1) {
  cat("Running Section 1: Super-spreader proportion × herd size\n")
  
  sizes  <- seq(from = 0, to = 750, length.out = 21)[-1]
  xi_vec <- seq(from = 0, to = 1, length.out = 21)
  
  param_grid <- expand_grid(size = sizes, infType = "spillover", SS_prop = xi_vec)
  
  results <- foreach(row = iter(param_grid, by = "row"),
                     .packages = c("tidyverse", "deSolve"), .combine = "rbind") %dopar% {
                       run_wl_sim(size = row$size, infType = row$infType, years = 1, reps = 500,
                                  pct = 0, prop_superSpreader = row$SS_prop,
                                  save_plots = FALSE, save_runs = FALSE) |>
                         filter(tstep == max(tstep)) |>
                         mutate(prev = `Total Infected` / N, size = row$size, prop_ss = row$SS_prop) |>
                         group_by(size, prop_ss) |>
                         summarise(mean_prev = mean(prev))
                     }
  
  write.csv(results, "data/ss_herdsize_meanprev.csv", row.names = FALSE)
}

###############################################################
## SECTION 2: Fadeout vs. initial exposure
###############################################################
if (section_id == 2) {
  cat("Running Section 2: Fadeout vs. initial exposure\n")
  
  sizes   <- seq(from = 0, to = 750, length.out = 21)[-1]
  initNum <- 0:25
  
  param_grid <- expand_grid(size = sizes, infType = "seeded", numExp = initNum)
  
  results <- foreach(row = iter(param_grid, by = "row"),
                     .packages = c("tidyverse", "deSolve"), .combine = "rbind") %dopar% {
                       run_wl_sim(size = row$size, infType = row$infType, years = 20, reps = 500,
                                  save_plots = FALSE, pct = row$numExp, save_runs = FALSE) |>
                         mutate(size = row$size, num_exp = row$numExp) |>
                         group_by(size, num_exp) |>
                         summarise(fadeout = mean(fadeout))
                     }
  
  write.csv(results, "data/herdsize_initnum_fadeout.csv", row.names = FALSE)
}

###############################################################
## SECTION 3: Herd size × harvest intensity
###############################################################
if (section_id == 3) {
  cat("Running Section 3: Herd size × harvest intensity\n")
  
  sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
  hunt_pct <- seq(from = 0, to = 1, length.out = 21)
  
  param_grid <- expand_grid(size = sizes, infType = "seeded", hunt = hunt_pct) |>
    mutate(overrides = purrr::map(hunt, ~list(eta_hunt = .x)))
  
  results <- foreach(row = iter(param_grid, by = "row"),
                     .packages = c("tidyverse", "deSolve"), .combine = "rbind") %dopar% {
                       df <- run_wl_sim(size = row$size, infType = row$infType, years = 20, reps = 500,
                                        pct = 0.025 * row$size, prop_superSpreader = 0,
                                        save_plots = FALSE, save_runs = FALSE,
                                        overrides = row$overrides[[1]]) %>%
                         mutate(size = row$size, pct_hunt = row$overrides[[1]]$eta_hunt,
                                prev_ratio = `Hunt Prevalence`)
                       
                       prev_sum <- df %>%
                         filter(quarter == 4) %>%
                         group_by(size, pct_hunt) %>%
                         summarise(prev_ratio = mean(prev_ratio, na.rm = TRUE), .groups = "drop")
                       
                       fade_sum <- df %>%
                         group_by(rep, size, pct_hunt) %>%
                         slice_max(tstep, with_ties = FALSE) %>%
                         ungroup() %>%
                         group_by(size, pct_hunt) %>%
                         summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                       
                       left_join(prev_sum, fade_sum, by = c("size", "pct_hunt"))
                     }
  
  results <- results |> mutate(fadeout = ifelse(fadeout == 0, NA, fadeout))
  write.csv(results, "data/hunt_herdsize.csv", row.names = FALSE)
}

###############################################################
## SECTION 4: Reduction in p1 × reduction in p2
###############################################################
if (section_id == 4) {
  cat("Running Section 4: Reduction in beta × reduction in p2\n")
  
  beta_baseline <- 0.05*0.5
  baseline_p2 <- list(p2_q1 = 0.2, p2_q2 = 0.3, p2_q3 = 0.00, p2_q4 = 0.15)
  scale_p1 <- function(r) list(deer_deer_contact_rate = (1 - r) * beta_baseline)
  scale_p2 <- function(r) {
    f <- 1 - r
    list(p2_q1 = f * baseline_p2$p2_q1, p2_q2 = f * baseline_p2$p2_q2,
         p2_q3 = f * baseline_p2$p2_q3, p2_q4 = f * baseline_p2$p2_q4)
  }
  
  p1_red_grid <- seq(0, 1, length.out = 21)
  farm_red_grid <- seq(0, 1, length.out = 21)
  hunt_levels <- c(0.2, 0.5, 0.95)
  size_K <- 250
  scenarios <- "spillover"
  
  param_grid <- expand_grid(p1_red = p1_red_grid, farm_red = farm_red_grid,
                            hunt = hunt_levels, infType = scenarios) %>%
    mutate(overrides = pmap(list(p1_red, farm_red, hunt),
                            function(p1r, red2, hnt) c(
                              scale_p1(p1r), scale_p2(red2), list(eta_hunt = hnt)
                            )))
  
  results <- foreach(row = iter(param_grid, by = "row"),
                     .packages = c("tidyverse", "deSolve"),
                     .combine = "rbind") %dopar% {
                       df <- run_wl_sim(size = size_K, infType = row$infType, years = 1, reps = 500,
                                        pct = 0, prop_superSpreader = 0, save_plots = FALSE,
                                        save_runs = FALSE, overrides = row$overrides[[1]]) %>%
                         mutate(p1_red = row$p1_red, farm_red = row$farm_red,
                                pct_hunt = row$overrides[[1]]$eta_hunt, scenario = row$infType)
                       
                       prev_sum <- df %>%
                         group_by(rep, p1_red, farm_red, pct_hunt, scenario) %>%
                         slice_max(tstep, with_ties = FALSE) %>%
                         ungroup() %>%
                         mutate(prev = `Total Infected` / N) %>%
                         group_by(p1_red, farm_red, pct_hunt, scenario) %>%
                         summarise(mean_prev = mean(prev, na.rm = TRUE), .groups = "drop")
                       
                       fade_sum <- df %>%
                         group_by(rep, p1_red, farm_red, pct_hunt, scenario) %>%
                         slice_max(tstep, with_ties = FALSE) %>%
                         ungroup() %>%
                         group_by(p1_red, farm_red, pct_hunt, scenario) %>%
                         summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                       
                       left_join(prev_sum, fade_sum, by = c("p1_red", "farm_red", "pct_hunt", "scenario"))
                     }
  
  results <- results |>
    mutate(fadeout = ifelse(fadeout == 0, NA_real_, fadeout))
  write.csv(results, "data/beta_p2_hunt.csv", row.names = FALSE)
}

###############################################################
## SECTION 5: Sensitivity analysis (LHS)
###############################################################
if (section_id == 5) {
  cat("Running Section 5: Sensitivity analysis (LHS)\n")
  ##-----------------------------------------------------------
  ## Sensitivity analysis
  ##-----------------------------------------------------------
  
  ######################
  ## Setup Parameters ##
  ######################
  # rm(list = ls())
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
  lambda_factor = 1.2
  save_pars = T
  save_runs = F
  save_plots = F
  save_plots_final = T
  recalc = T
  averaged = T
  ##parameters excluded from analysis: 
  # omega - birth pulse timing, verbose, integral type
  ######################
  
  ###############################
  ## Generate or load LHS pars ##
  ###############################
  # see if lh is in workspace or sensitivity_analysis folder; otherwise, generate new LHS set
  lh <- LHS_ParSets(k_lim = size_range, 
                    scenario = infection,
                    infected = pct,
                    pct = T,
                    SS_prop = prop_superSpreader,
                    verbose = 0,
                    Num_LHS_sets = N_LHS_sets,
                    file.path = 'data/',
                    file.out = save_pars,
                    file.name = par_file,
                    seed_q = fix_q,
                    seed = def_seed)
  
  #########################################
  ## run model and generate summary file ##
  #########################################
  # run with LHS parameters and generate summary file
  
  # Use to ID which param sets result in HUGE values for lambda
  get_lambda_quick <- function(pars, years = 1, reps = 10) {
    initial_state <- data.frame(
      S_0 = pars["S_0"], E1_0 = pars["E1_0"], I_0 = pars["I_0"],
      SuperS_0 = pars["SuperS_0"], SuperE1_0 = pars["SuperE1_0"], SuperI_0 = pars["SuperI_0"]
    )
    
    population_parameters <- data.frame(
      K = pars["K"], eta_hunt = pars["eta_hunt"], eta_nat = pars["eta_nat"],
      theta = pars["theta"], gamma = pars["gamma"], alpha_max = pars["alpha_max"],
      ksi = pars["ksi"], omega = pars["omega"], s = pars["s"], alpha = 0
    )
    
    disease_parameters <- data.frame(
      beta = pars["beta"], area = 1, p1 = pars["p1"],
      p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], p2_q4 = pars["p2_q4"],
      phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], sigma1_rate = pars["sigma1_rate"]
    )
    
    param_vals <- cbind(population_parameters, disease_parameters)
    
    # Short run, low verbosity, small n_reps
    sto_out <- wildlife_model(
      n_reps = reps, parameters = param_vals, initial_state = initial_state,
      nyears = years, verbose = -20, type = 'c', seed_quarter = 1, batch_name = "foo"
    )
    lambda_out <- getLambda_vec(sto_out, type = 'max')
    return(lambda_out)
  }
  
  
  library(doParallel)
  n.cores <- floor(detectCores() * 3/4)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  lambda_screen <- foreach(i = 1:nrow(lh), .combine = 'rbind') %dopar% {
    pars <- lh[i, ]
    lambda_out <- tryCatch({
      get_lambda_quick(pars, years = 4, reps = 10)
    }, error = function(e) NA_real_)
    data.frame(i = i, lambda = lambda_out)
  }
  
  
  
  exclude <- lambda_screen |> filter(lambda >= 1e4) |> pull(i) |> unique()
  
  if(def_seed == 43 & N_LHS_sets == 1000 & dim(lh)[1] == N_LHS_sets){
    lh <- lh[-c(exclude),] # lh[-c(363,367),] #remove these 2 entries due to excessive runtime issues -- getting massive vals for lambda O(10^15)
    write.csv(lh, paste0("data/",par_file,".csv"))
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
    disease_parameters<-data.frame(beta=pars["beta"], area=1, p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
    
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
    # if(save_runs){
    #   write.csv(sto_out, file = paste0('./sensitivity_analysis/sens_runs/sens_run_',infType, '_', i, '.csv'))
    # }
    
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
  
  
  
  TOC <- Sys.time()
  TOC-TIC
  
  # scale parameter values
  bTB_wl_scaled <- as.data.frame(R %>% mutate_at(vars('K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                                                      'alpha_max', 'ksi', 'omega', 's',
                                                      'beta', "p1", 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 
                                                      'phi', 'sigma1_mean', 'sigma1_rate', 'start_q'),
                                                 list(~scale(as.vector(.)))))
  
  write.csv(R, file = paste0("./data/LHS_summary_", infType, ".csv" ), row.names = F)
  write.csv(bTB_wl_scaled, file = paste0("./data/LHS_scaled_summary_", infType, ".csv" ), row.names = F)
}

if (section_id == 6) {
  # ---------------- grids ----------------
  sizes   <- c(250)                      # facet columns
  xi_vec  <- seq(0, 1, length.out = 31)          # x-axis (prop_superSpreader = ξ)
  phi_vec <- seq(1, 50,   length.out = 31)          # y-axis (SS farm-contact factor = φ)
  
  param_grid <- tidyr::expand_grid(
    size    = sizes,
    infType = "spillover",
    SS_prop = xi_vec,
    phi     = phi_vec
  ) %>%
    mutate(overrides = purrr::map(phi, ~list(phi = .x)))
  
  # -------------- parallel backend --------------
  n.cores <- max(1, floor(parallel::detectCores() * 3/4))
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  # -------------- run simulations --------------
  results <- foreach(row = iter(param_grid, by = "row"),
                     .packages = c("tidyverse", "deSolve"),
                     .combine  = "rbind") %dopar% {
                       gc()
                       
                       df <- run_wl_sim(
                         size               = row$size,
                         infType            = row$infType,
                         years              = 1,
                         reps               = 500,
                         pct                = 0,                  # spillover snapshot (no initial infection)
                         prop_superSpreader = row$SS_prop,        # ξ
                         save_plots         = FALSE,
                         save_runs = FALSE,
                         overrides          = row$overrides[[1]]  # φ
                       ) %>%
                         mutate(size = row$size,
                                prop_ss = row$SS_prop,
                                phi     = row$overrides[[1]]$phi)
                       
                       # Final-timestep per rep
                       final_by_rep <- df %>%
                         group_by(rep, size, prop_ss, phi) %>%
                         slice_max(tstep, with_ties = FALSE) %>%
                         ungroup()
                       
                       # Summaries
                       prev_sum <- final_by_rep %>%
                         mutate(prev = `Total Infected` / N) %>%
                         group_by(size, prop_ss, phi) %>%
                         summarise(mean_prev = mean(prev, na.rm = TRUE), .groups = "drop")
                       
                       fade_sum <- final_by_rep %>%
                         group_by(size, prop_ss, phi) %>%
                         summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                       
                       rm(df, final_by_rep); gc()
                       left_join(prev_sum, fade_sum, by = c("size", "prop_ss", "phi"))
                     }
  write.csv(results, "data/xi_phi_meanprev.csv")
}

stopCluster(cl)
cat("Section", section_id, "complete at:", Sys.time(), "\n")
