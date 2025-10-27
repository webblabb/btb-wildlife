
## -----------------------------------------------------------
## Hunting (eta_hunt) x Farm-contact reduction
## Scenarios: seeded, spillover
## Metrics: p(fadeout), mean infection prevalence (final tstep)
## -----------------------------------------------------------

rm(list = ls()); gc()
source("bTBwl_func.R")

library(tidyverse)
library(doParallel)
library(deSolve)   # if run_wl_sim needs it

# ---------- parallel backend ----------
n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# ---------- baseline farm-contact probabilities ----------
baseline_p2 <- list(
  p2_q1 = 0.01,
  p2_q2 = 0.015,
  p2_q3 = 0.00,
  p2_q4 = 0.0075
)

scaled_p2 <- function(reduction) {
  f <- 1 - reduction
  list(
    p2_q1 = f * baseline_p2$p2_q1,
    p2_q2 = f * baseline_p2$p2_q2,
    p2_q3 = f * baseline_p2$p2_q3,
    p2_q4 = f * baseline_p2$p2_q4
  )
}

# ---------- grids ----------
year_grid  <- 1 # c(1, 5, 10)   # deer_deer_contact_rate (p1)
hunt_grid     <- seq(0, 1, length.out = 11)   # eta_hunt (Q4 harvest proportion)
farm_red_grid <- seq(0, 1, length.out = 11) # c(0.00, 0.50, 1.00)          # 0%, 50%, 100% reduction
scenarios     <- c("spillover")

# Fix herd size to make a 2D surface per facet
size_K  <- 250
pctSeed <- 0 # 0.025 * size_K

param_grid <- tidyr::expand_grid(
  years  = year_grid,
  hunt     = hunt_grid,
  farm_red = farm_red_grid,
  infType  = scenarios
) %>%
  mutate(
    overrides = purrr::pmap(
      list(hunt, farm_red),
      function(hnt, red) c(
        list(
          eta_hunt               = hnt
        ),
        scaled_p2(red)
      )
    )
  )

# ---------- run experiments ----------
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse", "deSolve"),
                   .combine  = "rbind") %dopar% {
                     
                     # # Scenario-specific settings
                     # if (row$infType == "seeded") {
                     #   years <- 20; reps <- 300; pct <- pctSeed
                     # } else { # spillover
                     #   years <- 1;  reps <- 250; pct <- 0
                     # }
                     
                     df <- run_wl_sim(
                       size                = size_K,
                       infType             = row$infType,
                       years               = row$years,
                       reps                = 250,
                       pct                 = 0,
                       prop_superSpreader  = 0,
                       save_plots          = FALSE,
                       overrides           = row$overrides[[1]]
                     ) %>%
                       mutate(
                         years  = row$years,
                         pct_hunt = row$overrides[[1]]$eta_hunt,
                         farm_red = row$farm_red,
                         scenario = row$infType
                       )
                     
                     # --- mean infection prevalence (final timestep per rep) ---
                     prev_sum <- df %>%
                       group_by(rep, years , pct_hunt, farm_red, scenario) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup() %>%
                       mutate(prev = `Total Infected` / N) %>%
                       group_by(years, pct_hunt, farm_red, scenario) %>%
                       summarise(mean_prev = mean(prev, na.rm = TRUE), .groups = "drop")
                     
                     # --- p(fadeout) at final timestep ---
                     fade_sum <- df %>%
                       group_by(rep, years, pct_hunt, farm_red, scenario) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup() %>%
                       group_by(years , pct_hunt, farm_red, scenario) %>%
                       summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                     
                     rm(df); gc()
                     left_join(prev_sum, fade_sum,
                               by = c("years","pct_hunt","farm_red","scenario"))
                   }

stopCluster(cl)

# Optional: show zeros as NA (gray) in fadeout map
results <- results %>%
  mutate(
    fadeout        = ifelse(fadeout == 0, NA_real_, fadeout),
    farm_red_label = factor(paste0(round(100 * farm_red), "% reduction"),
                            levels = farm_red_grid, "% reduction"),
    scenario_label = factor(scenario, levels = c("seeded","spillover"))
  )

# ---------- save data ----------
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "data/processed/contact_hunt_farmReduction_seeded_spillover.csv",
          row.names = FALSE)

# ---------- figures ----------
library(ggplot2)

# p(fadeout)
fade_prob_map <- ggplot(
  data = results,
  aes(x = farm_red, y = pct_hunt, fill = fadeout)
) +
  geom_tile() +
  scale_fill_distiller(
    palette   = "Spectral",
    direction = -1,
    na.value  = "gray21",
    breaks    = round(seq(0, max(results$fadeout, na.rm = TRUE), length.out = 4), 2),
    limits    = c(0, max(results$fadeout, na.rm = TRUE)),
    name      = expression(italic(p) * "(fadeout)")
  ) +
  labs(
    x = "% reduction in cattle and fomite contacts",
    y = expression("percent harvested (" * eta[h] * ")")
  ) +
  # facet_grid(rows = vars(scenario_label), cols = vars(farm_red_label)) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold")
  ) +
  coord_cartesian(expand = FALSE)

fade_prob_map

# mean infection prevalence
mean_prev_map <- ggplot(
  data = results,
  aes(x = farm_red, y = pct_hunt, fill = mean_prev)
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    breaks = round(seq(0, max(results$mean_prev, na.rm = TRUE), length.out = 5), 3),
    limits = c(0, max(results$mean_prev, na.rm = TRUE)),
    name   = "mean infection\nprevalence"
  ) +
  labs(
    x = "% reduction in cattle and fomite contacts",
    y = expression("percent harvested (" * eta[h] * ")")
  ) +
  # facet_grid(rows = vars(scenario_label), cols = vars(farm_red_label)) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold")
  ) +
  coord_cartesian(expand = FALSE)

mean_prev_map

# save figs
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/fade_prob_map_contact_hunt_farmReduction_byScenario.png",
       fade_prob_map, width = 12, height = 7.5, dpi = 300)
ggsave("figures/mean_prev_map_contact_hunt_farmReduction_byScenario.png",
       mean_prev_map, width = 12, height = 7.5, dpi = 300)

param_grid <- expand_grid(
  size = c(25, 250, 500, 750),
  infType = c("seeded", "spillover")
)

for (i in 1:nrow(param_grid)) {
  row <- param_grid[i,]
  
  run_wl_sim(size = row$size, infType = row$infType) 
}


##########################################################################################
##########################################################################################
#######################     Multiple Herd Sizes -- Discrete time    ######################
##########################################################################################
##########################################################################################
system("g++ -L/usr/lib/x86_64-linux-gnu tb_wildlife_freq_cont.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu tb_wildlife_freq_disc.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")

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
supp = T
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
    func = SEI_model_full_freq,
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
  
  
  data_sto <- sto_out |> 
    select("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I") |> 
    # mutate(across(
    #   c(S, SS_S, E1, SS_E1, I, SS_I),     # class columns to scale
    #   ~ .x / N,                           # divide by N for that row
    #   .names = "{.col}"              # optional: create new *_prop columns
    # )) |>
    pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  #assign data frames to list
  #then clear data from workspace
  
  models_res[[1]] <- as.data.frame(data)
  remove(data)
  models_res[[2]] <- as.data.frame(sto_out)
  remove(sto_out)
  models_res[[3]] <- as.data.frame(data_sto)
  remove(data_sto)
  
  ##########
  
  # data_scaled <- data %>%
  #   group_by(time) %>%                                  # group by time
  #   mutate(total_N = value[variable == "N"],            # extract N for that time
  #          value = value / total_N) %>%                 # divide all classes by N
  #   ungroup() %>%
  #   filter(variable %in% c("S", "sS", "E", "I")) 
  # 
  # # map variable names to clear class labels
  # rename_map <- c(
  #   "S"      = "Susceptible",
  #   "SS_S"   = "Susceptible_SS",
  #   "sS"   = "Susceptible_SS",
  #   "E1"     = "Exposed",
  #   "E"     = "Exposed",
  #   "SS_E1"  = "Exposed_SS",
  #   "I"      = "Infectious",
  #   "SS_I"   = "Infectious_SS",
  #   "N"      = "Total"
  # )
  # 
  # data_sto <- data_sto |> 
  #   filter(variable != "N") %>%
  #   mutate(variable = recode(variable, !!!rename_map))
  # 
  # data_scaled <- data_scaled |> 
  #   filter(variable != "N") %>%
  #   mutate(variable = recode(variable, !!!rename_map))
  # 
  # ggplot() + 
  #   geom_line(data = data_sto, aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) + 
  #   geom_line(data = data_scaled, aes(x = time, y=value), size = 1, show.legend = FALSE) +
  #   scale_fill_distiller(palette = "Spectral") +
  #   # scale_color_manual(values=SEIcols) +
  #   facet_grid(~variable) +
  #   theme_bw() +
  #   theme(text = element_text(size = 16),
  #         panel.grid = element_blank())
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
  
  
  # Nplot <- ggplot() + 
  #   geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) + 
  #   geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
  #   scale_color_manual(values=SEIcols[7]) +
  #   scale_y_continuous(limits = c(0,NA)) +
  #   scale_x_continuous(breaks=seq(0,years*12,1), limits = c(0, years*12))
  
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
  
  if(save_plots){
    setwd(pth)
    
    plot <- ggarrange(Splot, Eplot, Iplot,SSSplot, ncol=3, nrow=2)
    plot_scl <- ggarrange(Splot_scl, Eplot_scl,  Iplot_scl, SSSplot_scl, ncol=3, nrow=2)
    
    plot
    plot_scl
    
    ggsave(paste0(infType,"_classes_",size,".pdf"), plot, width = 10, height = 7, dpi = 300)
    ggsave(paste0(infType,"_classes_",size,"_scaled.pdf"), plot_scl, width = 10, height = 7, dpi = 300)
    
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
      
      
      plot <- ggarrange(Splot, Eplot, Iplot,SSSplot, ncol=3, nrow=2)
      
      plot
      
      ggsave(paste0(infType,"_classes_",size,"_supp.pdf"), plot, width = 10, height = 7, dpi = 300)
    }
    
    remove(plot, plot_scl)
  }
  remove(Splot, SSSplot, Eplot, ESSplot, Iplot, ISSplot, Nplot)
  ###########
  
  
  ##################
  ## save outputs ##
  ##################
  
  if(save_runs){
    save(models_res, file = paste0('Disc_runs/', name_out, size, '-', Sys.Date(), '.RData'))
  }
  remove(models_res) #clear from workspace
  ##################
  
}

stopCluster(cl)









library(dplyr)

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


plot <- ggarrange(Splot, Eplot, Iplot,SSSplot, ncol=3, nrow=2)





#-----------------------\

# #############################################
# ## Hunting × Super-Susceptible & Contact Maps
# #############################################
# rm(list = ls())
# setwd(dirname(getActiveDocumentContext()$path))
# source("bTBwl_func.R")
# system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_DTMC.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")
# 
# years <- 20
# seedQuarter <- 1
# verbose <- 0
# lambda_factor <- 1.2
# reps <- 500
# herd_size <- 250
# infType <- "seeded"
# type_of_integral <- 3
# pth <- "/home/webblab/Documents/Brandon/bTB_wildlife_code/"
# save_runs <- TRUE
# 
# # Parameter ranges
# harvest_vec <- seq(0, 0.5, length.out = 21)[-1]   # annual harvest rates
# xi_vec      <- seq(0, 0.5, length.out = 21)[-1]   # super susceptible proportion
# phi_vec     <- seq(1, 50, length.out = 20)    # contact scaling
# 
# # Parallel setup
# n.cores <- floor(detectCores() * 3/4)
# cl <- makeCluster(n.cores)
# registerDoParallel(cl)
# clusterEvalQ(cl, { library(tidyverse); library(deSolve) })
# 
# #######################################################
# ## (1) Harvest × Super-Susceptible proportion (ξ)
# #######################################################
# R_xi <- foreach(i = 1:length(harvest_vec), .combine = 'rbind', .inorder = TRUE) %:%
#   foreach(j = 1:length(xi_vec), .combine = 'rbind', .inorder = TRUE) %dopar% {
#     
#     h_val  <- harvest_vec[i]
#     xi_val <- xi_vec[j]
#     
#     pars <- parameter_set_wl(k = herd_size, scenario = infType,
#                              initial_exposed = 0.02 * herd_size,
#                              SS_prop = xi_val,
#                              start_quarter = seedQuarter,
#                              verbose = verbose)
#     
#     # init_state <- data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"],
#     #                          SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
#     # 
#     # pop_pars <- data.frame(K=pars["K"], eta_hunt=h_val, eta_nat=pars["eta_nat"],
#     #                        theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"],
#     #                        ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"])
#     # 
#     # dis_pars <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"],
#     #                        p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"],
#     #                        phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
#     # 
#     # param_vals <- data.frame(merge(pop_pars, dis_pars))
#     
#     #stochastic model
#     initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
#     population_parameters<-data.frame(K=pars["K"], eta_hunt=h_val, eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"])
#     disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
# 
#     param_vals<-data.frame(merge(population_parameters,disease_parameters))
#     
#     #get lambda
#     tic = Sys.time()
#     sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 3, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
#     toc = Sys.time()
#     print(toc-tic)
#     lambda_out <- getLambda_vec(data = sto_out, type = 'max')
#     remove(sto_out)
#     
#     #run with new lambda
#     tic = Sys.time()
#     sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
#     toc = Sys.time()
#     print(toc-tic)
#     remove(initial_state,population_parameters,disease_parameters,param_vals)
#     
#     # summarize across reps
#     fade_prob <- mean(sto_out$fadeout)
#     mean_ratio <- mean(sto_out$`Hunt Prevalence`, na.rm=TRUE)
#     
#     data.frame(harvest=h_val, xi=xi_val,
#                fadeout=fade_prob,
#                mean_ratio=mean_ratio)
#   }
# 
# #######################################################
# ## (2) Harvest × Contact scaling (ϕ)
# #######################################################
# R_phi <- foreach(i = 1:length(harvest_vec), .combine = 'rbind', .inorder = TRUE) %:%
#   foreach(j = 1:length(phi_vec), .combine = 'rbind', .inorder = TRUE) %dopar% {
#     
#     h_val  <- harvest_vec[i]
#     phi_val <- phi_vec[j]
#     
#     pars <- parameter_set_wl(k = herd_size, scenario = infType,
#                              initial_exposed = 0.02 * herd_size,
#                              SS_prop = 0.05,                # fixed ξ baseline
#                              start_quarter = seedQuarter,
#                              verbose = verbose)
#     
#     #stochastic model
#     initial_state<-data.frame(S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"], SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"])
#     population_parameters<-data.frame(K=pars["K"], eta_hunt=h_val, eta_nat=pars["eta_nat"], theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"], ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"])
#     disease_parameters<-data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"], p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=phi_val, sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
#     
#     param_vals<-data.frame(merge(population_parameters,disease_parameters))
#     
#     #get lambda
#     tic = Sys.time()
#     sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = 3, seed_quarter = as.integer(pars["start_q"]), verbose = -3, batch_name = "test", type = 'c')
#     toc = Sys.time()
#     print(toc-tic)
#     lambda_out <- getLambda_vec(data = sto_out, type = 'max')
#     remove(sto_out)
#     
#     #run with new lambda
#     tic = Sys.time()
#     sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = initial_state, nyears = years, seed_quarter = as.integer(pars["start_q"]), verbose = verbose, batch_name = "test", type = 'd', lambda = lambda_out*lambda_factor, integrate_type = type_of_integral)
#     toc = Sys.time()
#     print(toc-tic)
#     remove(initial_state,population_parameters,disease_parameters,param_vals)
#     
#     fade_prob <- mean(sto_out$fadeout)
#     mean_ratio <- mean(sto_out$`Hunt Prevalence`, na.rm=TRUE)
#     
#     data.frame(harvest=h_val, phi=phi_val,
#                fadeout=fade_prob,
#                mean_ratio=mean_ratio)
#   }
# 
# stopCluster(cl)
# 
# #######################################################
# ## ---- Plotting (like your prior heatmaps) ----
# #######################################################
# library(reshape2)
# library(ggpubr)
# 
# # Harvest × Super-Susceptible
# fade_xi_melt <- R_xi %>% select(harvest, xi, fadeout)
# ratio_xi_melt <- R_xi %>% select(harvest, xi, mean_ratio)
# 
# fade_map_xi <- ggplot(fade_xi_melt, aes(x=harvest, y=xi, fill=fadeout)) +
#   geom_tile(color=NA) +
#   scale_fill_distiller(palette="Spectral", direction=-1, name=expression(italic(p) * "(fadeout)")) +
#   labs(x="percent harvested", y="super susceptible proportion (ξ)") +
#   theme_bw(base_size=14) +
#   theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
#         panel.grid=element_blank(),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=16))
# 
# ratio_map_xi <- ggplot(ratio_xi_melt, aes(x=harvest, y=xi, fill=mean_ratio)) +
#   geom_tile(color=NA) +
#   scale_fill_distiller(palette="Spectral", direction=-1, name="mean estimated prevalence ratio") +
#   labs(x="percent harvested", y="super susceptible proportion (ξ)") +
#   theme_bw(base_size=14) +
#   theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
#         panel.grid=element_blank(),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=16))
# 
# ratio_map_xi
# fade_map_xi
# 
# # Harvest × Contact scaling
# fade_phi_melt <- R_phi %>% select(harvest, phi, fadeout)
# ratio_phi_melt <- R_phi %>% select(harvest, phi, mean_ratio)
# 
# fade_map_phi <- ggplot(fade_phi_melt, aes(x=harvest, y=phi, fill=fadeout)) +
#   geom_tile(color=NA) +
#   scale_fill_distiller(palette="Spectral", direction=-1, name=expression(italic(p) * "(fadeout)")) +
#   labs(x="percent harvested", y="contact scaling (ϕ)") +
#   theme_bw(base_size=14) +
#   theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
#         panel.grid=element_blank(),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=16))
# 
# ratio_map_phi <- ggplot(ratio_phi_melt, aes(x=harvest, y=phi, fill=mean_ratio)) +
#   geom_tile(color=NA) +
#   scale_fill_distiller(palette="Spectral", direction=-1, name="mean estimated prevalence ratio") +
#   labs(x="percent harvested", y="contact scaling (ϕ)") +
#   theme_bw(base_size=14) +
#   theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
#         panel.grid=element_blank(),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=16))
# 
# ratio_map_phi
# fade_map_phi
# 
# # Save
# ggsave("fade_prob_map_harvest_xi.png", fade_map_xi, width=7.5, height=6, dpi=600)
# ggsave("prev_ratio_map_harvest_xi.png", ratio_map_xi, width=7.5, height=6, dpi=600)
# ggsave("fade_prob_map_harvest_phi.png", fade_map_phi, width=7.5, height=6, dpi=600)
# ggsave("prev_ratio_map_harvest_phi.png", ratio_map_phi, width=7.5, height=6, dpi=600)
# 
# # Optional composite
# ggarrange(ratio_map_xi, fade_map_xi, ratio_map_phi, fade_map_phi,
#           ncol=2, nrow=2,
#           labels=c("(a)", "(b)", "(c)", "(d)"))


#############################################
## Hunting × Super-Susceptible & Contact Maps (by Herd Size)
#############################################
rm(list = ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("bTBwl_func.R")

# --- Base configuration ---
years <- 20
seedQuarter <- 1
verbose <- 0
lambda_factor <- 1.2
reps <- 200
infType <- "seeded"
type_of_integral <- 3
pth <- "/home/webblab/Documents/Brandon/bTB_wildlife_code/"
save_runs <- TRUE

# Parameter ranges
herd_sizes <- c(25, 250, 500, 750)
harvest_vec <- seq(0, 0.5, length.out = 21)[-1]   # harvest rate
xi_vec      <- seq(0, 0.5, length.out = 21)[-1]   # ξ proportion
phi_vec     <- seq(1, 50, length.out = 20)        # ϕ contact scaling

# Parallel setup
n.cores <- floor(detectCores() * 3/4)
cl <- makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { library(tidyverse); library(deSolve) })

###########################################
## ---------- Outer Herd Loop ------------
###########################################
for (herd_size in herd_sizes) {
  
  message(paste0("Running simulations for herd size = ", herd_size))
  
  #----------------------------------------
  # (1) Harvest × Super-Susceptible proportion (ξ)
  #----------------------------------------
  R_xi <- foreach(i = 1:length(harvest_vec), .combine = 'rbind', .inorder = TRUE) %:%
    foreach(j = 1:length(xi_vec), .combine = 'rbind', .inorder = TRUE) %dopar% {
      
      h_val  <- harvest_vec[i]
      xi_val <- xi_vec[j]
      
      pars <- parameter_set_wl(k = herd_size, scenario = infType,
                               initial_exposed = 0.02 * herd_size,
                               SS_prop = xi_val,
                               start_quarter = seedQuarter,
                               verbose = verbose)
      
      # --- stochastic model ---
      init_state <- data.frame(
        S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"],
        SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"]
      )
      
      pop_pars <- data.frame(K=pars["K"], eta_hunt=h_val, eta_nat=pars["eta_nat"],
                             theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"],
                             ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"])
      
      dis_pars <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"],
                             p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"],
                             phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
      
      param_vals <- data.frame(merge(pop_pars, dis_pars))
      
      # --- calibrate λ ---
      sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = init_state,
                                nyears = 3, seed_quarter = as.integer(pars["start_q"]),
                                verbose = -3, batch_name = "test", type = 'c')
      lambda_out <- getLambda_vec(data = sto_out, type = 'max')
      remove(sto_out)
      
      # --- main run ---
      sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = init_state,
                                nyears = years, seed_quarter = as.integer(pars["start_q"]),
                                verbose = verbose, batch_name = "test", type = 'd',
                                lambda = lambda_out * lambda_factor, integrate_type = type_of_integral)
      
      fade_prob  <- mean(sto_out$fadeout)
      mean_ratio <- mean(sto_out$`Hunt Prevalence`, na.rm = TRUE)
      
      data.frame(herd_size = herd_size,
                 harvest = h_val,
                 xi = xi_val,
                 fadeout = fade_prob,
                 mean_ratio = mean_ratio)
    }
  
  #----------------------------------------
  # (2) Harvest × Contact scaling (ϕ)
  #----------------------------------------
  R_phi <- foreach(i = 1:length(harvest_vec), .combine = 'rbind', .inorder = TRUE) %:%
    foreach(j = 1:length(phi_vec), .combine = 'rbind', .inorder = TRUE) %dopar% {
      
      h_val  <- harvest_vec[i]
      phi_val <- phi_vec[j]
      
      pars <- parameter_set_wl(k = herd_size, scenario = infType,
                               initial_exposed = 0.02 * herd_size,
                               SS_prop = 0.05,  # fixed ξ
                               start_quarter = seedQuarter,
                               verbose = verbose)
      
      init_state <- data.frame(
        S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"],
        SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"]
      )
      
      pop_pars <- data.frame(K=pars["K"], eta_hunt=h_val, eta_nat=pars["eta_nat"],
                             theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"],
                             ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"])
      
      dis_pars <- data.frame(beta=pars["beta"], area=pars["area"], p1=pars["p1"],
                             p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"],
                             phi=phi_val, sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
      
      param_vals <- data.frame(merge(pop_pars, dis_pars))
      
      sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = init_state,
                                nyears = 3, seed_quarter = as.integer(pars["start_q"]),
                                verbose = -3, batch_name = "test", type = 'c')
      lambda_out <- getLambda_vec(data = sto_out, type = 'max')
      remove(sto_out)
      
      sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = init_state,
                                nyears = years, seed_quarter = as.integer(pars["start_q"]),
                                verbose = verbose, batch_name = "test", type = 'd',
                                lambda = lambda_out * lambda_factor, integrate_type = type_of_integral)
      
      fade_prob  <- mean(sto_out$fadeout)
      mean_ratio <- mean(sto_out$`Hunt Prevalence`, na.rm = TRUE)
      
      data.frame(herd_size = herd_size,
                 harvest = h_val,
                 phi = phi_val,
                 fadeout = fade_prob,
                 mean_ratio = mean_ratio)
    }
  
  #----------------------------------------
  # ---- Plotting section ----
  #----------------------------------------
  library(ggpubr)
  fade_xi_melt  <- R_xi %>% select(harvest, xi, fadeout)
  ratio_xi_melt <- R_xi %>% select(harvest, xi, mean_ratio)
  fade_phi_melt <- R_phi %>% select(harvest, phi, fadeout)
  ratio_phi_melt<- R_phi %>% select(harvest, phi, mean_ratio)
  
  fade_map_xi <- ggplot(fade_xi_melt, aes(x=harvest, y=xi, fill=fadeout)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1, name=expression(italic(p) * "(fadeout)")) +
    labs(x="percent harvested", y=expression("proportion of super susceptibles ("*xi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16))
  
  ratio_map_xi <- ggplot(ratio_xi_melt, aes(x=harvest, y=xi, fill=mean_ratio)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1, name="mean estimated\nprevalence ratio") +
    labs(x="percent harvested", y=expression("proportion of super susceptibles ("*xi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16))
  
  fade_map_phi <- ggplot(fade_phi_melt, aes(x=harvest, y=phi, fill=fadeout)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1, name=expression(italic(p) * "(fadeout)")) +
    labs(x="percent harvested", y=expression("super susceptible contact scaling ("*phi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16))
  
  ratio_map_phi <- ggplot(ratio_phi_melt, aes(x=harvest, y=phi, fill=mean_ratio)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1, name="mean estimated\nprevalence ratio") +
    labs(x="percent harvested", y=expression("super susceptible contact scaling ("*phi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16))
  
  # Save individual and composite plots
  ggsave(paste0("fade_prob_map_xi_K", herd_size, ".png"), fade_map_xi, width=7.5, height=6, dpi=600)
  ggsave(paste0("prev_ratio_map_xi_K", herd_size, ".png"), ratio_map_xi, width=7.5, height=6, dpi=600)
  ggsave(paste0("fade_prob_map_phi_K", herd_size, ".png"), fade_map_phi, width=7.5, height=6, dpi=600)
  ggsave(paste0("prev_ratio_map_phi_K", herd_size, ".png"), ratio_map_phi, width=7.5, height=6, dpi=600)
  
  # composite <- ggarrange(ratio_map_xi, fade_map_xi, ratio_map_phi, fade_map_phi,
  #                        ncol=2, nrow=2,
  #                        labels=c("(a)", "(b)", "(c)", "(d)"))
  # ggsave(paste0("hunt_xi_phi_maps_K", herd_size, ".png"), composite, width=12, height=10, dpi=600)
  rm(R_phi,R_xi)
  gc()
}

stopCluster(cl)




















#############################################
## Hunting × Super-Susceptible & Contact Maps (param grid + shared legends)
#############################################
rm(list = ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("bTBwl_func.R")

# --- Base configuration ---
years <- 20
seedQuarter <- 1
verbose <- 0
lambda_factor <- 1.2
reps <- 200
infType <- "seeded"
type_of_integral <- 3
pth <- "/home/webblab/Documents/Brandon/bTB_wildlife_code/"
save_runs <- TRUE

# Parameter ranges
herd_sizes <- c(25, 250, 500, 750)
harvest_vec <- seq(0, 0.5, length.out = 21)[-1]   # harvest rate
xi_vec      <- seq(0, 0.5, length.out = 21)[-1]   # ξ proportion
phi_vec     <- seq(1, 50, length.out = 20)        # ϕ contact scaling

# Define parameter grid
param_grid <- expand.grid(herd_size = herd_sizes, grid_type = c("xi", "phi"))

# Parallel setup
n.cores <- floor(detectCores() * 3/4)
cl <- makeCluster(n.cores)
registerDoParallel(cl)
clusterEvalQ(cl, { library(tidyverse); library(deSolve) })

###########################################
## ---------- Run all grids --------------
###########################################
all_results <- list()

for (k in seq_len(nrow(param_grid))) {
  herd_size <- param_grid$herd_size[k]
  grid_type <- param_grid$grid_type[k]
  
  message(paste0("Running simulations for herd size = ", herd_size, " (", grid_type, " grid)"))
  
  if (grid_type == "xi") {
    range_primary <- xi_vec
    fixed_SS <- NA
  } else {
    range_primary <- phi_vec
    fixed_SS <- 0.05
  }
  
  R_out <- foreach(i = 1:length(harvest_vec), .combine = 'rbind', .inorder = TRUE) %:%
    foreach(j = 1:length(range_primary), .combine = 'rbind', .inorder = TRUE) %dopar% {
      
      h_val <- harvest_vec[i]
      var_val <- range_primary[j]
      
      pars <- parameter_set_wl(
        k = herd_size,
        scenario = infType,
        initial_exposed = 0.02 * herd_size,
        SS_prop = ifelse(grid_type == "xi", var_val, fixed_SS),
        start_quarter = seedQuarter,
        verbose = verbose
      )
      
      # --- stochastic model ---
      init_state <- data.frame(
        S_0=pars["S_0"], E1_0=pars["E1_0"], I_0=pars["I_0"],
        SuperS_0=pars["SuperS_0"], SuperE1_0=pars["SuperE1_0"], SuperI_0=pars["SuperI_0"]
      )
      
      pop_pars <- data.frame(
        K=pars["K"], eta_hunt=h_val, eta_nat=pars["eta_nat"],
        theta=pars["theta"], gamma=pars["gamma"], alpha_max=pars["alpha_max"],
        ksi=pars["ksi"], omega=pars["omega"], s=pars["s"], alpha=pars["alpha"]
      )
      
      dis_pars <- data.frame(
        beta=pars["beta"], area=pars["area"], p1=pars["p1"],
        p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"],
        phi=ifelse(grid_type == "phi", var_val, pars["phi"]),
        sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"]
      )
      
      param_vals <- data.frame(merge(pop_pars, dis_pars))
      
      # --- calibrate λ ---
      sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = init_state,
                                nyears = 3, seed_quarter = as.integer(pars["start_q"]),
                                verbose = -3, batch_name = "test", type = 'c')
      lambda_out <- getLambda_vec(data = sto_out, type = 'max')
      remove(sto_out)
      
      # --- main run ---
      sto_out <- wildlife_model(n_reps = reps, parameters = param_vals, initial_state = init_state,
                                nyears = years, seed_quarter = as.integer(pars["start_q"]),
                                verbose = verbose, batch_name = "test", type = 'd',
                                lambda = lambda_out * lambda_factor, integrate_type = type_of_integral)
      
      fade_prob  <- mean(sto_out$fadeout)
      mean_ratio <- mean(sto_out$`Hunt Prevalence`, na.rm = TRUE)
      
      if (grid_type == "xi") {
        data.frame(herd_size = herd_size, harvest = h_val, xi = var_val,
                   fadeout = fade_prob, mean_ratio = mean_ratio, grid = "xi")
      } else {
        data.frame(herd_size = herd_size, harvest = h_val, phi = var_val,
                   fadeout = fade_prob, mean_ratio = mean_ratio, grid = "phi")
      }
    }
  
  all_results[[paste0(grid_type, "_", herd_size)]] <- R_out
}

stopCluster(cl)

###########################################
## ---- Combine and normalize scales -----
###########################################
all_df <- bind_rows(all_results)
fade_max  <- max(all_df$fadeout, na.rm = TRUE)
ratio_max <- max(all_df$mean_ratio, na.rm = TRUE)

###########################################
## ---- Plotting (shared color scales) ----
###########################################
library(ggpubr)

for (herd_size in herd_sizes) {
  fade_xi <- all_df %>% filter(grid == "xi", herd_size == !!herd_size)
  fade_phi <- all_df %>% filter(grid == "phi", herd_size == !!herd_size)
  
  fade_map_xi <- ggplot(fade_xi, aes(x=harvest, y=xi, fill=fadeout)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1,
                         name=expression(italic(p) * "(fadeout)"),
                         limits=c(0, fade_max)) +
    labs(x="percent harvested", y=expression("proportion of super susceptibles ("*xi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank())
  
  ratio_map_xi <- ggplot(fade_xi, aes(x=harvest, y=xi, fill=mean_ratio)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1,
                         name="mean estimated\nprevalence ratio",
                         limits=c(0, ratio_max)) +
    labs(x="percent harvested", y=expression("proportion of super susceptibles ("*xi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank())
  
  fade_map_phi <- ggplot(fade_phi, aes(x=harvest, y=phi, fill=fadeout)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1,
                         name=expression(italic(p) * "(fadeout)"),
                         limits=c(0, fade_max)) +
    labs(x="percent harvested", y=expression("contact scaling ("*phi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank())
  
  ratio_map_phi <- ggplot(fade_phi, aes(x=harvest, y=phi, fill=mean_ratio)) +
    geom_tile(color=NA) +
    scale_fill_distiller(palette="Spectral", direction=-1,
                         name="mean estimated\nprevalence ratio",
                         limits=c(0, ratio_max)) +
    labs(x="percent harvested", y=expression("contact scaling ("*phi*")")) +
    theme_bw(base_size=14) +
    theme(panel.border=element_rect(color="black", fill=NA, size=0.8),
          panel.grid=element_blank())
  
  composite <- ggarrange(ratio_map_xi, fade_map_xi, ratio_map_phi, fade_map_phi,
                         ncol=2, nrow=2, common.legend=TRUE, legend="right",
                         labels=c("(a)", "(b)", "(c)", "(d)")) %>%
    annotate_figure(top=text_grob(paste("Herd size =", herd_size), face="bold", size=16))
  
  ggsave(paste0("hunt_xi_phi_maps_K", herd_size, "_shared.png"), composite, width=12, height=10, dpi=600)
}







# ===
# Vary susc. and contact under spillover conditions

suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(foreach)
  library(doParallel)
  library(scales)
})

dir.create("data",     showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",  showWarnings = FALSE, recursive = TRUE)

source("bTBwl_func.R")

## --- constants ---
years_run <- 20
reps_run  <- 200
sizes     <- seq(0, 750, length.out = 10)[-1]
hunts     <- seq(0, 0.5, length.out = 10)
infType   <- "seeded"

## choose 3 levels (low/med/high)
init_prop_fixed <- 0.05              # keep initial exposure fixed for these analyses
ss_levels  <- c(low = 0.00, med = 0.05, high = 0.20)
phi_levels <- c(low = 1,    med = 10,   high = 30)

## --- parallel ---
n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

## ==========================================================
## A) Vary super-susceptible proportion (ξ): low/med/high
## ==========================================================
param_ss <- tidyr::expand_grid(
  size   = sizes,
  hunt   = hunts,
  ss_cat = names(ss_levels)
) |>
  mutate(
    ss_val    = unname(ss_levels[ss_cat]),
    overrides = purrr::map(hunt, ~list(eta_hunt = .x))
  )

results_ss <- foreach(row = iter(param_ss, by = "row"),
                      .packages = c("tidyverse","deSolve"),
                      .combine = "rbind") %dopar% {
                        
                        df <- run_wl_sim(
                          size               = row$size,
                          infType            = infType,
                          years              = years_run,
                          reps               = reps_run,
                          pct                = init_prop_fixed * row$size,   # fixed initial exposure fraction
                          prop_superSpreader = row$ss_val,
                          save_plots         = FALSE,
                          save_runs          = FALSE,
                          overrides          = row$overrides[[1]]
                        ) %>%
                          mutate(
                            size       = row$size,
                            pct_hunt   = row$overrides[[1]]$eta_hunt,
                            ss_cat     = row$ss_cat,
                            prev_ratio = `Hunt Prevalence`
                          )
                        
                        prev_sum <- df %>%
                          filter(quarter == 4) %>%
                          group_by(size, pct_hunt, ss_cat) %>%
                          summarise(prev_ratio = mean(prev_ratio, na.rm = TRUE), .groups = "drop")
                        
                        fade_sum <- df %>%
                          group_by(rep, size, pct_hunt, ss_cat) %>%
                          slice_max(tstep, with_ties = FALSE) %>%
                          ungroup() %>%
                          group_by(size, pct_hunt, ss_cat) %>%
                          summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                        
                        rm(df); gc()
                        left_join(prev_sum, fade_sum, by = c("size","pct_hunt","ss_cat"))
                      }

write.csv(results_ss, file.path("data", "hunt_by_propSS.csv"), row.names = FALSE)

results_ss <- results_ss %>%
  mutate(ss_lab = forcats::fct_recode(ss_cat,
                                      "ξ = low"  = "low",
                                      "ξ = med"  = "med",
                                      "ξ = high" = "high"))

fade_map_ss <- ggplot(results_ss, aes(x = size, y = pct_hunt, fill = fadeout)) +
  geom_tile() +
  facet_wrap(~ ss_lab, nrow = 1) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       breaks = round(seq(0, max(results_ss$fadeout, na.rm = TRUE), length.out = 4), 2),
                       limits = c(0, max(results_ss$fadeout, na.rm = TRUE)),
                       name = expression(italic(p) * "(fadeout)")) +
  scale_x_continuous(breaks = c(25, 250, 500, 750)) +
  labs(x = "herd size", y = expression("percent harvested ("*eta[h]*")")) +
  theme_bw(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid = element_blank())

fade_map_ss

ggsave(file.path("figures", "fade_prob_map_by_propSS.png"),
       fade_map_ss, width = 12, height = 6, dpi = 300)

mean_map_ss <- ggplot(results_ss, aes(x = size, y = pct_hunt, fill = prev_ratio)) +
  geom_tile(color = NA) +
  facet_wrap(~ ss_lab, nrow = 1) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       breaks = round(seq(0, max(results_ss$prev_ratio, na.rm = TRUE), length.out = 5), 2),
                       limits = c(0, round(max(results_ss$prev_ratio, na.rm = TRUE), 2)),
                       name = "mean estimated\nprevalence ratio") +
  scale_x_continuous(breaks = c(25, 250, 500, 750)) +
  labs(x = "herd size", y = expression("percent harvested ("*eta[h]*")")) +
  theme_bw(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid = element_blank())

mean_map_ss

ggsave(file.path("figures", "prev_est_map_by_propSS.png"),
       mean_map_ss, width = 12, height = 6, dpi = 300)

## ==========================================================
## B) Vary contact scaling (ϕ): low/med/high
## ==========================================================
param_phi <- tidyr::expand_grid(
  size    = sizes,
  hunt    = hunts,
  phi_cat = names(phi_levels)
) |>
  mutate(
    phi_val   = unname(phi_levels[phi_cat]),
    overrides = purrr::map2(hunt, phi_val, ~list(eta_hunt = .x, phi = .y))
  )

results_phi <- foreach(row = iter(param_phi, by = "row"),
                       .packages = c("tidyverse","deSolve"),
                       .combine = "rbind") %dopar% {
                         
                         df <- run_wl_sim(
                           size               = row$size,
                           infType            = "spillover",
                           years              = years_run,
                           reps               = reps_run,
                           pct                = init_prop_fixed * row$size,
                           prop_superSpreader = 0.05,                 # hold ξ fixed here
                           save_plots         = FALSE,
                           save_runs          = FALSE,
                           overrides          = row$overrides[[1]]    # sets both eta_hunt and phi
                         ) %>%
                           mutate(
                             size       = row$size,
                             pct_hunt   = row$overrides[[1]]$eta_hunt,
                             phi_cat    = row$phi_cat,
                             prev_ratio = `Hunt Prevalence`
                           )
                         
                         prev_sum <- df %>%
                           filter(quarter == 4) %>%
                           group_by(size, pct_hunt, phi_cat) %>%
                           summarise(prev_ratio = mean(prev_ratio, na.rm = TRUE), .groups = "drop")
                         
                         fade_sum <- df %>%
                           group_by(rep, size, pct_hunt, phi_cat) %>%
                           slice_max(tstep, with_ties = FALSE) %>%
                           ungroup() %>%
                           group_by(size, pct_hunt, phi_cat) %>%
                           summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                         
                         rm(df); gc()
                         left_join(prev_sum, fade_sum, by = c("size","pct_hunt","phi_cat"))
                       }

write.csv(results_phi, file.path("data", "hunt_by_phi.csv"), row.names = FALSE)

results_phi <- results_phi %>%
  mutate(phi_lab = forcats::fct_recode(phi_cat,
                                       "ϕ = low"  = "low",
                                       "ϕ = med"  = "med",
                                       "ϕ = high" = "high"))

fade_map_phi <- ggplot(results_phi, aes(x = size, y = pct_hunt, fill = fadeout)) +
  geom_tile() +
  facet_wrap(~ phi_lab, nrow = 1) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       breaks = round(seq(0, max(results_phi$fadeout, na.rm = TRUE), length.out = 4), 2),
                       limits = c(0, max(results_phi$fadeout, na.rm = TRUE)),
                       name = expression(italic(p) * "(fadeout)")) +
  scale_x_continuous(breaks = c(25, 250, 500, 750)) +
  labs(x = "herd size", y = expression("percent harvested ("*eta[h]*")")) +
  theme_bw(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid = element_blank())

fade_map_phi

ggsave(file.path("figures", "fade_prob_map_by_phi.png"),
       fade_map_phi, width = 12, height = 6, dpi = 300)

mean_map_phi <- ggplot(results_phi, aes(x = size, y = pct_hunt, fill = prev_ratio)) +
  geom_tile(color = NA) +
  facet_wrap(~ phi_lab, nrow = 1) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       breaks = round(seq(0, max(results_phi$prev_ratio, na.rm = TRUE), length.out = 5), 2),
                       limits = c(0, round(max(results_phi$prev_ratio, na.rm = TRUE), 2)),
                       name = "mean estimated\nprevalence ratio") +
  scale_x_continuous(breaks = c(25, 250, 500, 750)) +
  labs(x = "herd size", y = expression("percent harvested ("*eta[h]*")")) +
  theme_bw(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid = element_blank())

mean_map_phi

ggsave(file.path("figures", "prev_est_map_by_phi.png"),
       mean_map_phi, width = 12, height = 6, dpi = 300)

stopCluster(cl)










## -----------------------------------------------------------
## Figure: hunting rate (eta_h) × % reduction in p1
## Columns: herd size (three sizes)
## Scenario: spillover (edit 'scenarios' to include "seeded" if desired)
## Metrics: p(fadeout), mean infection prevalence (final tstep)
## -----------------------------------------------------------

rm(list = ls()); gc()
source("bTBwl_func.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(doParallel)
  library(deSolve)   # if run_wl_sim uses it
  library(scales)
})

# ---------- parallel backend ----------
n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# ---------- baselines ----------
p1_baseline <- 0.05   # baseline deer–deer contact (p1); adjust to your baseline

# helper to scale p1 by a % reduction in [0,1]
scale_p1 <- function(reduction) {
  # reduction: 0 => baseline; 1 => zero deer–deer contact
  list(deer_deer_contact_rate = (1 - reduction) * p1_baseline)
}

# ---------- grids ----------
p1_red_grid <- seq(0, 1, length.out = 21)      # y-axis: % reduction in p1
hunt_grid   <- seq(0, 1, length.out = 21)      # x-axis: hunting rate (eta_h)

sizes       <- c(25, 250, 500)                # three herd sizes (facets)
scenarios   <- c("seeded", "spillover")                  # change to c("seeded","spillover") if wanted

param_grid <- tidyr::expand_grid(
  p1_red  = p1_red_grid,
  hunt    = hunt_grid,
  size    = sizes,
  infType = scenarios
) %>%
  mutate(
    overrides = purrr::pmap(
      list(p1_red, hunt),
      function(p1r, hnt) c(
        scale_p1(p1r),
        list(eta_hunt = hnt)
      )
    )
  )

# ---------- run experiments ----------
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse", "deSolve"),
                   .combine  = "rbind") %dopar% {
                     
                     df <- run_wl_sim(
                       size                = row$size,
                       infType             = row$infType,
                       years               = 1,
                       reps                = 300,
                       pct                 = 0.02,
                       prop_superSpreader  = 0,
                       save_plots          = FALSE,
                       overrides           = row$overrides[[1]]
                     ) %>%
                       mutate(
                         size     = row$size,
                         p1_red   = row$p1_red,
                         pct_hunt = row$overrides[[1]]$eta_hunt,
                         scenario = row$infType
                       )
                     
                     # Final-timestep summaries by (hunt, p1_red, size, scenario)
                     final_by_rep <- df %>%
                       group_by(rep, size, p1_red, pct_hunt, scenario) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup()
                     
                     prev_sum <- final_by_rep %>%
                       mutate(prev = `Total Infected` / N) %>%
                       group_by(size, p1_red, pct_hunt, scenario) %>%
                       summarise(mean_prev = mean(prev, na.rm = TRUE), .groups = "drop")
                     
                     fade_sum <- final_by_rep %>%
                       group_by(size, p1_red, pct_hunt, scenario) %>%
                       summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                     
                     rm(df, final_by_rep); gc()
                     left_join(prev_sum, fade_sum,
                               by = c("size","p1_red","pct_hunt","scenario"))
                   }

stopCluster(cl)

# ---------- labels & save ----------
results <- results %>%
  mutate(
    fadeout      = ifelse(fadeout == 0, NA_real_, fadeout),
    size_label   = factor(paste0("N==", size), levels = paste0("N==", sizes)),
    hunt_label   = pct_hunt,                     # keep numeric for x-axis
    p1_red_pct   = 100 * p1_red                  # y-axis in %
  )

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "data/processed/hunt_vs_p1Reduction_bySize_spillover.csv",
          row.names = FALSE)

# ---------- figures ----------
# p(fadeout): x = hunting rate (eta_h), y = % reduction in p1
fade_prob_map <- ggplot(
  data = results,
  aes(x = hunt_label, y = p1_red_pct, fill = fadeout)
) +
  geom_tile() +
  scale_fill_distiller(
    palette   = "Spectral",
    direction = -1,
    na.value  = "gray21",
    breaks    = round(seq(0, max(results$fadeout, na.rm = TRUE), length.out = 4), 2),
    limits    = c(0, max(results$fadeout, na.rm = TRUE)),
    name      = expression(italic(p) * "(fadeout)")
  ) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  labs(
    x = expression("hunting rate  ("*eta[h]*")"),
    y = "% reduction in deer–deer contact (p1)"
  ) +
  facet_grid(row = vars(scenario), cols = vars(size_label), labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

# mean infection prevalence
mean_prev_map <- ggplot(
  data = results,
  aes(x = hunt_label, y = p1_red_pct, fill = mean_prev)
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    breaks = round(seq(0, max(results$mean_prev, na.rm = TRUE), length.out = 5), 2),
    limits = c(0, round(max(results$mean_prev, na.rm = TRUE), 3)),
    name   = "mean infection\nprevalence"
  ) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  labs(
    x = expression("hunting rate  ("*eta[h]*")"),
    y = "% reduction in deer–deer contact (p1)"
  ) +
  facet_grid(row = vars(scenario), cols = vars(size_label), labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

# preview
fade_prob_map
mean_prev_map

# save figs
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/fade_prob_map_hunt_vs_p1Reduction_bySize_spillover.pdf",
       fade_prob_map, width = 12, height = 5, dpi = 300)
ggsave("figures/mean_prev_map_hunt_vs_p1Reduction_bySize_spillover.pdf",
       mean_prev_map, width = 12, height = 5, dpi = 300)


## -----------------------------------------------------------
## Figure: prop_ss (xi) × SS contact factor (phi)
## Facets: columns = herd size (N); rows = scenario (seeded/spillover)
## Metrics: p(fadeout), mean infection prevalence (final tstep)
## -----------------------------------------------------------

rm(list = ls()); gc()
source("bTBwl_func.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(doParallel)
  library(deSolve)
  library(scales)
})

# ---------- parallel backend ----------
n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# ---------- grids ----------
# Proportion super-susceptible (xi)
prop_ss_grid <- seq(0, 1, length.out = 21)     # tweak range/step as you like

# SS contact factor (phi) multiplies farm-contact risk for SS class
phi_grid     <- seq(1, 50, length.out = 21)       # tweak range/step as you like

# Facet sizes
sizes        <- c(50, 250, 500)

# Scenarios (uncomment "seeded" if you want both)
scenarios    <- c("spillover")  # or c("seeded","spillover")

param_grid <- tidyr::expand_grid(
  prop_ss = prop_ss_grid,
  phi     = phi_grid,
  size    = sizes,
  infType = scenarios
) %>%
  mutate(
    overrides = purrr::map(phi, ~list(phi = .x))  # override phi inside the model
  )

# ---------- run experiments ----------
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse", "deSolve"),
                   .combine  = "rbind") %dopar% {
                     
                     df <- run_wl_sim(
                       size                = row$size,
                       infType             = row$infType,
                       years               = 1,
                       reps                = 500,
                       pct                 = 0,
                       prop_superSpreader  = row$prop_ss,             # <- x-axis driver
                       save_plots          = FALSE,
                       overrides           = row$overrides[[1]]       # <- y-axis driver (phi)
                     ) %>%
                       mutate(
                         size     = row$size,
                         prop_ss  = row$prop_ss,
                         phi      = row$overrides[[1]]$phi,
                         scenario = row$infType
                       )
                     
                     # Final-timestep summaries by (prop_ss, phi, size, scenario)
                     final_by_rep <- df %>%
                       group_by(rep, size, prop_ss, phi, scenario) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup()
                     
                     prev_sum <- final_by_rep %>%
                       mutate(prev = `Total Infected` / N) %>%
                       group_by(size, prop_ss, phi, scenario) %>%
                       summarise(mean_prev = mean(prev, na.rm = TRUE), .groups = "drop")
                     
                     fade_sum <- final_by_rep %>%
                       group_by(size, prop_ss, phi, scenario) %>%
                       summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                     
                     rm(df, final_by_rep); gc()
                     left_join(prev_sum, fade_sum, by = c("size","prop_ss","phi","scenario"))
                   }

stopCluster(cl)

# ---------- labels & save ----------
results <- results %>%
  mutate(
    fadeout     = ifelse(fadeout == 0, NA_real_, fadeout),
    size_label  = factor(paste0("N==", size), levels = paste0("N==", sort(unique(size)))),
    # keep numeric axes for continuous scales
    xi_label    = prop_ss,
    phi_label   = phi,
    scenario    = factor(scenario, levels = c("seeded","spillover"))
  )

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "data/processed/propSS_vs_phi_bySize_byScenario.csv", row.names = FALSE)

# ---------- figures ----------
# p(fadeout): x = prop_ss (xi), y = phi
fade_prob_map <- ggplot(
  data = results,
  aes(x = xi_label, y = phi_label, fill = fadeout)
) +
  geom_tile() +
  scale_fill_distiller(
    palette   = "Spectral",
    direction = -1,
    na.value  = "gray21",
    breaks    = round(seq(0, max(results$fadeout, na.rm = TRUE), length.out = 4), 2),
    limits    = c(0, max(results$fadeout, na.rm = TRUE)),
    name      = expression(italic(p) * "(fadeout)")
  ) +
  scale_x_continuous(expand = c(0,0), breaks = pretty(results$xi_label, n = 5)) +
  scale_y_continuous(expand = c(0,0), breaks = pretty(results$phi_label, n = 6)) +
  labs(
    x = expression("proportion super susceptibles  ("*xi*")"),
    y = expression("SS farm-contact factor  ("*phi*")")
  ) +
  facet_grid(row = vars(scenario), cols = vars(size_label), labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

# mean infection prevalence
mean_prev_map <- ggplot(
  data = results,
  aes(x = xi_label, y = phi_label, fill = mean_prev)
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    breaks = round(seq(0, max(results$mean_prev, na.rm = TRUE), length.out = 5), 2),
    limits = c(0, round(max(results$mean_prev, na.rm = TRUE), 3)),
    name   = "mean infection\nprevalence"
  ) +
  scale_x_continuous(expand = c(0,0), breaks = pretty(results$xi_label, n = 5)) +
  scale_y_continuous(expand = c(0,0), breaks = pretty(results$phi_label, n = 6)) +
  labs(
    x = expression("proportion super susceptibles  ("*xi*")"),
    y = expression("SS farm-contact factor  ("*phi*")")
  ) +
  facet_grid(row = vars(scenario), cols = vars(size_label), labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

# preview
fade_prob_map
mean_prev_map

# save figs
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/fade_prob_map_propSS_vs_phi_bySize_byScenario.pdf",
       fade_prob_map, width = 12, height = 5, dpi = 300)
ggsave("figures/mean_prev_map_propSS_vs_phi_bySize_byScenario.pdf",
       mean_prev_map, width = 12, height = 5, dpi = 300)
