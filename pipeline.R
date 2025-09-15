####################
## Super Susceptible ##
####################
# Run time: variable for grid size
print("Super susceptible epidemic extent simulations by proportion and herd size")
rm(list = ls())
source(file = "R/bTBwl_func.R")

years <- 1
times <- seq(from = 0, to = years * 12, by = 12 / 365) # daily time steps
seedQuarter <- 1
pct <- 0
verbose <- 0
lambda_factor <- 1.2
reps <- 500
sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
range_vec <- seq(from = 0, to = 1, length.out = 21) # prop SS - epidemic extent
infType <- "spillover"
runtype <- "SS_pub_"
type_of_integral <- 3

# Create all combinations of sizes and SS proportion as a tibble
param_grid <- expand_grid(
  size = as.integer(sizes),
  SS_prop = range_vec
)

# Initialize result tibbles
ss_results <- tibble(
  size = integer(),
  SS_prop = double(),
  min_prev = double(),
  median_prev = double(),
  mean_prev = double(),
  max_prev = double()
)
ss_epi_extent <- tibble(
  size = integer(),
  SS_prop = double(),
  tstep = double(),
  pct_inf = double(),
  num_inf = double()
)

for (row in 1:nrow(param_grid)) {
  size <- param_grid$size[row]
  SS_p <- param_grid$SS_prop[row]
  
  pars <- parameter_set_wl(
    k = size,
    scenario = infType,
    initial_exposed = pct * size,
    SS_prop = SS_p,
    start_quarter = seedQuarter,
    verbose = verbose
  )
  initial_state <- tibble(
    S_0 = pars["S_0"], E1_0 = pars["E1_0"], I_0 = pars["I_0"],
    SuperS_0 = pars["SuperS_0"], SuperE1_0 = pars["SuperE1_0"], SuperI_0 = pars["SuperI_0"]
  )
  population_parameters <- tibble(
    K = pars["K"], eta_hunt = pars["eta_hunt"], eta_nat = pars["eta_nat"], theta = pars["theta"],
    gamma = pars["gamma"], alpha_max = pars["alpha_max"], ksi = pars["ksi"], omega = pars["omega"],
    s = pars["s"], alpha = pars["alpha"]
  )
  disease_parameters <- tibble(
    beta = pars["beta"], area = pars["area"], p1 = pars["p1"],
    p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], p2_q4 = pars["p2_q4"],
    phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], sigma1_rate = pars["sigma1_rate"]
  )
  param_vals <- left_join(population_parameters, disease_parameters, by = character())
  
  # Pre-simulation for lambda
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = as.data.frame(param_vals),
    initial_state = as.data.frame(initial_state),
    nyears = 3,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = -3,
    batch_name = "test",
    type = 'c'
  )
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  rm(sto_out)
  
  # Main simulation
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = as.data.frame(param_vals),
    initial_state = as.data.frame(initial_state),
    nyears = years,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = verbose,
    batch_name = "test",
    type = 'd',
    lambda = lambda_out * lambda_factor,
    integrate_type = type_of_integral
  )
  
  # Assign epidemic extent metrics
  sto_out_df <- as.data.frame(sto_out[, c("tstep", "N", "Total Infected")])
  sto_out_df <- sto_out_df %>%
    group_by(tstep) %>%
    summarise(
      N = mean(N),
      num_inf = mean(`Total Infected`)
    ) %>%
    mutate(pct_inf = num_inf / N)
  
  ss_epi_extent <- ss_epi_extent %>% bind_rows(
    sto_out_df %>% select(tstep, pct_inf, num_inf) %>%
      mutate(size = size, SS_prop = SS_p)
  )
  
  # Summary statistics: min, median, mean, max of prevalence
  pct_prev <- sto_out_df$pct_inf
  stats <- c(
    min_prev = min(pct_prev, na.rm = TRUE),
    median_prev = median(pct_prev, na.rm = TRUE),
    mean_prev = mean(pct_prev, na.rm = TRUE),
    max_prev = max(pct_prev, na.rm = TRUE)
  )
  ss_results <- ss_results %>% add_row(
    size = size,
    SS_prop = SS_p,
    min_prev = stats["min_prev"],
    median_prev = stats["median_prev"],
    mean_prev = stats["mean_prev"],
    max_prev = stats["max_prev"]
  )
  
  rm(sto_out, sto_out_df, stats)
  gc()
}

# POST-PROCESSING: tidy data.frames for plotting
ss_stats_data <- ss_results %>% 
  rename(
    `herd size` = size,
    `proportion super susceptible` = SS_prop,
    `min prevalence` = min_prev,
    `median prevalence` = median_prev,
    `mean prevalence` = mean_prev,
    `max prevalence` = max_prev
  )

ss_extent_data <- ss_epi_extent %>% 
  rename(
    `herd size` = size,
    `proportion super susceptible` = SS_prop,
    tstep = tstep,
    `percent infected` = pct_inf,
    `number infected` = num_inf
  )

# Save to disk
ss_stats_name <- "data/ss_plotStat.RData"
ss_extent_name <- "data/ss_plotData.RData"
if(!file.exists(ss_stats_name)) save(ss_stats_data, file = ss_stats_name)
if(!file.exists(ss_extent_name)) save(ss_extent_data, file = ss_extent_name)

####################
##    Fadeout     ##
####################
# Run time: ~8 min for 20x21 grid
print("P(fadeout) by number initially exposed and herd size simulations")
rm(list = ls())
source(file = "R/bTBwl_func.R")

years <- 20
times <- seq(from = 0, to = years * 12, by = 12 / 365)
seedQuarter <- 1
prop_superSpreader <- 0
verbose <- 0
reps <- 500
lambda_factor <- 1.2
sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
pct <- seq(from = 0, to = .5, length.out = 21)
infType <- "seeded"
runtype <- "fade_"
type_of_integral <- 3

# Create all combinations of sizes and pct as a tibble
param_grid <- expand_grid(
  size = as.integer(sizes),
  pct = pct
)

# Initialize result tibbles
fade_results <- tibble(
  size = integer(),
  pct = double(),
  fadeout_prob = double(),
  fadeout_time = double()
)

for (row in 1:nrow(param_grid)) {
  size <- param_grid$size[row]
  p <- param_grid$pct[row]
  pars <- parameter_set_wl(
    k = size,
    scenario = infType,
    initial_exposed = p * size,
    SS_prop = 0,
    start_quarter = seedQuarter,
    verbose = verbose
  )
  initial_state <- tibble(
    S_0 = pars["S_0"], E1_0 = pars["E1_0"], I_0 = pars["I_0"],
    SuperS_0 = pars["SuperS_0"], SuperE1_0 = pars["SuperE1_0"], SuperI_0 = pars["SuperI_0"]
  )
  population_parameters <- tibble(
    K = pars["K"], eta_hunt = pars["eta_hunt"], eta_nat = pars["eta_nat"], theta = pars["theta"],
    gamma = pars["gamma"], alpha_max = pars["alpha_max"], ksi = pars["ksi"], omega = pars["omega"],
    s = pars["s"], alpha = pars["alpha"]
  )
  disease_parameters <- tibble(
    beta = pars["beta"], area = pars["area"], p1 = pars["p1"],
    p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], p2_q4 = pars["p2_q4"],
    phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], sigma1_rate = pars["sigma1_rate"]
  )
  param_vals <- left_join(population_parameters, disease_parameters, by = character())
  
  # Pre-simulation for lambda
  sto_out <- wildlife_model(
    n_reps = 200,
    parameters = as.data.frame(param_vals),
    initial_state = as.data.frame(initial_state),
    nyears = 5,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = -3,
    batch_name = "test",
    type = 'c'
  )
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  rm(sto_out)
  
  # Main simulation
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = as.data.frame(param_vals),
    initial_state = as.data.frame(initial_state),
    nyears = years,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = verbose,
    batch_name = "test",
    type = 'd',
    lambda = lambda_out * lambda_factor,
    integrate_type = type_of_integral
  )
  
  # Fadeout probability: mean of final 'fadeout' for all reps
  fadeout_final <- sto_out %>% 
    filter(tstep == max(tstep)) %>% 
    summarise(fadeout_prob = mean(fadeout, na.rm = TRUE)) %>% 
    pull(fadeout_prob)
  
  # Fadeout time: mean of max fadeout time per rep, for reps that faded out
  fade_times <- sto_out %>% 
    group_by(rep) %>% 
    summarise(ft = max(`fadeout time`, na.rm = TRUE)) %>% 
    filter(ft != 0)
  FOT <- ifelse(nrow(fade_times) > 0, mean(fade_times$ft, na.rm = TRUE), 1e10)
  
  fade_results <- fade_results %>% 
    add_row(size = size, pct = p, fadeout_prob = fadeout_final, fadeout_time = FOT)
  
  rm(sto_out, fade_times)
  gc()
}

# POST-PROCESSING: tidy data.frames for plotting
fade_data_m <- fade_results %>% 
  rename(
    `herd size` = size,
    `proportion initially exposed` = pct,
    `fadeout probability` = fadeout_prob
  )
fadeTime_data_m <- fade_results %>% 
  rename(
    `herd size` = size,
    `proportion initially exposed` = pct,
    `fadeout time` = fadeout_time
  )
fadeTime_data_m <- fadeTime_data_m %>% mutate(`fadeout time` = if_else(`fadeout time` > 1e6, NA_real_, `fadeout time`))
fade_data_m <- fade_data_m %>% mutate(`fadeout probability` = na_if(`fadeout probability`, 0))

# Save to disk
fade_data_name <- "data/fade_PlotDat.RData"
fadeTime_data_name <- "data/fadeTime_PlotDat.RData"
if(!file.exists(fade_data_name)) save(fade_data_m, file = fade_data_name)
if(!file.exists(fadeTime_data_name)) save(fadeTime_data_m, file = fadeTime_data_name)

####################
## Hunting  Rates ##
####################
# Run time: variable for grid size
print("Hunting effect simulations by harvest rate and herd size")
rm(list = ls())
source(file = "R/bTBwl_func.R")

years <- 20
times <- seq(from = 0, to = years * 12, by = 12 / 365) # daily time steps
seedQuarter <- 1
prop_superSpreader <- 0.00
pct <- .025 # explore low medium and high prevalence c(.025, .1, .25)
case <- "low"
verbose <- 0
reps <- 300
sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
range_vec <- seq(from = 0, to = .5, length.out = 21) # harvest rate
infType <- "seeded"
runtype <- "hunt_lowPrev_"
type_of_integral <- 3
lambda_factor <- 1.2

# Create all combinations of sizes and harvest rates as a tibble
param_grid <- expand_grid(
  size = as.integer(sizes),
  harvest_rate = range_vec
)

# Initialize result tibbles
hunt_results <- tibble(
  size = integer(),
  harvest_rate = double(),
  mean_prev = double(),
  median_prev = double()
)
hunt_epi_extent <- tibble(
  size = integer(),
  harvest_rate = double(),
  tstep = double(),
  mean_prev = double(),
  median_prev = double()
)
fadeout_results <- tibble(
  size = integer(),
  harvest_rate = double(),
  fadeout_prob = double()
)

for (row in 1:nrow(param_grid)) {
  size <- param_grid$size[row]
  hr <- param_grid$harvest_rate[row]
  
  pars <- parameter_set_wl(
    k = size,
    scenario = infType,
    initial_exposed = pct * size,
    SS_prop = 0,
    start_quarter = seedQuarter,
    verbose = verbose
  )
  initial_state <- tibble(
    S_0 = pars["S_0"], E1_0 = pars["E1_0"], I_0 = pars["I_0"],
    SuperS_0 = pars["SuperS_0"], SuperE1_0 = pars["SuperE1_0"], SuperI_0 = pars["SuperI_0"]
  )
  population_parameters <- tibble(
    K = pars["K"], eta_hunt = hr, eta_nat = pars["eta_nat"], theta = pars["theta"],
    gamma = pars["gamma"], alpha_max = pars["alpha_max"], ksi = pars["ksi"], omega = pars["omega"],
    s = pars["s"], alpha = pars["alpha"]
  )
  disease_parameters <- tibble(
    beta = pars["beta"], area = pars["area"], p1 = pars["p1"],
    p2_q1 = pars["p2_q1"], p2_q2 = pars["p2_q2"], p2_q3 = pars["p2_q3"], p2_q4 = pars["p2_q4"],
    phi = pars["phi"], sigma1_mean = pars["sigma1_mean"], sigma1_rate = pars["sigma1_rate"]
  )
  param_vals <- left_join(population_parameters, disease_parameters, by = character())
  
  # Pre-simulation for lambda
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = as.data.frame(param_vals),
    initial_state = as.data.frame(initial_state),
    nyears = 3,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = -3,
    batch_name = "test",
    type = 'c'
  )
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  rm(sto_out)
  
  # Main simulation
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = as.data.frame(param_vals),
    initial_state = as.data.frame(initial_state),
    nyears = years,
    seed_quarter = as.integer(pars["start_q"]),
    verbose = verbose,
    batch_name = "test",
    type = 'd',
    lambda = lambda_out * lambda_factor,
    integrate_type = type_of_integral
  )
  
  sto_out_df <- as.data.frame(sto_out)
  # Only consider quarter == 4 (hunting)
  hunt_quarter <- sto_out_df %>% filter(quarter == 4)
  
  # Prevalence estimates
  prev_mean <- mean(hunt_quarter$`Hunt Prevalence`, na.rm = TRUE)
  prev_median <- median(hunt_quarter$`Hunt Prevalence`, na.rm = TRUE)
  hunt_results <- hunt_results %>% add_row(
    size = size,
    harvest_rate = hr,
    mean_prev = prev_mean,
    median_prev = prev_median
  )
  
  # Aggregate metrics by tstep
  mean_Estim <- hunt_quarter %>%
    group_by(tstep) %>%
    summarise(mean_prev = mean(`Hunt Prevalence`, na.rm = TRUE))
  median_Estim <- hunt_quarter %>%
    group_by(tstep) %>%
    summarise(median_prev = median(`Hunt Prevalence`, na.rm = TRUE))
  
  # Merge mean and median by tstep
  epi_extent <- left_join(mean_Estim, median_Estim, by = "tstep") %>%
    mutate(size = size, harvest_rate = hr)
  hunt_epi_extent <- bind_rows(hunt_epi_extent, epi_extent)
  
  # Fadeout probability: mean of final 'fadeout' for all reps
  sto_out_agg <- sto_out_df %>%
    group_by(tstep) %>%
    summarise(fadeout = mean(fadeout, na.rm = TRUE))
  fade_final <- sto_out_agg$fadeout[nrow(sto_out_agg)]
  fadeout_results <- fadeout_results %>% add_row(
    size = size,
    harvest_rate = hr,
    fadeout_prob = fade_final
  )
  
  rm(sto_out, sto_out_df, hunt_quarter, mean_Estim, median_Estim, epi_extent, sto_out_agg)
  gc()
}

# POST-PROCESSING: tidy data.frames for plotting
hunt_stats_data <- hunt_results %>% 
  rename(
    `herd size` = size,
    `percent harvested` = harvest_rate,
    `mean prevalence` = mean_prev,
    `median prevalence` = median_prev
  )

hunt_extent_data <- hunt_epi_extent %>% 
  rename(
    `herd size` = size,
    `percent harvested` = harvest_rate,
    tstep = tstep,
    `mean prevalence` = mean_prev,
    `median prevalence` = median_prev
  )

fadeout_matrix <- fadeout_results %>% 
  pivot_wider(names_from = harvest_rate, values_from = fadeout_prob) %>%
  select(-size) %>%
  as.matrix()
rownames(fadeout_matrix) <- sizes
colnames(fadeout_matrix) <- range_vec

hunt_data_name <- paste0('data/hunt_plotData_', case, '.RData')
hunt_stat_name <- paste0('data/hunt_plotStat_', case, '.RData')
hunt_fade_name <- paste0('data/hunt_plotFade_', case, '.RData')

# Save to disk
if(!file.exists(hunt_data_name)) save(hunt_extent_data, file = hunt_data_name)
if(!file.exists(hunt_stat_name)) save(hunt_stats_data, file = hunt_stat_name)
if(!file.exists(hunt_fade_name)) save(fadeout_matrix, file = hunt_fade_name)

# ------------------------------------------------------------------------------
# Sensitivity analysis

######################
## Setup Parameters ##
######################
rm(list = ls())
library(tidyverse)
source(file = "R/bTBwl_func.R")
N_LHS_sets <- 1000
size_range <- c(50, 1000)
years <- 4
prop_superSpreader <- 0.1
pct <- .02
verbose <- 0
reps <- 500
def_seed <- 43
infection <- "seeded"
infType <- "seeded_q1"
fix_q <- 1
runtype <- "LHS_"
type_of_integral <- 3
par_file <- paste0('bTBwl_LHSpars_', infType)
lambda_factor <- 1.2
save_pars <- TRUE
save_runs <- TRUE
recalc <- FALSE

# Directory setup for serialization
raw_dir <- "data/raw/"
summary_dir <- "data/"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

###############################
## Generate or load LHS pars ##
###############################
  lh <- LHS_ParSets(
    k_lim = size_range,
    scenario = infection,
    infected = pct,
    pct = TRUE,
    SS_prop = prop_superSpreader,
    verbose = 0,
    Num_LHS_sets = N_LHS_sets,
    file.path = raw_dir,
    file.out = save_pars,
    file.name = par_file,
    seed_q = fix_q,
    seed = def_seed
  )
if ("lh" %in% ls()) {
  # Found in workspace
} else if (paste0('bTBwl_LHSpars_', infType, '.csv') %in% list.files(raw_dir)) {
  lh <- read_csv(paste0(raw_dir, 'bTBwl_LHSpars_', infType, '.csv'))
} else {
  lh <- LHS_ParSets(
    k_lim = size_range,
    scenario = infection,
    infected = pct,
    pct = TRUE,
    SS_prop = prop_superSpreader,
    verbose = 0,
    Num_LHS_sets = N_LHS_sets,
    file.path = raw_dir,
    file.out = save_pars,
    file.name = par_file,
    seed_q = fix_q,
    seed = def_seed
  )
}
###############################

#########################################
## Run model and save raw outputs      ##
#########################################
if (def_seed == 43 && N_LHS_sets == 1000 && nrow(lh) == N_LHS_sets) {
  lh <- lh[-c(363, 367), ] # Remove problematic entries
}

TIC <- Sys.time()

# SERIAL version (no cluster setup, for reproducibility)
sens_results <- vector("list", nrow(lh))
for (i in seq_len(nrow(lh))) {
# for (i in 321:max(nrow(lh))) {
  pars <- lh[i, ]
  
  initial_state <- data.frame(
    S_0 = pars[["S_0"]], E1_0 = pars[["E1_0"]], I_0 = pars[["I_0"]],
    SuperS_0 = pars[["SuperS_0"]], SuperE1_0 = pars[["SuperE1_0"]], SuperI_0 = pars[["SuperI_0"]]
  )
  population_parameters <- data.frame(
    K = pars[["K"]], eta_hunt = pars[["eta_hunt"]], eta_nat = pars[["eta_nat"]], theta = pars[["theta"]],
    gamma = pars[["gamma"]], alpha_max = pars[["alpha_max"]], ksi = pars[["ksi"]], omega = pars[["omega"]],
    s = pars[["s"]], alpha = 0
  )
  disease_parameters <- data.frame(
    beta = pars[["beta"]], area = 1, p1 = 1,
    p2_q1 = pars[["p2_q1"]], p2_q2 = pars[["p2_q2"]], p2_q3 = pars[["p2_q3"]], p2_q4 = pars[["p2_q4"]],
    phi = pars[["phi"]], sigma1_mean = pars[["sigma1_mean"]], sigma1_rate = pars[["sigma1_rate"]]
  )
  param_vals <- merge(population_parameters, disease_parameters)
  
  # Stochastic lambda estimation
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = param_vals,
    initial_state = initial_state,
    nyears = years,
    seed_quarter = as.integer(pars[["start_q"]]),
    verbose = -3,
    batch_name = "test",
    type = 'c'
  )
  lambda_out <- getLambda_vec(data = sto_out, type = 'max')
  
  # Main simulation
  sto_out <- wildlife_model(
    n_reps = reps,
    parameters = param_vals,
    initial_state = initial_state,
    nyears = years,
    seed_quarter = as.integer(pars[["start_q"]]),
    verbose = -19,
    batch_name = "test",
    type = 'd',
    lambda = lambda_out * lambda_factor,
    integrate_type = type_of_integral
  )
  
  # if (save_runs) {
  #   write_csv(as_tibble(sto_out), paste0(raw_dir, "sens_run_", infType, "_", i, ".csv"))
  # }
  
  # Mean hunt prevalence
  prev_est <- sto_out %>%
    filter(quarter == 4) %>%
    group_by(rep) %>%
    summarise(HuntPrev = mean(`Hunt Prevalence`, na.rm = TRUE)) %>%
    pull(HuntPrev)
  
  # Only last time step
  sto_final <- sto_out %>%
    filter(tstep == max(tstep)) %>%
    mutate(`Hunt Prevalence` = prev_est) %>%
    select(N, `Total Infected`, fadeout, `fadeout time`, `Hunt Prevalence`)
  
  # Merge with parameter set
  sto_final <- bind_cols(sto_final, pars[rep(1, reps), ])
  sens_results[[i]] <- sto_final
  gc()
}
TOC <- Sys.time()
print(TOC - TIC)

# Final sensitivity analysis dataframe
R <- bind_rows(sens_results)
# Scale parameter columns for sensitivity summary
scale_cols <- c('K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 'alpha_max', 'ksi',
                'omega', 's', 'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi',
                'sigma1_mean', 'sigma1_rate', 'start_q')
R_scaled <- R %>%
  mutate(across(all_of(scale_cols), ~ as.vector(scale(.))))

# Save final sensitivity analysis dataframes
write_csv(R, paste0(summary_dir, "LHS_summary_", infType, ".csv"))
write_csv(R_scaled, paste0(summary_dir, "LHS_scaled_summary_", infType, ".csv"))
