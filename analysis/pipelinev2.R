# ================================
# Clear workspace and set working directory
# ================================
rm(list = ls()) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary functions
source(file = "bTBwl_func.R")

# Compile C++ code
system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_DTMC.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_CTMC.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")

# ===============================
# Fadeout runs simulation Settings

param_grid <- expand.grid(
  type = c("c", "d"),
  infType = c("seeded", "spillover"),
  years = c(3, 7, 10),
  pct = 0.02,
  prop_superSpreader = c(0, 0.05, 0.1),
  size = c(10, 50, 100, 250, 500, 750, 1000),
  reps = c(200, 500)
)

# Apply the filtering criteria
filtered_param_grid <- param_grid %>%
  filter(
    (type == "c" & years == 10 & infType == "seeded" & 
       prop_superSpreader == 0.05 & reps == 500 & 
       size %in% c(10, 50, 100, 250, 500)) |
      
      (type == "c" & years == 3 & infType == "spillover" & 
         prop_superSpreader == 0.0 & reps == 200 & 
         size %in% c(10, 50, 100, 250, 500, 750, 1000)) |
      
      (type == "d" & years == 10 & infType == "seeded" & 
         prop_superSpreader == 0.1 & reps == 500 & 
         size %in% c(10, 50, 100, 250, 500, 750))
  )

# Check the result
filtered_param_grid

# ================================
# Run Simulations in Parallel
initialize_environment()
cl <- initialize_cluster(detectCores()-2)

foreach(row = iter(filtered_param_grid, by = "row"), .combine = 'cbind', .inorder = TRUE) %dopar% {
  run_pipeline(
    type = row$type, 
    years = row$years, 
    infType = row$infType, 
    runtype = "pub_", 
    pct = row$pct,
    prop_superSpreader = row$prop_superSpreader, 
    reps = row$reps, 
    size = row$size,
    pth = getwd(),
    scaled_plots = T,
    gen_plots = T,
    save_runs = T
  )
}
stopCluster(cl)

# ================================
# Super susceptible simulations
initialize_environment()
cl <- initialize_cluster(detectCores()-2)

# Generate all (pct, size) combinations
sizes <- round(seq(from = 0, to = 750, length.out = 31)[-1], 0)
prop_ss <- round(seq(from = 0, to = 1, length.out = 21), 2)
param_grid <- expand.grid(prop_ss = prop_ss, size = sizes)

# Run Simulations in Parallel
ss_data <- foreach(row = iter(param_grid, by = "row"), .combine = 'rbind', .inorder = FALSE) %dopar% {
  model_out <- run_pipeline(
    type = "d", 
    years = 1, 
    infType = "spillover", 
    runtype = "super_susc_", 
    pct = 0,
    prop_superSpreader = row$prop_ss, 
    reps = 500, 
    size = row$size,
    pth = getwd(),
    scaled_plots = F,
    gen_plots = F,
    save_runs = T
  )
  
  # Process stochastic output
  data_sto <- model_out[[2]][,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I")] %>%
    pivot_longer(cols = c("N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I"), names_to = "variable", values_to = "value")
  
  # Aggregate stochastic output
  sto_out <- model_out[[2]][,c("tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", "fadeout")]
  sto_out <- aggregate(. ~ tstep, data = sto_out, FUN = mean)
  
  # Compute prevalence
  pct_prev <- sto_out[,"Total Infected"] / sto_out[,'N']
  
  # Store data output
  ss_data <- data.frame(
    tstep = sto_out[,"tstep"],
    Total_Infected = sto_out[,"Total Infected"],
    pct_prev = pct_prev,
    prop_ss = row$prop_ss, 
    size = row$size
  )
  
  rm(sto_out, data_sto, model_out)
  ss_data
}

stopCluster(cl)

# Compute statistics
ss_stats <- ss_data |>
  group_by(prop_ss, size) |>
  summarise(
    min = min(pct_prev, na.rm = TRUE),
    median = median(pct_prev, na.rm = TRUE),
    mean = mean(pct_prev, na.rm = TRUE),
    max = max(pct_prev, na.rm = TRUE),
    .groups = "drop"
  )

# Save data for future use
save(ss_stats, file = "data/ss_stats.RData")
write_csv(ss_stats, "data/ss_stats.csv")

save(ss_data, file = "data/ss_data.RData")
write_csv(ss_data, "data/ss_data.csv")
# ================================
# Fadeout by number simulations
initialize_environment()
cl <- initialize_cluster(detectCores()-2)

# Generate all (pct, size) combinations
sizes <- round(seq(from = 0, to = 750, length.out = 5)[-1], 0)
num <- round(seq(from = 0, to = 25, length.out = 4), 0)
param_grid <- expand.grid(num = num, size = sizes)

# Run Simulations in Parallel
fade_data <- foreach(row = iter(param_grid, by = "row"), .combine = 'rbind', .inorder = TRUE) %dopar% {
  model_out <- run_pipeline(
    type = "d", 
    years = 20, 
    infType = "seeded", 
    runtype = "fade_", 
    pct = row$num,
    prop_superSpreader = 0, 
    reps = 500, 
    size = row$size,
    pth = getwd(),
    scaled_plots = F,
    gen_plots = F,
    save_runs = F
  )
  
  sto_out <- as.data.frame(model_out[[2]][,c("rep", "tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", 'fadeout', 'fadeout time') ])
  data_over_ranges <- aggregate(data=subset(sto_out, select=-c(rep)), .~tstep, FUN = mean) #determine average trajectory
  A <- aggregate(data=subset(sto_out, select=c(rep, `fadeout time`)), .~rep, FUN = max)
  FOT <- mean( A$`fadeout time`[A$`fadeout time`!=0] )
  if(is.nan(FOT)){FOT=1e10}
  remove(sto_out, A)
  
  # save fadeout results to data matrix
  p_fade <- data_over_ranges$fadeout[length(data_over_ranges$fadeout)] # extract mean fade value from last averaged row
  
  results <- data.frame(herd_size = row$size,
                        number_initially_exposed = row$pct,
                        fadeout_probability = p_fade,
                        fadeout_time = FOT)
  
  return(results)
}

stopCluster(cl)

# Save data for future use
save(fade_data, file = "data/fade_data.RData")
write_csv(fade_data, "data/fade_data.csv")




#======
library(data.table)
library(foreach)
library(doParallel)
library(tidyr)
library(readr)

initialize_environment()
cl <- initialize_cluster(detectCores()*3/4)

# Generate all (prop_ss, size) combinations
sizes <- round(seq(from = 0, to = 750, length.out = 31)[-1], 0)
prop_ss <- round(seq(from = 0, to = 1, length.out = 21), 2)
param_grid <- expand.grid(prop_ss = prop_ss, size = sizes)
param_list <- split(param_grid, seq_len(nrow(param_grid)))  # Avoid using iter()

# Run simulations in parallel
ss_data <- foreach(params = param_list, .combine = 'rbind', .inorder = FALSE) %dopar% {

  model_out <- run_simulation(type = "d", 
                              years = 1, 
                              lambda_years = 3, 
                              infType = "spillover", 
                              runtype = "super_susc_", 
                              pct = 0, 
                              prop_superSpreader = params$prop_ss, 
                              reps = 500, 
                              size = params$size, 
                              pth = getwd(), 
                              data_out = "wide", 
                              solve_ode = FALSE,
                              save_runs = TRUE)
  
  # Aggregate stochastic output
  sto_out <- model_out[[1]][,c("tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", "fadeout")]
  sto_out <- aggregate(. ~ tstep, data = sto_out, FUN = mean)
  
  # Compute prevalence
  pct_prev <- sto_out[,"Total Infected"] / sto_out[,'N']
  
  # Store data output
  ss_data <- data.frame(
    tstep = sto_out[,"tstep"],
    Total_Infected = sto_out[,"Total Infected"],
    pct_prev = pct_prev,
    prop_ss = params$prop_ss, 
    size = params$size
  )
  
  rm(sto_out, model_out)
  gc()
  
  ss_data
}

stopCluster(cl)

# Compute statistics
ss_stats <- ss_data |>
  group_by(prop_ss, size) |>
  summarise(
    min = min(pct_prev, na.rm = TRUE),
    median = median(pct_prev, na.rm = TRUE),
    mean = mean(pct_prev, na.rm = TRUE),
    max = max(pct_prev, na.rm = TRUE),
    .groups = "drop"
  )






initialize_environment()
cl <- initialize_cluster(detectCores() * 3 / 4)

# Generate all (prop_ss, size) combinations
sizes <- round(seq(from = 0, to = 750, length.out = 31)[-1], 0)
prop_ss <- round(seq(from = 0, to = 1, length.out = 21), 2)

# # Generate all (prop_ss, size) combinations
# sizes <- round(seq(from = 0, to = 750, by = 25)[-1], 0)
# prop_ss <- round(seq(from = 0, to = 1, by = 0.05), 2)
# 
# sizes <- sizes[sizes>375]

batch_sizes <- split(sizes, ceiling(seq_along(sizes)/5))  # Split into smaller groups
for (batch in batch_sizes) {
  # Run simulations in parallel
  ss_data <- foreach(size = batch, .combine = 'rbind', .inorder = FALSE) %dopar% {
    
    # Nested foreach loop for each combination of prop_ss and size
    data <- foreach(prop = prop_ss, .combine = 'rbind', .inorder = FALSE) %dopar% {
      
      model_out <- run_simulation(
        type = "d", 
        years = 1, 
        lambda_years = 3, 
        infType = "spillover", 
        runtype = "super_susc_", 
        pct = 0, 
        prop_superSpreader = prop, 
        reps = 500, 
        size = size, 
        pth = getwd(), 
        data_out = "wide", 
        solve_ode = FALSE,
        save_runs = FALSE
      )
      
      # Aggregate stochastic output
      sto_out <- model_out[[1]][, c("tstep", "N", "S", "SS_S", "E1", "SS_E1", "I", "SS_I", "Total Infected", "fadeout")]
      sto_out <- aggregate(. ~ tstep, data = sto_out, FUN = mean)
      
      # Compute prevalence
      pct_prev <- sto_out[, "Total Infected"] / sto_out[, 'N']
      
      # Store data output
      ss_data <- data.frame(
        tstep = sto_out[, "tstep"],
        Total_Infected = sto_out[, "Total Infected"],
        pct_prev = pct_prev,
        prop_ss = prop, 
        size = size
      )
      
      # save(ss_data, file = paste0('data/', "ss_spillover_", prop, '-', size, '-', Sys.Date(), '.csv'))
      
      rm(sto_out, model_out)
      gc()
      
      return(ss_data)
    }
    
    save(data, file = paste0('data/', "ss_spillover_", size, '-', Sys.Date(), '.csv'))
    
    return(data)
  }
  saveRDS(ss_data, file = paste0('data/ss_spillover_batch_', batch[1], '.rds'))
  rm(ss_data)
  gc()
}
stopCluster(cl)


