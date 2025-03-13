# ================================
# Clear workspace and set working directory
# ================================
rm(list = ls()) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary functions
source(file = "bTBwl_func.R")

# ===============================
# Fadeout runs simulation Settings
initialize_environment()
cl <- initialize_cluster(detectCores()-2)

param_grid <- expand.grid(type = c("c","d"),
                          infType = c("seeded","spillover"),
                          years = c(3, 10, 7),
                          pct = 0.02, 
                          prop_superSpreader = c(0, 0.05,0.1),
                          size = c(10, 50, 100, 250, 500, 750, 1000)) |>
  filter(!(type == "c" & size > 500)) |> # subset the runs a bit so we don't have ones that last forever
  filter(!(years %in% c(7,10) & infType == "spillover")) |>
  mutate(runtype = "pub")
  

# ================================
# Run Simulations in Parallel

foreach(row = iter(param_grid, by = "row"), .combine = 'cbind', .inorder = TRUE) %dopar% {
  run_pipeline(
    type = row$type, 
    years = row$years, 
    infType = row$infType, 
    runtype = row$runtype, 
    pct = row$pct,
    prop_superSpreader = row$prop_superSpreader, 
    reps = 500, 
    size = row$size,
    pth = getwd(),
    scaled_plots = T,
    gen_plots = T,
    save_runs = T
  )
}
stopCluster(cl)

# ================================
# Fadeout runs simulation Settings
initialize_environment()
cl <- initialize_cluster(detectCores()-2)

# Generate all (pct, size) combinations
sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
pct_seq <- seq(from = 0, to = 0.5, length.out = 21)
param_grid <- expand.grid(pct = pct_seq, size = sizes)

# ================================
# Run Simulations in Parallel

results_list <- foreach(row = iter(param_grid, by = "row"), .combine = 'cbind', .inorder = TRUE) %dopar% {
  run_pipeline(
    type = "d", 
    years = 20, 
    infType = "seeded", 
    runtype = "fade_", 
    pct = row$pct,
    prop_superSpreader = 0, 
    reps = 500, 
    size = row$size,
    pth = getwd(),
    scaled_plots = F,
    gen_plots = F,
    save_runs = T
  )
}
stopCluster(cl)
