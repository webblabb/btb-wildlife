# ================================
rm(list = ls()) #clear workspace for memory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../R/bTBwl_func.R")

run_pipeline <- function(type, years, infType, runtype, pct, prop_superSpreader, reps, sizes, pth, scaled_plots) {
  initialize_environment()
  parameters <- setup_parameters(years, infType, runtype, pct, prop_superSpreader, reps, sizes, pth)
  create_directories()
  cl <- initialize_cluster(parameters$n.cores)
  
  foreach(i = 1:length(parameters$sizes)) %dopar% {
    size <- as.integer(parameters$sizes[i])
    models_res <- run_simulation(size, parameters, type)
    plots <- generate_plots(models_res = models_res, parameters = parameters, 
                            size = size, type = type, scaled = scaled_plots)
    save_results(models_res, parameters, size)
  }
  
  stopCluster(cl)
}

# Run Continuous time simulations
run_pipeline(type = "c", years = 10, infType = "seeded", runtype = "stdPres_",
             pct = 0.02, prop_superSpreader = 0.05, reps = 500, sizes = c(10, 50, 100, 250, 500),
             pth = getwd(), scaled_plots = T)

# Run Pre-simulations for Lambda values
run_pipeline(type = "c", years = 3, infType = "spillover", runtype = "stdLambda_", pct = 0.02,
             prop_superSpreader = 0.0, reps = 200, sizes = c(10, 50, 100, 250, 500, 750, 1000),
             pth = getwd(), scaled_plots = T)

# Run Discrete time simulations
run_pipeline(type = "d", years = 7, infType = "seeded", runtype = "stdD_pub_", pct = 0.02, 
             prop_superSpreader = 0.1, reps = 500, sizes = c(10, 50, 100, 250, 500, 750), pth = getwd(), 
             scaled_plots = T)

# ===============
# Fadeout runs
# ================================
rm(list = ls()) # Clear workspace for memory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")

# Define simulation settings
years <- 20
sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
pct_seq <- seq(from = 0, to = 0.5, length.out = 21)
infType <- "seeded"
runtype <- "fade_"
type_of_integral <- 3
prop_superSpreader <- 0
verbose <- 0
reps <- 500
lambda_factor <- 1.2
pth <- getwd()
save_plots <- FALSE

initialize_environment()
cl <- initialize_cluster()

# Run pipeline for each pct value
results_list <- foreach(pct = pct_seq, .combine = 'cbind', .inorder = TRUE) %dopar% {
  run_pipeline(
    type = "d", 
    years = years, 
    infType = infType, 
    runtype = runtype, 
    pct = pct, 
    prop_superSpreader = prop_superSpreader, 
    reps = reps, 
    sizes = sizes, 
    pth = pth, 
    scaled_plots = save_plots
  )
}

stopCluster(cl)
