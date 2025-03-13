# ================================
rm(list = ls()) #clear workspace for memory
setwd(dirname(getActiveDocumentContext()$path))
source(file = "bTBwl_func.R")

run_pipeline <- function(type, years, infType, runtype, pct, prop_superSpreader, reps, sizes, pth, scaled_plots) {
  initialize_environment()
  parameters <- setup_parameters(years, infType, runtype, pct, prop_superSpreader, reps, sizes, pth)
  create_directories()
  cl <- initialize_cluster(parameters$n.cores)
  
  foreach(i = 1:length(parameters$sizes)) %dopar% {
    size <- as.integer(parameters$sizes[i])
    models_res <- run_simulation(size, parameters, type)
    plots <- generate_plots(model_res = models_res, parameters = parameters, 
                            size = size, type = type, scaled = scaled_plots)
    save_results(models_res, parameters, size)
  }
  
  stopCluster(cl)
}

# # Run Continuous time simulations
# run_pipeline(type = "c", years = 10, infType = "seeded", runtype = "stdPres_", 
#              pct = 0.02, prop_superSpreader = 0.05, reps = 500, sizes = c(10, 50, 100, 250, 500), 
#              pth = getwd(), scaled_plots = T)
# 
# # Run Pre-simulations for Lambda values
# run_pipeline(type = "c", years = 3, infType = "spillover", runtype = "stdLambda_", pct = 0.02, 
#              prop_superSpreader = 0.0, reps = 200, sizes = c(10, 50, 100, 250, 500, 750, 1000), 
#              pth = getwd(), scaled_plots = T)

# Run Discrete time simulations
run_pipeline(type = "d", years = 7, infType = "seeded", runtype = "stdD_pub_", pct = 0.02, 
             prop_superSpreader = 0.1, reps = 500, sizes = 100, pth = getwd(), scaled_plots = T)
