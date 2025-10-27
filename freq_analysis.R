###############################################################
## Compile and run code ##
###############################################################
system("g++ -L/usr/lib/x86_64-linux-gnu tb_wildlife_freq_cont.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu tb_wildlife_freq_disc.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")

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

dir.create("data",      showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw",  showWarnings = FALSE, recursive = TRUE)
dir.create("figures/supp",   showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Main text figures

## -----------------------------------------------------------
## Figure 1: size x infType

param_grid <- expand_grid(
  size = c(250),
  infType = c("seeded", "spillover")
)

for (i in 1:nrow(param_grid)) {
  row <- param_grid[i,]
  
  run_wl_sim(size = row$size, 
             infType = row$infType, 
             save_plots = TRUE,
             freq = T,
             supp = T) 
}

## -----------------------------------------------------------
## Figure 2: super-spreader prop x herd size
## -----------------------------------------------------------

## --- setup parallel backend ---
n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

sizes  <- seq(from = 0, to = 750, length.out = 21)[-1]
xi_vec <- seq(from = 0, to = 1, length.out = 11)

param_grid <- expand_grid(
  size    = sizes,
  infType = "spillover",
  SS_prop = xi_vec
)

## --- parallel loop over rows ---
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse", "deSolve"), .combine = "rbind" ) %dopar% {
                     
                     gc()
                     
                     run_wl_sim(
                       size               = row$size,
                       infType            = row$infType,
                       years              = 1,
                       reps               = 500,
                       pct                = 0,         # no initial infection
                       prop_superSpreader = row$SS_prop,
                       save_plots         = FALSE
                     ) |> 
                       filter(tstep == max(tstep)) |>
                       mutate(prev = `Total Infected`/N, size = row$size, prop_ss = row$SS_prop) |>
                       group_by(size, prop_ss) |>
                       summarise(mean_prev = mean(prev))
                     
                   }
stopCluster(cl)

mean_map <- ggplot(
  data = results,
  aes(
    x = size,
    y = prop_ss,
    fill = mean_prev
  )
) +
  geom_tile(color = NA) +  # remove borders between tiles
  scale_fill_distiller(palette = "Spectral",
                       breaks = round(seq(0,max(results$mean_prev), length.out = 4),3),   # manually specify desired breaks
                       limits = c(0, max(results$mean_prev))) +
  scale_x_continuous(breaks = c(25, 250, 500, 750, 1000)) +
  labs(
    fill = str_wrap("mean infection prevalence", 0),
    x = "herd size",
    y = expression("proportion of super susceptibles ("*xi*")")
  ) +
  theme_bw(base_size = 14) +   # clean white theme
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8), # black box around plot
    panel.grid = element_blank(),       # remove grid lines
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, 'cm'),
    legend.text = element_text(size = 11)
  ) +
  coord_cartesian(expand = FALSE)  # ensures tiles touch axes

mean_map

ggsave("mean_prevalence_map_freq.pdf", mean_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## Figure 3: fadeout. herd size x number initially exposed 
## -----------------------------------------------------------

n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

sizes   <- seq(from = 0, to = 750, length.out = 11)[-1]
initNum <- 0:10

param_grid <- tidyr::expand_grid(
  size    = sizes,
  infType = "seeded",
  numExp  = initNum
)

## --- parallel loop over rows ---
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse", "deSolve"), .combine = "rbind" ) %dopar% {
                     
                     run_wl_sim(
                       size       = row$size,
                       infType    = row$infType,
                       years      = 20,
                       reps       = 500,
                       save_plots = FALSE,
                       pct        = row$numExp
                     ) |> 
                       mutate(size = row$size, num_exp = row$numExp) |>
                       group_by(size, num_exp) |>
                       summarise(fadeout = mean(fadeout))
                   }
stopCluster(cl)

write.csv(results, "data/fadeout.csv")

fade_prob_map <- ggplot(
  data = results,
  aes(
    x = size,
    y = num_exp,
    fill = fadeout
  )
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,                # reverse if you want red = high, blue = low
    na.value = "gray21",
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    name = expression(italic(p) * "(fadeout)")
  ) +
  scale_x_continuous(breaks = c(25, 250, 500, 750, 1000)) +
  labs(
    x = "herd size",
    y = "number initially exposed"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, 'cm'),
    legend.text = element_text(size = 11),
  ) +
  coord_cartesian(expand = FALSE)

fade_prob_map

ggsave("figures/fade_prob_map_num_freq.pdf", fade_prob_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## Figure 5 (a and b): herd size x harvest intensity
## -----------------------------------------------------------

n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

sizes <- seq(from = 0, to = 750, length.out = 11)[-1]
hunt_pct <- seq(from = 0, to = .75, length.out = 11)

param_grid <- expand_grid(
  size    = sizes,
  infType = "seeded",
  hunt = hunt_pct
) |>
  mutate(overrides = purrr::map(hunt, ~list(eta_hunt = .x)))

## --- parallel loop over rows ---
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse", "deSolve"), .combine = "rbind" ) %dopar% {
                     
                     df <- run_wl_sim(
                       size       = row$size,
                       infType    = row$infType,
                       years      = 20,
                       reps       = 300,
                       pct        = 0.025*row$size,
                       prop_superSpreader = 0,
                       save_plots = FALSE,
                       overrides  = row$overrides[[1]]
                     ) %>%
                       mutate(
                         size      = row$size,
                         pct_hunt  = row$overrides[[1]]$eta_hunt,
                         prev_ratio = `Hunt Prevalence`
                       )
                     
                     prev_sum <- df %>%
                       filter(quarter == 4) %>%
                       group_by(size, pct_hunt) %>%
                       summarise(prev_ratio = mean(prev_ratio, na.rm = TRUE), .groups = "drop")
                     
                     # fadeout only at FINAL tstep per replicate
                     fade_sum <- df %>%
                       group_by(rep, size, pct_hunt) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup() %>%
                       group_by(size, pct_hunt) %>%
                       summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                     
                     rm(df)
                     gc()
                     
                     left_join(prev_sum, fade_sum, by = c("size", "pct_hunt"))
                   }
stopCluster(cl)

write.csv(results, "data/hunt_freq.csv")

# ---------
fade_prob_map_hunt <- ggplot(
  data = results,
  aes(
    x = size,
    y = pct_hunt,
    fill = fadeout
  )
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    na.value = "gray21",
    breaks = round(seq(0, max(results$fadeout, na.rm = TRUE), length.out = 4),2),
    limits = c(0, max(results$fadeout, na.rm = TRUE)),
    name = expression(italic(p) * "(fadeout)")
  ) +
  scale_x_continuous(breaks = c(25, 250, 500, 750, 1000)) +
  labs(
    x = "herd size",
    y = expression("percent harvested ("*eta[h]*")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, "cm"),
    legend.text = element_text(size = 11)
  ) +
  coord_cartesian(expand = FALSE)

fade_prob_map_hunt

ggsave("fade_prob_map_hunt_freq.pdf", fade_prob_map, width = 7.5, height = 6, dpi = 300)

mean_map <- ggplot(
  data = results,
  aes(
    x = size,
    y = pct_hunt,
    fill = prev_ratio
  )
) +
  geom_tile(color = NA) +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    breaks = round(seq(0, max(results$prev_ratio, na.rm = TRUE), length.out = 5), 1),
    limits = round(c(0, max(results$prev_ratio, na.rm = TRUE)),2),
    name = "mean estimated\nprevalence ratio"
  ) +
  scale_x_continuous(breaks = c(25, 250, 500, 750, 1000)) +
  labs(
    x = "herd size",
    y = expression("percent harvested ("*eta[h]*")")
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, "cm"),
    legend.text = element_text(size = 11)
  ) +
  coord_cartesian(expand = FALSE)

mean_map

ggsave("prev_est_map_hunt_freq.pdf", mean_map, width = 7.5, height = 6, dpi = 300)

##-----------------------------------------------------------
## Sensitivity analysis
##-----------------------------------------------------------

######################
## Setup Parameters ##
######################
rm(list = ls())
library(rstudioapi);
setwd(dirname(getActiveDocumentContext()$path))
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
pth = "/home/webblab/Documents/Brandon/bTB_wildlife_code/" #path to main directory - must contain model .exe files
lambda_factor = 1.2
save_pars = T
save_runs = T
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
if(recalc){
  print('calculating parameter sets')
  lh <- LHS_ParSets(k_lim = size_range, 
                    scenario = infection,
                    infected = pct,
                    pct = T,
                    SS_prop = prop_superSpreader,
                    verbose = 0,
                    Num_LHS_sets = N_LHS_sets,
                    file.path = '/home/webblab/Documents/Brandon/bTB_wildlife_code/sensitivity_analysis/parameters/',
                    file.out = save_pars,
                    file.name = par_file,
                    seed_q = fix_q,
                    seed = def_seed)
}
if( ("lh" %in% ls()) ){
  print('found in workspace')
}else if( paste0('bTBwl_LHSpars_', infType, '_summary.csv') %in% list.files("./sensitivity_analysis/parameters/") ){
  print('loading from file')
  lh <- read.csv(file = paste0('./sensitivity_analysis/parameters/bTBwl_LHSpars_', infType, '.csv'))
}else{
  print('recalculating parameter sets')
  lh <- LHS_ParSets(k_lim = size_range, 
                    scenario = infection,
                    infected = pct,
                    pct = T,
                    SS_prop = prop_superSpreader,
                    verbose = 0,
                    Num_LHS_sets = N_LHS_sets,
                    file.path = '/home/webblab/Documents/Brandon/bTB_wildlife_code/sensitivity_analysis/parameters/',
                    file.out = save_pars,
                    file.name = par_file,
                    seed_q = fix_q,
                    seed = def_seed)
}
###############################

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

stopCluster(cl)

exclude <- lambda_screen |> filter(lambda >= 1e4) |> pull(i) |> unique()

if(def_seed == 43 & N_LHS_sets == 1000 & dim(lh)[1] == N_LHS_sets){
  lh <- lh[-c(exclude),] # lh[-c(363,367),] #remove these 2 entries due to excessive runtime issues -- getting massive vals for lambda O(10^15)
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
  #disease_parameters<-data.frame(beta=pars["beta"], area=1, p1=1, p2_q1=pars["p2_q1"], p2_q2=pars["p2_q2"], p2_q3=pars["p2_q3"], p2_q4=pars["p2_q4"], phi=pars["phi"], sigma1_mean=pars["sigma1_mean"], sigma1_rate=pars["sigma1_rate"])
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
  if(save_runs){
    write.csv(sto_out, file = paste0('./sensitivity_analysis/sens_runs/sens_run_',infType, '_', i, '.csv'))
  }
  
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

stopCluster(cl)

TOC <- Sys.time()
TOC-TIC

c('K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 'alpha_max', 'ksi', 
  'omega', 's', 'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi', 
  'sigma1_mean', 'sigma1_rate', 'start_q')

#R <- cbind(R, lh[rep(seq_len(nrow(lh)), each = reps), ])
# scale parameter values
bTB_wl_scaled <- as.data.frame(R %>% mutate_at(vars('K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                                                    'alpha_max', 'ksi', 'omega', 's',
                                                    'beta', "p1", 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 
                                                    'phi', 'sigma1_mean', 'sigma1_rate', 'start_q'),
                                               list(~scale(as.vector(.)))))


infType <- paste0(infType, "_freq")
# R <- read.csv("./sensitivity_analysis/summary_files/LHS_summary_seeded_q1.csv")
# save summary file and rescaled summary
write.csv(R, file = paste0("./data/raw/LHS_summary_", infType, ".csv" ), row.names = F)
remove(R)
write.csv(bTB_wl_scaled, file = paste0("./data/raw/LHS_scaled_summary_", infType, ".csv" ), row.names = F)
#########################################

#####################
## Monotonic Plots ##
#####################
# load in data if not in workspace
# define infType variable and averaged 
# infType = 'seeded_q1'; remove(bTB_wl_scaled)
averaged = T
setwd(pth)

if( ("bTB_wl_scaled" %in% ls()) ){
  print('found in workspace')
}else if( paste0('LHS_scaled_summary_', infType, '.csv') %in% list.files("./sensitivity_analysis/summary_files/") ){
  print('loading from file')
  bTB_wl_scaled <- read.csv(file = paste0('./data/raw/LHS_scaled_summary_', infType, '.csv'))
}
mono_plots = T
bTB_wl_scaled <- bTB_wl_scaled[,!(names(bTB_wl_scaled) %in% c('N', "pars", "S_0", "E1_0", "SuperS_0", "SuperE1_0", "I_0", "SuperI_0"))]
names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max', 'ksi', 'omega', 's', 
                          'beta', "p1", 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi', 'sigma1_mean', 'sigma1_rate', 
                          'start_q')

if(averaged){
  bTB_wl_agg <- aggregate.data.frame(x = bTB_wl_scaled, by = list(rep(1:dim(lh)[1], times=1, each=reps)), FUN = mean)
  bTB_wl_agg <- bTB_wl_agg[,names(bTB_wl_agg) %in% names(bTB_wl_scaled)]
  bTB_wl_scaled <- bTB_wl_agg
  infType <- paste0(infType, '_avg')
  write.csv(bTB_wl_scaled, file = paste0("./data/LHS_scaled_summary_", infType, ".csv" ), row.names = F)
}

#####################

setwd('./sensitivity_analysis/monotonic_plots/')

#--- Ensure consistent column names (do this BEFORE anything else) ---#
names(bTB_wl_scaled) <- make.names(names(bTB_wl_scaled))
# Check it worked
print(names(bTB_wl_scaled)[1:10])
# Should show: "Total.Infected", "fadeout", "fadeout.time", "Hunt.Prevalence", etc.

#--- Variables ---#
responses <- c("Total.Infected", "fadeout", "fadeout.time", "Hunt.Prevalence")

parameters <- c(
  "K","eta_hunt","eta_nat","theta","gamma","alpha_max","ksi","omega","s",
  "beta","p1","p2_q1","p2_q2","p2_q3","p2_q4","phi","sigma1_mean","sigma1_rate","start_q"
)

#--- Helper function ---#
plot_monotonicity <- function(param, infType) {
  if (!param %in% names(bTB_wl_scaled)) return(invisible())
  
  plots <- map(
    responses,
    \(resp) {
      if (!resp %in% names(bTB_wl_scaled)) return(NULL)
      ggplot(bTB_wl_scaled, aes_string(x = param, y = resp)) +
        geom_point(alpha = 0.6, size = 1.2) +
        # geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = "#2171b5") +
        labs(x = param, y = resp) +
        theme_bw(base_size = 20)
    }
  ) |> compact()
  
  if (length(plots) > 0) {
    outfile <- paste0("monotonic_", param, "_", infType, ".jpeg")
    jpeg(outfile, width = 840, height = 840, units = "px", res = 100)
    print(ggarrange(plotlist = plots, nrow = 2, ncol = 2))
    dev.off()
  }
}

setwd("./sensitivity_analysis/monotonic_plots/")
if (!dir.exists(getwd())) dir.create(getwd(), recursive = TRUE)

if (exists("mono_plots") && isTRUE(mono_plots)) {
  walk(parameters, plot_monotonicity, infType = infType)
}

setwd(pth)
#####################

names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max', 'ksi', 'omega', 's', 
                          'beta', "p1", 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'phi', 'sigma1_mean', 'sigma1_rate', 
                          'start_q')

# --------------
# PRCC
blue<-c('#2171b5')
responses <- c("Total Infected", "fadeout", "fadeout time", "Hunt Prevalence")
parameters <- c("K","eta_hunt","eta_nat","theta","gamma","alpha_max",
                "ksi","omega","s","beta","p1","phi","sigma1_mean","sigma1_rate")

# Run PRCCs for all responses
prcc_df <- map_dfr(responses, ~ run_prcc(
  data = bTB_wl_scaled,
  response = .x,
  parameters = parameters
))

param_order_prcc <- prcc_df %>%
  filter(response == "Total Infected") %>%
  arrange((est)) %>%  
  pull(var)

prcc_df <- prcc_df %>%
  mutate(var = factor(var, levels = param_order_prcc),
         est = ifelse(p.value >= 0.05, NA, est),
         response = factor(response, levels = c("Total Infected", "fadeout", "fadeout time", "Hunt Prevalence")),
         response = recode(response,
                           "Total Infected" = "Total Infected",
                           "fadeout" = "fadeout",
                           "fadeout time" = "fadeout time",
                           "Hunt Prevalence" = "Hunt Prevalence"))

sens <- ggplot(prcc_df, aes(y = est, x = var)) + 
  geom_col(fill = blue) +
  coord_flip() +
  facet_grid(~response) + 
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_x_discrete(labels = function(x) TeX(param_labels[x])) +
  labs(y = "proportional effect", x = NULL) +
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

sens

ggsave("sens_plot_prcc_freq.pdf", sens, width = 14, height = 6, dpi = 300)

# --------------
# Regression
responses <- c("`Total Infected`", "fadeout", "`fadeout time`", "`Hunt Prevalence`")

# Run LMs for all responses
lm_df <- map_dfr(responses, ~ run_lm(
  data = bTB_wl_scaled,
  response = .x,
  parameters = parameters
))

lm_df <- lm_df %>%
  mutate(parameter = factor(parameter, levels = param_order_prcc),
         estimate = ifelse(p_value >= 0.05, NA, estimate),
         response = factor(response, levels = c("`Total Infected`", "fadeout", "`fadeout time`", "`Hunt Prevalence`")),
         response = recode(response,
                           "`Total Infected`" = "Total Infected",
                           "fadeout" = "fadeout",
                           "`fadeout time`" = "fadeout time",
                           "`Hunt Prevalence`" = "Hunt Prevalence"))

reg <- ggplot(lm_df, aes(y = estimate, x = parameter)) + 
  geom_col(fill = blue) +
  coord_flip() +
  facet_grid(~response) + 
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_x_discrete(labels = function(x) TeX(param_labels[x])) +
  labs(y = "proportional effect", x = NULL) +
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

reg

ggsave("sens_plot_reg_freq.pdf", reg, width = 14, height = 6, dpi = 300)


# --------------
# Regression with interactions
responses <- c("`Total Infected`", "fadeout", "`fadeout time`", "`Hunt Prevalence`")

# Run LMs for all responses
lm_df <- map_dfr(responses, ~ run_lm_int(
  data = bTB_wl_scaled,
  response = .x,
  parameters = parameters
))

lm_df <- lm_df %>%
  filter(p_value < 0.01 & abs(estimate) > 0.001) %>%
  mutate(response = factor(response, levels = c("`Total Infected`", "fadeout", "`fadeout time`", "`Hunt Prevalence`")),
         response = recode(response,
                           "`Total Infected`" = "Total Infected",
                           "fadeout" = "fadeout",
                           "`fadeout time`" = "fadeout time",
                           "`Hunt Prevalence`" = "Hunt Prevalence"))

# Run LMs for all responses
lm_df <- map_dfr(responses, ~ run_lm_int(
  data = bTB_wl_scaled,
  response = .x,
  parameters = parameters
))

# Identify parameters that meet significance threshold in any response
sig_params <- lm_df %>%
  filter(p_value < 0.01, abs(estimate) > 0.01) %>%
  pull(parameter) %>%
  unique()

# Keep only rows for those parameters (across *all* responses)
lm_df <- lm_df %>%
  filter(parameter %in% sig_params) %>%
  mutate(response = factor(response, levels = c("`Total Infected`", "fadeout", "`fadeout time`", "`Hunt Prevalence`")),
         response = recode(response,
                           "`Total Infected`" = "Total Infected",
                           "fadeout" = "fadeout",
                           "`fadeout time`" = "fadeout time",
                           "`Hunt Prevalence`" = "Hunt Prevalence"),
         estimate = ifelse(p_value >= 0.05, 0, estimate))

# Determine unified parameter order based on Total Infected
param_order <- lm_df %>%
  filter(response == "Total Infected") %>%
  arrange(desc(estimate)) %>%
  pull(parameter)

# Apply the unified order to all facets
lm_df <- lm_df %>%
  mutate(parameter = factor(parameter, levels = param_order))

# produce a label vector in the same order as your factor levels
label_vec <- setNames(sapply(levels(lm_df$parameter), make_label),
                      levels(lm_df$parameter))

# plot
reg <- ggplot(lm_df, aes(y = estimate, x = parameter)) +
  geom_col(fill = blue) +
  coord_flip() +
  facet_grid(~response) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_x_discrete(labels = label_vec) +
  labs(y = "proportional effect", x = NULL) +
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

reg

ggsave("sens_plot_reg_int_freq.pdf", reg, width = 17, height = 10, dpi = 300)