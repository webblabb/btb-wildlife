###############################################################
## Compile and run code ##
###############################################################

system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_CTMC.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu bTB_wildlifeModel_DTMC.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")

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
# ------------------------------------------------------------

## -----------------------------------------------------------
## Figure 1: size x infType

param_grid <- expand_grid(
  size = c(250),
  infType = c("seeded", "spillover")
)

for (i in 1:nrow(param_grid)) {
  row <- param_grid[i,]
  
  run_wl_sim(size = row$size, infType = row$infType, save_plots = TRUE, ode = T) 
}

## -----------------------------------------------------------
## Figure 2: super-spreader prop x herd size

sizes  <- seq(from = 0, to = 750, length.out = 21)[-1]
xi_vec <- seq(from = 0, to = 1, length.out = 21)

param_grid <- expand_grid(
  size    = sizes,
  infType = "spillover",
  SS_prop = xi_vec
)

n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

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

write.csv(results, "data/ss_herdsize_meanprev.csv")

results <- read.csv("data/ss_herdsize_meanprev.csv")

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
                       breaks = round(seq(0,max(results$mean_prev), length.out = 4),1),   # manually specify desired breaks
                       limits = c(0, 1)) +
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

ggsave("figures/mean_prevalence_map.pdf", mean_map, width = 7.5, height = 6, dpi = 300)

# --------------
results <- read.csv("data/xi_phi_meanprev.csv") |>
  mutate(fadeout = ifelse(fadeout == 0, NA, fadeout))

fade_prob_map <- ggplot(
  data = results,
  aes(
    x = phi,
    y = prop_ss,
    fill = fadeout
  )
) +
  geom_tile(color = NA) +
  scale_fill_distiller(palette = "Spectral",
                       breaks = round(seq(0,max(results$mean_prev), length.out = 4),1),   # manually specify desired breaks
                       limits = c(0, 1)) +
  scale_x_continuous(breaks = round(seq(1,50,length.out = 5),0)) +
  labs(
    x = expression("super susceptible contact scaling ("*phi*")"),
    y = expression("proportion of super susceptibles ("*xi*")"),
    fill = expression(italic(p) * "(fadeout)")
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

ggsave("figures/strong_rare_fade_map.pdf", fade_prob_map, width = 7.5, height = 6, dpi = 300)


mean_sv_map <- ggplot(
  data = results,
  aes(
    x = phi,
    y = prop_ss,
    fill = mean_prev
  )
) +
  geom_tile(color = NA) +  # remove borders between tiles
  scale_fill_distiller(palette = "Spectral",
                       breaks = round(seq(0,max(results$mean_prev), length.out = 4),1),   # manually specify desired breaks
                       limits = c(0, 1)) +
  scale_x_continuous(breaks = round(seq(1,50,length.out = 5),0)) +
  labs(
    x = expression("super susceptible contact scaling ("*phi*")"),
    y = expression("proportion of super susceptibles ("*xi*")"),
    fill = "mean\ninfection\nprevalence"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, 'cm')
  ) +
  coord_cartesian(expand = FALSE)

mean_sv_map

ggsave("figures/strong_rare_ss_map.pdf", mean_sv_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## Figure 3: fadeout. herd size x number initially exposed 

sizes   <- seq(from = 0, to = 750, length.out = 21)[-1]
initNum <- 0:25

param_grid <- tidyr::expand_grid(
  size    = sizes,
  infType = "seeded",
  numExp  = initNum
)

n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

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

write.csv(results, "data/herdsize_initnum_fadeout.csv")

results <- read.csv("data/herdsize_initnum_fadeout.csv") |>
  mutate(fadeout = ifelse(fadeout == 0, NA, fadeout))

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

ggsave("figures/fade_prob_map_num.pdf", fade_prob_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## Figure 5 (a and b): herd size x harvest intensity
## -----------------------------------------------------------

n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
hunt_pct <- seq(from = 0, to = 1, length.out = 21)

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

results <- results |>
  mutate(fadeout = ifelse(fadeout == 0, NA, fadeout))

write.csv(results, "data/hunt_herdsize.csv")

results <- read.csv("data/hunt_herdsize.csv") |>
  mutate(fadeout = ifelse(fadeout == 0, NA, fadeout))
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
  breaks = round(seq(0, max(results$fadeout, na.rm = TRUE), length.out = 4),1),
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

ggsave("figures/fade_prob_map_hunt.pdf", fade_prob_map_hunt, width = 7.5, height = 6, dpi = 300)

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
    breaks = round(seq(0, max(results$prev_ratio, na.rm = TRUE), length.out = 5), 2),
    limits = round(c(0, max(results$prev_ratio, na.rm = TRUE)),4),
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

ggsave("figures/prev_est_map_hunt.pdf", mean_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## % reduction in p1  x  % reduction in farm-contact (p2’s)
## Columns: hunting level (eta_hunt = 0, 0.5, 1)
## Rows: scenarios (seeded, spillover)
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

# ---------- baselines ----------
p1_baseline <- 0.05   # your fixed baseline p1
baseline_p2 <- list(
  p2_q1 = 0.01,
  p2_q2 = 0.015,
  p2_q3 = 0.00,
  p2_q4 = 0.0075
)

# Build override lists from % reductions r in [0,1]
scale_p1 <- function(reduction) {
  # 0 => baseline; 1 => zero deer-deer contact
  list(deer_deer_contact_rate = (1 - reduction) * p1_baseline)
}

scale_p2 <- function(reduction) {
  f <- 1 - reduction
  list(
    p2_q1 = f * baseline_p2$p2_q1,
    p2_q2 = f * baseline_p2$p2_q2,
    p2_q3 = f * baseline_p2$p2_q3,
    p2_q4 = f * baseline_p2$p2_q4
  )
}

p1_red_grid    <- seq(0, 1, length.out = 21)   # % reduction of p1 (0..1)
farm_red_grid  <- seq(0, 1, length.out = 21)   # % reduction of farm-contact (0..1)
hunt_levels    <- c(0.2, 0.5, 0.95)               # fixed columns
scenarios      <- c("spillover")

# Fixed herd size for 2D surfaces
size_K  <- 250
# pctSeed <- 0.025 * size_K

param_grid <- tidyr::expand_grid(
  p1_red   = p1_red_grid,
  farm_red = farm_red_grid,
  hunt     = hunt_levels,
  infType  = scenarios
) %>%
  mutate(
    overrides = purrr::pmap(
      list(p1_red, farm_red, hunt),
      function(p1r, red2, hnt) c(
        scale_p1(p1r),                # scaled deer-deer contact
        scale_p2(red2),               # scaled farm-contact probs
        list(eta_hunt = hnt)          # hunting per facet col
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
                     #   years <- 1;  reps <- 500; pct <- 0
                     # }
                     
                     df <- run_wl_sim(
                       size                = size_K,
                       infType             = row$infType,
                       years               = 1,
                       reps                = 250,
                       pct                 = 0,
                       prop_superSpreader  = 0,
                       save_plots          = FALSE,
                       overrides           = row$overrides[[1]]
                     ) %>%
                       mutate(
                         p1_red    = row$p1_red,
                         farm_red  = row$farm_red,
                         pct_hunt  = row$overrides[[1]]$eta_hunt,
                         scenario  = row$infType
                       )
                     
                     # Mean infection prevalence (final timestep per rep)
                     prev_sum <- df %>%
                       group_by(rep, p1_red, farm_red, pct_hunt, scenario) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup() %>%
                       mutate(prev = `Total Infected` / N) %>%
                       group_by(p1_red, farm_red, pct_hunt, scenario) %>%
                       summarise(mean_prev = mean(prev, na.rm = TRUE), .groups = "drop")
                     
                     # p(fadeout) at final timestep
                     fade_sum <- df %>%
                       group_by(rep, p1_red, farm_red, pct_hunt, scenario) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup() %>%
                       group_by(p1_red, farm_red, pct_hunt, scenario) %>%
                       summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                     
                     rm(df); gc()
                     left_join(prev_sum, fade_sum,
                               by = c("p1_red","farm_red","pct_hunt","scenario"))
                   }

stopCluster(cl)

write.csv(results, "data/p1_p2_hunt.csv", row.names = FALSE)

results <- read.csv("data/beta_p2_hunt.csv") 

p1_red_grid    <- seq(0, 1, length.out = 21)   # % reduction of p1 (0..1)
farm_red_grid  <- seq(0, 1, length.out = 21)   # % reduction of farm-contact (0..1)
hunt_levels    <- c(0.2, 0.5, 0.95)               # fixed columns
scenarios      <- c("spillover")

results <- results %>%
  mutate(
    fadeout        = ifelse(fadeout == 0, NA_real_, fadeout),
    # make LaTeX/plotmath style labels: eta[h]==0, eta[h]==0.5, etc.
    hunt_label = sprintf("eta[h]==%s", number(pct_hunt, accuracy = 0.01, trim = TRUE)),
    hunt_label = factor(
      hunt_label,
      levels = sprintf("eta[h]==%s", number(hunt_levels, accuracy = 0.01, trim = TRUE))
    ),
    scenario_label = factor(scenario, levels = c("seeded","spillover"),
                            labels = c("Seeded", "Spillover")),
    p1_red_pct     = 100 * p1_red,
    farm_red_pct   = 100 * farm_red
  )


# p(fadeout): x = % reduction in p1, y = % reduction in farm-contact
fade_prob_map <- ggplot(
  data = results,
  aes(x = p1_red_pct, y = farm_red_pct, fill = fadeout)
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
  scale_x_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  labs(
    x = expression("% reduction in deer–deer transmission ("*beta[wild]*")"), # "% reduction in deer–deer contact (p1)",
    y = expression("% reduction in cattle/fomite transmission ("*beta[fq]*")")
  ) +
  facet_grid( cols = vars(hunt_label),labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

fade_prob_map

# mean infection prevalence
mean_prev_map <- ggplot(
  data = results,
  aes(x = p1_red_pct, y = farm_red_pct, fill = mean_prev)
) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    breaks = round(seq(0, max(results$mean_prev, na.rm = TRUE), length.out = 5), 2),
    limits = c(0, round(max(results$mean_prev, na.rm = TRUE),1)),
    name   = "mean infection\nprevalence"
  ) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  labs(
    x = expression("% reduction in deer–deer transmission ("*beta[wild]*")"), # "% reduction in deer–deer contact (p1)",
    y = expression("% reduction in cattle/fomite contacts ("*beta[fq]*")")
  ) +
  facet_grid( cols = vars(hunt_label),labeller = label_parsed) +
  theme_bw(base_size = 14) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid       = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 16),
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.key.width = unit(0.75, "cm"),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )

mean_prev_map

ggsave("figures/fade_prob_map_beta_p2_hunt.pdf", fade_prob_map, width = 12, height = 5, dpi = 300)
ggsave("figures/mean_prev_map_beta_p2_hunt.pdf", mean_prev_map, width = 12, height = 5, dpi = 300)

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

bTB_wl_scaled <- read.csv("data/LHS_scaled_summary_seeded_q1.csv")
lh <- read.csv("data/bTBwl_LHSpars_seeded_q1.csv")
#####################
## Monotonic Plots ##
#####################
# load in data if not in workspace
# define infType variable and averaged 
# infType = 'seeded_q1'; remove(bTB_wl_scaled)
averaged = T

if( ("bTB_wl_scaled" %in% ls()) ){
  print('found in workspace')
}else if( paste0('LHS_scaled_summary_', infType, '.csv') %in% list.files("./sensitivity_analysis/summary_files/") ){
  print('loading from file')
  bTB_wl_scaled <- read.csv(file = paste0('./data/raw/LHS_scaled_summary_', infType, '.csv'))
}
mono_plots = T
bTB_wl_scaled <- bTB_wl_scaled[,!(names(bTB_wl_scaled) %in% c('N', "pars", "p1", "phi","ksi", "S_0", "E1_0", "SuperS_0", "SuperE1_0", "I_0", "SuperI_0"))]
names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max', 'omega', 's', 
                          'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'sigma1_mean', 'sigma1_rate', 
                          'start_q')

reps = 500
if(averaged){
  bTB_wl_agg <- aggregate.data.frame(x = bTB_wl_scaled, by = list(rep(1:dim(lh)[1], times=1, each=reps)), FUN = mean)
  bTB_wl_agg <- bTB_wl_agg[,names(bTB_wl_agg) %in% names(bTB_wl_scaled)]
  bTB_wl_scaled <- bTB_wl_agg
  infType = "seeded_q1"
  infType <- paste0(infType, '_avg')
  write.csv(bTB_wl_scaled, file = paste0("./data/LHS_scaled_summary_", infType, ".csv" ), row.names = F)
}

#####################

setwd('./figures/supp/')

#--- Ensure consistent column names (do this BEFORE anything else) ---#
names(bTB_wl_scaled) <- make.names(names(bTB_wl_scaled))
# Check it worked
print(names(bTB_wl_scaled)[1:10])
# Should show: "Total.Infected", "fadeout", "fadeout.time", "Hunt.Prevalence", etc.

#--- Variables ---#
responses <- c("Total.Infected", "fadeout", "fadeout.time", "Hunt.Prevalence")

parameters <- c(
  "K","eta_hunt","eta_nat","theta","gamma","alpha_max","omega","s",
  "beta","p2_q1","p2_q2","p2_q3","p2_q4","sigma1_mean","sigma1_rate","start_q"
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

if (exists("mono_plots") && isTRUE(mono_plots)) {
  walk(parameters, plot_monotonicity, infType = infType)
}

#####################

names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max', 'omega', 's', 
                          'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'sigma1_mean', 'sigma1_rate', 
                          'start_q')

# --------------
# PRCC
setwd("../..")
blue<-c('#2171b5')
responses <- c("Total Infected", "fadeout", "fadeout time", "Hunt Prevalence")
parameters <- c("K","eta_hunt","eta_nat","theta","gamma","alpha_max",
                "omega","s","beta","sigma1_mean","sigma1_rate")

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
  theme_bw(base_size = 20) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18)
  )

sens

ggsave("figures/sens_plot_prcc.pdf", sens, width = 14, height = 7.5, dpi = 300)

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
  theme_bw(base_size = 20) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18)
  )

reg

ggsave("figures/sens_plot_reg.pdf", reg, width = 14, height = 7.5, dpi = 300)


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
  arrange((estimate)) %>%
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
  theme_bw(base_size = 20) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18)
  )

reg

ggsave("figures/sens_plot_reg_int.pdf", reg, width = 18, height = 10, dpi = 300)

## Supplementary figures
## -----------------------------------------------------------
## Figure S1 (supplement)
## -----------------------------------------------------------

sizes   <- seq(0, 750, length.out = 31)[-1]
pct_vec <- seq(0, 0.5, length.out = 21)

param_grid <- expand_grid(
  size    = sizes,
  infType = "seeded",
  pct     = pct_vec
) |>
  mutate(overrides = list(list(SS_prop = 0)))

for (i in seq_len(nrow(param_grid))) {
  row <- param_grid[i, ]
  run_wl_sim(
    size       = row$size,
    infType    = row$infType,
    years      = 20,
    reps       = 500,
    save_plots = FALSE,
    pct        = row$pct,
    overrides  = row$overrides[[1]]
  )
}


###############################################################
## Herd size fixed (250); vary initial exposure proportion × harvest
###############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(scales)
})

dir.create("data",     showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",  showWarnings = FALSE, recursive = TRUE)

source("bTBwl_func.R")

## ---- settings ----
size_fixed   <- 250
infType      <- "seeded"
years_run    <- 20
reps_run     <- 200
hunt_pct     <- seq(0, 0.5, length.out = 10)
init_prop    <- seq(0, 0.25, length.out = 10)  # proportion of herd initially exposed
ss_prop      <- 0.05                              # super-susceptible proportion held fixed

## ---- parallel ----
n.cores <- max(1, floor(parallel::detectCores() * 3/4))
cl <- makeCluster(n.cores)
registerDoParallel(cl)

param_grid <- tidyr::expand_grid(
  size      = seq(0, 750, length.out = 10)[-1],
  infType   = "seeded",
  init_prop = c(0.05, 0.1, 0.25),
  hunt      = seq(0, 0.5, length.out = 10)
) |>
  mutate(overrides = purrr::map(hunt, ~list(eta_hunt = .x)))

## ---- sims ----
results <- foreach(row = iter(param_grid, by = "row"),
                   .packages = c("tidyverse","deSolve"),
                   .combine = "rbind") %dopar% {
                     
                     df <- run_wl_sim(
                       size                = row$size,
                       infType             = row$infType,
                       years               = years_run,
                       reps                = reps_run,
                       pct                 = row$init_prop * row$size,  # number initially exposed
                       prop_superSpreader  = 0.05,
                       save_plots          = FALSE,
                       save_runs = FALSE,
                       overrides           = row$overrides[[1]]
                     ) %>%
                       mutate(
                         size       = row$size,
                         pct_hunt   = row$overrides[[1]]$eta_hunt,
                         init_prop  = row$init_prop,
                         prev_ratio = `Hunt Prevalence`
                       )
                     
                     prev_sum <- df %>%
                       filter(quarter == 4) %>%
                       group_by(size, pct_hunt, init_prop) %>%
                       summarise(prev_ratio = mean(prev_ratio, na.rm = TRUE), .groups = "drop")
                     
                     fade_sum <- df %>%
                       group_by(rep, size, pct_hunt, init_prop) %>%
                       slice_max(tstep, with_ties = FALSE) %>%
                       ungroup() %>%
                       group_by(size, pct_hunt, init_prop) %>%
                       summarise(fadeout = mean(fadeout, na.rm = TRUE), .groups = "drop")
                     
                     rm(df); gc()
                     left_join(prev_sum, fade_sum, by = c("size","pct_hunt","init_prop"))
                   }

stopCluster(cl)

write.csv(results, file.path("data", "hunt_initprop.csv"), row.names = FALSE)

## ---- plots ----
fade_map <- ggplot(results, aes(x = size, y = pct_hunt, fill = fadeout)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    na.value = "gray21",
    breaks = round(seq(0, max(results$fadeout, na.rm = TRUE), length.out = 4), 2),
    limits = c(0, max(results$fadeout, na.rm = TRUE)),
    name = expression(italic(p) * "(fadeout)")
  ) +
  facet_grid(~init_prop) +
  scale_x_continuous(breaks = seq(min(init_prop), max(init_prop), by = 0.05)) +
  labs(
    x = "herd size",
    y = expression("percent harvested ("*eta[h]*")"),
    title = "Fadeout probability by initial exposure and harvest (K = 250)"
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

fade_map

ggsave(file.path("figures", "fade_prob_map_initprop_K250.pdf"),
       fade_map, width = 7.5, height = 6, dpi = 300)

mean_map <- ggplot(results, aes(x = size, y = pct_hunt, fill = prev_ratio)) +
  geom_tile(color = NA) +
  scale_fill_distiller(
    palette = "Spectral",
    direction = -1,
    breaks = round(seq(0, max(results$prev_ratio, na.rm = TRUE), length.out = 5), 2),
    limits = c(0, round(max(results$prev_ratio, na.rm = TRUE), 2)),
    name = "mean estimated\nprevalence ratio"
  ) +
  facet_grid(~init_prop) +
  scale_x_continuous(breaks = seq(min(init_prop), max(init_prop), by = 0.05)) +
  labs(
    x = "herd size",
    y = expression("percent harvested ("*eta[h]*")"),
    title = "Mean estimated prevalence ratio (quarter 4) (K = 250)"
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

ggsave(file.path("figures", "prev_est_map_initprop_K250.pdf"),
       mean_map, width = 7.5, height = 6, dpi = 300)
