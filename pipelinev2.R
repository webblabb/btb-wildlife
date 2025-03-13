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
                          years = 10,
                          pct = 0.02, 
                          prop_superSpreader = c(0, 0.05,0.1),
                          size = c(10, 50, 100, 250, 500, 750, 1000),
                          reps = 500) |>
  filter(!(type == "c" & size > 500)) |> # subset the runs a bit so we don't have ones that last forever
  filter(!(type == "c" & prop_superSpreader > 0 & infType == "spillover")) |> # subset the runs a bit so we don't have ones that last forever
  filter(!(type == "d" & prop_superSpreader != 0.1 & infType == "spillover")) |> # subset the runs a bit so we don't have ones that last forever
  filter(!(type == "c" & prop_superSpreader != 0 & infType == "spillover")) |> # subset the runs a bit so we don't have ones that last forever
  mutate(reps = case_when(type == 'c' & infType == "spillover" ~ 200,
                          TRUE ~ reps),
         years = case_when(type == 'c' & infType == "spillover" ~ 3,
                           TRUE ~ years))
  

# ================================
# Run Simulations in Parallel

foreach(row = iter(param_grid, by = "row"), .combine = 'cbind', .inorder = TRUE) %dopar% {
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

# ===============================
# Plot fadeout simulations
# Process results
sizes <- seq(from = 0, to = 750, length.out = 31)[-1]
pct_seq <- seq(from = 0, to = 0.5, length.out = 21)

fade_data <- matrix(data = results_list[1,], nrow = length(sizes), ncol = length(pct_seq), byrow = TRUE)
fadeTime_data <- matrix(data = results_list[2,], nrow = length(sizes), ncol = length(pct_seq), byrow = TRUE)

fade_data_m <- as.data.frame(fade_data) %>%
  rownames_to_column(var = "herd_size") %>%
  pivot_longer(cols = -herd_size, names_to = "proportion_initially_exposed", values_to = "fadeout_probability") %>%
  mutate(across(c(herd_size, proportion_initially_exposed), as.numeric))

fadeTime_data_m <- as.data.frame(fadeTime_data) %>%
  rownames_to_column(var = "herd_size") %>%
  pivot_longer(cols = -herd_size, names_to = "proportion_initially_exposed", values_to = "fadeout_time") %>%
  mutate(across(c(herd_size, proportion_initially_exposed), as.numeric),
         fadeout_time = ifelse(fadeout_time > 1e6, NA, fadeout_time))

# Save data for future use
save(fade_data_m, file = "fade_PlotDat.RData")
save(fadeTime_data_m, file = "fadeTime_PlotDat.RData")

# Generate plots
fade_prob_map <- ggplot(data = fade_data_m, aes(x = herd_size, y = proportion_initially_exposed, fill = fadeout_probability)) + 
  geom_tile() +
  labs(fill = "Fadeout Probability") +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21') +
  ggtitle("Fadeout probability with fixed proportion initially exposed")

fade_time_map <- ggplot(data = fadeTime_data_m, aes(x = herd_size, y = proportion_initially_exposed, fill = fadeout_time)) + 
  geom_tile() +
  labs(fill = "Fadeout Time (months)") +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21') +
  ggtitle("Fadeout time with fixed proportion initially exposed")

fade_prob <- ggplot(data = fade_data_m, aes(x = proportion_initially_exposed, y = fadeout_probability, group = herd_size, color = herd_size)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral") +
  ggtitle("Fadeout Probability vs Proportion Initially Exposed")

fade_time <- ggplot(data = fadeTime_data_m, aes(x = proportion_initially_exposed, y = fadeout_time, group = herd_size, color = herd_size)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral") +
  ggtitle("Fadeout Time vs Proportion Initially Exposed")


ggsave(filename = paste0("results/Fadeout_maps_", Sys.Date(), ".jpeg"), plot = ggarrange(fade_prob_map, fade_time_map, ncol = 2, nrow = 1))
ggsave(filename = paste0("results/FadeoutProb_maps_", Sys.Date(), ".jpeg"), plot = fade_prob_map)
ggsave(filename = paste0("results/FadeoutTime_maps_", Sys.Date(), ".jpeg"), plot = fade_time_map)
ggsave(filename = paste0("results/Fadeout_plots_", Sys.Date(), ".jpeg"), plot = ggarrange(fade_prob, fade_time, ncol = 2, nrow = 1))


setwd(getwd())