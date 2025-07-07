# ================================
# Fadeout runs simulation Settings
initialize_environment()
cl <- initialize_cluster(detectCores()-2)

# Generate all (pct, size) combinations
sizes <- round(seq(from = 0, to = 750, length.out = 5)[-1], 0)
pct_seq <- round(seq(from = 0, to = 0.5, length.out = 4), 2)
param_grid <- expand.grid(pct = pct_seq, size = sizes)

# ================================
# Run Simulations in Parallel

fade_data <- foreach(row = iter(param_grid, by = "row"), .combine = 'rbind', .inorder = TRUE) %dopar% {
  model_out <- run_pipeline(
    type = "d", 
    years = 20, 
    infType = "seeded", 
    runtype = "fade_", 
    pct = row$pct * row$size,
    prop_superSpreader = 0, 
    reps = 500, 
    size = row$size,
    pth = getwd(),
    scaled_plots = F,
    gen_plots = F,
    save_runs = T
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
                        proportion_initially_exposed = row$pct,
                        fadeout_probability = p_fade,
                        fadeout_time = FOT)
  
  return(results)
}

stopCluster(cl)

# ===============================
# Plot fadeout simulations
# Process results

# Save data for future use
save(fade_data, file = "data/fade_data.RData")

load("data/fade_data.RData")

# Generate plots
fade_prob_map <- ggplot(data = fade_data, aes(x = herd_size, y = proportion_initially_exposed, fill = fadeout_probability)) + 
  geom_tile() +
  labs(fill = "Fadeout Probability") +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21') +
  ggtitle("Fadeout probability with fixed proportion initially exposed")

fade_time_map <- ggplot(data = fade_data, aes(x = herd_size, y = proportion_initially_exposed, fill = fadeout_time)) + 
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
