# ================================
# Clear workspace and set working directory
# ================================
rm(list = ls()) 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary functions
source(file = "bTBwl_func.R")

# ===============================
# Discrete model trajectories

load("data/pub_seeded_0.02-0.1-250-2025-03-13.RData")

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


Nplot <- ggplot() +
  geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) +
  geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
  scale_color_manual(values=SEIcols[7]) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_x_continuous(breaks=seq(0, 7*12,1), limits = c(0, 7*12))

Splot <- ggplot() +
  geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) +
  geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[6]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))

Eplot <- ggplot() +
  geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) +
  geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[4]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))

Iplot <- ggplot() +
  geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) +
  geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[2]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))

SSSplot <- ggplot() +
  geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) +
  geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[5]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))


print(ggarrange(Splot, Eplot,  Iplot, SSSplot, ncol=3, nrow=2))
plot <- ggarrange(Splot, Eplot, Iplot, SSSplot, ncol = 3, nrow = 2)
  
ggsave(filename = paste0("results/fig2a.jpeg"), plot = plot, width = 10, height = 5, dpi = 300, units = "in")
  
load("data/pub_spillover_0.02-0.1-100-2025-03-13.RData")

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


Nplot <- ggplot() +
  geom_line(data = sto_by_class[["N"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .9) +
  geom_line(data = det_by_class[["N"]], aes(x = time, y=value), size = 1) +
  scale_color_manual(values=SEIcols[7]) +
  scale_y_continuous(limits = c(0,NA)) +
  scale_x_continuous(breaks=seq(0, 7*12,1), limits = c(0, 7*12))

Splot <- ggplot() +
  geom_line(data = sto_by_class[["S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) +
  geom_line(data = det_by_class[["S"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[6]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = 'Susceptible' ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))

Eplot <- ggplot() +
  geom_line(data = sto_by_class[["E1"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) +
  geom_line(data = det_by_class[["E"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[4]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = 'Exposed' ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))

Iplot <- ggplot() +
  geom_line(data = sto_by_class[["I"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .7, show.legend = FALSE) +
  geom_line(data = det_by_class[["I"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[2]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = 'Infectious' ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))

SSSplot <- ggplot() +
  geom_line(data = sto_by_class[["SS_S"]], aes(x = tstep, y=value, color=variable, group = rep), size = 1, alpha = .6, show.legend = FALSE) +
  geom_line(data = det_by_class[["sS"]], aes(x = time, y=value), size = 1, show.legend = FALSE) +
  scale_color_manual(values=SEIcols[5]) +
  theme(text = element_text(size = 16)) +
  scale_y_continuous(name = expression('Susceptible'['SS']) ,limits = c(0,NA)) +
  scale_x_continuous(name = 'time (months)', breaks=seq(0, 7*12,12), limits = c(0, 7*12))


print(ggarrange(Splot, Eplot,  Iplot, SSSplot, ncol=3, nrow=2))
plot <- ggarrange(Splot, Eplot, Iplot, SSSplot, ncol = 3, nrow = 2) 
ggsave(filename = paste0("results/fig2b.jpeg"), plot = plot, width = 10, height = 5, dpi = 300, units = "in")

# ======================================================
# Fadeout maps
load("data/fade_data.RData")
fade_data <- read_csv("data/fade_data.csv")

# Generate plots
fade_prob_map <- ggplot(data = fade_data, aes(x = herd_size, y = number_initially_exposed, fill = fadeout_probability)) + 
  geom_tile() +
  labs(fill = "Fadeout Probability") +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21') +
  ggtitle("Fadeout probability with fixed proportion initially exposed")

fade_time_map <- ggplot(data = fade_data, aes(x = herd_size, y = number_initially_exposed, fill = fadeout_time)) + 
  geom_tile() +
  labs(fill = "Fadeout Time (months)") +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21') +
  ggtitle("Fadeout time with fixed proportion initially exposed")

fade_prob <- ggplot(data = fade_data_m, aes(x = number_initially_exposed, y = fadeout_probability, group = herd_size, color = herd_size)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral") +
  ggtitle("Fadeout Probability vs Proportion Initially Exposed")

fade_time <- ggplot(data = fadeTime_data_m, aes(x = number_initially_exposed, y = fadeout_time, group = herd_size, color = herd_size)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral") +
  ggtitle("Fadeout Time vs Proportion Initially Exposed")


ggsave(filename = paste0("results/Fadeout_maps_", Sys.Date(), ".jpeg"), plot = ggarrange(fade_prob_map, fade_time_map, ncol = 2, nrow = 1))
ggsave(filename = paste0("results/FadeoutProb_maps_", Sys.Date(), ".jpeg"), plot = fade_prob_map)
ggsave(filename = paste0("results/FadeoutTime_maps_", Sys.Date(), ".jpeg"), plot = fade_time_map)
ggsave(filename = paste0("results/Fadeout_plots_", Sys.Date(), ".jpeg"), plot = ggarrange(fade_prob, fade_time, ncol = 2, nrow = 1))
