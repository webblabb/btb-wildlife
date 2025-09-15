# ------------------------------------------------------------------------------
# Run all simulations
# Note that this can take awhile
source("pipeline.R")

# ------------------------------------------------------------------------------
# Super susceptibles plot

#load data from plot directory if not in workspace and also in directory
if(!("ss_data_full" %in% ls()) && !("ss_stat_full" %in% ls()) ){
  if(!exists(ss_data_name) && !exists(ss_stat_name)){
    load(file = ss_data_name)
    load(file = ss_stat_name)
  }
}
#####################

################
# metric plots #
################
min_mat <- matrix(data = ss_stat_full$min, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(min_mat) <- sizes
colnames(min_mat) <- range_vec

median_mat <- matrix(data = ss_stat_full$median, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(median_mat) <- sizes
colnames(median_mat) <- range_vec

mean_mat <- matrix(data = ss_stat_full$mean, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(mean_mat) <- sizes
colnames(mean_mat) <- range_vec

max_mat <- matrix(data = ss_stat_full$max, nrow = length(sizes), ncol = length(range_vec), byrow = T)
rownames(max_mat) <- sizes
colnames(max_mat) <- range_vec

melty_min_mat <- melt(min_mat)
melty_median_mat <- melt(median_mat)
melty_mean_mat <- melt(mean_mat)
melty_max_mat <- melt(max_mat)

names(melty_min_mat) <- c('herd size','proportion of super susceptibles','min. infection prevalence')
names(melty_median_mat) <- c('herd size','proportion of super susceptibles','median infection prevalence')
names(melty_mean_mat) <- c('herd size','proportion of super susceptibles','mean infection prevalence')
names(melty_max_mat) <- c('herd size','proportion of super susceptibles','max. infection prevalence')

min_map <- ggplot(data = melty_min_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`min. infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_median_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral")

median_map <- ggplot(data = melty_median_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`median infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_median_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral")
median_map <- annotate_figure(median_map, top = text_grob("Median infection prevalence\nwith constant force of infection", 
                                                          color = "black", face = "bold", size = 18))
median_map

mean_map <- ggplot(data = melty_mean_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`mean infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_mean_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral")
mean_map <- annotate_figure(mean_map, top = text_grob("Mean infection prevalence with constant force of infection", 
                                                      color = "black", face = "bold", size = 18))
mean_map

max_map <- ggplot(data = melty_max_mat, aes(x=`herd size`, y=`proportion of super susceptibles`, fill=`max. infection prevalence`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(melty_median_mat[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral")
################

################
## Regression ##
################

#ss proportion: No interaction regression 
setwd(paste0(pth,'Disc_runs/ss_runs/regression/'))
single.model.ss_prev = as.formula( `mean inf. prevalence` ~ `proportion of super susceptibles` + `herd size`)
lm.ss_prev = lm(single.model.ss_prev, data = as.data.frame(scale(melty_mean_mat)))
write.csv(lm.ss_prev, file=paste0("lmOutput_mean_ssPrev_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.ss_prev)
tab_model(lm.ss_prev, file = paste0("reg_ssPrev_", infType, Sys.Date(), ".doc"))
setwd(pth)
################

# ------------------------------------------------------------------------------
# Fadeout plot
#load data from plot directory if not in workspace and also in directory
if(!("fade_Ndata_m" %in% ls()) && !("fadeTime_Ndata_m" %in% ls()) ){
  if(!exists(fade_data_name) && !exists(fadeTime_data_name)){
    load(file = fade_data_name)
    load(file = fadeTime_data_name)
  }
}

fade_Ndata_m$`fadeout probability`[fade_Ndata_m$`fadeout probability` == 0] <- NA

fade_prob_map <- ggplot(data = fade_Ndata_m, aes(x=`herd size`, y=`number initially exposed`, fill=`fadeout probability`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fade_Ndata_m[3]), 1)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_prob_map <- annotate_figure(fade_prob_map, top = text_grob("Fadeout probability with fixed number initially exposed", 
                                                                color = "black", face = "bold", size = 18))
fade_prob_map

fade_time_map <- ggplot(data = fadeTime_Ndata_m, aes(x=`herd size`, y=`number initially exposed`, fill=`fadeout time`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fadeTime_Ndata_m[3]), 1)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  labs(fill='fadeout\ntime (months)') +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_time_map <- annotate_figure(fade_time_map, top = text_grob("Fadeout time with fixed number initially exposed", 
                                                                color = "black", face = "bold", size = 18))
fade_time_map

fade_prob <- ggplot(data = fade_Ndata_m, aes(x=`number initially exposed`, y=`fadeout probability`, group=`herd size`, color=`herd size`)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral")

fade_time <- ggplot(data = fadeTime_Ndata_m, aes(x=`number initially exposed`, y=`fadeout time`, group=`herd size`, color=`herd size`)) + 
  geom_line() +
  geom_point() +
  scale_colour_distiller(palette = "Spectral")


if(save_plots_final){
  jpeg(filename = paste0("Fadeout_maps", '_num-', Sys.Date(), ".jpeg"))
  print(ggarrange(fade_prob_map, fade_time_map, ncol=2, nrow=1))
  dev.off()
  
  jpeg(filename = paste0("FadeoutProb_maps", '_num-', Sys.Date(), ".jpeg"))
  print(fade_prob_map)
  dev.off()
  
  jpeg(filename = paste0("FadeoutTime_maps", '_num-', Sys.Date(), ".jpeg"))
  print(fade_time_map)
  dev.off()
  
  jpeg(filename = paste0("Fadeout_plots", '_num-', Sys.Date(), ".jpeg"))
  print(ggarrange(fade_prob, fade_time, ncol=2, nrow=1))
  dev.off()
}
setwd(pth)
####################

################
## Regression ##
################

setwd(paste0(pth,'Disc_runs/fadeout_runs/regression/'))
single.model.fadeNum = as.formula( `fadeout probability` ~ `number initially infected` + `herd size`)
lm.fadeNum = lm(single.model.fadeNum, data = as.data.frame(scale(fade_Ndata_m)))
write.csv(lm.fadeNum, file=paste0("lmOutput_FadeNum_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.fadeNum)
tab_model(lm.fadeNum, file = paste0("reg_FadeNum_", infType, Sys.Date(), ".doc"))
setwd(pth)


# ------------------------------------------------------------------------------
# Hunt plot 1
#####################
## Prev. Est. maps ##
#####################
mean_mat_melt <- hunt_stats_data |>
  select(`herd size`, `percent harvested`, `mean estimated prevalence ratio` = `mean prevalence`)

# median_mat_melt <- melt(median_mat)
# names(median_mat_melt) <- c('herd size', 'percent harvested', 'median estimated prevalence ratio')
# mean_mat_melt <- melt(mean_mat)
# names(mean_mat_melt) <- c('herd size', 'percent harvested', 'mean estimated prevalence ratio')
# 
# median_map <- ggplot(data = median_mat_melt, aes(x=`herd size`, y=`percent harvested`, fill=`median estimated prevalence ratio`)) + 
#   geom_tile() +
#   labs(fill=str_wrap(names(median_mat_melt[3]), 18)) +
#   theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
#   scale_fill_distiller(palette = "Spectral")
# median_map <- annotate_figure(median_map, top = text_grob("Median estimated prevalence ratio", 
#                                                           color = "black", face = "bold", size = 18))
# median_map

mean_map <- ggplot(data = mean_mat_melt, aes(x=`herd size`, y=`percent harvested`, fill=`mean estimated prevalence ratio`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(mean_mat_melt[3]), 18)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 12)) +
  scale_fill_distiller(palette = "Spectral")
mean_map <- annotate_figure(mean_map, top = text_grob("Mean estimated prevalence ratio", 
                                                      color = "black", face = "bold", size = 18))
mean_map

#####################

################
## Regression ##
################
#mean hunt ratio: No interaction regression 
single.model.mean_est = as.formula( `mean estimated prevalence ratio` ~ `percent harvested` + `herd size`)
lm.mean_est = lm(single.model.mean_est, data = as.data.frame(scale(mean_mat_melt)))
write.csv(lm.mean_est, file=paste0("lmOutputMeanEst_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.mean_est)
sjPlot::tab_model(lm.mean_est, file = paste0("reg_huntMean_", infType, Sys.Date(), ".doc"))
setwd(pth)

# ------------------------------------------------------------------------------
# Fadeout plot 2
# mean_mat <- matrix(data = prevEst_full$mean, nrow = length(sizes), ncol = length(range_vec), byrow = T)
# rownames(mean_mat) <- sizes
# colnames(mean_mat) <- range_vec
# 
# median_mat <- matrix(data = prevEst_full$median, nrow = length(sizes), ncol = length(range_vec), byrow = T)
# rownames(median_mat) <- sizes
# colnames(median_mat) <- range_vec
#####################

######################
## Fadeout heatmaps ##
######################
fade_data_full <- fadeout_matrix
fade_data_melt <- melt(fade_data_full)
names(fade_data_melt) <- c('herd size', 'percent harvested', 'fadeout probability')
fade_data_melt$`fadeout probability`[fade_data_melt$`fadeout probability` == 0] <- NA

fade_prob_map <- ggplot(data = fade_data_melt, aes(x=`herd size`, y=`percent harvested`, fill=`fadeout probability`)) + 
  geom_tile() +
  labs(fill=str_wrap(names(fade_data_melt[3]), 0)) +
  theme(legend.key.width = unit(.75, 'cm'), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  scale_fill_distiller(palette = "Spectral", na.value = 'gray21')
fade_prob_map <- annotate_figure(fade_prob_map, top = text_grob("Fadeout probability across different hunter harvest rates", 
                                                                color = "black", face = "bold", size = 18))
fade_prob_map

######################

################
## Regression ##
################

#hunt fadeout: No interaction regression 
setwd(paste0(pth,'Disc_runs/hunt_runs/regression/'))
single.model.hunt_fade = as.formula( `fadeout probability` ~ `percent harvested` + `herd size`)
lm.hunt_fade = lm(single.model.hunt_fade, data = as.data.frame(scale(fade_data_melt)))
write.csv(lm.hunt_fade, file=paste0("lmOutput_huntFade_", infType, Sys.Date(), ".csv"), row.names = F)
summary(lm.hunt_fade)
sjPlot::tab_model(lm.hunt_fade, file = paste0("reg_huntFade_", infType, Sys.Date(), ".doc"))
setwd(pth)
################

# -------------------------------------------------------
# Sensitivity analysis plot
bTB_wl_scaled<-R_scaled

####################
## Logistic Model ##
####################
library(sjPlot)
library(sjmisc)
library(sjlabelled)

single.model = as.formula( paste( 'RESPONSE~', paste(parameters, collapse = '+') ) )
p_vals_lm <- data.frame(parameters = parameter_names)
#TruPrev: No interaction regression 
single.model.TruPrev = update.formula(single.model, `Total Infected` ~ .)
lm.TruPrev.PRCC = lm(single.model.TruPrev, data = bTB_wl_scaled)
write.csv(tidy(lm.TruPrev.PRCC), file=paste0("lmOutputTruPrev_", infType, ".csv"), row.names = F)
summary(lm.TruPrev.PRCC)
tab_model(lm.TruPrev.PRCC)

values.in.PRCC <- lm.TruPrev.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC)) #regression results scaled by the largest abs regression value, premsize here

lmTruPrevScale<-data.frame(values.in.PRCC.scaled, row.names = parameters)

p_vals_lm$TruPrev <- summary(lm.TruPrev.PRCC)$coefficients[-1,4]


#Fadeout: No interaction regression 
single.model.Fadeout = update.formula(single.model, fadeout ~ .)
lm.Fadeout.PRCC = lm(single.model.Fadeout, data = bTB_wl_scaled)
write.csv(lm.Fadeout.PRCC, file=paste0("lmOutputFadeout_", infType, ".csv"), row.names = F)
summary(lm.Fadeout.PRCC)
tab_model(lm.Fadeout.PRCC)

values.in.PRCC <- lm.Fadeout.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC))  #regression results scaled by the largest abs regression value

lmFadeoutScale<-data.frame(values.in.PRCC.scaled,row.names = parameters)

p_vals_lm$Fadeout <- summary(lm.Fadeout.PRCC)$coefficients[-1,4]

#Fadeout Time: No interaction regression 
single.model.FadeoutTime = update.formula(single.model, `fadeout time` ~ .)
lm.FadeoutTime.PRCC = lm(single.model.FadeoutTime, data = bTB_wl_scaled)
write.csv(lm.FadeoutTime.PRCC, file=paste0("lmOutputFadeoutTime_", infType, ".csv"), row.names = F)
summary(lm.FadeoutTime.PRCC)
tab_model(lm.FadeoutTime.PRCC)

values.in.PRCC <- lm.FadeoutTime.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC))  #regression results scaled by the largest abs regression value

lmFadeoutTimeScale<-data.frame(values.in.PRCC.scaled,row.names = parameters)

p_vals_lm$FadeoutTime <- summary(lm.FadeoutTime.PRCC)$coefficients[-1,4]

#Hunt Prevalence: No interaction regression 
single.model.HuntPrev = update.formula(single.model, `Hunt Prevalence` ~ .)
lm.HuntPrev.PRCC = lm(single.model.HuntPrev, data = bTB_wl_scaled)
write.csv(lm.HuntPrev.PRCC, file=paste0("lmOutputHuntPrev_", infType, ".csv"), row.names = F)
summary(lm.HuntPrev.PRCC)
tab_model(lm.HuntPrev.PRCC)

values.in.PRCC <- lm.HuntPrev.PRCC$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled<-values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC))  #regression results scaled by the largest abs regression value

lmHuntPrevScale<-data.frame(values.in.PRCC.scaled,row.names = parameters)
p_vals_lm$HuntPrev <- summary(lm.HuntPrev.PRCC)$coefficients[-1,4]

#Plot 
lmTruPrevScale$parameters <- parameter_names
lmFadeoutTimeScale$parameters <- parameter_names
lmHuntPrevScale$parameters <- parameter_names
lmFadeoutScale$parameters <- parameter_names

#ordering values for plot
lmTruPrevScale <-  lmTruPrevScale[order(match(lmTruPrevScale$parameters,       prcc.TruPrev$parameters)),]
lmFadeoutTimeScale <-lmFadeoutTimeScale[order(match(lmFadeoutTimeScale$parameters,   prcc.TruPrev$parameters)),]
lmHuntPrevScale <-  lmHuntPrevScale[order(match(lmHuntPrevScale$parameters,       prcc.TruPrev$parameters)),]
lmFadeoutScale <- lmFadeoutScale[order(match(lmFadeoutScale$parameters,     prcc.TruPrev$parameters)),]

p_vals_lm <- p_vals_lm[order(match(p_vals_lm$parameters,     prcc.TruPrev$parameters)),]


lmScale <- cbind(TruPrev=lmTruPrevScale, FadeoutTime=lmFadeoutTimeScale$values.in.PRCC.scaled, HuntPrev=lmHuntPrevScale$values.in.PRCC.scaled, Fadeout=lmFadeoutScale$values.in.PRCC.scaled)
# uncomment the line below if you want lm plots in ascending order -- i.e. probably not identical to PRCC ordering
#lmScale <- lmScale[order(lmScale$TruPrev.values.in.PRCC.scaled),]
lmScale <- lmScale[,c("TruPrev.values.in.PRCC.scaled", "TruPrev.parameters", "Fadeout", "FadeoutTime", "HuntPrev")]
lmScale.long <- melt(lmScale,id=c("TruPrev.parameters"))

p_vals_lm.long <- melt(p_vals_lm,id=c("parameters"))

colnames(lmScale.long)[1] <- "parameters"
lmScale.long$parameters2 <- factor(lmScale.long$parameters,levels=lmScale$TruPrev.parameters)
lmScale.long$time_ordered <- factor(lmScale.long$variable,levels=c('TruPrev.values.in.PRCC.scaled','Fadeout','FadeoutTime','HuntPrev'))

#Renaming facets
variable_names2 <- list(
  "TruPrev.values.in.PRCC.scaled" = "True Prevalence",
  "Fadeout" = "Fadeout",
  "FadeoutTime" = "Fadeout Time" ,
  "HuntPrev" = "Observed Prevalence"
)

variable_labeller <- function(variable,value){
  return(variable_names2[value])
}


if(significant){
  lmScale.long$value[p_vals_lm.long$value >= .05] <- NA
}

#plot single model results
reg <- ggplot() + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14)) +
  geom_col(data = lmScale.long, aes(y = value, x = parameters2, fill = blue)) + 
  scale_y_discrete(limits=c(-1, -.5, 0, .5, 1)) +
  coord_flip() + 
  facet_grid(~time_ordered, labeller= variable_labeller) +
  scale_fill_manual (values = blue) + 
  theme(panel.spacing.y = unit(1.5, "lines")) + 
  theme(strip.text = element_text(size = 12))
reg <- annotate_figure(reg, top = text_grob("Regression Analysis", 
                                            color = "black", face = "bold", size = 18))
reg

jpeg(filename = paste0("figures/SingleModel_Trans_", infType, ".jpeg"), width = 1440, height = 840, units = 'px', res = 100)
reg
dev.off()


####################