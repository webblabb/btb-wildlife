# ============================================================
# Setup
# ============================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(grid)
library(broom)
library(purrr)
library(sensitivity)   # for pcc()

# ------------------------------------------------------------
# Inputs you already have (verify/adjust if names changed)
# ------------------------------------------------------------
# bTB_wl_scaled is assumed already cleaned/renamed like you posted.
# Response + parameter sets:
responses  <- c('Total Infected','fadeout','fadeout time','Hunt Prevalence')
parameters <- setdiff(names(bTB_wl_scaled), responses)

# Pretty labels (mapped to your parameters order)
parameter_names <- c(
  'carrying.capacity','hunt.mortality','base.mortality','dens.dep.asymetry',
  'dens.dep.mortality','max.birth.rate','proportion.SS','birth.timing',
  'birth.duration','transmission.rate','SS.contact.factor',
  'latency.mean','latency.rate','start.quarter'
)

# If your parameters exclude 'start_q' then drop its pretty label:
if (!"start_q" %in% parameters) {
  parameter_names <- parameter_names[parameter_names != 'start.quarter']
}

# Blue for bars
blue <- '#2171b5'

# ============================================================
# 1) PRCC utilities (sensitivity::pcc)
# ============================================================

# Compute PRCC + bootstrap CIs for one response
.compute_prcc <- function(dat, params, resp_col, nboot = 500, write_stub = NULL, infType = NULL) {
  # Build X (predictors) and y (response)
  Xdf <- dat[, params, drop = FALSE]
  X   <- data.matrix(Xdf)
  y   <- dat[[resp_col]]; if (is.data.frame(y)) y <- y[[1]]; y <- as.numeric(y)
  
  # Drop non-finite
  ok <- is.finite(y) & apply(X, 1, function(r) all(is.finite(r)))
  X <- X[ok, , drop = FALSE]; y <- y[ok]
  
  # Make safe column names for pcc()
  orig  <- colnames(X)
  safe  <- make.names(orig, unique = TRUE)
  colnames(X) <- safe
  
  # Run PRCC
  fit <- pcc(X = as.data.frame(X), y = y, rank = TRUE, nboot = nboot)
  
  out <- as.data.frame(fit$PRCC)
  out$parameters_safe <- rownames(out); rownames(out) <- NULL
  out$parameters <- orig[match(out$parameters_safe, safe)]
  
  out <- out %>%
    transmute(
      parameters,
      est    = original,
      min_ci = `min. c.i.`,
      max_ci = `max. c.i.`
    )
  
  # Optional CSV
  if (!is.null(write_stub) && !is.null(infType)) {
    write.csv(out, paste0(write_stub, "_", infType, ".csv"), row.names = FALSE)
  }
  out
}

# Run PRCC for all responses; returns long & wide forms
.run_all_prcc <- function(dat, params, responses, infType, nboot = 500) {
  lst <- map(responses, ~ .compute_prcc(dat, params, .x, nboot = nboot,
                                        write_stub = paste0("transSensitivityPRCC_", make.names(.x)),
                                        infType = infType))
  names(lst) <- responses
  
  # Map pretty labels immediately (assumes same length/order as params)
  # Build a lookup from original param order to pretty names
  param_lookup <- tibble(parameters = params, pretty = parameter_names[seq_along(params)])
  lst <- map(lst, ~ .x %>% left_join(param_lookup, by = "parameters") %>%
               mutate(parameters = pretty) %>%
               select(-pretty))
  
  # Wide by parameter
  wide <- lst[[1]] %>% select(parameters, TruPrev = est)
  wide <- wide %>% left_join(lst[['fadeout']]        %>% select(parameters, Fadeout = est), by = "parameters")
  wide <- wide %>% left_join(lst[['fadeout time']]   %>% select(parameters, FadeoutTime = est), by = "parameters")
  wide <- wide %>% left_join(lst[['Hunt Prevalence']]%>% select(parameters, HuntPrev = est), by = "parameters")
  
  # Attach CIs (optional use in plotting/filtering)
  ci_long <- bind_rows(
    lst[[1]]              %>% mutate(metric = 'TruPrev'),
    lst[['fadeout']]      %>% mutate(metric = 'Fadeout'),
    lst[['fadeout time']] %>% mutate(metric = 'FadeoutTime'),
    lst[['Hunt Prevalence']] %>% mutate(metric = 'HuntPrev')
  )
  
  # Order by TruPrev (NA last)
  wide <- wide %>% arrange(is.na(TruPrev), TruPrev)
  
  long <- wide %>%
    pivot_longer(cols = c(TruPrev, Fadeout, FadeoutTime, HuntPrev),
                 names_to = "metric", values_to = "value")
  
  long$parameters2 <- factor(long$parameters, levels = wide$parameters)
  long$metric <- factor(long$metric,
                        levels = c("TruPrev","Fadeout","FadeoutTime","HuntPrev"),
                        labels = c("True Prevalence","Fadeout","Fadeout Time","Observed Prevalence"))
  
  list(wide = wide, long = long, ci_long = ci_long)
}

# ============================================================
# 2) Regression utilities (lm), scaled coefs + p-values
# ============================================================

.fit_lm_model <- function(data, params, resp_col, write_stub = NULL, infType = NULL) {
  form <- as.formula(paste(resp_col, "~", paste(params, collapse = " + ")))
  model <- lm(form, data = data)
  td <- tidy(model)
  if (!is.null(write_stub) && !is.null(infType)) {
    write.csv(td, paste0(write_stub, "_", infType, ".csv"), row.names = FALSE)
  }
  # Coefs without intercept
  coefs <- coef(model)[-1]
  coefs_scaled <- coefs / max(abs(coefs), na.rm = TRUE)
  tibble(
    parameters = params,
    scaled_est = as.numeric(coefs_scaled),
    pval = summary(model)$coefficients[-1, 4]
  )
}

.run_all_lm <- function(data, params, responses, prcc_order_params, infType) {
  lm_long <- map_dfr(responses, function(resp) {
    .fit_lm_model(data, params, resp,
                  write_stub = paste0("lmOutput_", make.names(resp)),
                  infType = infType) %>%
      mutate(response = resp)
  })
  
  # pretty label map
  param_lookup <- tibble(parameters = params, pretty = parameter_names[seq_along(params)])
  lm_long <- lm_long %>%
    left_join(param_lookup, by = "parameters") %>%
    mutate(parameters = pretty) %>%
    select(-pretty)
  
  # Order parameters by PRCC TruPrev ordering
  lm_long$parameters <- factor(lm_long$parameters, levels = prcc_order_params)
  lm_long$response <- factor(lm_long$response,
                             levels = responses,
                             labels = c("True Prevalence","Fadeout","Fadeout Time","Observed Prevalence"))
  lm_long
}

# ============================================================
# 3) Plotters
# ============================================================

.plot_prcc <- function(prcc_long, ci_long = NULL, mask_cross0 = TRUE) {
  # If masking by CI cross 0, set NA where CI spans 0
  data_plot <- prcc_long
  if (mask_cross0 && !is.null(ci_long)) {
    # prepare a CI table aligned to long
    ci_tbl <- ci_long %>%
      mutate(metric = recode(metric,
                             TruPrev = "True Prevalence",
                             Fadeout = "Fadeout",
                             FadeoutTime = "Fadeout Time",
                             HuntPrev = "Observed Prevalence")) %>%
      select(parameters, metric, min_ci, max_ci)
    data_plot <- data_plot %>%
      left_join(ci_tbl, by = c("parameters" = "parameters", "metric" = "metric")) %>%
      mutate(value = ifelse(!is.na(min_ci) & !is.na(max_ci) &
                              ((min_ci <= 0 & max_ci >= 0)), NA, value))
  }
  
  ggplot(data_plot, aes(x = parameters2, y = value, fill = I(blue))) +
    geom_col(na.rm = TRUE) +
    coord_flip() +
    facet_grid(~ metric) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      strip.text = element_text(size = 12),
      panel.spacing.y = unit(1.5, "lines")
    )
}

.plot_lm <- function(lm_long, mask_p = 0.05) {
  data_plot <- lm_long %>%
    mutate(scaled_est = ifelse(!is.na(pval) & pval >= mask_p, NA, scaled_est))
  
  ggplot(data_plot, aes(x = parameters, y = scaled_est, fill = I(blue))) +
    geom_col(na.rm = TRUE) +
    coord_flip() +
    facet_grid(~ response) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      strip.text = element_text(size = 12),
      panel.spacing.y = unit(1.5, "lines")
    )
}

# ============================================================
# 4) Run everything
# ============================================================

# PRCC
prcc_out <- .run_all_prcc(bTB_wl_scaled, parameters, responses, infType, nboot = 0)
# Regression (ordered by PRCC TruPrev)
lm_long <- .run_all_lm(bTB_wl_scaled, parameters, responses,
                       prcc_order_params = prcc_out$wide$parameters, infType = infType)

# Plots
sens_plot <- .plot_prcc(prcc_out$long, prcc_out$ci_long, mask_cross0 = TRUE)
sens_plot <- annotate_figure(sens_plot, top = text_grob("PRCC Analysis", face = "bold", size = 18))

reg_plot  <- .plot_lm(lm_long, mask_p = 0.05)
reg_plot  <- annotate_figure(reg_plot,  top = text_grob("Regression Analysis", face = "bold", size = 18))

# Combined (PRCC on top, Regression on bottom)
combined <- ggarrange(
  sens_plot, reg_plot,
  ncol = 1, nrow = 2,
  heights = c(1, 1.1)  # regression slightly taller
)

# ============================================================
# 5) Save figures
# ============================================================
dir.create("figures", showWarnings = FALSE)

jpeg(filename = paste0("figures/PRCC_Trans_", infType, ".jpeg"),
     width = 1440, height = 840, units = 'px', res = 100)
print(sens_plot)
dev.off()

jpeg(filename = paste0("figures/SingleModel_Trans_", infType, ".jpeg"),
     width = 1440, height = 840, units = 'px', res = 100)
print(reg_plot)
dev.off()

jpeg(filename = paste0("figures/Combined_PRCC_Reg_", infType, ".jpeg"),
     width = 1440, height = 1200, units = 'px', res = 110)
print(combined)
dev.off()





















# --- Setup ---
library(sensitivity)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(grid)

# --- Helper: PRCC with bootstrap CIs (no p-values), no filtering here ---
.compute_prcc <- function(dat, params, resp_col, file_stub, infType, nboot = 100) {
  # Build X, y
  Xdf <- dat[, params, drop = FALSE]
  X   <- data.matrix(Xdf)                 # coerce to numeric matrix
  y   <- dat[[resp_col]]
  if (is.data.frame(y)) y <- y[[1]]
  if (is.list(y))       y <- unlist(y, recursive = FALSE, use.names = FALSE)
  y <- as.numeric(y)
  
  # Drop NA/Inf rows
  ok <- is.finite(y) & apply(X, 1, function(r) all(is.finite(r)))
  X <- X[ok, , drop = FALSE]
  y <- y[ok]
  
  # Sanitize predictor names for formula safety inside pcc()
  orig_names <- colnames(X)                    # original param codes
  safe_names <- make.names(orig_names, unique = TRUE)
  colnames(X) <- safe_names
  
  # Run PRCC with bootstrap CIs
  fit <- pcc(X = as.data.frame(X), y = y, rank = TRUE, nboot = nboot)
  
  # Extract & map back to original names
  out <- as.data.frame(fit$PRCC)
  out$parameters_safe <- rownames(out)
  rownames(out) <- NULL
  
  # Map safe -> original codes
  out$parameters <- orig_names[match(out$parameters_safe, safe_names)]
  
  out <- out |>
    dplyr::transmute(
      parameters,                # original param code
      est   = original,
      # min_ci = `min. c.i.`,
      # max_ci = `max. c.i.`
    )
  
  write.csv(out, file = paste0(file_stub, "_", infType, ".csv"), row.names = FALSE)
  out
}

# --- Run PRCCs for each response (unfiltered) ---
prcc.TruPrev_raw     <- .compute_prcc(dat = bTB_wl_scaled, params = parameters, resp_col = response[1], "transSensitivityPRCC_truPrev", infType, nboot = 0)
prcc.Fadeout_raw     <- .compute_prcc(bTB_wl_scaled, parameters, response[2], "transSensitivityPRCC_fadeout", infType, nboot = 0)
prcc.FadeoutTime_raw <- .compute_prcc(bTB_wl_scaled, parameters, response[3], "transSensitivityPRCC_fadeoutTime", infType, nboot = 0)
prcc.HuntPrev_raw    <- .compute_prcc(bTB_wl_scaled, parameters, response[4], "transSensitivityPRCC_obsPrev", infType, nboot = 0)

# --- Rename parameters to pretty labels FIRST ---
parameter_names <- c(
  'carrying.capacity','hunt.mortality','base.mortality','dens.dep.asymetry',
  'dens.dep.mortality','max.birth.rate','proportion.SS','birth.timing',
  'birth.duration','transmission.rate','SS.contact.factor',
  'latency.mean','latency.rate'
)

rename_after <- function(df, pretty_vec) {
  # assumes length(pretty_vec) == number of parameters and original order
  df$parameters <- pretty_vec
  df
}

prcc.TruPrev     <- rename_after(prcc.TruPrev_raw,     parameter_names)
prcc.Fadeout     <- rename_after(prcc.Fadeout_raw,     parameter_names)
prcc.FadeoutTime <- rename_after(prcc.FadeoutTime_raw, parameter_names)
prcc.HuntPrev    <- rename_after(prcc.HuntPrev_raw,    parameter_names)

# --- NOW filter rows whose CI crosses zero (per response, after renaming) ---
drop_cross0 <- function(df) {
  df %>% filter( (min_ci > 0 & max_ci > 0) | (min_ci < 0 & max_ci < 0) )
}

# prcc.TruPrev     <- drop_cross0(prcc.TruPrev)
# prcc.Fadeout     <- drop_cross0(prcc.Fadeout)
# prcc.FadeoutTime <- drop_cross0(prcc.FadeoutTime)
# prcc.HuntPrev    <- drop_cross0(prcc.HuntPrev)

# --- Build a wide table by parameter (left-join so lengths can differ) ---
wide <- prcc.TruPrev %>% select(parameters, TruPrev = est)
wide <- wide %>% left_join(prcc.Fadeout %>% select(parameters, Fadeout = est), by = "parameters")
wide <- wide %>% left_join(prcc.FadeoutTime %>% select(parameters, FadeoutTime = est), by = "parameters")
wide <- wide %>% left_join(prcc.HuntPrev %>% select(parameters, HuntPrev = est), by = "parameters")

# Order by TruPrev estimate (NA placed last)
wide <- wide %>% arrange(is.na(TruPrev), TruPrev)

# Long format for plotting; keep NA rows—they simply won’t draw a bar in that facet
prcc.trans.long <- wide %>%
  pivot_longer(cols = c(TruPrev, Fadeout, FadeoutTime, HuntPrev),
               names_to = "metric", values_to = "value")

# Factor for facet order and y-axis order (based on TruPrev)
prcc.trans.long$parameters2 <- factor(prcc.trans.long$parameters, levels = wide$parameters)
prcc.trans.long$metric <- factor(prcc.trans.long$metric,
                                 levels = c("TruPrev","Fadeout","FadeoutTime","HuntPrev"),
                                 labels = c("True Prevalence","Fadeout","Fadeout Time","Observed Prevalence"))

# --- Plot ---
blue <- c('#2171b5')

sens <- ggplot() +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    strip.text = element_text(size = 12),
    panel.spacing.y = unit(1.5, "lines")
  ) +
  geom_col(data = prcc.trans.long, aes(y = value, x = parameters2, fill = I(blue)), na.rm = TRUE) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
  labs(y = NULL, x = NULL) +
  coord_flip() +
  facet_grid(~ metric)

sens <- annotate_figure(
  sens,
  top = text_grob("PRCC Analysis", color = "black", face = "bold", size = 18)
)

print(sens)

jpeg(filename = paste0("figures/PRCC_Trans_", infType, ".jpeg"),
     width = 1440, height = 840, units = 'px', res = 100)
print(sens)
dev.off()


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
write.csv(tidy(lm.Fadeout.PRCC), file=paste0("lmOutputFadeout_", infType, ".csv"), row.names = F)
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
write.csv(tidy(lm.FadeoutTime.PRCC), file=paste0("lmOutputFadeoutTime_", infType, ".csv"), row.names = F)
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
write.csv(tidy(lm.HuntPrev.PRCC), file=paste0("lmOutputHuntPrev_", infType, ".csv"), row.names = F)
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

order <- wide$parameters

#ordering values for plot
lmTruPrevScale <-  lmTruPrevScale[order(match(lmTruPrevScale$parameters,       order)),]
lmFadeoutTimeScale <-lmFadeoutTimeScale[order(match(lmFadeoutTimeScale$parameters,   order)),]
lmHuntPrevScale <-  lmHuntPrevScale[order(match(lmHuntPrevScale$parameters,       order)),]
lmFadeoutScale <- lmFadeoutScale[order(match(lmFadeoutScale$parameters,    order)),]

p_vals_lm <- p_vals_lm[order(match(p_vals_lm$parameters,     order)),]


lmScale <- cbind(TruPrev=lmTruPrevScale, FadeoutTime=lmFadeoutTimeScale$values.in.PRCC.scaled, HuntPrev=lmHuntPrevScale$values.in.PRCC.scaled, Fadeout=lmFadeoutScale$values.in.PRCC.scaled)
# uncomment the line below if you want lm plots in ascending order -- i.e. probably not identical to PRCC ordering
#lmScale <- lmScale[order(lmScale$TruPrev.values.in.PRCC.scaled),]
lmScale <- lmScale[,c("TruPrev.values.in.PRCC.scaled", "TruPrev.parameters", "Fadeout", "FadeoutTime", "HuntPrev")]
lmScale.long <- melt(lmScale,id=c("TruPrev.parameters"))

p_vals_lm.long <- melt(p_vals_lm,id=c("parameters"))

colnames(lmScale.long)[1] <- "parameters"
lmScale.long$parameters2 <- factor(lmScale.long$parameters,levels= lmScale$TruPrev.parameters)
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

significant=T
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
