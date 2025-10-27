# ------------------------------------------------------------
# Main text figures
# ------------------------------------------------------------
system("g++ -L/usr/lib/x86_64-linux-gnu tb_wildlife_freq_cont.cpp -lgsl -lgslcblas -lm -o wl_model_CTMC.exe")
system("g++ -L/usr/lib/x86_64-linux-gnu tb_wildlife_freq_disc.cpp -lgsl -lgslcblas -lm -o wl_model_DTMC.exe")
## -----------------------------------------------------------
## Figure 1: size x infType

param_grid <- expand_grid(
  size = c(250),
  infType = c("seeded", "spillover")
)

for (i in 1:nrow(param_grid)) {
  row <- param_grid[i,]
  
  run_wl_sim(size = row$size, infType = row$infType, save_plots = TRUE, freq = T) 
}

## -----------------------------------------------------------
## Figure 2: super-spreader prop x herd size
results <- read.csv("data/ss_herdsize_meanprev_freq.csv")

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
                       breaks = round(seq(0,max(results$mean_prev), length.out = 4),4),   # manually specify desired breaks
                       limits = c(0, round(max(results$mean_prev),3))) +
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

ggsave("figures/supp/mean_prevalence_map_freq.pdf", mean_map, width = 7.5, height = 6, dpi = 300)

# --------------
results <- read.csv("data/xi_phi_meanprev_freq.csv") |>
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
                       breaks = round(seq(0,max(results$fadeout), length.out = 4),2),   # manually specify desired breaks
                       limits = c(0, max(results$fadeout))) +
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

ggsave("figures/supp/strong_rare_fade_map_freq.pdf", fade_prob_map, width = 7.5, height = 6, dpi = 300)

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
                       breaks = round(seq(0,max(results$mean_prev), length.out = 4),4),   # manually specify desired breaks
                       limits = c(0, max(results$mean_prev))) +
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

ggsave("figures/supp/strong_rare_ss_map_freq.pdf", mean_sv_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## Figure 3: fadeout. herd size x number initially exposed 

results <- read.csv("data/herdsize_initnum_fadeout_freq.csv") |>
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

ggsave("figures/supp/fade_prob_map_num_freq.pdf", fade_prob_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## Figure 5 (a and b): herd size x harvest intensity
## -----------------------------------------------------------

results <- read.csv("data/hunt_herdsize_freq.csv")|>
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

ggsave("figures/supp/fade_prob_map_hunt_freq.pdf", fade_prob_map_hunt, width = 7.5, height = 6, dpi = 300)

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
    limits = round(c(0, max(results$prev_ratio, na.rm = TRUE)),5),
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

ggsave("figures/supp/prev_est_map_hunt_freq.pdf", mean_map, width = 7.5, height = 6, dpi = 300)

## -----------------------------------------------------------
## % reduction in beta  x  % reduction in farm-contact (p2’s)
## 
## -----------------------------------------------------------

results <- read.csv("data/beta_p2_hunt_freq.csv") 

p1_red_grid    <- seq(0, 1, length.out = 21)   # % reduction of p1 (0..1)
farm_red_grid  <- seq(0, 1, length.out = 21)   # % reduction of farm-contact (0..1)
hunt_levels    <- c(0.2, 0.5, 0.95)               # fixed columns
scenarios      <- c("spillover")

results <- results %>%
  mutate(
    fadeout        = ifelse(fadeout == 0, NA_real_, fadeout),
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
    y = expression("% reduction in cattle/fomite contacts ("*p[q]*")")
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
    breaks = round(seq(0, max(results$mean_prev, na.rm = TRUE), length.out = 5), 5),
    limits = c(0, round(max(results$mean_prev, na.rm = TRUE),5)),
    name   = "mean infection\nprevalence"
  ) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,100,25)) +
  labs(
    x = expression("% reduction in deer–deer transmission ("*beta[wild]*")"), # "% reduction in deer–deer contact (p1)",
    y = expression("% reduction in cattle/fomite contacts ("*p[q]*")")
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

ggsave("figures/supp/fade_prob_map_beta_p2_hunt_freq.pdf", fade_prob_map, width = 12, height = 5, dpi = 300)
ggsave("figures/supp/mean_prev_map_beta_p2_hunt_freq.pdf", mean_prev_map, width = 12, height = 5, dpi = 300)

# ---------------------------
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
bTB_wl_scaled <- read.csv("data/LHS_scaled_summary_seeded_q1_freq.csv")

averaged = T
mono_plots = T
bTB_wl_scaled <- bTB_wl_scaled[,!(names(bTB_wl_scaled) %in% c("p1","ksi","phi",'N', "pars", "S_0", "E1_0", "SuperS_0", "SuperE1_0", "I_0", "SuperI_0"))]
names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max',  'omega', 's', 
                          'beta',  'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4',  'sigma1_mean', 'sigma1_rate', 
                          'start_q')
reps = 500
if(averaged){
  lh <- read.csv("data/bTBwl_LHSpars_seeded_q1_freq_freq.csv")
  infType = 'seeded_q1'
  bTB_wl_agg <- aggregate.data.frame(x = bTB_wl_scaled, by = list(rep(1:dim(lh)[1], times=1, each=reps)), FUN = mean)
  bTB_wl_agg <- bTB_wl_agg[,names(bTB_wl_agg) %in% names(bTB_wl_scaled)]
  bTB_wl_scaled <- bTB_wl_agg
  infType <- paste0(infType, '_avg')
  write.csv(bTB_wl_scaled, file = paste0("./data/LHS_scaled_summary_", infType, "_freq.csv" ), row.names = F)
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
    outfile <- paste0("monotonic_", param, "_", infType, "_freq.jpeg")
    jpeg(outfile, width = 840, height = 840, units = "px", res = 100)
    print(ggarrange(plotlist = plots, nrow = 2, ncol = 2))
    dev.off()
  }
}

if (exists("mono_plots") && isTRUE(mono_plots)) {
  walk(parameters, plot_monotonicity, infType = infType)
}

#####################
setwd("../..")
names(bTB_wl_scaled) <- c('Total Infected', 'fadeout', 'fadeout time', 'Hunt Prevalence', 
                          'K', 'eta_hunt', 'eta_nat', 'theta', 'gamma', 
                          'alpha_max', 'omega', 's', 
                          'beta', 'p2_q1', 'p2_q2', 'p2_q3', 'p2_q4', 'sigma1_mean', 'sigma1_rate', 
                          'start_q')

# --------------
# PRCC
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
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

sens

ggsave("figures/supp/sens_plot_prcc_freq.pdf", sens, width = 14, height = 6, dpi = 300)

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

ggsave("figures/supp/sens_plot_reg_freq.pdf", reg, width = 14, height = 6, dpi = 300)


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

ggsave("figures/supp/sens_plot_reg_int_freq.pdf", reg, width = 17, height = 10, dpi = 300)