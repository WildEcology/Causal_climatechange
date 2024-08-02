## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-07-21
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------



librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, 
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwich, rasterVis, broom.mixed, scales,RColorBrewer)

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


rwidat <- read_csv("Data/Data_processed/whitebark/PANEL_data_rwi.csv")
demdat <- read_csv("Data/Data_processed/whitebark/plot_demdat.csv")

paneldat <- rwidat %>% 
  left_join(demdat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Panel models
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe_mod = feols(value ~ tmax + ppt | plot_id_needle,
                data= paneldat)


fe_mod_lag = feols(value ~ tmax +  laggedtmax + ppt +  laggedprecip | plot_id_needle,
                    data = paneldat)

summary(fe_mod_lag)


fe_mod_lag5 = feols(value ~ lag(tmax, 5) + lag(ppt, 5)| plot_id_needle,
                    data= paneldat)
summary(fe_mod_lag5)


fe_lm_lag = feols(value ~ tmax +  laggedtmax + ppt +  laggedprecip + factor(plot_id_needle),
                   data = paneldat)

summary(fe_lm_lag)

fe_mod_lag_RS = feols(value ~ tmax:plot_id_needle +  laggedtmax:plot_id_needle + ppt:plot_id_needle +  laggedprecip:plot_id_needle | plot_id_needle,
                   data = paneldat)

summary(fe_mod_lag_RS)

fe_RS <- lm(value ~ (tmax +  ppt + laggedtmax + laggedprecip) * plot_id_needle,
                          data = paneldat)
summary(fe_RS)


tab_model(fe_mod, fe_mod_lag, fe_mod_lag_RS,show.aic = T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LMM models
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmm_mod = lmer(value ~ tmax +  ppt +  (1|plot_id_needle), data = paneldat)

lmm_mod_lag = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + (1|plot_id_needle),
                   data = paneldat)

lmm_mod_lag_Rslopes = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + 
                             (1 + tmax | plot_id_needle)+ (1 + laggedtmax | plot_id_needle) + 
                             (1 + ppt | plot_id_needle)+ (1 + laggedprecip | plot_id_needle),
                   data = paneldat)

#lm_mod_elev = lmer(value ~ tmax +  ppt + elevation + (1|plot_id_needle), data = paneldat)
lmm_mod_elev = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + elevation + (1|plot_id_needle), data = paneldat)

#lmm_mod_slope = lmer(value ~ tmax +  ppt + slope + (1|plot_id_needle), data = paneldat)
lmm_mod_slope = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + slope + (1|plot_id_needle), data = paneldat)

#lmm_mod_elevslope = lmer(value ~ tmax +  ppt + slope + elevation+ (1|plot_id_needle), data = paneldat)
lmm_mod_elevslope = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + slope + elevation+ aspect+(1|plot_id_needle), data = paneldat)

# lmm_mod_elevslope_Rslopes = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + slope + elevation+ aspect +
#                                    (1 + tmax | plot_id_needle)+ (1 + laggedtmax | plot_id_needle) + 
#                                    (1 + ppt | plot_id_needle)+ (1 + laggedprecip | plot_id_needle), data = paneldat)

lmm_mod_elevslope_Rslopes = lmer(value ~ tmax +  laggedtmax + ppt +  laggedprecip + slope + elevation+ aspect+
                                   (1 + tmax | plot_id_needle)+ (1 + laggedtmax | plot_id_needle) + 
                                   (1 + ppt | plot_id_needle)+ (1 + laggedprecip | plot_id_needle), data = paneldat)

#lmm_mod_elevslope_Rslopes = lmer(value ~ tmax*plot_id_needle +  laggedtmax*plot_id_needle + ppt*plot_id_needle +  laggedprecip*plot_id_needle+ slope + elevation+ aspect+(1|plot_id_needle), data = paneldat)
# lm_mod_lag_square = lmer(value ~ tmax + I(tmax^2) + laggedtmax+ I(laggedtmax^2) + 
#                             ppt + I(ppt^2) + laggedprecip + 
#                             I(laggedprecip^2)+  (1|plot_id_needle), data = paneldat)


tab_model(fe_mod, fe_mod_lag, lmm_mod, lmm_mod_lag, lmm_mod_lag_Rslopes, lmm_mod_elevslope, lmm_mod_elevslope_Rslopes, show.aic = T, digits = 4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                           Function that extracts estimates and CIs
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Extract fixed effects estimates and confidence intervals
# model_list <- list(fe_mod=fe_mod, fe_mod_square = fe_mod_square, fe_mod_lag = fe_mod_lag, 
#                fe_mod_lag_square = fe_mod_lag_square,
#                lmm_mod = lmm_mod, lmm_mod_square = lmm_mod_square, 
#                lmm_mod_lag = lmm_mod_lag, lmm_mod_lag_square = lmm_mod_lag_square)

# model_list <- list(fe_mod=fe_mod, fe_mod_lag = fe_mod_lag,fe_lag_RS = fe_RS,
#                lmm_mod = lmm_mod, lmm_mod_lag = lmm_mod_lag,
#                lmm_mod_elev = lmm_mod_elev, lmm_mod_slope=lmm_mod_slope, lmm_mod_elevslope = lmm_mod_elevslope,
#                lmm_mod_elevslope_Rslopes = lmm_mod_elevslope_Rslopes, lmm_mod_lag_Rslopes = lmm_mod_lag_Rslopes)

model_list <- list(fe=fe_mod, fe_lag = fe_mod_lag,fe_lm_lag = fe_lm_lag, fe_lag_RS = fe_RS,
                   lmm= lmm_mod, lmm_lag = lmm_mod_lag,lmm_lag_Rslopes = lmm_mod_lag_Rslopes,
                   lmm_elev = lmm_mod_elev, lmm_slope=lmm_mod_slope, lmm_ElevSlope = lmm_mod_elevslope,
                   lmm_lag_ES_Rslopes = lmm_mod_elevslope_Rslopes)

extract_model_info <- function(model_list) {
  
  # Initialize an empty dataframe to store results
  results <- data.frame(
    model = character(),
    term = character(),
    estimate = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each model in the list
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    model_summary <- tidy(model, conf.int = TRUE)
    
    # Add the model name to the dataframe
    model_summary$model <- model_name
    
    # Select relevant columns and bind to results
    results <- rbind(results, model_summary[, c("model", "term", "estimate", "conf.low", "conf.high")])
  }
  
  return(results)
}

model_info_df <- extract_model_info(model_list)


modeldat <- model_info_df %>% 
  na.omit() %>% 
  filter(term!="(Intercept)") %>% 
  mutate(modeltype = ifelse(str_detect(model, "lmm"), "LMM", "FE Model"))
  

# # Plot using ggplot2
# modeldat %>% 
#   filter(term == "tmax") %>% 
#   ggplot(aes(y = estimate, x = model, ymin = conf.low, ymax = conf.high, color = modeltype)) +
#     geom_pointrange() +
#     labs(
#        y = "Temperature coefficient",
#        x = "Model")
#   #ylim(0, .025)+
#   #geom_hline(yintercept = 0, linetype="dashed")
# 
# 
# modeldat %>% 
#   filter(term == "laggedtmax") %>% 
#   ggplot(aes(y = estimate, x = model, ymin = conf.low, ymax = conf.high, color = modeltype)) +
#   geom_hline(yintercept = 0, linetype="dashed")+
#   geom_pointrange() +
#   labs(
#        y = "Lagged Tmax",
#        x = "Estimate")
#   #theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust = .35))
#   
# modeldat %>% 
#   filter(term == "tmax") %>% 
#   ggplot(aes(y = estimate, x = model, ymin = conf.low, ymax = conf.high, color = modeltype)) +
#   #geom_hline(yintercept = 0, linetype="dashed")+
#   geom_pointrange() +
#   labs(
#     y = "Tmax",
#     x = "Estimate")
#   #theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust = .35))
#   #ylim(0, .06)
# 
# 
# modeldat %>% 
#   filter(term == "ppt") %>% 
#   ggplot(aes(y = estimate, x = model, ymin = conf.low, ymax = conf.high, color = modeltype)) +
#   #geom_hline(yintercept = 0, linetype="dashed")+
#   geom_pointrange() +
#   labs(
#     y = "Precip coefficient",
#     x = "Model")
# 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                  MODEL COMPARISON
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


compare <- modeldat %>% 
  filter(term %in%c("tmax", "laggedtmax")) %>%
  filter(model%in%c("fe_lag" ,  "lmm_lag")) %>% 
  mutate(term_rename = ifelse(term == "tmax", "Temperature", "Temperature lagged")) %>% 
  mutate(Model= ifelse(model == "fe_lag", "FE panel model", "LMM")) %>% 
  ggplot(aes(y = estimate, x = term_rename, ymin = conf.low, ymax = conf.high, color = Model)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_pointrange(position = position_dodge(width = .5)) +
  scale_color_manual(values = c("cyan4","orchid4"))+
  labs(
    y = "Coefficient",
    x = "Explanatory variable") +
  theme(
    axis.title.x = element_text(vjust = -1),  # Increase spacing for the x-axis label
    axis.text.x = element_text(margin = margin(t = 5))  # Increase spacing for the x-axis text
  )
compare

compareppt <- modeldat %>% 
  filter(term %in%c("ppt", "laggedprecip")) %>%
  filter(model%in%c("fe_lag" ,  "lmm_lag")) %>% 
  mutate(term_rename = ifelse(term == "ppt", "Precipitation",
                              ifelse(term == "laggedprecip", "Precipitation lagged", NA))) %>% 
  mutate(Model= ifelse(model == "fe_lag", "FE panel model", "LMM")) %>% 
  ggplot(aes(y = estimate, x = term_rename, ymin = conf.low, ymax = conf.high, color = Model)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_pointrange(position = position_dodge(width = .5)) +
  scale_color_manual(values = c("cyan4","orchid4"))+
  labs(
    y = "Coefficient",
    x = "Explanatory variable") +
  theme(
    axis.title.x = element_text(vjust = -1),  # Increase spacing for the x-axis label
    axis.text.x = element_text(margin = margin(t = 5))  # Increase spacing for the x-axis text
  )


compare + compareppt +  plot_layout(guides = "collect") & plot_annotation(tag_levels = "A")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   ROBUSTNESS CHECKS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


robusttmax <- modeldat %>% 
  filter(term %in%c("tmax", "laggedtmax")) %>%
  #filter(model%in%c("fe_mod_lag" ,  "glm_mod_lag")) %>% 
  mutate(term_rename = ifelse(term == "tmax", "Temperature", "Temperature lagged")) %>% 
  #mutate(Model= ifelse(model == "fe", "FE panel model", "GLMM")) %>% 
  ggplot(aes(y = estimate, x = term_rename, ymin = conf.low, ymax = conf.high, color = model)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_pointrange(position = position_dodge(width = .5)) +
  scale_color_manual(values = c("cyan4","gold1", "orchid4",  "royalblue4", "red4", "skyblue4", "darkgrey", "black", "lightblue", "darkgreen", "orange")) +
  labs(
    y = "Coefficient",
    x = "Explanatory variable") +
  theme(
    axis.title.x = element_text(vjust = -1),  # Increase spacing for the x-axis label
    axis.text.x = element_text(margin = margin(t = 5))  # Increase spacing for the x-axis text
  )
robusttmax






robustppt <- modeldat %>% 
  filter(term %in%c("ppt", "laggedprecip")) %>%
  mutate(term_rename = ifelse(term == "ppt", "Precipitation",
                               ifelse(term == "laggedprecip", "Precipitation lagged", NA))) %>% 
  #mutate(Model= ifelse(model == "fe_mod_lag", "FE panel model", "GLMM")) %>% 
  ggplot(aes(y = estimate, x = term_rename, ymin = conf.low, ymax = conf.high, color = model)) +
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_pointrange(position = position_dodge(width = .5)) +
  scale_color_manual(values = c("cyan4","gold1", "orchid4",  "royalblue4", "red4", "skyblue4", "darkgrey", "black", "lightblue", "darkgreen", "orange")) +
  labs(
    y = "Coefficient",
    x = "Explanatory variable") +
  theme(
    axis.title.x = element_text(vjust = -1),  # Increase spacing for the x-axis label
    axis.text.x = element_text(margin = margin(t = 5))  # Increase spacing for the x-axis text
  )
robustppt

robusttmax + robustppt +  plot_layout(guides = "collect") & plot_annotation(tag_levels = "A")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   comparing elevation correlations
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tmaxelev <- paneldat %>% 
  ggplot(aes(x=elevation, y=tmax, color = "tmax"))+
  geom_point()+
  geom_smooth(se=F)+
  ylab("Temperature")+
  xlab("Elevation (m)")+
  stat_cor(method = "pearson")+
  scale_color_manual (values = "darkred")

tmaxelev

pptelev <- paneldat %>% 
  ggplot(aes(x=elevation, y=ppt, color = "ppt"))+
  geom_point()+
  geom_smooth(se=F)+
  ylab("Precipitation")+
  xlab("Elevation (m)")+
  scale_color_manual (values = "darkblue")+
  stat_cor(method = "pearson")

pptelev


tmaxelev + pptelev


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   RWI and climate figures
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## RWI FIGURE
rwidat <- paneldat %>% 
  group_by(plot_id_needle, year) %>% 
  summarize(rwi = mean(value, na.omit=T))

meanrwi <- paneldat %>% 
  group_by(year) %>% 
  summarize(meanrwi = mean(value, na.omit=T))

combrwi <- rwidat %>% 
  left_join(meanrwi)


# Function to darken colors
darken_color <- function(color, factor = 0.3) {
  rgb_val <- col2rgb(color)
  darkened_rgb <- rgb_val * (1 - factor)
  darkened_rgb <- pmax(pmin(darkened_rgb, 255), 0) # Ensure values are within 0-255
  rgb(darkened_rgb[1], darkened_rgb[2], darkened_rgb[3], maxColorValue = 255)
}

# Generate a blue-green palette with 27 colors from RColorBrewer
base_palette <- brewer.pal(n = 9, name = "Blues")

# Extend the palette to 27 colors by interpolating between the base colors
extended_palette <- colorRampPalette(base_palette)(27)

# Darken the colors using the custom darken_color function
darkened_palette <- sapply(extended_palette, darken_color, factor = 0.1)

theme_set(
  theme_bw(base_size = 17)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

combrwi %>% 
  ggplot(aes(x = year, y=rwi, color=plot_id_needle))+
  geom_point(size=.2)+
  geom_line(linewidth=.1)+
  scale_color_manual(values = extended_palette) +
  geom_line(aes(x=year, y=meanrwi), color = "black", size = 1)+
  #scale_color_manual(values="darkgrey")+
  theme(legend.position = "none")+
  ylab("RWI")+
  xlab("Year")+
  scale_x_continuous(breaks = seq(1901, 2021, by = 10))




