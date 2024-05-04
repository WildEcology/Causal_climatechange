## ---------------------------
##
## Script name: Tree-ring Figures
##
## Author: Dr. Joan Dudney
##
## Date Created: 2022-03-09
##
## Copyright (c) Joan Dudney, 2022
## Email: jdudney@berkeley.edu
##
## ---------------------------
##
## Notes: this code creates figures 3 & 4
## as well as the majority of supplemental figures
##
## ---------------------------

librarian::shelf(sjPlot, ggeffects, patchwork, tidyverse, 
                 lme4, plotrix, ggpubr, mgcv, nlme, fixest, plotrix, egg, ggpmisc,
                 mvtnorm, clubSandwichm, furr)


theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# parallelizing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_cores <- availableCores() - 6
future::plan(multisession, workers = n_cores)

my_seed <- 5597

n_mc <- 10000

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  READING IN AND CLEANING DATA   
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## cleaning tree-ring cores
cores <- read_csv("Data/Data_cleaned/Whitebark/wbp_rwi_v4.csv")
plotnames <- read_csv("Data/Data_cleaned/Whitebark/plot_names_normalized_v2.csv") 
treelook <- read_csv("Data/Data_cleaned/Whitebark/wbp_tree_lookup.csv")[-1]

## taking out duplicated cores
dup.cores <-c("FS12120A","FS1510A","FS1530A","FSD2B04B","FSD2B16B","FSD4B07B",
              "FSD4B16A","FSD5C18A","SE04B12B","SE04B28A","SE18A21A","SE24A11A",
              "SE24A27A","SE4309A","SE4316A","SE4319B","SE4329A","SE67A28B","SE67C12B",
              "YO03A06B","YO03A11B","YO03A13A","YO03A16A","YO03A24B","YO03A30B","YO25A29B")

mm <- match(dup.cores, treelook$COREID_ALT)
dup.cores[!is.na(mm)]<- as.character(treelook$COREID_2[na.omit(mm)])

plots <- plotnames%>%
  select(plot_id_needle, plot_label, plot_label_2)%>%
  rename(plot=plot_label_2)

cores_labels <- treelook %>% 
  select(COREID_2, COREID, plot_label_2, plot_id_needle, Tree_Num, lat,long, DBH_cm, Ht_m, Canopy_Pos) %>%   
  filter(!COREID_2 %in% dup.cores) %>% 
  filter(COREID_2!="site_01_06") %>% ## removing wpbr tree
  rename_all(tolower)

## combining datasets
cores_cleaner <- cores %>% 
  pivot_longer(-"...1", names_to = "COREID_2") %>% 
  rename_all(tolower)%>%
  rename(year="...1") %>% 
  left_join(cores_labels) %>% 
  filter(coreid_2!="site_01_06") 

## prism data
prism <- read_csv("Data/Data_cleaned/Whitebark/prism_plots_1900.csv")

# without-summer water year (Oct 1-June 30)
prism_growing <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpd=mean(vpdmax, na.rm=T), tmax = mean(tmax, na.rm=T),
            ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() 

summer_vals <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month)) %>% 
  #summermonths=ifelse(month%in%c(6:9), "summer")) %>% 
  filter(month%in%c(6:9)) %>%
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, year) %>% 
  summarize(vpdsummer=mean(vpdmax, na.rm=T), tmaxsummer = mean(tmax, na.rm=T),
            pptsummer = sum(ppt, na.rm=T))%>%
  ungroup() 

lagged_growing <- prism_growing %>% 
  filter(year< 2018) %>% 
  mutate(year = year +1) %>% 
  rename(laggedtmax = tmax, laggedprecip = ppt) %>% 
  select(-vpd)

## combining datasets
treenum <- cores_cleaner %>% 
  distinct(plot_id_needle, tree_num) %>% 
  mutate(tree_id=seq(1:771))

# full dataframe
clim_core1 <- cores_cleaner %>%
  #filter(year>1900) %>%
  filter(!is.na(value)) %>% 
  left_join(prism_growing) %>% 
  left_join(lagged_growing) %>%
  left_join(summer_vals) %>% 
  left_join(treenum)

## removing outliers
clim_core <- clim_core1 %>% 
  filter(value<2.202722) %>% 
  filter(value>0.2425516) %>% 
  mutate(tree_id = paste0(tree_num, plot_id_needle))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sensitivity
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datquants <- clim_core %>% 
  mutate(thresh = ifelse(tmax<=3.82, "0-15",
                         ifelse(tmax>3.82 & tmax<=4.76, "16-30",
                                ifelse(tmax>4.76 & tmax <= 5.45, "31-45",
                                       ifelse(tmax > 5.45 & tmax <= 6.07, "31-45",
                                              ifelse(tmax > 6.07 & tmax <= 7.47, "61-85", "86-100"))))))
  

threshcode = c("0-15","16-30","31-45","31-45","61-85", "86-100")

est_sens <- function(threshcode){

  clim_df <- datquants %>%
              filter(thresh == threshcode) %>%
              select(tmax,laggedtmax, ppt, plot_id_needle, tree_id, value, year)

  mod_clim <- lm(value~tmax+ppt + laggedtmax + year + plot_id_needle , data = clim_df)
  mod_dat <- data.frame(thresh = threshcode,
                        tmax = as.numeric(mod_clim$coefficients[2]),
                        ppt = as.numeric(mod_clim$coefficients[3]),
                        laggedtmax = as.numeric(mod_clim$coefficients[4]))

  return(mod_dat)
}


df_mod_res=data.frame()

for (i in threshcode) {
  code = i
  est_coef=est_sens(code)
  df_mod_res=rbind(df_mod_res, est_coef)

}

summary(mod_clim) 

sens_dat_fig = df_mod_res

currentppt <- sens_dat_fig %>% 
  ggplot(aes(x= factor(thresh), y = tmax, group = 1))+
  geom_point()+
  geom_line() +
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Max. temp (°C) quantiles")+
  ylab("Temperature sensitivity")

pptfig <- sens_dat_fig %>% 
  ggplot(aes(x= factor(thresh), y = ppt, group = 1))+
  geom_line(color = "darkgrey")+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Max. temp (°C) quantiles")+
  ylab("Precipitation sensitivity")

laggedtmaxfig <- sens_dat_fig %>% 
  ggplot(aes(x= factor(thresh), y = laggedtmax, group = 1))+
  geom_line(color = "darkgrey")+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Max. temp (°C) quantiles")+
  ylab("Lagged temperature sensitivity")

pptfig + currentppt + laggedtmaxfig + plot_annotation(tag_levels="A") 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# looking at sensitivity across time
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dat_timeperiod <- clim_core %>% 
  mutate(decade_label = paste0((year %/% 10) * 10, "s")) %>% 
  mutate(year5_label = paste0((year - 1) %/% 5 * 5 + 5)) %>% 
  filter(year5_label>1900)



time_sens_summer <- function(label){
  
  clim_data <- dat_timeperiod %>%
    filter(year5_label == label) %>%
    select(tmax,laggedtmax, pptsummer,ppt, plot_id_needle, tree_id, value, year)
  
  mod_clim <- lm(value~ tmax + pptsummer +laggedtmax+ ppt+year + plot_id_needle , data = clim_data)
  mod_dat <- data.frame(timeperiod = label,
                        tmax = as.numeric(mod_clim$coefficients[2]),
                        pptsummer = as.numeric(mod_clim$coefficients[3]),
                        laggedtmax = as.numeric(mod_clim$coefficients[4]),
                        ppt = as.numeric(mod_clim$coefficients[5]))
  
  return(mod_dat)
}

labelnew <- unique(dat_timeperiod$year5_label)
newlabel <- labelnew

df_mod_time=data.frame()

for (i in newlabel) {
  code = i
  est_coef=time_sens_summer(code)
  df_mod_time=rbind(df_mod_time, est_coef)
  
}

time_dat_fig = df_mod_time


timetemp <- time_dat_fig %>% 
  ggplot(aes(x= factor(timeperiod), y = tmax, group = 1))+
  geom_point()+
  geom_line() +
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Timeperiod")+
  ylab("Temperature sensitivity")

ppttime <- time_dat_fig %>% 
  ggplot(aes(x= factor(timeperiod), y = pptsummer, group = 1))+
  geom_line(color = "darkgrey")+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Timeperiod")+
  ylab("Precipitation sensitivity")

laggedtime <- time_dat_fig %>% 
  ggplot(aes(x= factor(timeperiod), y = laggedtmax, group = 1))+
  geom_line(color = "darkgrey")+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  xlab("Max. temp (°C)")+
  ylab("Lagged temperature sensitivity")

timetemp+ppttime+laggedtime



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sensitivity across years and quantiles
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# dat_timeperiod <- clim_core %>% 
#   mutate(decade_label = paste0((year %/% 10) * 10, "s")) %>% 
#   mutate(year5_label = paste0((year - 1) %/% 5 * 5 + 5)) %>% 
#   filter(year5_label>1900) %>% 
#   mutate(thresh = ifelse(tmax<=3.82, "0-15",
#                          ifelse(tmax>3.82 & tmax<=4.76, "16-30",
#                                 ifelse(tmax>4.76 & tmax <= 5.45, "31-45",
#                                        ifelse(tmax > 5.45 & tmax <= 6.07, "31-45",
#                                               ifelse(tmax > 6.07 & tmax <= 7.47, "61-85", "86-100")))))) %>% 
#   mutate(tmaxfifths = ifelse(tmax<=4.19, "20th",
#                              ifelse(tmax>4.19 & tmax<=5.24, "21-40th",
#                                     ifelse(tmax>5.24 & tmax <= 6.07, "41-60",
#                                            ifelse(tmax > 6.07 & tmax <= 7.12, "61-80", ">80th")))))
# 
# 
# threshcode = unique(dat_timeperiod$tmaxfifths)
# newlabel <- unique(dat_timeperiod$decade_label)
# 
# #time_thresh <- function(label, threshcode) {
# df_mod_time=data.frame()
# 
# for (label in newlabel) {
#   for (thresh in threshcode) {
#   
#     filtered_data <- tryCatch({
#       filter_result <- dat_timeperiod %>%
#         filter(decade_label==label) %>% 
#         #filter(tmaxfifths == thresh)
#       filter_result  # Return the filtered data
#     }, error = function(e) {
#       NULL  # If an error occurs, return NULL
#     })
#     
#     # Check if filtered_data is NULL or empty
#     if (is.null(filtered_data) || nrow(filtered_data) == 0) {
#       # If filtered_data is NULL or empty, skip to the next iteration
#       next
#     }
#     
#     # Fit a simple linear model with the filtered dataset
#     mod_clim <- lm(value~ tmax + pptsummer + year + plot_id_needle , data = filtered_data)
#     mod_dat <- data.frame(timeperiod = label,
#                           thresh = threshcode,
#                           tmax = as.numeric(mod_clim$coefficients[2]),
#                           ppt = as.numeric(mod_clim$coefficients[3]))
#     
#     return(mod_dat)
#     df_mod_time=rbind(df_mod_time, mod_dat)
#   }
# } 
#   
# 
# df_mod_time
# 
# for (label in newlabel) {
#   for (thresh in threshcode) {
#   
#     est_coef=time_thresh(label, thresh)
#     df_mod_time=rbind(df_mod_time, est_coef)
#   
#   }
# }
# 
# time_dat_fig = df_mod_time
# 
# filtered_data <- tryCatch({
#   filter_result <- larger_dataset %>%
#     filter(.data[[variable]] == some_condition)
#   filter_result  # Return the filtered data
# }, error = function(e) {
#   NULL  # If an error occurs, return NULL
# })
# 
# # Check if filtered_data is NULL
# if (is.null(filtered_data)) {
#   # If filtered_data is NULL, skip to the next iteration
#   next
# }
# 
# timetemp <- time_dat_fig %>% 
#   ggplot(aes(x= factor(timeperiod), y = tmax, group = 1))+
#   geom_point()+
#   geom_line() +
#   geom_hline(yintercept = 0, linetype="dashed", color="grey")+
#   xlab("Timeperiod")+
#   ylab("Temperature sensitivity")
# 
# ppttime <- time_dat_fig %>% 
#   ggplot(aes(x= factor(timeperiod), y = ppt, group = 1))+
#   geom_line(color = "darkgrey")+
#   geom_point()+
#   geom_hline(yintercept = 0, linetype="dashed", color="grey")+
#   xlab("Timeperiod")+
#   ylab("Precipitation sensitivity")
# 
# 
# # Creating a larger dataset with random data for demonstration purposes
# set.seed(123)  # for reproducibility
# 
# # Define the number of rows and columns for the dataset
# num_rows <- 100
# num_cols <- 5
# 
# # Generate random data for the dataset
# larger_dataset <- data.frame(matrix(rnorm(num_rows * num_cols), ncol = num_cols))
# 
# # Assign meaningful column names
# colnames(larger_dataset) <- c("var1", "var2", "var3", "var4", "var5")
# 
# # Adding a response variable for linear model
# larger_dataset$response_variable <- rnorm(num_rows)
# 
# 
# # Define the variables you want to iterate through
# variables_to_iterate <- colnames(larger_dataset)
# 
# # Iterate through the variables
# for (variable in variables_to_iterate) {
#   # Attempt to filter the dataset
#   filtered_data <- tryCatch({
#     filter_result <- larger_dataset %>%
#       filter(.data[[variable]] == some_condition)
#     filter_result  # Return the filtered data
#   }, error = function(e) {
#     NULL  # If an error occurs, return NULL
#   })
#   
#   # Check if filtered_data is NULL
#   if (is.null(filtered_data)) {
#     # If filtered_data is NULL, skip to the next iteration
#     next
#   }
#   
#   # If filtered_data is not NULL, continue with further processing
#   # Fit a simple linear model with the filtered dataset
#   model <- lm(response_variable ~ ., data = filtered_data)
#   
#   # Print the summary of the model
#   print(summary(model))
# }
# 
# 
# # Set seed for reproducibility
# set.seed(123)
# 
# # Create a dataset with multiple variables
# num_rows <- 100
# num_vars <- 5
# 
# larger_dataset <- data.frame(matrix(rnorm(num_rows * num_vars), ncol = num_vars))
# colnames(larger_dataset) <- paste0("var", 1:num_vars)
# 
# # Adding a response variable for linear model
# larger_dataset$response_variable <- rnorm(num_rows)
# 
# # Define the variables you want to iterate through
# variables_to_iterate <- colnames(larger_dataset)[1:num_vars]
# 
# # Iterate through the variables
# for (variable in variables_to_iterate) {
#   # Attempt to filter the dataset
#   filtered_data <- tryCatch({
#     filter_result <- larger_dataset %>%
#       filter(.data[[variable]] > 0)  # Example filter condition (change as needed)
#     filter_result  # Return the filtered data
#   }, error = function(e) {
#     NULL  # If an error occurs, return NULL
#   })
#   
#   # Check if filtered_data is NULL or empty
#   if (is.null(filtered_data) || nrow(filtered_data) == 0) {
#     # If filtered_data is NULL or empty, skip to the next iteration
#     next
#   }
#   
#   # If filtered_data is not NULL or empty, continue with further processing
#   # Fit a simple linear model with the filtered dataset
#   model <- lm(response_variable ~ ., data = filtered_data)
#   
#   # Print the summary of the model
#   cat("Linear Model Summary for Filtered Data (Variable:", variable, ")\n")
#   print(summary(model))
#   cat("\n")
# }
# 
# arger_dataset$response_variable <- rnorm(num_rows)
# 
# # Define the variables you want to iterate through
# variables_to_iterate <- colnames(larger_dataset)[1:num_vars]
# 
# # Iterate through the variables
# for (variable1 in variables_to_iterate) {
#   for (variable2 in variables_to_iterate) {
#     # Skip if both variables are the same
#     
#     # Attempt to filter the dataset
#     filtered_data <- tryCatch({
#       filter_result <- larger_dataset %>%
#         filter(.data[[variable1]] > 0 & .data[[variable2]] < 0)  # Example filter condition (change as needed)
#       filter_result  # Return the filtered data
#     }, error = function(e) {
#       NULL  # If an error occurs, return NULL
#     })
#     
#     # Check if filtered_data is NULL or empty
#     if (is.null(filtered_data) || nrow(filtered_data) == 0) {
#       # If filtered_data is NULL or empty, skip to the next iteration
#       next
#     }
#     
#     # If filtered_data is not NULL or empty, continue with further processing
#     # Fit a simple linear model with the filtered dataset
#     model <- lm(response_variable ~ ., data = filtered_data)
#     
#     # Print the summary of the model
#     cat("Linear Model Summary for Filtered Data (Variables:", variable1, "&", variable2, ")\n")
#     print(summary(model))
#     cat("\n")
#   }
# }
