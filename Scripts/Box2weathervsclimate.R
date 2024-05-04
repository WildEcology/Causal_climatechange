## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-04-05
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
                 mvtnorm, clubSandwich)


theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select



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
  ungroup() %>% 
  mutate(ppt2=ppt^2, tmax2=tmax^2)


recent10years <- prism_growing %>% 
  filter(year>2010)
recent10years$Years = as.factor(10)

earlygrowing <- prism_growing %>% 
  filter(year>1970&year<1990)
earlygrowing$Years = as.factor(20)
  
earliest <- prism_growing %>% 
  filter(year>1940&year<1946)

earliest$Years = as.factor(5)

newclim <- prism_growing %>% 
  mutate(decade = ifelse(year>2010, "2010era20yrs",
                         ifelse(year>1970&year<200, "1970era30yrs",
                                ifelse(year>1940&year<1950, "1940era10yrs", "other"))))


mergclim <- prism_growing %>%
  mutate(Years = as.factor(100)) %>% 
  full_join(earlygrowing) %>% 
  full_join(earliest) %>% 
  full_join(recent10years) %>% 
  mutate(Years = factor(Years, levels = c(5, 10, 20, 100))) %>% 
  group_by(year, Years) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)

mergclim %>% 
  ggplot(aes(x=tmax, color = Years, fill = Years))+
  geom_density(alpha = .5)+
  scale_fill_manual(values = c( "#58a4b0" , "#373f51" , "#FDAF7B","#824D74"))+
  scale_color_manual(values = c( "#58a4b0" , "#373f51" , "#FDAF7B","#824D74"))+
  xlim(2, 10)+
  xlab("Temperature (°C)")+
  ylab("Density")
  


  

colors = c("#a9bcd0", "#58a4b0" , "#373f51" , "#BE7B72","#824D74", "#FDAF7B")
