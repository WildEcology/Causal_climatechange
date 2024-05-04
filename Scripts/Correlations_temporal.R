# ---------------------------
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
                 mvtnorm, clubSandwich)


theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select

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

lagged_growing <- prism_growing %>% 
  filter(year< 2018) %>% 
  mutate(year = year +1) %>% 
  rename(laggedtmax = tmax, laggedprecip = ppt) %>% 
  select(-vpd)

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

fall_vals <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month)) %>% 
  #summermonths=ifelse(month%in%c(6:9), "summer")) %>% 
  filter(month%in%c(8:9)) %>%
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, year) %>% 
  summarize(vpdfall=mean(vpdmax, na.rm=T), tmaxfall = mean(tmax, na.rm=T),
            pptfall = sum(ppt, na.rm=T))%>%
  ungroup()

winter_vals <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month)) %>% 
  #summermonths=ifelse(month%in%c(6:9), "summer")) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  filter(!month%in%c(6:10)) %>% 
  summarize(vpdwint=mean(vpdmax, na.rm=T), tmaxwint = mean(tmax, na.rm=T),
            pptwint = sum(ppt, na.rm=T))%>%
  ungroup() %>% 
  rename(year=growing)


## combining datasets
treenum <- cores_cleaner %>% 
  distinct(plot_id_needle, tree_num) %>% 
  mutate(tree_id=seq(1:771))

# full dataframe
clim_core1 <- cores_cleaner %>%
  filter(year>1899) %>%
  filter(!is.na(value)) %>% 
  left_join(prism_growing) %>% 
  left_join(summer_vals) %>% 
  left_join(winter_vals) %>% 
  left_join(treenum) %>% 
  left_join(lagged_growing) %>% 
  left_join(fall_vals) 

## removing outliers
clim_core <- clim_core1 %>% 
  filter(value<2.202722) %>% 
  filter(value>0.2425516) %>% 
  mutate(tree_id = paste0(tree_num, plot_id_needle))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# correlations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

corrdat <- clim_core %>% 
  # mutate(tmax15 = ifelse(tmax<=3.82, "15th",
  #                        ifelse(tmax>3.82 & tmax<=4.76, "16-30th",
  #                               ifelse(tmax>4.76 & tmax <= 5.45, "31-45",
  #                                      ifelse(tmax > 5.45 & tmax <= 6.07, "46-60",
  #                                             ifelse(tmax > 6.07 & tmax <= 7.47, "61-85th", ">85th")))))) %>%
  # mutate(tmaxquants = ifelse(tmax<=4.90, "33th",
  #                            ifelse(tmax>4.90 & tmax<=6.35, "34-66th",
  #                                   ifelse(tmax>6.35, ">67", "NA")))) %>% 
  # mutate(tmaxfifths = ifelse(tmax<=4.19, "20th",
  #                            ifelse(tmax>4.19 & tmax<=5.24, "21-40th",
  #                                   ifelse(tmax>5.24 & tmax <= 6.07, "41-60",
  #                                          ifelse(tmax > 6.07 & tmax <= 7.12, "61-80", ">80th"))))) %>% 
  mutate(decade_label = paste0((year %/% 10) * 10)) %>% 
  group_by(decade_label, year) %>% 
  reframe(corppt = cor(ppt, value), cortmax = cor(tmax, value),
          corsummerppt = cor(pptsummer, value), corsummertmax = cor(tmaxsummer, value),
          corwintppt = cor(pptwint, value), corwinttmax = cor(tmaxwint, value),
          corfallppt = cor(pptfall, value), corfalltmax = cor(tmaxfall, value))


corrdat %>% 
  ggplot(aes(x=year, y=corfallppt, group = decade_label, col = decade_label))+
  geom_point()+
  geom_line()+
  geom_smooth(method="lm")+
  geom_hline(yintercept = 0, linetype="dashed")

corrdat %>% 
  ggplot(aes(x=year, y=corwintppt, group = decade_label, col = decade_label))+
  geom_point()+
  geom_line()+
  geom_smooth(method="lm")+
  geom_hline(yintercept = 0, linetype="dashed")

corrdat %>% 
  ggplot(aes(x=year, y=corppt))+
  geom_point()+
  geom_smooth()



# fe_mod = feols(value ~ tmax + ppt  | plot_id_needle + year,
#                 data = paneldat)
# summary(fe_mod)

paneldatcheck <- clim_core %>% 
  select(tmax, ppt, value, tree_id, plot_id_needle, year, laggedprecip, laggedtmax, lat, long, pptsummer, tmaxsummer) %>% 
  na.omit()


fe_modcheck = feNmlm(value ~ tmax + (tmax^2) +laggedtmax+ (laggedtmax^2) + pptsummer + (pptsummer^2)   | plot_id_needle + year,
                data = paneldatcheck, family = "gaussian")

fe_modcheck = feNmlm(value ~ tmax + (tmax^2) +  ppt + (ppt^2)  | plot_id_needle + year,
                     data = paneldatcheck, family = "gaussian")

summary(fe_modcheck, vcov = "conley")

fe_modnew = feNmlm(value ~ tmax + (tmax^2) + ppt + (ppt^2) | plot_id_needle + year,
                   data = clim_core, family = "gaussian")

summary(fe_modnew, vcov = "conley")


newdat <- paneldat
newdat$pred =predict(fe_mod)

ggplot(newdat, aes(x=ppt, y=pred))+
  geom_smooth(span=1, method="loess")

ggplot(paneldatcheck, aes(x = laggedtmax, y = value))+
  geom_smooth



