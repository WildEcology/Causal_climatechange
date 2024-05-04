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
  mutate(year = year + 1) %>% 
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
  left_join(treenum)

## removing outliers
clim_core <- clim_core1 %>% 
  filter(value<2.202722) %>% 
  filter(value>0.2425516) %>% 
  mutate(tree_id = paste0(tree_num, plot_id_needle))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FE panel model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paneldat <- clim_core %>% 
  select(tmax,  ppt, value, tree_id, laggedtmax, laggedprecip,plot_id_needle, year, lat, long) %>% 
  na.omit()

fe_mod = feNmlm(value ~ tmax + (tmax^2) + ppt + (ppt^2) | plot_id_needle + year,
                 data = paneldat, family = "gaussian")
# 
# summary(fe_mod, vcov = "conley")

# fe_mod = feNmlm(value ~ tmax + (tmax^2) +laggedtmax+ (laggedtmax^2) + ppt + (ppt^2)  | plot_id_needle + year,
#                      data = paneldat, family = "gaussian")

summary(fe_mod,vcov = "conley")

tab_model(fe_mod)

newdat <- paneldat
newdat$pred =predict(fe_mod)

ggplot(newdat, aes(x=laggedtmax, y=pred))+
  geom_smooth()


## checking if this is the same just adding in a squared term (do we need year fixed effects? year + clim are perfectly correlated)
# femod1 = feols(value ~ tmax + tmax2 + ppt + ppt2 | plot_id_needle + year,
#           data = paneldat)
# summary(femod1)
# 
# tab_model(femod1) 

## how to add conley standard errors
## can you also account for two-way clustering
## distributed lag model? Add lagged effects?

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Temperature and RWI figure
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## why choose 1980? Likely this is when
meant <- mean(prism_growing$tmax)

sumrydat <- prism_growing %>% 
  group_by(year) %>% 
  summarise(m.tmax = mean(tmax)) 

before80 <- sumrydat %>% 
  filter(year<1980)

after80 <- sumrydat %>% 
  filter(year>1979)

mean(after80$m.tmax)
mean(before80$m.tmax)

mean(after80$m.tmax) - mean(before80$m.tmax)


temp_fig=ggplot(sumrydat, aes(x=year, y=m.tmax, color="m.tmax"))+
  scale_color_manual(values = "#412234",
                     labels="",
                     name="")+
  geom_point(alpha=.9)+
  geom_line(alpha=.9)+
  geom_hline(yintercept = meant, 
             colour = "#412234", linetype=2)+
  geom_rect(aes(xmin=1980, xmax=2020, ymin=6.05552, ymax=6.1),
            linetype=0, alpha=.01,fill="red")+
  geom_rect(aes(xmin=1901, xmax=1980, ymin=5.653887, ymax=5.6),
            linetype=0, alpha=.01,fill="blue")+
  ylab("Max. temperature (°C)")+
  guides(color="none")+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 8))

temp_fig


rwi_sum <- paneldat %>% 
  group_by(year) %>% 
  summarize(m.rwi=mean(value))

meanrwi <- mean(rwi_sum$m.rwi)

before80rwi <- rwi_sum %>% 
  filter(year<1980)

after80rwi <- rwi_sum %>% 
  filter(year>1979)

mean(after80rwi$m.rwi)
mean(before80rwi$m.rwi)

t.test(after80rwi$m.rwi, before80rwi$m.rwi)

rwi_fig=ggplot(rwi_sum, aes(x=year, y=m.rwi, color="blue"))+
  scale_color_manual(values = "#2a9d8f",
                     labels="",
                     name="")+
  geom_point(alpha=.9)+
  geom_line(alpha=.9)+
  geom_rect(aes(xmin=1980, xmax=2018, ymin=0.9950582, ymax=0.9950582+.003),
            linetype=0, alpha=.01,fill="red")+
  geom_rect(aes(xmin=1901, xmax=1980, ymin=0.9893391, ymax=0.9893391+.003),
            linetype=0, alpha=.01,fill="blue")+
  geom_hline(yintercept = meanrwi, 
             colour = "#412234", linetype=2)+
  ylab("RWI")+
  guides(color="none")+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 8))

rwi_fig

temp_fig/rwi_fig +
  plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MC simulation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feVcov <- vcov(fe_mod, vcov = "conley")
coef_vector = fe_mod$coefficients

draw = rmvnorm(n = 100000, mean = coef_vector, sigma = feVcov)


## replacing precip and temp (do we need to do both if temp is the only sig? 
## probably not)

# bef80rwi <- paneldat %>% 
#   filter(year<1980)
#
# bef80rwi <- paneldat %>% 
#   filter(year<1980)
# 
# aft80 <- paneldat %>% 
#   filter(year>=1980)
# 
# mean(aft80$tmax)-
# mean(bef80rwi$tmax)
# 
# aft80rwi <- paneldat %>% 
#   filter(year>=2000) %>% 
#   mutate(c.tmax = tmax - 0.8442259) %>% 
#   select(-tmax) %>% 
#   rename(tmax=c.tmax) %>%
#   full_join(bef80rwi)

# aft80rwi <- paneldat %>% 
#   filter(year>1979) %>% 
#   mutate(c.tmax = tmax - 0.4016339) %>% 
#   select(-tmax) %>% 
#   rename(tmax=c.tmax) %>%
#   full_join(bef80rwi)

## maybe estimate the mean over all years and then SD of the difference and then use
## this to randomly draw the change in a MC simulation?


## creating a counterfactual scenario

## 
bef00 <- paneldat %>%
  filter(year<2000) %>% 
  select(plot_id_needle, year, tmax, laggedtmax, tree_id, laggedtmax)

aft00 <- paneldat %>%
  filter(year>=2000)

lin_counter <- lm(tmax~year+plot_id_needle, dat = bef00)
lin_lagged <- lm(laggedtmax~year+plot_id_needle, dat = bef00)

summary(lin_counter)
summary(lin_lagged)


##==============================================================================================
##   NOTE:  USE MONTE CARLO TO ESTIMATE ERROR                 
##==============================================================================================
predcounter <- ggpredict(lin_counter, terms = c("year[2000:2018]", "plot_id_needle"))
predcounter_lagged <- ggpredict(lin_lagged, terms = c("year[2000:2018]", "plot_id_needle"))

predscenario <- predcounter %>% 
  as.data.frame() %>% 
  rename(tmax = predicted, year = x, plot_id_needle = group) %>% 
  select(tmax, year, plot_id_needle) %>% 
  mutate(type = "predicted")

predscenariolagged <- predcounter_lagged %>% 
  as.data.frame() %>% 
  rename(laggedtmax = predicted, year = x, plot_id_needle = group) %>% 
  select(laggedtmax, year, plot_id_needle) %>% 
  mutate(type = "predicted")

historicclim <- bef00 %>% 
  select(tmax,laggedtmax, year, plot_id_needle) %>% 
  mutate(type = "historic") %>% 
  distinct()

counterfactual <- predscenario %>%
  left_join(predscenariolagged) %>% 
  full_join(historicclim)

counterfactual %>% 
  ggplot(aes(x=year, y=laggedtmax))+
  geom_point()+
  geom_line()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Comparing predicted to historic tmax
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## comparing predicted scenario to observed
comparepred <- aft00 %>%
  select(tmax, laggedtmax, year, plot_id_needle) %>% 
  mutate(type = "actual") %>% 
  full_join(predscenario)
  
comparepred %>% 
  ggplot(aes(x = year, y=tmax, fill = type, color=type))+
  geom_point()+
  geom_smooth()


#comparing means of actual and predicted
meanpred <- comparepred %>% 
  group_by(type) %>% 
  summarize(m = mean(tmax), l = mean(laggedtmax))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# merging counterfactual to create counterfactual dataset
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counterdat_all <- paneldat %>%
  select(-tmax, -laggedtmax) %>%
  left_join(counterfactual, by = c("plot_id_needle", "year"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# estimating RWI in response to counter and actual
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## two datasets for comparison

newpaneldat <- paneldat

dat_thresh=paneldat%>%
  distinct(plot_id_needle, year, tmax) %>%
  select(tmax) %>% 
  reframe(qs = quantile(tmax, probs = c(0.05, .15, .30, .45, .60, .85)))

dat_thresh

quartiles <- paneldat %>% 
  distinct(plot_id_needle, tmax, year) %>% 
  mutate(tmaxquants = ifelse(tmax<=4.90, "33th",
                             ifelse(tmax>4.90 & tmax<=6.35, "34-66th",
                                    ifelse(tmax>6.35, ">67", "NA")))) %>% 
  mutate(tmaxfifths = ifelse(tmax<=4.19, "20th",
                           ifelse(tmax>4.19 & tmax<=5.24, "21-40th",
                                  ifelse(tmax>5.24 & tmax <= 6.07, "41-60",
                                         ifelse(tmax > 6.07 & tmax <= 7.12, "61-80", ">80th"))))) %>% 
  mutate(tmax15 = ifelse(tmax<=3.82, "0-15",
                             ifelse(tmax>3.82 & tmax<=4.76, "16-30",
                                    ifelse(tmax>4.76 & tmax <= 5.45, "31-45",
                                           ifelse(tmax > 5.45 & tmax <= 6.07, "46-60",
                                                  ifelse(tmax > 6.07 & tmax <= 7.47, "61-85", "86-100")))))) %>% 
  mutate(tmaxmost = ifelse(tmax < 2.68, "5th",
                            ifelse(tmax >= 2.68 & tmax <=3.82, "6-15th",
                                 ifelse(tmax>3.82 & tmax<=4.76, "16-30th",
                                        ifelse(tmax>4.76 & tmax <= 5.45, "31-45",
                                               ifelse(tmax > 5.45 & tmax <= 6.07, "46-60",
                                                      ifelse(tmax > 6.07 & tmax <= 7.47, "61-85th", ">85th")))))))

test_thresh <- paneldat %>% 
  left_join(quartiles) %>% 
  mutate(timeperiod = ifelse(year>2000, "after", 'before')) %>% 
  group_by(tmaxmost, timeperiod) %>% 
  reframe(meanrwi = mean(value, na.rm=T))


mcsimdat <- paneldat %>% 
  left_join(quartiles)

df_monte=data.frame()

for (i in 1:1000){
 
  ##run the monte carlo simulation
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  mcsimdat$vals_real = predict(modified_fe_mod, newdata = paneldat)
  mcsimdat$vals_count = predict(modified_fe_mod, newdata = counterdat_all)
  
 
  mcsimdat$diff = mcsimdat$vals_real - mcsimdat$vals_count
  #mcsimdat$it = paste(i)
  
  dat_summary <- mcsimdat %>%
    group_by(tmax15) %>%
    summarise(m.diff = mean(diff), meantmax = mean(tmax)) %>%
    mutate(it = paste(i))

  tot=data.frame(dat_summary)
  df_monte=rbind(df_monte, tot)
}

df_monte_mean=data.frame()

for (i in 1:1000){
  
  ##run the monte carlo simulation
  d = draw[i,]
  modified_fe_mod = fe_mod
  modified_fe_mod$coefficients = d
  mcsimdat$vals_real = predict(modified_fe_mod, newdata = paneldat)
  mcsimdat$vals_count = predict(modified_fe_mod, newdata = counterdat_all)
  
  
  mcsimdat$diff = mcsimdat$vals_real - mcsimdat$vals_count
  #mcsimdat$it = paste(i)
  
  dat_summary <- mcsimdat %>%
    summarise(m.diff = mean(diff), meantmax = mean(tmax)) %>%
    mutate(it = paste(i))
  
  tot=data.frame(dat_summary)
  df_monte_mean=rbind(df_monte_mean, tot)
}

newdat_sim_mean <- df_monte_mean
mean(newdat_sim_mean$m.diff)
quantile(newdat_sim_mean$m.diff, probs = c(.05, .95))

newdat_sim <- df_monte 
mean(newdat_sim$m.diff)
quantile(newdat_sim$m.diff, probs = c(.05, .95))

newdat_sim %>% 
  ggplot(aes(x=tmax15, y = m.diff))+
  geom_boxplot()

quantilestat <- newdat_sim %>% 
  group_by(tmax15) %>% 
  reframe(quantslower = quantile(m.diff, probs = c(.05)),
          quantshigher = quantile(m.diff, probs = c(.95)),
          meandiff = mean(m.diff))

quantilestat %>% 
  ggplot(aes(x = tmax15, y=meandiff*100, fill = "darkgreen", color = "black"))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  geom_bar(stat="identity",position="dodge", width = 0.5, alpha=.7)+
  geom_errorbar(aes(ymin = quantslower*100, ymax = quantshigher*100), alpha=.3, 
                width=.2,position = position_dodge(0.5),size=.5 )+
  scale_color_manual(values="black")+
  scale_fill_manual(values = "#58a4b0")+
  xlab("Quantiles")+
  guides(fill=F, color = F)+
  ylab("Δ RWI")+
  ylim(-1.1, .9)+
  scale_x_discrete(limits =c("0-15" = "0-15" , 
                             "16-30" = "16-30",
                             "31-45"="31-45" ,
                             "46-60" = "46-60",
                             "61-85" = "61-85",
                             "86-100" =  "86-100"))




