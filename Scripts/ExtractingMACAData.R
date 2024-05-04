## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-04-21
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

## next steps -- read in all future scenarios
## extract values for whitebark pine range



# Packages
librarian::shelf(tidyverse, terra, sgsR, tmap, tmaptools, sf,  ncdf4, CFtime)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# reading in whitebark pine plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## reading in plot lat/longs
plotdat <- read_csv("Data/Data_raw/plot_names_normalized_v2.csv")

latlong <- plotdat %>% 
  select(lat, long, plot_id_needle) %>% 
  na.omit()

colnames(latlong)[c(2)] = c("lon")

## project
latlons <- terra::vect(latlong, geom=c('lon','lat'), crs="EPSG:4326")
plot(latlons, col="blue")

prjcrs = "+proj=utm +zone=11 +datum=WGS84"

cat(crs(latlons), "\n")
crs(latlons, proj=TRUE)

projlatlon <- project(latlons, prjcrs)
plot(projlatlon)
crs(projlatlon, proj=TRUE)

latlons_sf <- sf::st_as_sf(projlatlon)
plot(latlons_sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extracting MACA data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## setwd to maca file location
wdir <- "/Users/treelife/Documents/Papers/Causal inference paper/Causal_climatechange/Data/MACA_data/"
setwd(wdir)




## reading in dem data
pptdat <- rast("pptRCP4.5_21.nc")
crs(pptdat) <- "EPSG:4326"
projppt <- project(pptdat, prjcrs)
hist(projppt$precipitation_4)
plot(projppt[[1]])

## checking that both datasets match
cat(crs(projppt), "\n")
crs(projppt, proj=T)
crs(projlatlon, proj = T)

## checking to see if the points and rasters overlap correctly
tmap_mode("view")

firstfig <- tm_shape(latlons_sf)+
    tm_dots()
secondfig <- tm_shape(projppt[[2]])+
  tm_raster()

firstfig + secondfig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extracting the dates for each netcdf file to become column labels
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc_data <- nc_open("pptRCP8.5_1.nc")
summary(nc_data)

t <- ncatt_get(nc_data, "time", attname="units")

cf <- CFtime(nc_data$dim$time$units, nc_data$dim$time$calendar, nc_data$dim$time$vals)
dates <- CFtimestamp(cf)
range(dates)

nc_close(nc_data) 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extracting precip and temperature values
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pptvalues <- terra::extract(projppt, projlatlon, df=T)
colnames(pptvalues)[c(2:1081)] = dates

plotnames = projlatlon$plot_id_needle
pptvalues$plot_id_needle <- plotnames

longppt <- pptvalues %>% 
  select(-ID) %>% 
  pivot_longer(-plot_id_needle) %>% 
  separate(name, c("year", "month", "day")) %>% 
  mutate(year = as.numeric(year), month = as.numeric(month), day = as.numeric(day))

hist(longppt$value)
range(longppt$year)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimating lagged and current precip values
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maca_growing_ppt <- longppt %>%
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  group_by(plot_id_needle, growing) %>% 
  summarize(ppt=sum(value, na.rm=T))%>%
  rename(year=growing, macappt = ppt) %>% 
  filter(year>2010 & year < 2100) %>% ## 2010 and 2100 are not full years
  ungroup() 

## get lagged precip and temp
maca_ppt_lagged <- longppt %>%
  mutate(growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September to be consistent
  group_by(plot_id_needle, growing) %>%
  filter(growing < 2099) %>% ## creating a variable with lagged year
  summarize(ppt=sum(value, na.rm=T))%>%
  rename(lagged_year=growing, laggedmacappt = ppt) %>% 
  filter(lagged_year>2010) %>% ## 2010 is not a full year
  mutate(year = lagged_year +1) %>% 
  ungroup() 


## merging lagged precip and current precip

macappt <- maca_growing_ppt %>% 
  left_join(maca_ppt_lagged)






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OLD CODE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# # Create some example spatial points
# points <- st_as_sf(projlatlon)
# 
# # Get the bounding box of the points
# bbox <- st_bbox(points)
# 
# # Convert the bounding box to a polygon
# bbox_polygon <- st_as_sfc(st_bbox(points))
# bbox_vect <- vect(bbox_polygon)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Making sure the data layers are correct; visualization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# nc_data <- nc_open("pptRCP8.5_1.nc")
# summary(nc_data)
# print(nc_data)
# 
# lon <- ncvar_get(nc_data, "lon")
# lon[lon > 180] <- lon[lon > 180] - 360
# 
# lat <- ncvar_get(nc_data, "lat", verbose = F)
# t <- ncatt_get(nc_data, "time", attname="units")
# 
# head(lon) 
# 
# precip.array <- ncvar_get(nc_data, "precipitation") # store the data in a 3-dimensional array
# dim(precip.array) 
# 
# fillvalue <- ncatt_get(nc_data, "precipitation", "_FillValue")
# fillvalue
# 
# cf <- CFtime(nc_data$dim$time$units, nc_data$dim$time$calendar, nc_data$dim$time$vals)
# dates <- CFtimestamp(cf)
# range(dates)
# 
# nc_close(nc_data) 
# 



# precip.array[precip.array == fillvalue$value] <- NA
# precip.slice <- precip.array[, , 1]
# dim(precip.slice)
# 
# r <- raster(t(precip.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs="EPSG:4326")
# r <- flip(r, direction='y')
# plot(r)
# 
# tmap_mode("view")
# 
# firstfig <- tm_shape(latlons_sf)+
#     tm_dots()
# secondfig <- tm_shape(r)+
#   tm_raster()
# 
# bothfig <- firstfig +
#   tm_shape(r)+
#   tm_raster(style= "pretty",
#             title="Precip", alpha = 0.1)
# bothfig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Working with the data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# time <- ncvar_get(nc_data, "time", attname="units")
# head(time)
# 
# tunits <- ncatt_get(nc_data, "time", "units")
# 
# head(tunits)
# nt <- dim(time,)
# 
# attributes(nc_data$var)
# attributes(nc_data$dim)
# 
# lat <- ncvar_get(nc_data, "lat")
# nlat <- dim(lat) #to check it matches the metadata: 23
# lon <- ncvar_get(nc_data, "lon")
# nlon <- dim(lon) #to check, should be 24
# 
# # Check your lat lon dimensions match the information in the metadata we explored before:
# print(c(nlon, nlat))
# 
# # Read the latitude and longitude variables
# lon <- ncvar_get(nc_data, "lon")
# lat <- ncvar_get(nc_data, "lat", verbose = F)
# t <- ncvar_get(nc_data, "time")
# 
# ## collating precip data into array
# precip.array <- ncvar_get(nc_data, "precipitation") # store the data in a 3-dimensional array
# dim(precip.array)
# 
# ## what is the na value
# fillvalue <- ncatt_get(nc_data, "precipitation", "_FillValue")
# fillvalue
# 
# ## netcdf file
# nc_close(nc)
# 
# ## fill -9999 with NA
# precip.array[precip.array == fillvalue$value] <- NA
# 
# ## having a look
# precip.slice <- precip.array[, , 1]
# dim(precip.slice)
# 
# 
# dates <- CFtimestamp(cf)
# 
# dim(dates)
# range(dates)
# 
# lonlattime <- as.matrix(expand.grid(lon,lat,dates))
# vec_long <- as.vector(precip.array)
# length(vec_long)
# 
# all_obs <- data.frame(cbind(lonlattime, vec_long))
# ## note: do not use head() here; takes too long
# dim(all_obs)
# 
# colnames(all_obs) <- c("lon", "lat","date","ppt_mm")
# 
# ## removing NAs and getting lat/long to match with whitebark pine plots
# rmNa_allobs <- all_obs %>% 
#   na.omit() %>%
#   mutate(long = as.numeric(lon), lat = as.numeric(lat),
#          ppt_mm = as.numeric(ppt_mm)) %>% 
#   mutate(lat = round(lat, digits = 4), long = round(long, digits = 4))
# 
# range(rmNa_allobs$lat)
# 
# 
# ## not a huge number removed
# dim(rmNa_allobs)
# 
# sliced = rmNa_allobs %>% 
#   filter(lon == -120.8972 & lat == 35.3546)
# 
# sliced %>%
#   slice(999:1016) %>% 
#   ggplot(aes(x=date, y=ppt_mm, group = 1))+
#   geom_point()+
#   geom_line()
# 
# ## merging with the dataset
# head(latlong)
# 
# colnames(latlong)
# 
# alldat_plot <- latlong %>% 
#   right_join(rmNa_allobs)
# 
# check <- rmNa_allobs %>% 
#   filter(lat == 35.3546)
# 
# sliced <- rmNa_allobs %>% 
#   slice(1:500)
# 
# plot(sliced$ppt_mm)
# 
