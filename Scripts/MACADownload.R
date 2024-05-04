## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-04-20
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes:
## precip is in mm
## Fill value = -9999
## coordinates: time lat lon 
## coordinate_system: WGS84,EPSG:4326
## ---------------------------

## packages

librarian::shelf(HelpersMG, ncdf4, tidyverse, chron, lattice, RColorBrewer, terra, magrittr)

select=dplyr::select

theme_set(
  theme_bw(base_size = 11)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


  
#downloaded VPD data from https://climate.northwestknowledge.net/MACA/data_portal.php


## reading in URLs
urls <- read_csv("Data/MACA_files/MACA_URLs.csv")

## setting the working directory where files will be downloaded
wdir <- "/Users/treelife/Documents/Papers/Causal inference paper/Causal_climatechange/Data/MACA_data/"
setwd(wdir)

sliceurls <- urls %>% 
  filter(Var!=tmn) 

dim(sliceurls)

## downloading through the for loop took about 50 minutes

# downloading the data
# for (i in 1:80) {
# 
#   # Extract filename
#   files <- sliceurls %>%
#     slice(i) %>%
#     mutate(name = paste0(Var, Scenario, "_", i, ".nc"))
# 
#   filename <- files$name
# 
#   url <- files$File
# 
#   # Construct full path to save the file
#   save_path <- file.path(wdir, filename)
# 
#   # Download the file
#   download.file(url, destfile = save_path, mode = "wb")
# 
# }

t1 = terra::rast("pptRCP8.5_1.nc")
names(t1)
t1[names(t1)[2]]





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OLD CODE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# #Open the NetCDF file
# nc <- nc_open("pptRCP8.5_1.nc")
# summary(nc)
# print(nc)
# 
# # Read the latitude and longitude variables
# lon <- ncvar_get(nc, "lon")
# lat <- ncvar_get(nc, "lat", verbose = F)
# t <- ncvar_get(nc, "time")
# 
# ## collating precip data into array
# precip.array <- ncvar_get(nc, "precipitation") # store the data in a 3-dimensional array
# dim(precip.array) 
# 
# ## what is the na value
# fillvalue <- ncatt_get(nc, "precipitation", "_FillValue")
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
# 
# dim(precip.slice)
# 
# r <- rast(t(precip.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), 
#           ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# 
# 
# # Create a terra SpatRaster object for lat and lon
# r <- rast(lat, lon)
# 
# # Create an empty list to store values for each time step
# values_list <- list()
# 
# # Loop through each time step
# #for (i in 1:nc$nvars) {
#   
#   # Read data for the i-th variable (assuming it contains the desired variable)
#   data <- ncvar_get(nc, nc$var[i]$name)
#   
#   # Extract values at points
#   # Assuming 'points' is a SpatialPoints or a SpatialPointsDataFrame object
#   values <- extract(r, points, at = data)
#   
#   # Add extracted values to the list
#   values_list[[i]] <- values
# }
# 
# # Close the NetCDF file
# nc_close(nc)
# 
# # Convert the list of values to a data frame
# values_df <- do.call(cbind, values_list)
# 
# # Add meaningful column names if needed
# # colnames(values_df) <- c("value1", "value2", ...)
# 
# # Display or further process the extracted values
# print(values_df)
# 
# filenames <- list.files(pattern='*.nc',full.names=TRUE)
# soils <- raster::brick(filenames, lvar=4, level=1)
# 
# ##reading inn rpc4.5
# ncfname <- filenames 
# dname <- "vpd"  # note: tmp means temperature (not temporary)
# 
# ncin <- nc_open(ncfname)
# print(ncin)
# 
# 
# plotdata=read_csv("latlongs_july15.csv")
# plotdata$lon=plotdata$long+360
# plotdata=plotdata%>%dplyr::select(-long)
# 
# vpd <- brick("2016_2020.nc")
# 
# xy <- plotdata[,c(3,1)] #Column 1 is longitude and column 2 is latitude
# xy
# spts <- SpatialPoints(xy, proj4string=CRS("+proj=longlat +datum=WGS84"))
# #Extract data by spatial point
# vpd2 <- extract(vpd, spts)
# vpd3=as.data.frame(vpd2)
# vpd3$plot=1:154
# 
# mean_vpd=vpd3%>%
#   pivot_longer(-plot)%>%
#   group_by(plot)%>%
#   summarize(vpd_kpa=mean(value))%>%
#   mutate(vpd_hpa=vpd_kpa*10, plot=1:154)##convert kpa to hpa
# 
# 
# #write_csv(mean_vpd, "maca2016_2020.csv")
# 
# ##calculating change in vpd
# mac_2020=read_csv("maca data/maca2016_2020.csv")
# mac_2020=mac_2020%>%select(plot, vpd_hpa)%>%set_colnames(c("plot", "vpd2020"))
# 
# mac_2056=read_csv("maca data/maca2056_2060.csv")
# mac_2056=mac_2056%>%select(plot, vpd_hpa)%>%set_colnames(c("plot", "vpd2056"))
# 
# mac_2099=read_csv("maca data/future2099.csv")
# mac_2099=mac_2099%>%select(plot, vpd_hpa)%>%set_colnames(c("plot", "vpd2099"))
# 
# 
# macdat=mac_2020%>%
#   left_join(mac_2056)%>%
#   pivot_longer(-plot)
# 
# ggplot(macdat, aes(x=value, fill=name))+
#   geom_density(alpha=.5)+
#   scale_fill_manual(values=c("#69b3a2", "#495DB2"), 
#                     labels=c("2016_2020", "2056-2060"),
#                     name="RCP4.5")+
#   ylab("Density")+
#   scale_x_continuous(name="Mean predicted VPD (hPa)", breaks = scales::pretty_breaks(n = 10))
# 
# 
# macdatdiff=mac_2020%>%
#   left_join(mac_2056)%>%
#   mutate(perdiff=((vpd2056-vpd2020)/vpd2020)*100)%>%
#   summarize(mean(perdiff))

##csv for mean per diff
#write_csv(macdatdiff, "maca_plot_perdiff.csv")

##8.2% increase in vpd by 2056

# # get longitude and latitude
# lon <- ncvar_get(ncin,"lon")
# nlon <- dim(lon)
# 
# head(lon)
# dim(lon)
# max(lon)
# min(lon)
# 
# lat <- ncvar_get(ncin,"lat")
# nlat <- dim(lat)
# head(lat)
# max(lat)
# min(lat)
# 
# print(c(nlon,nlat))
# 
# # get time
# time <- ncvar_get(ncin,"time")
# time
# 
# tunits <- ncatt_get(ncin,"time","units")
# nt <- dim(time)
# nt
# tunits
# 
# 
# tmp_array <- ncvar_get(ncin,dname)
# dlname <- ncatt_get(ncin,dname,"long_name")
# dunits <- ncatt_get(ncin,dname,"units")
# fillvalue <- ncatt_get(ncin,dname,"_FillValue")
# dim(tmp_array)
# 
# 
# title <- ncatt_get(ncin,0,"title")
# institution <- ncatt_get(ncin,0,"institution")
# datasource <- ncatt_get(ncin,0,"source")
# references <- ncatt_get(ncin,0,"references")
# history <- ncatt_get(ncin,0,"history")
# Conventions <- ncatt_get(ncin,0,"Conventions")
# 
# m <-60
# tmp_slice <- tmp_array[,,m]
# image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))
# 
# help(expand.grid)
# grid <- expand.grid(lon=lon, lat=lat)
# #cutpts <- c(-50,-40,-30,-20,-10,0,10,20,30,40,50)
# levelplot(tmp_slice ~ lon * lat, data=grid, cuts=11, pretty=T, 
#           col.regions=(rev(brewer.pal(10,"RdBu"))))
# 
# # create dataframe -- reshape data
# # matrix (nlon*nlat rows by 2 cols) of lons and lats
# lonlat <- as.matrix(expand.grid(lon,lat))
# dim(lonlat)
# 
# 
# 
# # vector of `vpd` values
# tmp_vec <- as.vector(tmp_slice)
# length(tmp_vec)
# 
# # create dataframe and add names
# tmp_df01 <- data.frame(cbind(lonlat, tmp_slice))
# names(tmp_df01) <- c("lon","lat",paste(dname,as.character(m), sep="_"))
# head(na.omit(tmp_df01), 10)
# 
# plotdata=read_csv("latlongs_july15.csv")
# plotdata$lon=plotdata$long+360
# plotdata=plotdata%>%dplyr::select(-long)
# 
# 
# 
# 
# 
# 
# ggplot(mean_vpd, aes(x=vpd_hpa))+
#   geom_density()
# 
# colnames(mean_vpd)="vpd"
# plotmeans=data.frame(mean_vpd, vpd3$plot)
# 
# 
# temp3 <- t(temp2) #transpose raster object
# colnames(temp2) <- col.name[1,] #It would be better if you have the location names corresponding to the points
# head(temp3)
# 
# col.name=t(plotdata[2])

