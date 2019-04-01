#load required packages
library(raster)
library(data.table)
library(rgdal)
library(rgeos)
library(sp)
library(dplyr)
library(tidyr)
library(gsubfn)
library(purrr)
library(RCurl)
library(naniar)
library(sf)

##################################
####Bathymetry and derivitives####
##################################

#Bathymetry

#read in raster file with GEBCO bathymetry data @ 30s resolution (~1 km) for NWA
NWA_bathy_gebco <- raster("Data/Rasters/Bathymetry/GEBCO_30sec.tif")
NWA_bathy_gebco <- crop(NWA_bathy_gebco, extent(-70.9, -42, 41, 58.6))
values(NWA_bathy_gebco)[values(NWA_bathy_gebco) > 0] <- NA
writeRaster(NWA_bathy_gebco,'Data/Rasters/Bathymetry/NWA_GEBCO_30sec.tif', overwrite = T)

#Gulf

#crop to Gulf region study area and reproject/resample to UTM 20, Datum: NAD83, 1 km resolution

Gulf_SA <- readOGR("Data/Shapefiles/Gulf_RVsurveyAgg.shp")
r1000 <- raster(extent(Gulf_SA),res = 1000, crs = proj4string(Gulf_SA))
Gulf_bathy1km <- projectRaster(NWA_bathy_gebco,r1000, method = 'bilinear', filename = "Data/Rasters/Gulf_bathy1km.tif", overwrite = T)
Gulf_bathy <- aggregate(Gulf_bathy1km, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/Gulf_bathy.tif", overwrite = T)
#Gulf_bathy <- raster("Data/Rasters/Gulf_bathy.tif")
#Gulf_bathy <- crop(Gulf_bathy, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_bathy.tif", overwrite = T)

#derive other terrain characteristics from bathymetric data
Gulf_Slope <- terrain(Gulf_bathy1km, opt = 'slope', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/Gulf_slope.tif', overwrite = T) #slope (degrees using Horn algorithm)
Gulf_Slope <- aggregate(Gulf_Slope, fact = 4, fun = mean, na.rm = T, filename = 'Data/Rasters/Gulf_slope.tif', overwrite = T)
#Gulf_Slope <- raster("Data/Rasters/Gulf_Slope.tif")
#Gulf_Slope <- crop(Gulf_Slope, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_Slope.tif", overwrite = T)

Gulf_Aspect <- terrain(Gulf_bathy1km, opt = 'aspect', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/Gulf_aspect.tif', overwrite = T) #Aspect (degrees from N)
Gulf_Aspect <- aggregate(Gulf_Aspect, fact = 4, fun = mean, na.rm = T, filename = 'Data/Rasters/Gulf_aspect.tif', overwrite = T)
#Gulf_Aspect <- raster("Data/Rasters/Gulf_Aspect.tif")
#Gulf_Aspect <- crop(Gulf_Aspect, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_Aspect.tif", overwrite = T)

#calculate various broad and fine scale bathymetric position indices

#matrix defining window/scale of BPI
f1 <- matrix(1, nrow = 3, ncol = 3) #8 neareast neighbours (NN), radius of 1 grid cell, 1 km
f2 <- matrix(1, nrow = 5, ncol = 5) #24 NN, radius of 2 grid cells, 2 km
f5 <- matrix(1, nrow = 11, ncol = 11) #120 NN, radius of 5 grid cells,  5 km
f10 <- matrix(1, nrow = 21, ncol = 21) #440 NN, radius of 10 grid cells, 10 km
f20 <- matrix(1, nrow = 41, ncol = 41) #1681 NN, radius of 20 grid cells, 20 km

#focal function taking difference between depth in cell and mean of depth in window of given radius

Gulf_bpi1 <- focal(Gulf_bathy1km, w = f1, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Gulf_bpi1.tif", overwrite = T)
Gulf_bpi2 <- focal(Gulf_bathy1km, w = f2, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Gulf_bpi2.tif", overwrite = T)
Gulf_bpi5 <- focal(Gulf_bathy1km, w = f5, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Gulf_bpi5.tif", overwrite = T)
Gulf_bpi10 <- focal(Gulf_bathy1km, w = f10, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Gulf_bpi10.tif", overwrite = T)
Gulf_bpi20 <- focal(Gulf_bathy1km, w = f20, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Gulf_bpi20.tif", overwrite = T)

#aggregate bpi indices at 1000m res to 4km

Gulf_bpi1 <- aggregate(Gulf_bpi1, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/Gulf_bpi1.tif", overwrite = T)
Gulf_bpi2 <- aggregate(Gulf_bpi2, fact = 4, fun = mean, na.rm = T,filename = "Data/Rasters/Gulf_bpi2.tif", overwrite = T)
Gulf_bpi5 <- aggregate(Gulf_bpi5, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/Gulf_bpi5.tif", overwrite = T)
Gulf_bpi10 <- aggregate(Gulf_bpi10, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/Gulf_bpi10.tif", overwrite = T)
Gulf_bpi20 <- aggregate(Gulf_bpi20, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/Gulf_bpi20.tif", overwrite = T)

# Gulf_bpi1 <- raster("Data/Rasters/Gulf_bpi1.tif")
# Gulf_bpi1 <- crop(Gulf_bpi1, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_bpi1.tif", overwrite = T)
# Gulf_bpi2 <- raster("Data/Rasters/Gulf_bpi2.tif")
# Gulf_bpi2 <- crop(Gulf_bpi2, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_bpi2.tif", overwrite = T)
# Gulf_bpi5 <- raster("Data/Rasters/Gulf_bpi5.tif")
# Gulf_bpi5 <- crop(Gulf_bpi5, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_bpi5.tif", overwrite = T)
# Gulf_bpi10 <- raster("Data/Rasters/Gulf_bpi10.tif")
# Gulf_bpi10 <- crop(Gulf_bpi10, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_bpi10.tif", overwrite = T)
# Gulf_bpi20 <- raster("Data/Rasters/Gulf_bpi20.tif")
# Gulf_bpi20 <- crop(Gulf_bpi20, extent(280798,720798,5059640,5447640), filename = "Data/Rasters/Gulf_bpi20.tif", overwrite = T)

###########################
####Chl a @ sea surface####
###########################

#Read in extended NWA chl layers, crop and resample
#to Newfoundland region @ 1km res, aggregate to 4 km and write raster
#read temp _> reproject temp -> aggregate reprojected

Gulf_SA <- readOGR("Data/Shapefiles/Gulf_RVsurveyAgg.shp")
r1000 <- raster(extent(Gulf_SA),res = 1000, crs = proj4string(Gulf_SA))
Rpth <- "Data/Rasters/"

chl_files <- list.files(path = "Data/Rasters/", full.names = F, pattern = 'NWA.+chl')
Gulf_chl  <- chl_files %>% 
  map(~ raster(paste(Rpth, .x, sep = ''))) %>% 
  map(~ projectRaster(.x, r1000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','Gulf',names(.x)),'.tif',sep = ''), overwrite = T)) %>%
  map(~ aggregate(.x, fact = 4, fun = mean, filename = paste(Rpth,names(.x),'.tif',sep = ''), overwrite = T))

# chl_list <- list.files('Data/Rasters/', full.names = T, pattern  = 'Gulf.+chl')
# chl_stack <- stack(chl_list)
# chl_stack <- crop(chl_stack, extent(280798,720798,5059640,5447640))
# chl_names <- paste(Rpth, names(chl_stack),'.tif', sep = '')
# writeRaster(chl_stack, filename = chl_names, bylayer = T, overwrite = T)

###########################
####Oceanographic Data ####
###########################

#Functions and objects required to convert pt shapefiles of environmental
#variable summaries to raster format (4km resolution)

source("Code/VoronoiResample.R")

Spth <- "Data/Shapefiles/"
Rpth <- "Data/Rasters/"
extentGulf <- extent(-66.5,-59.5,45,49.75)
Gulf_SA <- readOGR("Data/Shapefiles/Gulf_RVsurveyAgg.shp")

###mixed layer depth#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions
# and objects called by 'mixed layer depth' chunk)


#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

source("Code/Voronoipolygons2.R") # voronoi (thiessen) polygon function
input.file <- "mld_NWA_pts.shp"

NWA_pts <- readOGR(paste(Spth,input.file,sep='')) #read in pt shapefile
summary_stats <- colnames(NWA_pts@data)
Gulf_pts <- crop(NWA_pts, extentGulf) #crop pts to generous extent for Maritimes
Gulf_pts <- spTransform(Gulf_pts,proj4string(Gulf_SA)) #reproject to UTM 20
vor_Gulf <- voronoipolygons(Gulf_pts) # create voronoi polygons for each point
proj4string(vor_Gulf) <- proj4string(Gulf_SA)

# create fishnet grid of 4km X 4km polygons spanning Gulf study area
r.template <- raster(extent(Gulf_SA), res = 4000, crs = proj4string(Gulf_SA))
fishnet <- rasterToPolygons(r.template)
fishnet$layer <- c(1:length(fishnet$layer))

# convert voronoi and fishnet polygon df's to sf objects
vor_Gulf <- st_as_sf(vor_Gulf)
fishnet <- st_as_sf(fishnet)

# spatial join of fishnet grid with voronoi polygons
#group by fishnet polygon id and aggregate
Gulf_4km <- st_join(fishnet, vor_Gulf) %>% 
  group_by(layer) %>% 
  summarise_at(summary_stats, mean, na.rm = T) %>% 
  as(.,'Spatial')

#rasterize fishnet grid & save rasters

Gulf_mn_ann_mld <- rasterize(Gulf_4km, r.template, field = Gulf_4km@data$mn_mld_ann, filename = "Data/Rasters/Gulf_mn_ann_mld.tif", overwrite = T)
Gulf_mn_sum_mld <- rasterize(Gulf_4km, r.template, field = Gulf_4km@data$mn_mld_sum, filename = "Data/Rasters/Gulf_mn_sum_mld.tif", overwrite = T)

####Bottom Stress####  

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Stress' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

VoronoiResample(Spth,Rpth,"BtmStr_NWA_pts.shp",extentGulf,Gulf_SA,"Gulf",4000,"BtmStr")

####SST#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'SST' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

VoronoiResample(Spth,Rpth,"SST_NWA_pts.shp",extentGulf,Gulf_SA,"Gulf",4000,"SST")

####Bottom Temperature####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Temperature' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

VoronoiResample(Spth,Rpth,"BT_NWA_pts.shp",extentGulf,Gulf_SA,"Gulf",4000,"BT")

####Bottom Salinity####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Salinity' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

VoronoiResample(Spth,Rpth,"BSal_NWA_pts.shp",extentGulf,Gulf_SA,"Gulf",4000,"BSal")

####VelEW####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelEW' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

VoronoiResample(Spth,Rpth,"VelEW_NWA_pts.shp",extentGulf,Gulf_SA,"Gulf",4000,"VelEW")

####VelNS####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelNS' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for Gulf region

VoronoiResample(Spth,Rpth,"VelNS_NWA_pts.shp",extentGulf,Gulf_SA,"Gulf",4000,"VelNS")

############################
####Primary Productivity####
############################

#Read in extended NWA PP layers,
#Crop to Newfoundland, reproject to UTM 21 and resample to 4km resolution

Gulf_SA <- readOGR("Data/Shapefiles/Gulf_RVsurveyAgg.shp")
r4000 <- raster(extent(Gulf_SA),res = 4000, crs = proj4string(Gulf_SA))
Rpth <- "Data/Rasters/"

PP_files <- list.files(path = Rpth, full.names = F, pattern = 'NWA.+PP')

Gulf_PP <- PP_files %>% 
  map(~ raster(paste(Rpth,.x, sep = ''))) %>%
  map(~ projectRaster(.x, r4000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','Gulf',names(.x)),'.tif',sep = ''), overwrite = T))

####DO,pH,Nutrients####

Gulf_SA <- readOGR("Data/Shapefiles/Gulf_RVsurveyAgg.shp") #Gulf Study Area shapefile
r4000 <- raster(extent(Gulf_SA),res = 4000, crs = proj4string(Gulf_SA)) # template raster in UTM 21N with 4 km resolution
Rpth <- "Data/Rasters/" #Raster file path
NWA.ext <- extent(-70.9, -42, 41, 58.6) #Extent object for Northwest Atl.
LandBorders <- readOGR("Data/Shapefiles/Gulf_landborders.shp") #Shapefile w/ provincial boundaries
LandBuffer <- gBuffer(LandBorders, width = 5000) #create 5 km buffer around land borders
writeOGR(as(LandBuffer,'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'Gulf_LandBuffer_5km', driver = 'ESRI Shapefile')
library(sdmpredictors)
sdm_datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)
sdm_layers <- list_layers(sdm_datasets, terrestrial=FALSE, marine=TRUE, monthly=TRUE,
                  version=NULL)
BO_list <- list('BO_dissox','BO_nitrate','BO_ph', 'BO_phosphate','BO_silicate') 
names(BO_list) <- c('DO','nitrate','pH', 'phosphate','silicate')
#download environmental predictor layers from Bio-Oracle database, crop to study area,
#disaggregate to smaller cell size, resample and reproject to match other raster layers,
#mask cells within 5 km radius of land, write to file
BO_layers <- BO_list %>% 
  map(~ load_layers(.x)) %>%
  map(~ crop(.x, NWA.ext)) %>%  
  map(~ disaggregate(.x, fact = 2, method = 'bilinear', na.rm = T)) %>%  
  map(~ projectRaster(.x, r4000, method = 'bilinear', na.rm = T)) %>% 
  map(~ mask(.x, LandBuffer, inverse = T, filename = paste(Rpth,'Gulf_',gsub('BO_','',names(.x)),'.tif',sep = ''), overwrite = T))
