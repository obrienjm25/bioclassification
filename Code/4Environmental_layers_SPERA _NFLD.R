#load required packages
library(marmap)
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
library(corrplot)

##################################
####Bathymetry and derivitives####
##################################

#Bathymetry

#read in raster file with CHS bathymetry data @ 15s resolution (~0.5 km) for NWA
NWA_bathy <- raster("Data/Rasters/Bathymetry/NWA_bathymetry_15sec_CHS.tif")

#Newfoundland

#crop to Newfoundland region study area and reproject/resample to UTM 21, Datum: NAD83, 500 m resolution

NL_SA <- readOGR("Data/Shapefiles/NL_RVsurveyAgg.shp")
r500 <- raster(extent(NL_SA),res = 500, crs = proj4string(NL_SA))
NL_bathy500m <- projectRaster(NWA_bathy,r500, method = 'bilinear', filename = "Data/Rasters/NL_bathy500m.tif", overwrite = T)
NL_bathy <- aggregate(NL_bathy500m, fact = 8, fun = mean, filename = "Data/Rasters/NL_bathy.tif", overwrite = T)

#derive other terrain characteristics from bathymetric data
NL_Slope <- terrain(NL_bathy500m, opt = 'slope', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/NL_slope.tif', overwrite = T) #slope (degrees using Horn algorithm)
NL_Slope <- aggregate(NL_Slope, fact = 8, fun = mean, filename = 'Data/Rasters/NL_slope.tif', overwrite = T)
NL_Aspect <- terrain(NL_bathy500m, opt = 'aspect', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/NL_aspect.tif', overwrite = T) #Aspect (degrees from N)
NL_Aspect <- aggregate(NL_Aspect, fact = 8, fun = mean, filename = 'Data/Rasters/NL_aspect.tif', overwrite = T)

#calculate various broad and fine scale bathymetric position indices

#matrix defining window/scale of BPI
f0.5 <- matrix(1, nrow = 3, ncol = 3) #8 neareast neighbours (NN), radius of 1 grid cell, 0.5km
f1 <- matrix(1, nrow = 5, ncol = 5) #24 NN, radius of 2 grid cells, 1 km
f2 <- matrix(1, nrow = 9, ncol = 9) #80 NN, radius of 4 grid cells, ~ km
f5 <- matrix(1, nrow = 21, ncol = 21) #440 NN, radius of 10 grid cells, 5 km
f10 <- matrix(1, nrow = 41, ncol = 41) #1681 NN, radius of 20 grid cells, 10 km
f20 <- matrix(1, nrow = 81, ncol = 81) #6561 NN, radius of 40 grid cells, 20 km

#focal function taking difference between depth in cell and mean of depth in window of given radius
NL_bpi0.5 <- focal(NL_bathy500m, w = f0.5, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/NL_bpi0.5.tif", overwrite = T)
NL_bpi1 <- focal(NL_bathy500m, w = f1, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/NL_bpi1.tif", overwrite = T)
NL_bpi2 <- focal(NL_bathy500m, w = f2, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/NL_bpi2.tif", overwrite = T)
NL_bpi5 <- focal(NL_bathy500m, w = f5, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/NL_bpi5.tif", overwrite = T)
NL_bpi10 <- focal(NL_bathy500m, w = f10, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/NL_bpi10.tif", overwrite = T)
NL_bpi20 <- focal(NL_bathy500m, w = f20, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/NL_bpi20.tif", overwrite = T)

#aggregate bpi indices at 500m res to 4km
NL_bpi0.5 <- aggregate(NL_bpi0.5, fact = 8, fun = mean, filename = "Data/Rasters/NL_bpi0.5.tif", overwrite = T)
NL_bpi1 <- aggregate(NL_bpi1, fact = 8, fun = mean, filename = "Data/Rasters/NL_bpi1.tif", overwrite = T)
NL_bpi2 <- aggregate(NL_bpi2, fact = 8, fun = mean, filename = "Data/Rasters/NL_bpi2.tif", overwrite = T)
NL_bpi5 <- aggregate(NL_bpi5, fact = 8, fun = mean, filename = "Data/Rasters/NL_bpi5.tif", overwrite = T)
NL_bpi10 <- aggregate(NL_bpi10, fact = 8, fun = mean, filename = "Data/Rasters/NL_bpi10.tif", overwrite = T)
NL_bpi20 <- aggregate(NL_bpi20, fact = 8, fun = mean, filename = "Data/Rasters/NL_bpi20.tif", overwrite = T)

###########################
####Chl a @ sea surface####
###########################

#Read in extended NWA chl layers, crop and resample
#to Newfoundland region @ 1km res, aggregate to 4 km and write raster
#read temp _> reproject temp -> aggregate reprojected

NL_SA <- readOGR("Data/Shapefiles/NL_RVsurveyAgg.shp")
r1000 <- raster(extent(NL_SA),res = 1000, crs = proj4string(NL_SA))
Rpth <- "Data/Rasters/"

chl_files <- list.files(path = "Data/Rasters/", full.names = F, pattern = 'NWA.+chl')
NL_chl  <- chl_files %>% 
  map(~ raster(paste(Rpth, .x, sep = ''))) %>% 
  map(~ projectRaster(.x, r1000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','NL',names(.x)),'.tif',sep = ''), overwrite = T)) %>%
  map(~ aggregate(.x, fact = 4, fun = mean, filename = paste(Rpth,names(.x),'.tif',sep = ''), overwrite = T))
  
###########################
####Oceanographic Data ####
###########################

#Functions and objects required to convert pt shapefiles of environmental
#variable summaries to raster format (4km resolution)

source("Code/VoronoiResample.R")

Spth <- "Data/Shapefiles/"
Rpth <- "Data/Rasters/"
extentNL <- extent(-63,-42,42,58)
NL_SA <- readOGR("Data/Shapefiles/NL_RVsurveyAgg.shp")

###mixed layer depth#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions
# and objects called by 'mixed layer depth' chunk)


#Convert Pt shapefile for NWA to raster with 4 km res for NL region

source("Code/Voronoipolygons2.R") # voronoi (thiessen) polygon function
input.file <- "mld_NWA_pts.shp"

NWA_pts <- readOGR(paste(Spth,input.file,sep='')) #read in pt shapefile
summary_stats <- colnames(NWA_pts@data)
NL_pts <- crop(NWA_pts, extentNL) #crop pts to generous extent for Maritimes
NL_pts <- spTransform(NL_pts,proj4string(NL_SA)) #reproject to UTM 21
vor_NL <- voronoipolygons(NL_pts) # create voronoi polygons for each point
proj4string(vor_NL) <- proj4string(NL_SA)

# create fishnet grid of 4km X 4km polygons spanning Maritimes study area
r.template <- raster(extent(NL_SA), res = 4000, crs = proj4string(NL_SA))
fishnet <- rasterToPolygons(r.template)
fishnet$layer <- c(1:length(fishnet$layer))

# convert voronoi and fishnet polygon df's to sf objects
vor_NL <- st_as_sf(vor_NL)
fishnet <- st_as_sf(fishnet)

# spatial join of fishnet grid with voronoi polygons
#group by fishnet polygon id and aggregate
NL_4km <- st_join(fishnet, vor_NL) %>% 
  group_by(layer) %>% 
  summarise_at(summary_stats, mean, na.rm = T) %>% 
  as(.,'Spatial')

#rasterize fishnet grid & save rasters

NL_mn_ann_mld <- rasterize(NL_4km, r.template, field = NL_4km@data$mn_mld_ann, filename = "Data/Rasters/NL_mn_ann_mld.tif", overwrite = T)
NL_mn_sum_mld <- rasterize(NL_4km, r.template, field = NL_4km@data$mn_mld_sum, filename = "Data/Rasters/NL_mn_sum_mld.tif", overwrite = T)

####Bottom Stress####  

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Stress' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for NL region

VoronoiResample(Spth,Rpth,"BtmStr_NWA_pts.shp",extentNL,NL_SA,"NL",4000,"BtmStr")

####SST#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'SST' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for NL region

VoronoiResample(Spth,Rpth,"SST_NWA_pts.shp",extentNL,NL_SA,"NL",4000,"SST")

####Bottom Temperature####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Temperature' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for NL region

VoronoiResample(Spth,Rpth,"BT_NWA_pts.shp",extentNL,NL_SA,"NL",4000,"BT")

####Bottom Salinity####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Salinity' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for NL region

VoronoiResample(Spth,Rpth,"BSal_NWA_pts.shp",extentNL,NL_SA,"NL",4000,"BSal")

####VelEW####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelEW' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for NL region

VoronoiResample(Spth,Rpth,"VelEW_NWA_pts.shp",extentNL,NL_SA,"NL",4000,"VelEW")

####VelNS####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelNS' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for NL region

VoronoiResample(Spth,Rpth,"VelNS_NWA_pts.shp",extentNL,NL_SA,"NL",4000,"VelNS")

############################
####Primary Productivity####
############################

#Read in extended NWA PP layers,
#Crop to Newfoundland, reproject to UTM 21 and resample to 4km resolution

NL_SA <- readOGR("Data/Shapefiles/NL_RVsurveyAgg.shp")
r4000 <- raster(extent(NL_SA),res = 4000, crs = proj4string(NL_SA))
Rpth <- "Data/Rasters/"

PP_files <- list.files(path = Rpth, full.names = F, pattern = 'NWA.+PP')

NL_PP <- PP_files %>% 
  map(~ raster(paste(Rpth,.x, sep = ''))) %>%
  map(~ projectRaster(.x, r4000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','NL',names(.x)),'.tif',sep = ''), overwrite = T))

####DO,pH,Nutrients####

NL_SA <- readOGR("Data/Shapefiles/NL_RVsurveyAgg.shp") #NL Study Area shapefile
r4000 <- raster(extent(NL_SA),res = 4000, crs = proj4string(NL_SA)) # template raster in UTM 21N with 4 km resolution
Rpth <- "Data/Rasters/" #Raster file path
NWA.ext <- extent(-70.9, -42, 41, 58.6) #Extent object for Northwest Atl.
LandBorders <- readOGR("Data/Shapefiles/NL_landborders.shp") #Shapefile w/ provincial $ St. Pierre $ Miquelon boundaries
LandBuffer <- gBuffer(LandBorders, width = 5000) #create 5 km buffer around land borders
writeOGR(as(LandBuffer,'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'NL_LandBuffer_5km', driver = 'ESRI Shapefile')
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
  map(~ mask(.x, LandBuffer, inverse = T, filename = paste(Rpth,'NL_',gsub('BO_','',names(.x)),'.tif',sep = ''), overwrite = T))
