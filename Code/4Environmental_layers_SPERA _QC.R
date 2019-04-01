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

#QC

#crop to QC region study area and reproject/resample to UTM 20, Datum: NAD83, 1 km resolution

QC_SA <- readOGR("Data/Shapefiles/QC_StudyArea.shp")
r1000 <- raster(extent(QC_SA),res = 1000, crs = proj4string(QC_SA))
QC_bathy1km <- projectRaster(NWA_bathy_gebco,r1000, method = 'bilinear', filename = "Data/Rasters/QC_bathy1km.tif", overwrite = T)
QC_bathy <- aggregate(QC_bathy1km, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/QC_bathy.tif", overwrite = T)
QC_bathy <- crop(QC_bathy, extent(22174.67,1026175,5192424,5820424), filename = "Data/Rasters/QC_bathy.tif", overwrite = T) #align extent with other rasters

#derive other terrain characteristics from bathymetric data
QC_Slope <- terrain(QC_bathy1km, opt = 'slope', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/QC_slope.tif', overwrite = T) #slope (degrees using Horn algorithm)
QC_Slope <- aggregate(QC_Slope, fact = 4, fun = mean, na.rm = T, filename = 'Data/Rasters/QC_slope.tif', overwrite = T)
QC_Slope <- crop(QC_Slope, extent(22174.67,1026175,5192424,5820424), filename = "Data/Rasters/QC_Slope.tif", overwrite = T) #align extent with other rasters

QC_Aspect <- terrain(QC_bathy1km, opt = 'aspect', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/QC_aspect.tif', overwrite = T) #Aspect (degrees from N)
QC_Aspect <- aggregate(QC_Aspect, fact = 4, fun = mean, na.rm = T, filename = 'Data/Rasters/QC_aspect.tif', overwrite = T)
QC_Aspect <- crop(QC_Aspect, extent(22174.67,1026175,5192424,5820424), filename = "Data/Rasters/QC_Aspect.tif", overwrite = T) #align extent with other rasters

#calculate various broad and fine scale bathymetric position indices

#matrix defining window/scale of BPI
f1 <- matrix(1, nrow = 3, ncol = 3) #8 neareast neighbours (NN), radius of 1 grid cell, 1 km
f2 <- matrix(1, nrow = 5, ncol = 5) #24 NN, radius of 2 grid cells, 2 km
f5 <- matrix(1, nrow = 11, ncol = 11) #120 NN, radius of 5 grid cells,  5 km
f10 <- matrix(1, nrow = 21, ncol = 21) #440 NN, radius of 10 grid cells, 10 km
f20 <- matrix(1, nrow = 41, ncol = 41) #1681 NN, radius of 20 grid cells, 20 km

#focal function taking difference between depth in cell and mean of depth in window of given radius

QC_bpi1 <- focal(QC_bathy1km, w = f1, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/QC_bpi1.tif", overwrite = T)
QC_bpi2 <- focal(QC_bathy1km, w = f2, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/QC_bpi2.tif", overwrite = T)
QC_bpi5 <- focal(QC_bathy1km, w = f5, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/QC_bpi5.tif", overwrite = T)
QC_bpi10 <- focal(QC_bathy1km, w = f10, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/QC_bpi10.tif", overwrite = T)
QC_bpi20 <- focal(QC_bathy1km, w = f20, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/QC_bpi20.tif", overwrite = T)

#aggregate bpi indices at 1000m res to 4km

QC_bpi1 <- aggregate(QC_bpi1, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/QC_bpi1.tif", overwrite = T)
QC_bpi2 <- aggregate(QC_bpi2, fact = 4, fun = mean, na.rm = T,filename = "Data/Rasters/QC_bpi2.tif", overwrite = T)
QC_bpi5 <- aggregate(QC_bpi5, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/QC_bpi5.tif", overwrite = T)
QC_bpi10 <- aggregate(QC_bpi10, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/QC_bpi10.tif", overwrite = T)
QC_bpi20 <- aggregate(QC_bpi20, fact = 4, fun = mean, na.rm = T, filename = "Data/Rasters/QC_bpi20.tif", overwrite = T)

#align extent with other rasters
bpi_list <- list.files('Data/Rasters/', full.names = T, pattern  = 'QC.+bpi')
bpi_stack <- stack(bpi_list)
bpi_stack <- crop(bpi_stack, extent(22174.67,1026175,5192424,5820424))
bpi_names <- paste(Rpth, names(bpi_stack),'.tif', sep = '')
writeRaster(bpi_stack, filename = bpi_names, bylayer = T, overwrite = T)


###########################
####Chl a @ sea surface####
###########################

#Read in extended NWA chl layers, crop and resample
#to Quebec region @ 1km res, aggregate to 4 km and write raster
#read temp _> reproject temp -> aggregate reprojected

QC_SA <- readOGR("Data/Shapefiles/QC_StudyArea.shp")
r1000 <- raster(extent(QC_SA),res = 1000, crs = proj4string(QC_SA))
Rpth <- "Data/Rasters/"

chl_files <- list.files(path = "Data/Rasters/", full.names = F, pattern = 'NWA.+chl')
QC_chl  <- chl_files %>% 
  map(~ raster(paste(Rpth, .x, sep = ''))) %>% 
  map(~ projectRaster(.x, r1000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','QC',names(.x)),'.tif',sep = ''), overwrite = T)) %>%
  map(~ aggregate(.x, fact = 4, fun = mean, filename = paste(Rpth,names(.x),'.tif',sep = ''), overwrite = T))

#align extent with other rasters
chl_list <- list.files('Data/Rasters/', full.names = T, pattern  = 'QC.+chl')
chl_stack <- stack(chl_list)
chl_stack <- crop(chl_stack, extent(22174.67,1026175,5192424,5820424))
chl_names <- paste(Rpth, names(chl_stack),'.tif', sep = '')
writeRaster(chl_stack, filename = chl_names, bylayer = T, overwrite = T)

###########################
####Oceanographic Data ####
###########################

#Functions and objects required to convert pt shapefiles of environmental
#variable summaries to raster format (4km resolution)

source("Code/VoronoiResample.R")

Spth <- "Data/Shapefiles/"
Rpth <- "Data/Rasters/"
extentQC <- extent(-70,-55,46,52.75)
QC_SA <- readOGR("Data/Shapefiles/QC_StudyArea.shp")

###mixed layer depth#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions
# and objects called by 'mixed layer depth' chunk)


#Convert Pt shapefile for NWA to raster with 4 km res for QC region

source("Code/Voronoipolygons2.R") # voronoi (thiessen) polygon function
input.file <- "mld_NWA_pts.shp"

NWA_pts <- readOGR(paste(Spth,input.file,sep='')) #read in pt shapefile
summary_stats <- colnames(NWA_pts@data)
QC_pts <- crop(NWA_pts, extentQC) #crop pts to generous extent for Quebec
QC_pts <- spTransform(QC_pts,proj4string(QC_SA)) #reproject to UTM 20
vor_QC <- voronoipolygons(QC_pts) # create voronoi polygons for each point
proj4string(vor_QC) <- proj4string(QC_SA)

# create fishnet grid of 4km X 4km polygons spanning QC study area
r.template <- raster(extent(QC_SA), res = 4000, crs = proj4string(QC_SA))
fishnet <- rasterToPolygons(r.template)
fishnet$layer <- c(1:length(fishnet$layer))

# convert voronoi and fishnet polygon df's to sf objects
vor_QC <- st_as_sf(vor_QC)
fishnet <- st_as_sf(fishnet)

# spatial join of fishnet grid with voronoi polygons
#group by fishnet polygon id and aggregate
QC_4km <- st_join(fishnet, vor_QC) %>% 
  group_by(layer) %>% 
  summarise_at(summary_stats, mean, na.rm = T) %>% 
  as(.,'Spatial')

#rasterize fishnet grid & save rasters

QC_mn_ann_mld <- rasterize(QC_4km, r.template, field = QC_4km@data$mn_mld_ann, filename = "Data/Rasters/QC_mn_ann_mld.tif", overwrite = T)
QC_mn_sum_mld <- rasterize(QC_4km, r.template, field = QC_4km@data$mn_mld_sum, filename = "Data/Rasters/QC_mn_sum_mld.tif", overwrite = T)

####Bottom Stress####  

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Stress' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for QC region

VoronoiResample(Spth,Rpth,"BtmStr_NWA_pts.shp",extentQC,QC_SA,"QC",4000,"BtmStr")

####SST#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'SST' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for QC region

VoronoiResample(Spth,Rpth,"SST_NWA_pts.shp",extentQC,QC_SA,"QC",4000,"SST")

####Bottom Temperature####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Temperature' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for QC region

VoronoiResample(Spth,Rpth,"BT_NWA_pts.shp",extentQC,QC_SA,"QC",4000,"BT")

####Bottom Salinity####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Salinity' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for QC region

VoronoiResample(Spth,Rpth,"BSal_NWA_pts.shp",extentQC,QC_SA,"QC",4000,"BSal")

####VelEW####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelEW' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for QC region

VoronoiResample(Spth,Rpth,"VelEW_NWA_pts.shp",extentQC,QC_SA,"QC",4000,"VelEW")

####VelNS####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelNS' chunk)

#Convert Pt shapefile for NWA to raster with 4 km res for QC region

VoronoiResample(Spth,Rpth,"VelNS_NWA_pts.shp",extentQC,QC_SA,"QC",4000,"VelNS")

############################
####Primary Productivity####
############################

#Read in extended NWA PP layers,
#Crop to Quebec region extent, reproject to UTM 20 and resample to 4km resolution

QC_SA <- readOGR("Data/Shapefiles/QC_StudyArea.shp")
r4000 <- raster(extent(QC_SA),res = 4000, crs = proj4string(QC_SA))
Rpth <- "Data/Rasters/"

PP_files <- list.files(path = Rpth, full.names = F, pattern = 'NWA.+PP')

QC_PP <- PP_files %>% 
  map(~ raster(paste(Rpth,.x, sep = ''))) %>%
  map(~ projectRaster(.x, r4000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','QC',names(.x)),'.tif',sep = ''), overwrite = T))

####DO,pH,Nutrients####

QC_SA <- readOGR("Data/Shapefiles/QC_StudyArea.shp") #QC Study Area shapefile
r4000 <- raster(extent(QC_SA),res = 4000, crs = proj4string(QC_SA)) # template raster in UTM 20N with 4 km resolution
Rpth <- "Data/Rasters/" #Raster file path
NWA.ext <- extent(-70.9, -42, 41, 58.6) #Extent object for Northwest Atl.
LandBorders <- readOGR("Data/Shapefiles/QC_landborders.shp") #Shapefile w/ provincial boundaries
LandBuffer <- gBuffer(LandBorders, width = 5000) #create 5 km buffer around land borders
writeOGR(as(LandBuffer,'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'QC_LandBuffer_5km', driver = 'ESRI Shapefile')
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
  map(~ mask(.x, LandBuffer, inverse = T, filename = paste(Rpth,'QC_',gsub('BO_','',names(.x)),'.tif',sep = ''), overwrite = T))
