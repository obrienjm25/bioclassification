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
CHS15s <- raster("Data/Rasters/Bathymetry/CHS_bathymetry_15sec.tif")
NWA_bathy <- calc(CHS15s, fun = function(x){x * -1}, "Data/Rasters/Bathymetry/NWA_bathymetry_15sec_CHS.tif", overwrite = T)

#Maritimes

#crop to Maritimes region study area and reproject/resample to UTM 20, Datum: NAD83, 500 m resolution

MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")
r500 <- raster(extent(MaritimesSA),res = 500, crs = proj4string(MaritimesSA))
Mar_bathy500m <- projectRaster(NWA_bathy,r500, method = 'bilinear', filename = "Data/Rasters/Mar_bathy500m.tif", overwrite = T)
Mar_bathy <- aggregate(Mar_bathy500m, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bathy.tif", overwrite = T)

#derive other terrain characteristics from bathymetric data
Mar_Slope <- terrain(Mar_bathy500m, opt = 'slope', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/Mar_slope.tif', overwrite = T) #slope (degrees using Horn algorithm)
Mar_Slope <- aggregate(Mar_Slope, fact = 8, fun = mean, filename = 'Data/Rasters/Mar_slope.tif', overwrite = T)
Mar_Aspect <- terrain(Mar_bathy500m, opt = 'aspect', neighbors = 8, unit = 'degrees', filename = 'Data/Rasters/Mar_aspect.tif', overwrite = T) #Aspect (degrees from N)
Mar_Aspect <- aggregate(Mar_Aspect, fact = 8, fun = mean, filename = 'Data/Rasters/Mar_aspect.tif', overwrite = T)

#calculate various broad and fine scale bathymetric position indices

#matrix defining window/scale of BPI
f0.5 <- matrix(1, nrow = 3, ncol = 3) #8 neareast neighbours (NN), radius of 1 grid cell, 0.5km
f1 <- matrix(1, nrow = 5, ncol = 5) #24 NN, radius of 2 grid cells, 1 km
f2 <- matrix(1, nrow = 9, ncol = 9) #80 NN, radius of 4 grid cells, ~ km
f5 <- matrix(1, nrow = 21, ncol = 21) #440 NN, radius of 10 grid cells, 5 km
f10 <- matrix(1, nrow = 41, ncol = 41) #1681 NN, radius of 20 grid cells, 10 km
f20 <- matrix(1, nrow = 81, ncol = 81) #6561 NN, radius of 40 grid cells, 20 km

#focal function taking difference between depth in cell and mean of depth in window of given radius
Mar_bpi0.5 <- focal(Mar_bathy500m, w = f0.5, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Mar_bpi0.5.tif", overwrite = T)
Mar_bpi1 <- focal(Mar_bathy500m, w = f1, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Mar_bpi1.tif", overwrite = T)
Mar_bpi2 <- focal(Mar_bathy500m, w = f2, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Mar_bpi2.tif", overwrite = T)
Mar_bpi5 <- focal(Mar_bathy500m, w = f5, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Mar_bpi5.tif", overwrite = T)
Mar_bpi10 <- focal(Mar_bathy500m, w = f10, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Mar_bpi10.tif", overwrite = T)
Mar_bpi20 <- focal(Mar_bathy500m, w = f20, fun = function(x,...) x[(length(x)+1)/2] - mean(x[-(length(x)+1)/2], na.rm = T), pad = TRUE, padValue = NA, filename = "Data/Rasters/Mar_bpi20.tif", overwrite = T)

#aggregate bpi indices at 500m res to 4km
Mar_bpi0.5 <- aggregate(Mar_bpi0.5, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bpi0.5.tif", overwrite = T)
Mar_bpi1 <- aggregate(Mar_bpi1, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bpi1.tif", overwrite = T)
Mar_bpi2 <- aggregate(Mar_bpi2, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bpi2.tif", overwrite = T)
Mar_bpi5 <- aggregate(Mar_bpi5, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bpi5.tif", overwrite = T)
Mar_bpi10 <- aggregate(Mar_bpi10, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bpi10.tif", overwrite = T)
Mar_bpi20 <- aggregate(Mar_bpi20, fact = 8, fun = mean, filename = "Data/Rasters/Mar_bpi20.tif", overwrite = T)

###########################
####Chl a @ sea surface####
###########################

#Create vectors of indexing variables
year.MODIS <- as.character(2007:2011)
year.viirs <- as.character(2012:2016)
month <- c('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')
winter <- month[1:3]
spring <- month[4:6]
summer <- month[7:9]
fall <- month[10:12]

#create empty RasterStacks for each summary variable
mn_win_chl <- stack()
mn_spr_chl <- stack()
mn_sum_chl <- stack()
mn_fal_chl <- stack()
mn_ann_chl <- stack()
max_ann_chl <- stack()
min_ann_chl <- stack()
ran_ann_chl <- stack()

chl.path <- "U:/MPA_group/Data/chl_a_2007_2016/" #path to subfolder on R drive for chl_a data

#Loop through years covered reliably by AquaModis and Viirs respectively
#get names of all files in modis and viirs ftp folders
#Download all composites for that year to subfolder on R drive
#read into R as RasterStack
#replace all -1's with NA, and cap chl a values @ 7 ug/mL
#compute seasonal and annual summary statistics from stack layers for each gridcell
#add summary layers for each year to respective RasterStack for each summary variable
#write to file
#remove composite geotiff files from R: drive subfolder

for (year in c(year.MODIS,year.viirs)){
  
  if(year %in% year.MODIS){
    ftp <- "ftp://ftp.dfo-mpo.gc.ca/bometrics/modis/extNA_geotiff/"
  } else ftp <- "ftp://ftp.dfo-mpo.gc.ca/bometrics/viirs/extNA_geotiff/"
  
  allfiles <- getURL(ftp, dirlistonly=T)
  allfiles <- strsplit(allfiles, "\r*\n")[[1]]
  filenames <- allfiles[grepl(year,allfiles)]
  
  for (i in filenames){
    download.file(paste(ftp,i,sep = ""), paste(chl.path,i,sep=""), method = 'auto', mode = 'wb', quiet = T)
  }
  
  list.path <- list.files(path = chl.path, full.names = T,
                        pattern = "\\.tif$")
  
  s <- stack(list.path, quick = T)
  NAvalue(s) <- -1
  
  mn_win_chl <- stack(mn_win_chl, mean(s[[grep(paste(winter, collapse = "|"), s@layers)]], na.rm = T))
  mn_spr_chl <- stack(mn_spr_chl, mean(s[[grep(paste(spring, collapse = "|"), s@layers)]], na.rm = T))
  mn_sum_chl <- stack(mn_sum_chl, mean(s[[grep(paste(summer, collapse = "|"), s@layers)]], na.rm = T))
  mn_fal_chl <- stack(mn_fal_chl, mean(s[[grep(paste(fall, collapse = "|"), s@layers)]], na.rm = T))
  mn_ann_chl <- stack(mn_ann_chl, mean(s, na.rm = T))
  max_ann_chl <- stack(max_ann_chl, max(s, na.rm = T))
  min_ann_chl <- stack(min_ann_chl, min(s, na.rm = T))
  ran_ann_chl <- stack(ran_ann_chl, max(s, na.rm = T) - min(s, na.rm = T))

  rm(s)

  for (i in filenames){
    file.remove(paste(chl.path,i,sep=''))
  }
}

if(max(cellStats(max_ann_chl, max)) >= 7000){
  print("WARNING: suspiciously high chl a values")
} else print("HURRAY!!!")

#create layers that take average of summary variables across years using
#map function to write raster file for NWA, crop and resample
#to maritimes @ 1km res, aggregate to 4 km and write raster
#write temp _> reproject temp -> aggregate reprojected

MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")
r1000 <- raster(extent(MaritimesSA),res = 1000, crs = proj4string(MaritimesSA))
Rpth <- "Data/Rasters/"

chl_varlist <- c('mn_win_chl', 'mn_spr_chl', 'mn_sum_chl', 'mn_fal_chl', 'mn_ann_chl', 'max_ann_chl', 'min_ann_chl', 'ran_ann_chl')

chl_varlist %>% 
  map(~ writeRaster(mean(get(.x), na.rm = T), filename = paste(Rpth,'NWA_avg_',.x,'.tif',sep = ''), overwrite = T)) %>%
  map(~ projectRaster(.x, r1000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','Mar',names(.x)),'.tif',sep = ''), overwrite = T)) %>%
  map(~ aggregate(.x, fact = 4, fun = mean, filename = paste(Rpth,names(.x),'.tif',sep = ''), overwrite = T))
  
###########################
####Oceanographic Data ####
###########################

#grouping variables
months <- c(paste("M", c(1:12), sep = "_"))
layers <- c("Bottom_UVTS","BottomStress","MLD","Surface_UVTS")
vars <- c('VelEW', 'VelNS','BT','BSal','BtmStr','mld','SST')
year <- as.character(c(2007:2015))
winter <- months[1:3]
spring <- months[4:6]
summer <- months[7:9]
fall <- months[10:12]

#functions to be used to calculate annual summaries
rowMin <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  min = min(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11))
  return(min)}

rowMax <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  max = max(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11))
  return(max)}

rowMax2 <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  max = max(abs(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11)))
  return(max)} #for variables where +/- indicate direction

#Objects and functions required to convert pt shapefiles of environmental
#variable summaries to raster format (4km resolution)

source("Code/VoronoiResample.R")

Spth <- "Data/Shapefiles/"
Rpth <- "Data/Rasters/"
extentMar <- extent(-68,-55,41,48)
MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")

#create list, each element containing a list of textfile names for each environmental variable of interest
#Data provided as monthly predictions for pts on irregular grid extending from SW NS halfway up Labrador (2007-2015) 

bnam.files <- list()

for(i in 1:length(layers)){
  bnam.files[[i]] <- list.files(path = "Data/Zeliang_OceanographicData/", full.names = T,
                                pattern = layers[i])
  names(bnam.files)[i] <- layers[i]
}

####mixed layer depth#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'mixed layer depth' chunk)

#read in data for all years, append into one large table, reshape to months as cols,
#compute annual and seasonal summary variables for each point from month cols
#Group by year to get average for each variable across years
#write to pt shapefile

filename <- "mld_NWA_pts"
varname <- "mld"
df.long <- bnam.files[['MLD']] %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20..)", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T)))
 
df.long <- data.table::rbindlist(df.long) %>% data.frame()#14063868 observations of 5 variables 
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>% 
  separate(.,lonlat,c("x","y"), sep = "/") %>%
  mutate(mn_mld_ann = rowMeans(select(.,starts_with("M_"))),mn_mld_sum = rowMeans(select(.,one_of(summer)))) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

source("Code/Voronoipolygons2.R") # voronoi (thiessen) polygon function
input.file <- "mld_NWA_pts.shp"

NWA_pts <- readOGR(paste(Spth,input.file,sep='')) #read in pt shapefile
summary_stats <- colnames(NWA_pts@data)
Mar_pts <- crop(NWA_pts, extentMar) #crop pts to generous extent for Maritimes
Mar_pts <- spTransform(Mar_pts,proj4string(MaritimesSA)) #reproject to UTM 20
vor_Mar <- voronoipolygons(Mar_pts) # create voronoi polygons for each point
proj4string(vor_Mar) <- proj4string(MaritimesSA)

# create fishnet grid of 4km X 4km polygons spanning Maritimes study area
r.template <- raster(extent(MaritimesSA), res = 4000, crs = proj4string(MaritimesSA))
fishnet <- rasterToPolygons(r.template)
fishnet$layer <- c(1:length(fishnet$layer))

# convert voronoi and fishnet polygon df's to sf objects
vor_Mar <- st_as_sf(vor_Mar)
fishnet <- st_as_sf(fishnet)

# spatial join of fishnet grid with voronoi polygons
#group by fishnet polygon id and aggregate
Mar_4km <- st_join(fishnet, vor_Mar) %>% 
  group_by(layer) %>% 
  summarise_at(summary_stats, mean, na.rm = T) %>% 
  as(.,'Spatial')

#rasterize fishnet grid & save rasters

Mar_mn_ann_mld <- rasterize(Mar_4km, r.template, field = Mar_4km@data$mn_mld_ann, filename = "Data/Rasters/Mar_mn_ann_mld.tif", overwrite = T)
Mar_mn_sum_mld <- rasterize(Mar_4km, r.template, field = Mar_4km@data$mn_mld_sum, filename = "Data/Rasters/Mar_mn_sum_mld.tif", overwrite = T)

####Bottom Stress####  

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Stress' chunk)

filename <- "BtmStr_NWA_pts"
varname <- "BtmStr"
df.long <- bnam.files[['BottomStress']] %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20..)", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T)))

df.long <- data.table::rbindlist(df.long) %>% data.frame()#14063868 observations of 5 variables 
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(mn_ann = rowMeans(select(.,starts_with("M_")))) %>%   
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

VoronoiResample(Spth,Rpth,"BtmStr_NWA_pts.shp",extentMar,MaritimesSA,"Mar",4000,"BtmStr")

####SST#### 

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'SST' chunk)

filename <- "SST_NWA_pts"
varname <- "SST"
df.long <- bnam.files[['Surface_UVTS']] %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20..)", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T)))

df.long <- data.table::rbindlist(df.long) %>% data.frame() %>%   
  select(-c('V3','V4','V6'))#14063868 observations of 5 variables 
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(mn_ann = rowMeans(select(.,starts_with("M_"))),mn_win = rowMeans(select(.,one_of(winter))),
         mn_sum = rowMeans(select(.,one_of(summer))),mn_spr = rowMeans(select(.,one_of(spring))),
         mn_fall = rowMeans(select(.,one_of(fall))),max_ann = pmap_dbl(.,rowMax),
         min_ann = pmap_dbl(.,rowMin),range_ann = max_ann - min_ann) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

VoronoiResample(Spth,Rpth,"SST_NWA_pts.shp",extentMar,MaritimesSA,"Mar",4000,"SST")

####Bottom Temperature####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Temperature' chunk)

filename <- "BT_NWA_pts"
varname <- "BT"
df.long <- bnam.files[['Bottom_UVTS']] %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20..)", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T)))

df.bottom <- data.table::rbindlist(df.long) %>% data.frame()
df.long <- df.bottom %>% select(-c('V3','V4','V6'))#14063868 observations of 5 variables 
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(mn_ann = rowMeans(select(.,starts_with("M_"))),mn_win = rowMeans(select(.,one_of(winter))),
         mn_sum = rowMeans(select(.,one_of(summer))),mn_spr = rowMeans(select(.,one_of(spring))),
         mn_fall = rowMeans(select(.,one_of(fall))),max_ann = pmap_dbl(.,rowMax),
         min_ann = pmap_dbl(.,rowMin),range_ann = max_ann - min_ann) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

VoronoiResample(Spth,Rpth,"BT_NWA_pts.shp",extentMar,MaritimesSA,"Mar",4000,"BT")

####Bottom Salinity####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'Bottom Salinity' chunk)

filename <- "BSal_NWA_pts"
varname <- "BSal"

df.long <- df.bottom %>% select(-c('V3','V4','V5'))#14063868 observations of 5 variables  
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(mn_ann = rowMeans(select(.,starts_with("M_"))),mn_win = rowMeans(select(.,one_of(winter))),
         mn_sum = rowMeans(select(.,one_of(summer))),mn_spr = rowMeans(select(.,one_of(spring))),
         mn_fall = rowMeans(select(.,one_of(fall))),max_ann = pmap_dbl(.,rowMax),
         min_ann = pmap_dbl(.,rowMin),range_ann = max_ann - min_ann) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

VoronoiResample(Spth,Rpth,"BSal_NWA_pts.shp",extentMar,MaritimesSA,"Mar",4000,"BSal")

####VelEW####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelEW' chunk)

filename <- "VelEW_NWA_pts"
varname <- "VelEW"

df.long <- df.bottom %>% select(-c('V4','V5','V6'))#14063868 observations of 5 variables  
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>% 
  separate(.,lonlat,c("x","y"), sep = "/") %>%
  mutate(mn_ann = rowMeans(select(.,starts_with("M_"))),max_ann = pmap_dbl(.,rowMax2)) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

VoronoiResample(Spth,Rpth,"VelEW_NWA_pts.shp",extentMar,MaritimesSA,"Mar",4000,"VelEW")

####VelNS####

#README!! 
#(Run code in 'Oceanographic Data' chunk first for functions, grouping variables,
# and objects called by 'VelNS' chunk)

filename <- "VelNS_NWA_pts"
varname <- "VelNS"

df.long <- df.bottom %>% select(-c('V3','V5','V6'))#14063868 observations of 5 variables  
names(df.long) <- c("x","y",varname,"year","month")
df.long[which(df.long[3]==0),3] <- NA
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>% 
  separate(.,lonlat,c("x","y"), sep = "/") %>%
  mutate(mn_ann = rowMeans(select(.,starts_with("M_"))),max_ann = pmap_dbl(.,rowMax2)) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame()

coordinates(df.wide) <- ~x+y
proj4string(df.wide) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
writeOGR(df.wide, dsn = "Data/Shapefiles", layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

#Pt shapefile to raster with 4 km res

VoronoiResample(Spth,Rpth,"VelNS_NWA_pts.shp",extentMar,MaritimesSA,"Mar",4000,"VelNS")

############################
####Primary Productivity####
############################

#Create vectors of indexing variables
year.MODIS <- as.character(2007:2011)
year.viirs <- as.character(2012:2016)
spr_sum <- c('04_','05_','06_','07_','08_','09_')


#create empty RasterStacks for each summary variable
mn_ann_PP <- stack()
mn_sprsum_PP <- stack()
max_ann_PP <- stack()
max_sprsum_PP <- stack()

PP.path <- "Data/Rasters/PP/" #path to subfolder on R drive for PP data

#Loop through years covered reliably by AquaModis and Viirs respectively
#get names of all files in modis and viirs ftp folders
#Download all monthly composites for that year to subfolder on R drive
#read into R as RasterStack
#missing value flag is -9
#compute seasonal and annual summary statistics from stack layers for each gridcell
#add summary layers for each year to respective RasterStack for each summary stat
#remove composite geotiff files from R: drive subfolder
#write to file

for (year in c(year.MODIS,year.viirs)){
  
  if(year %in% year.MODIS){
    ftp <- "ftp://ftp.dfo-mpo.gc.ca/bometrics/PP/2002-2014_monthly_geotiffs/"
  } else ftp <- "ftp://ftp.dfo-mpo.gc.ca/bometrics/PP/2012-2016_monthly_geotiffs/"
  
  allfiles <- getURL(ftp, dirlistonly=T)
  allfiles <- strsplit(allfiles, "\r*\n")[[1]]
  filenames <- allfiles[grepl(year,allfiles)]
  
  for (i in filenames){
    download.file(paste(ftp,i,sep = ""), paste(PP.path,i,sep=""), method = 'auto', mode = 'wb', quiet = T)
  }
  
  list.path <- list.files(path = PP.path, full.names = T,
                          pattern = "\\.tif$")
  
  s <- stack(list.path, quick = T)
  NAvalue(s) <- -9
  
  mn_ann_PP <- stack(mn_ann_PP, mean(s, na.rm = T))
  mn_sprsum_PP <- stack(mn_sprsum_PP, mean(s[[grep(paste(spr_sum, collapse = "|"), s@layers)]], na.rm = T))
  max_ann_PP <- stack(max_ann_PP, max(s, na.rm = T))
  max_sprsum_PP <- stack(max_sprsum_PP, max(s[[grep(paste(spr_sum, collapse = "|"), s@layers)]], na.rm = T))
  
  rm(s)
  
  for (i in filenames){
    file.remove(paste(PP.path,i,sep=''))
  }
}

#Take average of summary stats across years for NWA
#Crop to Maritimes, reproject to UTM 20 and resample to 4km resolution

MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")
r4000 <- raster(extent(MaritimesSA),res = 4000, crs = proj4string(MaritimesSA))
Rpth <- "Data/Rasters/"

PP_varlist <- c('mn_ann_PP', 'mn_sprsum_PP', 'max_ann_PP', 'max_sprsum_PP')

PP_varlist %>% 
  map(~ writeRaster(mean(get(.x), na.rm = T), filename = paste(Rpth,'NWA_avg_',.x,'.tif',sep = ''), overwrite = T)) %>%
  map(~ projectRaster(.x, r4000, method = 'bilinear', filename = paste(Rpth,gsub('NWA','Mar',names(.x)),'.tif',sep = ''), overwrite = T))

#####################################
####NRCan Habitat Template Layers####
#####################################

####Scope for Growth####

#read xyz file
SFG.pts <- fread('Data/NRCan_HabitatTemplate/ScopeForGrowthXYZ.txt', skip = 4, col.names = c('x','y','SFG.scaled'), data.table = F)
coordinates(SFG.pts) <- ~ x+y # create SPDF from xyz file
#set CRS
proj4string(SFG.pts) <- CRS('+proj=tmerc +datum=WGS84 +units=m +no_defs +ellps=WGS84 +lon_0=-63 +lat_0=0 +k_0=0.9996 +x_0=500000 +y_0=0')
#convert SPDF to raster
r1 <- raster(res = 500, crs = proj4string(SFG.pts), ext = extent(SFG.pts) + 500)
SFG.raster <- rasterize(SFG.pts,r1, field = 'SFG.scaled') 
SFG.raster <- aggregate(SFG.raster, fact = 8) #aggregate to 4 km resolution
#Change projection to UTM Zone 20N & align with other environmental layers
template.raster <- raster('Data/Rasters/Mar_mn_ann_BT.tif') 
SFG.raster <- projectRaster(from = SFG.raster, to = template.raster, method = 'bilinear', filename = 'Data/Rasters/Mar_ScopeforGrowth.tif', overwrite = T)
compareRaster(SFG.raster,template.raster) # double check alignment

####Natural Disturbance####

#read xyz file
NatDist.pts <- fread('Data/NRCan_HabitatTemplate/DisturbanceXYZ.txt', skip = 4, col.names = c('x','y','ND.scaled'), data.table = F)
coordinates(NatDist.pts) <- ~ x+y # create SPDF from xyz file
#set CRS
proj4string(NatDist.pts) <- CRS('+proj=tmerc +datum=WGS84 +units=m +no_defs +ellps=WGS84 +lon_0=-63 +lat_0=0 +k_0=0.9996 +x_0=500000 +y_0=0')
#convert SPDF to raster
r1 <- raster(res = 500, crs = proj4string(NatDist.pts), ext = extent(NatDist.pts) + 500)
NatDist.raster <- rasterize(NatDist.pts,r1, field = 'ND.scaled') 
NatDist.raster <- aggregate(NatDist.raster, fact = 8) #aggregate to 4 km resolution
#Change projection to UTM Zone 20N & align with other environmental layers
template.raster <- raster('Data/Rasters/Mar_mn_ann_BT.tif') 
NatDist.raster <- projectRaster(from = NatDist.raster, to = template.raster, method = 'bilinear', filename = 'Data/Rasters/Mar_Disturbance.tif')
compareRaster(NatDist.raster,template.raster) # double check alignment

####Grain Size####

grain <- raster('Data/Rasters/grainsizemm/prj.adf') #500m res, UTM20, ellps=GRS80
grain2 <- raster:: aggregate(grain, fact = 8) #aggregate to 4 km resolution

#align extent and origin with other environmental rasters
template.raster <- raster('Data/Rasters/Mar_mn_ann_BT.tif') 
grain2 <- resample(grain2,template.raster, method = 'bilinear', filename = 'Data/Rasters/Mar_GrainSize_mm.tif', overwrite = T)
compareRaster(grain2, template.raster)# double check alignment


####DO,pH,Nutrients####

MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp") #Maritimes Study Area shapefile
r4000 <- raster(extent(MaritimesSA),res = 4000, crs = proj4string(MaritimesSA)) # template raster in UTM 20N with 4 km resolution
Rpth <- "Data/Rasters/" #Raster file path
NWA.ext <- extent(-70.9, -47.3, 41, 58.6) #Extent object for Northwest Atl.
LandMaritimes <- readOGR("Data/Shapefiles/LandBordersMaritimes.shp") #Shapefile w/ NS, NB, and ME landborders
canada <- raster::getData("GADM", country = "CAN", level = 1)
QC.NFLD.PEI <- c("Newfoundland and Labrador", "Prince Edward Island", "Qu\u{e9}bec")
canada <- canada[canada$NAME_1 %in% QC.NFLD.PEI,] #spPolygons df w/ NFLD, QC, and PEI borders
canada <- spTransform(canada, CRS('+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
canada <- crop(canada, extent(MaritimesSA))
Land <- rbind(LandMaritimes, canada) #join land borders dataframes
LandBuffer <- gBuffer(Land, width = 5000) #create 5 km buffer around land borders
writeOGR(Land, dsn = 'Data/Shapefiles', layer = 'Maritimes_prov_borders', driver = 'ESRI Shapefile')
writeOGR(as(LandBuffer,'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'Maritimes_LandBuffer_5km', driver = 'ESRI Shapefile')

library(sdmpredictors)
sdm_datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)
sdm_layers <- list_layers(sdm_datasets, terrestrial=FALSE, marine=TRUE, monthly=TRUE,
                  version=NULL)
BO_list <- list('BO_dissox','BO_nitrate','BO_ph', 'BO_phosphate','BO_silicate') 
names(BO_list) <- c('BO_dissox','BO_nitrate','BO_ph', 'BO_phosphate','BO_silicate')
#download environmental predictor layers from Bio-Oracle database, crop to study area,
#disaggregate to smaller cell size, resample and reproject to match other raster layers,
#mask cells within 5 km radius of land, write to file
BO_layers <- BO_list %>% 
  map(~ load_layers(.x)) %>%
  map(~ crop(.x, NWA.ext)) %>%  
  map(~ disaggregate(.x, fact = 2, method = 'bilinear', na.rm = T)) %>%  
  map(~ projectRaster(.x, r4000, method = 'bilinear', na.rm = T)) %>% 
  map(~ mask(.x, LandBuffer, inverse = T, filename = paste(Rpth,'Mar_',names(.x),'.tif',sep = ''), overwrite = T))