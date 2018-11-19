##Code to create shapefiles outlining study areas within DFO regions and adjacent land borders

setwd("E:/Biological Classification")

library(rgeos)
library(rgdal)
library(raster)
library(sp)

####Maritimes Region ####

#Get land boundaries of provinces/states bordering maritimes region

canada <- getData("GADM", country = "CAN", level = 1)
USA <- getData("GADM", country = "USA", level = 1)
CanUS <- bind(canada, USA)
ProvStat <- c("Nova Scotia", "New Brunswick", "Maine")
LandBordersMaritimes <- CanUS[CanUS$NAME_1 %in% ProvStat,]
LandBordersMaritimes <- spTransform(LandBordersMaritimes, CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
writeOGR(LandBordersMaritimes, dsn = "Data/Shapefiles", layer = "LandBordersMaritimes", driver = "ESRI Shapefile")

#Get shapefiles to create polygon of gridded region
MaritimeRegion <- readOGR("Data/Shapefiles/MaritimesPlanningArea.shp") #maritime planning region
MaritimeProj <- spTransform(MaritimeRegion,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
strata <- readOGR("Data/Shapefiles/MaritimesRegionStrataBoundaries.shp") #RV survey strata boundaries
strataProj <- spTransform(strata,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#Create polygon of Maritime region restricted offshore by limits of RV survey strata boundaries
MaritimesStudyArea <- aggregate(cover(MaritimeProj, strataProj))
MaritimesStudyArea <- SpatialPolygonsDataFrame(MaritimesStudyArea, data.frame(Name = "Maritimes RV survey area", Area = NA))
writeOGR(MaritimesStudyArea, dsn = "Data/Shapefiles", layer = "MaritimesStudyArea", driver = "ESRI Shapefile")
