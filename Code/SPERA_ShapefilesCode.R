##Code to create shapefiles outlining study areas within DFO regions and adjacent land borders

library(rgeos)
library(rgdal)
library(raster)
library(sp)
library(maptools)

####Maritimes Region ####

#Get land boundaries of provinces/states bordering maritimes region

canada <- getData("GADM", country = "CAN", level = 1)
USA <- getData("GADM", country = "USA", level = 1)
CanUS <- bind(canada, USA)
ProvStat <- c("Nova Scotia", "New Brunswick", "Maine", "Qu\u{e9}bec")
LandBordersMaritimes <- CanUS[CanUS$NAME_1 %in% ProvStat,]
LandBordersMaritimes <- spTransform(LandBordersMaritimes, CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
writeOGR(LandBordersMaritimes, dsn = "Data/Shapefiles", layer = "LandBordersMaritimes", driver = "ESRI Shapefile", overwrite_layer = T)
LandBorders_cropped <- crop(LandBordersMaritimes, extent(c(-147090.6,1029071,4796509,5329882)))
writeOGR(LandBorders_cropped, dsn = "Data/Shapefiles", layer = "Maritimes_prov_borders", driver = "ESRI Shapefile", overwrite_layer = T)


#Get shapefiles to create polygon of gridded region
MaritimeRegion <- readOGR("Data/Shapefiles/MaritimesPlanningArea.shp") #maritime planning region
MaritimeProj <- spTransform(MaritimeRegion,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
strata <- readOGR("Data/Shapefiles/MaritimesRegionStrataBoundaries.shp") #RV survey strata boundaries
strataProj <- spTransform(strata,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#Create polygon of Maritime region restricted offshore by limits of RV survey strata boundaries
MaritimesStudyArea <- aggregate(cover(MaritimeProj, strataProj))
MaritimesStudyArea <- SpatialPolygonsDataFrame(MaritimesStudyArea, data.frame(Name = "Maritimes RV survey area", Area = NA))
writeOGR(MaritimesStudyArea, dsn = "Data/Shapefiles", layer = "MaritimesStudyArea", driver = "ESRI Shapefile")

#Get land boundaries of provinces bordering Newfoundland region
canada <- raster::getData("GADM", country = 'CAN', level = 1)
spm <- raster::getData("GADM", country = 'SPM', level = 1)
prov <- c("Newfoundland and Labrador", "Prince Edward Island", "Qu\u{e9}bec", "Nova Scotia", "New Brunswick")
NLregion <- canada[canada$NAME_1 %in% prov,]
NLregion <- rbind(NLregion,spm)
NLregion <- spTransform(NLregion, CRS('+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
RVstrataNLlonlat <- readOGR('Data/Shapefiles/NF_SamplingStrata_20140514.shp')
RVstrataNL <- spTransform(RVstrataNLlonlat, CRS('+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
NLregion <- crop(NLregion, (extent(RVstrataNL) + 50000))
writeOGR(NLregion, dsn = 'Data/Shapefiles', layer = 'NL_landborders', driver = 'ESRI Shapefile')

#Get land boundaries of provinces bordering Gulf region and create shapefile for Gulf study area
canada <- raster::getData("GADM", country = 'CAN', level = 1)
prov <- c("Newfoundland and Labrador", "Prince Edward Island", "Qu\u{e9}bec", "Nova Scotia", "New Brunswick")
Gulfregion <- canada[canada$NAME_1 %in% prov,]
Gulfregion <- spTransform(Gulfregion, CRS('+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
load('Data/Shapefiles/SGSL strata borders.rda')
RVstrataGulf <- PolySet2SpatialPolygons(strat.S)
RVstrataGulf <- spTransform(RVstrataGulf, CRS('+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
RVstrataGulf <- as(RVstrataGulf, 'SpatialPolygonsDataFrame')
RVstrataGulf@data$ID <- 1:length(RVstrataGulf)
Gulfregion <- crop(Gulfregion, (extent(RVstrataGulf) + 50000))
writeOGR(Gulfregion, dsn = 'Data/Shapefiles', layer = 'Gulf_landborders', driver = 'ESRI Shapefile')
RVstrataGulfAgg <- gBuffer(RVstrataGulf, byid = F, width = 0)
writeOGR(as(RVstrataGulfAgg, 'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'Gulf_RVsurveyAgg', driver = 'ESRI Shapefile')

#Get land boundaries of provinces bordering QC region and create shapefile for QC study area
canada <- raster::getData("GADM", country = 'CAN', level = 1)
USA <- raster::getData('GADM', country = 'USA', level = 1)
spm <- raster::getData("GADM", country = 'SPM', level = 1)
prov_stat <- c("Maine", "Newfoundland and Labrador", "Prince Edward Island", "Qu\u{e9}bec", "Nova Scotia", "New Brunswick")
QCregion <- bind(canada,USA)
QCregion <- QCregion[QCregion$NAME_1 %in% prov_stat,]
QCregion <- bind(QCregion, spm)
QCregion <- spTransform(QCregion, CRS('+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
load('Data/Shapefiles/NGSL strata borders.rda')
RVstrataQC <- PolySet2SpatialPolygons(strat.N)
RVstrataQC <- spTransform(RVstrataQC, CRS('+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
RVstrataQC <- as(RVstrataQC, 'SpatialPolygonsDataFrame')
RVstrataQC@data$ID <- 1:length(RVstrataQC)
QCregion <- crop(QCregion, (extent(RVstrataQC) + 50000))
#writeOGR(QCregion, dsn = 'Data/Shapefiles', layer = 'QC_landborders', driver = 'ESRI Shapefile')
#writeOGR(as(RVstrataQC, 'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'QC_RVsurvey', driver = 'ESRI Shapefile')
