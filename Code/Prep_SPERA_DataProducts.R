###########################################################################################
##############UPDATE ABIOTIC LAYERS USED IN SPERA BIOCLASSIFICATION######################## 
##############TO BE MADE AVAILABLE AS DELIVERABLE DATA PRODUCTS############################
###########################################################################################

library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(dplyr)
library(tidyr)
library(purrr)

####GULF####

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Gulf', full.names = T) #get list of raster files
Gulf_SA <- readOGR('Data/Shapefiles/Gulf_RVsurveyAgg.shp') # Gulf_SA (RV survey boundaries)
LandBuffer <- readOGR('Data/Shapefiles/Gulf_LandBuffer_5km.shp') # 5km buffer around land points
env_pred <- stack(raster.list) %>% #stack rasters in list
  `names<-`(.,gsub('Gulf_','',names(.))) %>% #remove 'Gulf_' from layer names
  mask(., Gulf_SA) %>% #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) #apply 5 km land buffer

writeRaster(env_pred, 'U:/MPA_group/Data/SPERA_DataProducts/sGOSL', suffix = names(env_pred), 
  overwrite = T, bylayer = T, format = 'GTiff') #write multi-layer object into separate files


####MARITIMES####

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Mar', full.names = T) #get list of raster files
RV_strata <- readOGR('Data/Shapefiles/MaritimesRegionStrataAgg.shp') #Maritimes RV survey boundaries
LandBuffer <- readOGR('Data/Shapefiles/Maritimes_LandBuffer_5km.shp') # 5km buffer around land points
env_pred <- stack(raster.list) %>%  #stack rasters in list
  `names<-`(.,gsub('Mar_','',names(.))) %>% #remove 'Mar_' from layer names
  mask(., RV_strata) %>%  #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) #also mask cells overlapped by 5km land buffer

writeRaster(env_pred, 'U:/MPA_group/Data/SPERA_DataProducts/Maritimes', suffix = names(env_pred), 
            overwrite = T, bylayer = T, format = 'GTiff') #write multi-layer object into separate files

####Newfoundland####

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'NL', full.names = T) #get list of raster files
NL_SA <- readOGR('Data/Shapefiles/NL_RVsurveyAgg.shp') # NL_SA (RV survey boundaries)
LandBuffer <- readOGR('Data/Shapefiles/NL_LandBuffer_5km.shp') # 5km buffer around land points
env_pred <- stack(raster.list) %>% #stack rasters in list
  `names<-`(.,gsub('NL_','',names(.))) %>% #remove 'NL_' from layer names
  mask(., NL_SA) %>% #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) #also mask cells overlapped by 5km land buffer

writeRaster(env_pred, 'U:/MPA_group/Data/SPERA_DataProducts/NL', suffix = names(env_pred), 
            overwrite = T, bylayer = T, format = 'GTiff') #write multi-layer object into separate files

####QUEBEC####

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'QC', full.names = T) #get list of raster files
QC_SA <- readOGR('Data/Shapefiles/QC_StudyArea.shp') # QC_SA (RV survey boundaries)
LandBuffer <- readOGR('Data/Shapefiles/QC_LandBuffer_5km.shp') # 5km buffer around land points
QC_strata_1_5 <- readOGR('Data/Shapefiles/QC_Str1_5_agg.shp') #strata without survey data
env_pred <- stack(raster.list) %>%  #stack rasters in list
  `names<-`(.,gsub('QC_','',names(.))) %>% #remove 'QC_' from layer names
  mask(., QC_SA) %>% #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) %>% #also mask cells overlapped by 5km land buffer
  mask(., QC_strata_1_5, inverse = T) #also mask cells overlapped by unsampled strata

writeRaster(env_pred, 'U:/MPA_group/Data/SPERA_DataProducts/QC', suffix = names(env_pred), 
            overwrite = T, bylayer = T, format = 'GTiff') #write multi-layer object into separate files

#Edit XML files for Quebec 

library(XML)

#Get list of XML files
XML.list <- grep('QC_(?!PredClust).+\\.xml', list.files(path = 'U:/MPA_group/Data/SPERA_DataProducts/',full.names = T), perl=T, value=T)

maxnode <- "//metadata/eainfo/detailed/attr/attrdomv/rdom/rdommax" #location of node with range max
minnode <- "//metadata/eainfo/detailed/attr/attrdomv/rdom/rdommin" #location of node with range min
XML.trees <- map(XML.list, ~ xmlTreeParse(., useInternalNodes = T)) #read in XML trees for each QC raster layer

#replace range domain values in XML trees with min and max value in raster layers
map2(XML.trees,names(env_pred), ~ `xmlValue<-`(.x[[minnode]],value = as.character(cellStats(env_pred[[.y]],min, na.rm = T))))
map2(XML.trees,names(env_pred), ~ `xmlValue<-`(.x[[maxnode]],value = as.character(cellStats(env_pred[[.y]],max, na.rm = T))))

#write edited XML files
map2(XML.trees, names(env_pred), ~ saveXML(.x,file = paste0('U:/MPA_group/Data/SPERA_DataProducts/QC_',.y,'.tif.xml')))
