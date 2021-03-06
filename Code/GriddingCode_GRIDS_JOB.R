## Code to make a grid layer to group survey point observations

source("Code/GridFilter2.R")

library(rgeos)
library(rgdal)
library(raster)
library(sp)

####Maritimes####

MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")
LandMaritimes <- readOGR("Data/Shapefiles/LandBordersMaritimes.shp")

grids <- c(75000,70000,65000,60000,55000,50000,
           45000,40000,35000,30000,25000,20000,
           15000,10000,9000,8000,7000,6000,5000,
           4000,3000)

Griddata <- list()

for (i in 1:length(grids)){
  
  print(paste("working on grid ~ ",grids[i]/1000,"km resolution",sep=" "))
  
  tempgrid <- GridFilter2(MaritimesSA, resol=grids[i], LandMaritimes)
  
  Griddata[[i]] <- tempgrid
  plot(tempgrid,main=paste(grids[i]/1000,"km resolution",sep=" "))
  
  rm(tempgrid)
}

rm(list=setdiff(ls(), "Griddata"))

save.image("Data/MaritimeGridsWholeCell.RData")

####Newfoundland####

study_area <- readOGR('Data/Shapefiles/NL_RVsurveyAgg.shp')
Land <- readOGR("Data/Shapefiles/NL_landborders.shp")

grids <- c(75000,70000,65000,60000,55000,50000,
           45000,40000,35000,30000,25000,20000,
           15000,10000,9000,8000,7000,6000,5000,
           4000,3000)

Griddata <- list()

for (i in 1:length(grids)){
  
  print(paste("working on grid ~ ",grids[i]/1000,"km resolution",sep=" "))
  
  tempgrid <- GridFilter2(study_area, resol=grids[i], Land)
  
  Griddata[[i]] <- tempgrid
  
  rm(tempgrid)
}

rm(list=setdiff(ls(), "Griddata"))

save.image("Data/NLGridsWholeCell.RData")

####Gulf####

study_area <- readOGR('Data/Shapefiles/Gulf_RVsurveyAgg.shp')
Land <- readOGR("Data/Shapefiles/Gulf_landborders.shp")

grids <- c(75000,70000,65000,60000,55000,50000,
           45000,40000,35000,30000,25000,20000,
           15000,10000,9000,8000,7000,6000,5000,
           4000,3000)

Griddata <- list()

for (i in 1:length(grids)){
  
  print(paste("working on grid ~ ",grids[i]/1000,"km resolution",sep=" "))
  
  tempgrid <- GridFilter2(study_area, resol=grids[i], Land)
  
  Griddata[[i]] <- tempgrid
  
  rm(tempgrid)
}

rm(list=setdiff(ls(), "Griddata"))

save.image("Data/GulfGridsWholeCell.RData")

####Quebec####

study_area <- readOGR('Data/Shapefiles/QC_StudyArea.shp')
Land <- readOGR("Data/Shapefiles/QC_landborders.shp")

grids <- c(75000,70000,65000,60000,55000,50000,
           45000,40000,35000,30000,25000,20000,
           15000,10000,9000,8000,7000,6000,5000,
           4000,3000)

Griddata <- list()

for (i in 1:length(grids)){
  
  print(paste("working on grid ~ ",grids[i]/1000,"km resolution",sep=" "))
  
  tempgrid <- GridFilter2(study_area, resol=grids[i], Land)
  
  Griddata[[i]] <- tempgrid
  
  rm(tempgrid)
}

rm(list=setdiff(ls(), "Griddata"))

save.image("Data/QCGridsWholeCell.RData")

