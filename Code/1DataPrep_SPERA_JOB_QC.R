########################################
# Prepare data for dendrogram and obtain species groupings for benthic assemblages 
# Katie Gale (katie.gale@dfo-mpo.gc.ca) _ Created 2015
# Formatted and uploaded 2017 
#
# As published in 
# Rubidge, E., Gale, K.S.P., Curtis, J.M.R., McClelland, E., Feyrer, L., Bodtker, K., and Robb, C.
# 2016. Methodology of the Pacific Marine Ecological Classification System and its Application
# to the Northern and Southern Shelf Bioregions. DFO Can. Sci. Advis. Sec. Res. Doc.
# 2016/035. xi + 124 p.
#

## Adapted code for Biological classification analysis - Eastern Canadian groundfish - Stanley, Heaslip, Jeffery, O'Brien
#########################################

pkgs <- list("maptools", "plyr", "rgdal", "sp", "reshape", "reshape2", "tidyr", "dplyr", "raster", "rgeos", "ggplot2", "ggmap", "spatstat")
invisible(lapply(pkgs, library, character.only = T))

######
# 1. Select the datasources to be included to use in this analysis 
#In our case we will partition the survey regions. We might also look into grouping seasons
######


##
#Import data 
#This data file contains all records of all taxa for all habitats and years from DFO's groundfish 
#biological surveys in QC Region
#the first step of this file is paring it down to the dataset we want
QC <- readRDS('Data/QC_invertsAdded.rds')

#Examine distribution of entries among years and seasons
table(QC$year, QC$month) #all surveys past 2007 done b/w Jul-Sep (mainly Aug)

#Examine dataframe for duplicates
duplicates <- duplicated(QC)
table(duplicates) #68781 FALSE, 53 TRUE (duplicates)
table(QC[duplicates, 'year']) #duplicates from 2013-2017, all Bryozoa
QC <- distinct(QC) #remove duplicates
QC$ID <- paste(QC$year, QC$month, QC$day, QC$species,QC$DD_lat, QC$DD_lon, QC$spec_code, 
               QC$vessel, QC$survey,QC$nbpc,QC$stn,QC$stratum, sep = "_") #unique ID for each set/species
length(unique(QC$ID)) #63452
length(QC$ID) #68781 ## 5327 observations with repeated entries for each set

length(unique(QC$common_name)) # 249 ## equal to number of unique taxa
length(unique(QC$species)) #249

#Examine which unique ID's are repeated, but not true duplicates
Freq_UniqueID <- data.frame(table(QC$ID))
names(Freq_UniqueID) <- c("ID", "Freq")
Freq_UniqueID$ID <- as.character(Freq_UniqueID$ID)
Freq_UniqueID <- Freq_UniqueID[order(Freq_UniqueID$Freq, decreasing = T),]
Duplicates <- Freq_UniqueID[which(Freq_UniqueID$Freq > 1),]
DuplicateYears <- Duplicates %>% separate(., ID, c("Year", "Month", "Day", "Species", "latitude", "longitude",
                                                   'spec_code','vessel','survey','nbpc','station','stratum'), sep="_") %>%
  dplyr::select(.,Year)
DuplicateYears <- sort(unique(DuplicateYears$Year)) #Duplicates in every year from 2007-2017
DuplicateYears 
DuplicateTaxa <- Duplicates %>% separate(., ID, c("Year", "Month", "Day", "Species", "latitude", "longitude",
                                                   'spec_code','vessel','survey','nbpc','station','stratum'), sep="_") %>%
  dplyr::select(.,Species)
DuplicateTaxa <- sort(unique(DuplicateTaxa$Species))
DuplicateTaxa #5 taxa that are repeated, possibly because catch is split into size fractions
unique(Duplicates$Freq) #Repeat IDs are repeated anywhere from 2-9 times

#Take closer look at repeat entries

repeats <- QC[which(QC$ID %in% Duplicates$ID),1:40] #Repeats for 5 taxa due to splitting catch into size fractions/shell types

#correct redundant categories in horizontal.position and vertical.position

QC$horizontal.position <- gsub("brackish;marine", "brackish; marine", QC$horizontal.position)
QC$vertical.position <- gsub("pelagic-oceanic", "pelagic,oceanic", QC$vertical.position)

###### 
# 2. Filter dataframe to include observations 2007
#and after, groundfish and inverts (exclude pelagic fish)
#####

BottomHabitats <- c('bathydemersal', 'benthic', 'benthopelagic', 'demersal', 'sessile', 'reef-associated', 'bathypelagic')

QC <- QC %>% 
  filter(QC$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) %>% 
  filter(common_name != 'Atlantic Herring') #64102 observations

length(unique(QC$common_name)) #231 unique taxa

#Unique ID that can be used to go from long to wide format, where each row is a sample station.
QC$ID <- paste(QC$year, QC$month, QC$day, QC$DD_lat, QC$DD_lon, sep = "_") 
                
length(unique(QC$ID)) #1923 unique sets

#Max longitude values negative
QC$DD_lon <- -1*QC$DD_lon

#cast to wide format

WideQC <- dcast(QC[,c("ID","species","abundance")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(QC[!duplicated(QC$ID),],ID,DD_lat,DD_lon,year,month,day,
                     stn,stratum,vertical.position),by="ID")%>%
  dplyr::select(.,ID,DD_lon,DD_lat,year,month,day,stn,stratum,vertical.position,unique(QC$species))%>%
  data.frame()

###biotic data to point shapefile - convert to coordinate system that matches the grid
coordinates(WideQC)<-~DD_lon+DD_lat #tell R what the coordates are
proj4string(WideQC)<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system
QCdata <- spTransform(WideQC,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# writeOGR(QCdata, dsn = "Data/Shapefiles", layer = 'WideQC_surveydata', driver = 'ESRI Shapefile',overwrite_layer = T)

#load polygon for QC study area
QC_SA <- readOGR('Data/Shapefiles/QC_StudyArea.shp')
QC_Str1_5 <- readOGR('Data/Shapefiles/QC_Str1_5_agg.shp')
QC_SA <- erase(QC_SA,QC_Str1_5)

#Calculate density of survey pts and avg distance to nearest neighbour in survey area

Survey.owin <- as.owin(QC_SA) #Need to change class of survey polygon & pts to use spatstat functions
Survey.pts <- ppp(coordinates(QCdata)[,1], coordinates(QCdata)[,2], window = Survey.owin)
intensity(Survey.pts) # ~1.44 pts per 100 sq. km. 1 pt per ~69 sq. km (8.3 x 8.3 km)
NearNeigh <- nndist(Survey.pts, k = 1)
stats::quantile(NearNeigh, probs = c(0.1,0.25,0.5,0.75,0.8,0.90,0.95, na.rm =TRUE)) 
#95% of surveys with NN within 8.1 km, median of 3.4 km
mean(NearNeigh) # ~3.7 km on average to nearest neighbour
#standard tow is 0.8 nautical miles (1.48 km)
#Grid size should be no smaller than 1 X 1 km and no larger than 9 X 9 km. Choose 4 x 4 km.

#join wide format survey dataframe with grid
load("Data/QCGridsWholeCell.RData")
QCGrid <- Griddata[[20]] # 4 km resolution 
GridID<- over(QCdata,sp::geometry(QCGrid))
joinedgrid<-cbind(QCdata@data,GridID ) #Bind the data 
head(joinedgrid)

plot(QC_SA,axes=T)
points(QCdata@coords[!is.na(GridID),1],QCdata@coords[!is.na(GridID),2],col="red",pch=19,cex=0.5)
points(QCdata@coords[is.na(GridID),1],QCdata@coords[is.na(GridID),2],col="black",pch=19,cex=0.5)

#each record got assigned to a gridcell - each gridcell has an ID ("GridID"). 
#if gridID=NA, then that point was outside of our study area
populated<-(unique(joinedgrid$GridID))[complete.cases(unique(joinedgrid$GridID))] #get list of populated grid cells

QCGrid@data$id <- 1:nrow(QCGrid@data)
length(QCGrid) #8349 grid cells initiallly 
subgrid<-QCGrid[populated,] #select only those grid cells that our data populates

#View it and see it makes sense
plot(QC_SA,axes=T)
plot(subgrid,add=T)
points(QCdata@coords[!is.na(GridID),1],QCdata@coords[!is.na(GridID),2],col="red",cex=0.1)
points(QCdata@coords[is.na(GridID),1],QCdata@coords[is.na(GridID),2],col="black",pch = 16,cex=0.5)
length(populated)
# 1580 grid cells are populated by all benthic 
#megainverts + fish species (~ 18.9% of cells in regional grid. This is underestimate because study area polygon includes
#survey strata for which there are no data provided)

#Write this shapefile for later use
writeOGR(subgrid, dsn = "Data/QCGrids", layer = "QC_4km_Grid", driver = 'ESRI Shapefile', overwrite_layer = T)

goodco <- joinedgrid[!is.na(joinedgrid$GridID),] #select records that were in a grid cell 
nrow(goodco) # 1882 sets were in the QC Study Area

# write.csv(goodco, "Data/QCGriddedData.csv")

## set up data to view how many observations are in each grid cell
NumGrid <- data.frame(table(goodco$GridID))
colnames(NumGrid) <- c("Grid","Frequency")
NumGrid$Grid <- as.numeric(as.character(NumGrid$Grid))
subgrid@data$freq <- NumGrid[match(subgrid@data$id,NumGrid$Grid),"Frequency"]
QCGrid@data$freq <- NumGrid[match(QCGrid@data$id, NumGrid$Grid), "Frequency"]
QCGrid@data$freq[is.na(QCGrid@data$freq)] = 0

####Examine allocation and distribution of sets among grids

#Histogram
hist(NumGrid$Frequency, main = 'Allocation of sets among grid cells', xlab = 'Sets per grid', col = 'grey90')
abline(v = mean(NumGrid$Frequency), col = 'red', lwd = 2)
abline(v = median(NumGrid$Frequency), col = 'blue', lwd = 2, lty = 2)

mean(NumGrid$Frequency) # 1.19
median(NumGrid$Frequency) # 1
max(NumGrid$Frequency) # 6
sd(NumGrid$Frequency) # 0.48

#fortify for ggplot mapping
FortData <- fortify(subgrid, region = "id")
FortData.df <- join(FortData, subgrid@data, by = "id")
FortQC <- fortify(QC_SA)

p1 <- ggplot() + 
  geom_polygon(data=FortQC,
               aes(long, lat, group = group),fill=NA,col="black")+
  geom_polygon(data = FortData.df, 
               aes(long, lat, group = group, fill = freq)) + 
  coord_equal()+
  scale_fill_gradient(low="yellow", high="red")+
  theme_bw()+
  labs(x="UTM X (Zone 21)",y="UTM Y (Zone 21)",fill="Number of sets");p1

# ggsave("Output/SetsByGrid_QC_4km_Analysis.png",p1)

#Spatial autocorrelation of survey sets

library("spdep", lib.loc="~/R/win-library/3.5")
#nb.QCGrid <- poly2nb(QCGrid, row.names = QCGrid$id)
nb.QCGrid <- knn2nb(knearneigh(coordinates(QCGrid), k = 8,RANN = FALSE))
nb.QCGrid.wts <- nb2listw(nb.QCGrid, style = "B")
moran(QCGrid$freq, nb.QCGrid.wts, n = length(nb.QCGrid.wts$neighbours), S0 = Szero(nb.QCGrid.wts))
moran.mc(QCGrid$freq, nb.QCGrid.wts, nsim = 99) # p = 0.01, I = 0.054

#############################################
# 3. melt and cast to get site x species matrix -- We saved time and did this up front. The joined data
# up at this point is already in wide format. We do however, need to group the observations by GridID and 
#ID the cells as Presence - Absence

IDvars <-  c("ID","year","month","day","stn","stratum","vertical.position","GridID")
Species <- setdiff(names(goodco),IDvars)

## Trim out barren sites -- those sites with only one species

SiteSummary <- goodco%>%select_(.,.dots=c("GridID",Species))%>%group_by(GridID)%>%
  summarise_all(funs(sum(.)))%>%ungroup()%>%data.frame()

Indval <- SiteSummary$GridID
SpeciesSiteSummary <- SiteSummary[,Species]
SpeciesSiteSummary[SpeciesSiteSummary>0]=1
SpeciesSiteSummary[is.na(SpeciesSiteSummary)] <- 0
unique(unlist(SpeciesSiteSummary)) ## just 0 and 1 so we can now add across rows to figure which Grids are barren (1 or 0 species)
speccount <- rowSums(SpeciesSiteSummary)
range(speccount) ## 5 - 81
table(speccount) ## 0 barren sites

BarrenSites <- Indval[speccount<2]
IncludedSites <- setdiff(Indval,BarrenSites)

##############################
# drop "Rare species"
############################

## Trim out species that are only observed in fewer than 1% of stations

dropThreshold <- (ceiling(nrow(SpeciesSiteSummary)*0.01))

print(paste("Dropthreshold is",dropThreshold,"out of",nrow(SpeciesSiteSummary),"sites, or ~1%")) #16 out 1580 sites

SpeciesCountSummary <- colSums(SpeciesSiteSummary)

removedSpecies <- names(SpeciesCountSummary)[which(SpeciesCountSummary<dropThreshold)]
removedSpecies
IncludedSpecies <- setdiff(Species,removedSpecies)
length(IncludedSpecies) #180 unique species

### Trim the data ----------

Output <- SiteSummary[SiteSummary$GridID%in%IncludedSites,IncludedSpecies]
rownames(Output) <- SiteSummary[SiteSummary$GridID%in%IncludedSites,"GridID"] 
Output[Output>0] <- 1 # convert finally to precence (1) and absence (0)

##Save the data
write.csv(Output, "Data/ClusterData4km_QC.csv", row.names=T) #1580 obs. of 180 species

#Move onto cluster analysis
