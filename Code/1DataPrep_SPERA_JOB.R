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

#subset by these data types
#Load in data from assembled database
#fourregions_merged<-read.csv("Lydia_Stevens/Spreadsheets/fourregions_merged_130318.csv",stringsAsFactors = F) #792999 observations
#length(unique(fourregions_merged$common_name)) #146 unique taxa
#Maritimes2 <- fourregions_merged[fourregions_merged$region=="MARITIME",] #147367 observations
#length(unique(Maritimes$common_name)) #134 unique taxa
#rm(fourregions_merged);gc()
#save.image("Data/MaritimesData.RData")


##
#Import data 
#This data file contains all records of all taxa for all habitats and years from DFO's groundfish biological surveys in Maritime Region
#the first step of this file is paring it down to the dataset we want
load("Data/MaritimesData.RData")


#Examine dataframe for duplicates
#Maritimes$ID <- paste(Maritimes$year, Maritimes$set, Maritimes$month, Maritimes$day, Maritimes$common_name, sep = "_") #unique ID for each set/species
#length(unique(Maritimes$ID)) #182125
#length(Maritimes$ID) #182389 ## 264 Duplicate entries exist

#length(unique(Maritimes$common_name)) # 171 ## equal to number of unique species
#length(unique(Maritimes$spec_code)) #167
#length(unique(Maritimes$spec)) #196
#Fewer common_names/species than original species code. Perhaps some duplicates exist because unique taxa were 
#incorrectly reclassified

#Examine which unique ID's are duplicated
#Freq_UniqueID <- data.frame(table(Maritimes$ID))
#names(Freq_UniqueID) <- c("ID", "Freq")
#Freq_UniqueID$ID <- as.character(Freq_UniqueID$ID)
#Freq_UniqueID <- Freq_UniqueID[order(Freq_UniqueID$Freq, decreasing = T),]
#Duplicates <- Freq_UniqueID[which(Freq_UniqueID$Freq > 1),]
#DuplicateNames <- Duplicates %>% separate(., ID, c("Year", "Set", "Month", "Day", "common_name"), sep="_") %>%
#  dplyr::select(.,common_name)
#DuplicateNames <- sort(unique(DuplicateNames$common_name)) #26 Duplicated species


#Correct species misclassifications

Maritimes$species[Maritimes$spec == "340"] <- "Aspidophoroides monopterygius"
Maritimes$common_name[Maritimes$spec == "340"] <- "Alligatorfish"
Maritimes$species[Maritimes$spec == "2527"] <- "Hyas araneus"
Maritimes$common_name[Maritimes$spec == "2527"] <- "Toad crab"
Maritimes$species[Maritimes$spec == "601"] <- "Simenchelys parasitica"
Maritimes$common_name[Maritimes$spec == "601"] <- "Snubnose eel"
Maritimes$species[Maritimes$spec == "4212"] <- "Buccinum scalariforme"
Maritimes$common_name[Maritimes$spec == "4212"] <- "Silky buccinum"
Maritimes$species[Maritimes$spec == "503"] <- "Liparis atlanticus"
Maritimes$common_name[Maritimes$spec == "503"] <- "Atlantic snailfish"

#correct redundant categories in horizontal.position and vertical.position

Maritimes$horizontal.position <- gsub("brackish;marine", "brackish; marine", Maritimes$horizontal.position)
Maritimes$vertical.position <- gsub("pelagic-oceanic", "pelagic,oceanic", Maritimes$vertical.position)

###### 
# 2. Filter dataframe to include observations 2007 and after, groundfish and inverts
#####

BottomHabitats <- c('bathydemersal', 'benthic', 'benthopelagic', 'demersal', 'sessile', 'reef-associated', 'bathypelagic')

MaritimesSub <- Maritimes %>% 
  filter(Maritimes$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) #51309 observations

length(unique(MaritimesSub$common_name)) #160 unique taxa

#Unique ID that can be used to go from long to wide format, where each row is a sample sation.
MaritimesSub$ID <- paste(MaritimesSub$year,MaritimesSub$month,MaritimesSub$day, MaritimesSub$set, sep="_")
length(unique(MaritimesSub$ID)) #3412 unique sets

#cast to wide format
WideMaritimes_number <- dcast(MaritimesSub[,c("ID","species","totno")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(MaritimesSub[!duplicated(MaritimesSub$ID),],ID,latitude,longitude,year,month,day,
                     set,strat,season,vertical.position),by="ID")%>%
  dplyr::select(.,ID,longitude,latitude,year,month,day,set,strat,season,vertical.position,unique(MaritimesSub$species))%>%
  data.frame()

WideMaritimes_wgt <- dcast(MaritimesSub[,c("ID","species","totwgt")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(MaritimesSub[!duplicated(MaritimesSub$ID),],ID,latitude,longitude,year,month,day,
                     set,strat,season,vertical.position),by="ID")%>%
  dplyr::select(.,ID,longitude,latitude,year,month,day,set,strat,season,vertical.position,unique(MaritimesSub$species))%>%
  data.frame()


###biotic data to point shapefile - convert to coordinate system that matches the grid
coordinates(WideMaritimes_wgt)<-~longitude+latitude #tell R what the coordates are
proj4string(WideMaritimes_wgt)<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system
Mardata <- spTransform(WideMaritimes_wgt,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
#writeSpatialShape(Mardata, "Data/WideMaritimes_wgt_ptshp.shp")

#Maritimes region planning area
MaritimeRegion <- readOGR("Data/MaritimesPlanningRegion/MaritimesPlanningArea.shp")
MaritimeProj <- spTransform(MaritimeRegion,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#Calculate density of survey pts and avg distance to nearest neighbour in survey area

RVsurvey <- readOGR("Data/MaritimesPlanningRegion/MaritimesRegionStrataBoundaries.shp")
RVsurvey <- spTransform(RVsurvey, CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
SurveyBoundary <- aggregate(intersect(MaritimeProj, RVsurvey))
Survey.owin <- as.owin(SurveyBoundary) #Need to change class of survey polygon & pts to use spatstat functions
Survey.pts <- ppp(coordinates(Mardata)[,1], coordinates(Mardata)[,2], window = Survey.owin)
intensity(Survey.pts) # ~1.4 pts per 100 sq. km
NearNeigh <- nndist(Survey.pts, k = 1)
stats::quantile(NearNeigh, probs = c(0.5,0.75,0.8,0.90,0.95, na.rm =TRUE)) #95% of surveys with NN within 8.5 km
mean(NearNeigh) # ~3.8 km on average to nearest neighbour

#join wide format survey dataframe with grid
load("Data/MaritimeGridsWholeCell.RData")
MaritimeGrid <- Griddata[[20]] # 4 km resolution 
GridID<- over(Mardata,sp::geometry(MaritimeGrid))
joinedgrid<-cbind(Mardata@data,GridID ) #Bind the data 
head(joinedgrid)

plot(MaritimeGrid,axes=T)
points(Mardata@coords[!is.na(GridID),1],Mardata@coords[!is.na(GridID),2],col="red",pch=19,cex=0.5)
points(Mardata@coords[is.na(GridID),1],Mardata@coords[is.na(GridID),2],col="black",pch=19,cex=0.5)

#each record got assigned to a gridcell - each gridcell has an ID ("GridID"). 
#if gridID=NA, then that point was outside of our study area
populated<-(unique(joinedgrid$GridID))[complete.cases(unique(joinedgrid$GridID))] #get list of populated grid cells

MaritimeGrid@data$id <- 1:nrow(MaritimeGrid@data)
length(MaritimeGrid) #14780 grid cells initiallly 
subgrid<-MaritimeGrid[populated,] #select only those grid cells that our data populates

#View it and see it makes sense
plot(MaritimeProj,axes=T)
plot(subgrid,add=T)
points(Mardata@coords[!is.na(GridID),1],Mardata@coords[!is.na(GridID),2],col="red",pch=19,cex=0.5)
points(Mardata@coords[is.na(GridID),1],Mardata@coords[is.na(GridID),2],col="black",pch=19,cex=0.5)
length(populated) #   2529 grid cells are populated by all benthic megainverts + fish species (~ 60% of cells in regional grid)

#Write this shapefile for later use
#writeSpatialShape(subgrid,"Data/MaritimeGrids/Maritimes_4km_Grid.shp")

goodco <- joinedgrid[!is.na(joinedgrid$GridID),] #select records that were in a grid cell 
nrow(goodco) # 3018 sets were in the Maritime Study Area

#write.csv(goodco, "Data/MaritimesGriddedData.csv")

## set up data to view how many observations are in each grid cell
NumGrid <- data.frame(table(goodco$GridID))
colnames(NumGrid) <- c("Grid","Frequency")
NumGrid$Grid <- as.numeric(as.character(NumGrid$Grid))
subgrid@data$freq <- NumGrid[match(subgrid@data$id,NumGrid$Grid),"Frequency"]
MaritimeGrid@data$freq <- NumGrid[match(MaritimeGrid@data$id, NumGrid$Grid), "Frequency"]
MaritimeGrid@data$freq[is.na(MaritimeGrid@data$freq)] = 0

####Examine allocation and distribution of sets among grids

#Histogram
hist(NumGrid$Frequency, main = 'Allocation of sets among grid cells', xlab = 'Sets per grid', col = 'grey90')
abline(v = mean(NumGrid$Frequency), col = 'red', lwd = 2)
abline(v = median(NumGrid$Frequency), col = 'blue', lwd = 2, lty = 2)

mean(NumGrid$Frequency) # 1.19
median(NumGrid$Frequency) # 1
max(NumGrid$Frequency) # 5
sd(NumGrid$Frequency) # 0.47

#fortify for ggplot mapping
FortData <- fortify(subgrid, region = "id")
FortData.df <- join(FortData, subgrid@data, by = "id")
FortMaritimes <- fortify(MaritimeProj)

p1 <- ggplot() + 
  geom_polygon(data=FortMaritimes,
               aes(long, lat, group = group),fill=NA,col="black")+
  geom_polygon(data = FortData.df, 
               aes(long, lat, group = group, fill = freq)) + 
  coord_equal()+
  scale_fill_gradient(low="yellow", high="red")+
  theme_bw()+
  labs(x="UTM X (Zone 20)",y="UTM Y (Zone 20)",fill="Number of sets");p1

#ggsave("Output/SetsByGrid_Maritimes_20km_Analysis.png",p1)

#Spatial autocorrelation of survey sets

library("spdep", lib.loc="~/R/win-library/3.5")
#nb.MarGrid <- poly2nb(MaritimeGrid, row.names = MaritimeGrid$id)
nb.MarGrid <- knn2nb(knearneigh(coordinates(MaritimeGrid), k = 8,RANN = FALSE))
nb.MarGrid.wts <- nb2listw(nb.MarGrid, style = "B")
moran(MaritimeGrid$freq, nb.MarGrid.wts, n = length(nb.MarGrid.wts$neighbours), S0 = Szero(nb.MarGrid.wts))
moran.mc(MaritimeGrid$freq, nb.MarGrid.wts, nsim = 99)

#############################################
# 3. melt and cast to get site x species matrix -- We saved time and did this up front. The joined data
# up at this point is already in wide format. We do however, need to group the observations by GridID and 
#ID the cells as Presence - Absence

IDvars <-  c("ID","year","month","day","set","strat","season","vertical.position","GridID")
Species <- setdiff(names(goodco),IDvars)

## Trim out barren sites -- those sites with only one species

SiteSummary <- goodco%>%select_(.,.dots=c("GridID",Species))%>%group_by(GridID)%>%
  summarise_all(funs(sum(.)))%>%ungroup()%>%data.frame()

Indval <- SiteSummary$GridID
SpeciesSiteSummary <- SiteSummary[,Species]
SpeciesSiteSummary[SpeciesSiteSummary>0]=1
unique(unlist(SpeciesSiteSummary)) ## just 0 and 1 so we can now add across rows to figure which Grids are barren (1 or 0 species)
speccount <- rowSums(SpeciesSiteSummary)
range(speccount) ## 1 - 39
table(speccount) ## 2 barren sites

BarrenSites <- Indval[speccount<2]
IncludedSites <- setdiff(Indval,BarrenSites)

##############################
# drop "Rare species"
############################

## Trim out species that are only observed in fewer than 1% of stations

dropThreshold <- (ceiling(nrow(SpeciesSiteSummary)*0.01))

print(paste("Dropthreshold is",dropThreshold,"out of",nrow(SpeciesSiteSummary),"sites, or ~1%"))

SpeciesCountSummary <- colSums(SpeciesSiteSummary)

removedSpecies <- names(SpeciesCountSummary)[which(SpeciesCountSummary<dropThreshold)]
removedSpecies
IncludedSpecies <- setdiff(Species,removedSpecies)
length(IncludedSpecies) #114 unique species

### Trim the data ----------

Output <- SiteSummary[SiteSummary$GridID%in%IncludedSites,IncludedSpecies]
rownames(Output) <- SiteSummary[SiteSummary$GridID%in%IncludedSites,"GridID"] 
Output[Output>0] <- 1 # convert finally to precence (1) and absence (0)
Output <- dplyr::select(Output, -Clupea.harengus) #113 unique species
##Save the data
#write.csv(Output, "Data/ClusterData4km.csv", row.names=T)

#Move onto cluster analysis
