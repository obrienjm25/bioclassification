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
#biological surveys in Gulf Region from 2007-2017
#the first step of this file is paring it down to the dataset we want
gulf <- readRDS("Data/Gulf_invertsAdded.rds")

#Examine distribution of entries among years and seasons
table(gulf$year, gulf$month) #all surveys done August to October (mostly September)

#Examine dataframe for duplicates
duplicates <- duplicated(gulf)
table(duplicates) #45348 FALSE, 0 TRUE
gulf$ID <- paste(gulf$year, gulf$month, gulf$day, gulf$species,
                     gulf$set.number, gulf$stratum, sep = "_") #unique ID for each set/species
length(unique(gulf$ID)) #44473
length(gulf$ID) #45348 ## 875 Duplicate entries exist

length(unique(gulf$common_name)) # 187 ## equal to number of unique species
length(unique(gulf$species)) #187
length(unique(gulf$spec_code)) #302
#unique spec_code and unique taxa don't match up. Duplicates likely due to aggregation of original taxa list

#Examine which unique ID's are duplicated, which years, and which taxa
Freq_UniqueID <- data.frame(table(gulf$ID))
names(Freq_UniqueID) <- c("ID", "Freq")
Freq_UniqueID$ID <- as.character(Freq_UniqueID$ID)
Freq_UniqueID <- Freq_UniqueID[order(Freq_UniqueID$Freq, decreasing = T),]
Duplicates <- Freq_UniqueID[which(Freq_UniqueID$Freq > 1),]
DuplicateYears <- Duplicates %>% separate(., ID, c("Year", "Month", "Day", "Species", "Set", "Strat"), sep="_") %>%
  dplyr::select(.,Year)
DuplicateYears <- sort(unique(DuplicateYears$Year)) #Duplicates in all years (not systematic)
DuplicateSpecies <- Duplicates %>% separate(., ID, c("Year", "Month", "Day", "Species", "Set", "Strat"), sep="_") %>%
 dplyr::select(.,Species)
unique(DuplicateSpecies$Species) #all taxa that are regroupings of multiple taxa in original taxa list

#duplicates will disappear when the dataframe is re-organized to wide format

#correct redundant categories in horizontal.position and vertical.position

gulf$horizontal.position <- gsub("brackish;marine", "brackish; marine", gulf$horizontal.position)
gulf$vertical.position <- gsub("pelagic-oceanic", "pelagic,oceanic", gulf$vertical.position)

###### 
# 2. Filter dataframe to include observations 2007
#and after, groundfish and inverts (exclude pelagic fish)
#####

BottomHabitats <- c('bathydemersal', 'benthic', 'benthopelagic', 'demersal', 'sessile', 'reef-associated', 'bathypelagic')

gulf <- gulf %>% 
  filter(gulf$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) #42756 observations

length(unique(gulf$common_name)) #175 unique taxa

#Unique ID that can be used to go from long to wide format, where each row is a sample sation.
gulf$ID <- paste(gulf$year,gulf$month,gulf$day, gulf$set.number, gulf$stratum, sep="_")
length(unique(gulf$ID)) #1712 unique sets

#cast to wide format

Widegulf_wgt <- dcast(gulf[,c("ID","species","biomass")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(gulf[!duplicated(gulf$ID),],ID,latitude,longitude,year,month,day,
                     set.number,stratum,vertical.position),by="ID")%>%
  dplyr::select(.,ID,longitude,latitude,year,month,day,set.number,stratum,vertical.position,unique(gulf$species))%>%
  data.frame()

###biotic data to point shapefile - convert to coordinate system that matches the grid
coordinates(Widegulf_wgt)<-~longitude+latitude #tell R what the coordates are
proj4string(Widegulf_wgt)<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system
gulfdata <- spTransform(Widegulf_wgt,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
writeOGR(gulfdata, dsn = "Data/Shapefiles", layer = 'WideGulf_surveydata', 
         driver = 'ESRI Shapefile', overwrite_layer = T)

#load polygon for Gulf study area
Gulf_SA <- readOGR('Data/Shapefiles/Gulf_RVsurveyAgg.shp')

#Calculate density of survey pts and avg distance to nearest neighbour in survey area

Survey.owin <- as.owin(Gulf_SA) #Need to change class of survey polygon & pts to use spatstat functions
Survey.pts <- ppp(coordinates(gulfdata)[,1], coordinates(gulfdata)[,2], window = Survey.owin)
intensity(Survey.pts) # ~2.3 pts per 100 sq. km. 1 pt per ~43 sq. km (6.5 x 6.5 km)
NearNeigh <- nndist(Survey.pts, k = 1)
stats::quantile(NearNeigh, probs = c(0.1,0.25,0.5,0.75,0.8,0.90,0.95, na.rm =TRUE)) 
#95% of surveys with NN within 6.5 km, median of 2.7 km
mean(NearNeigh) # ~3 km on average to nearest neighbour
#standard tow is 1.75 nautical miles (3.24 km)
#Grid size should be no smaller than 3 X 3 km and no larger than 7 X 7 km

#join wide format survey dataframe with grid
load("Data/GulfGridsWholeCell.RData")
GulfGrid <- Griddata[[20]] # 4 km resolution 
GridID<- over(gulfdata,sp::geometry(GulfGrid))
joinedgrid<-cbind(gulfdata@data,GridID ) #Bind the data 
head(joinedgrid)

plot(Gulf_SA,axes=T)
points(gulfdata@coords[!is.na(GridID),1],gulfdata@coords[!is.na(GridID),2],col="red",pch=19,cex=0.5)
points(gulfdata@coords[is.na(GridID),1],gulfdata@coords[is.na(GridID),2],col="black",pch=19,cex=0.5)

#each record got assigned to a gridcell - each gridcell has an ID ("GridID"). 
#if gridID=NA, then that point was outside of our study area
populated<-(unique(joinedgrid$GridID))[complete.cases(unique(joinedgrid$GridID))] #get list of populated grid cells

GulfGrid@data$id <- 1:nrow(GulfGrid@data)
length(GulfGrid) #4554 grid cells initiallly 
subgrid<-GulfGrid[populated,] #select only those grid cells that our data populates

#View it and see it makes sense
plot(Gulf_SA,axes=T)
plot(subgrid,add=T)
points(gulfdata@coords[!is.na(GridID),1],gulfdata@coords[!is.na(GridID),2],col="red",cex=0.1)
points(gulfdata@coords[is.na(GridID),1],gulfdata@coords[is.na(GridID),2],col="black",pch = 16,cex=0.5)
length(populated) #   1326 grid cells are populated by all benthic megainverts + fish species (~ 29.1% of cells in regional grid)

#Write this shapefile for later use
writeOGR(subgrid, dsn = "Data/GulfGrids", layer = "Gulf_4km_Grid", driver = 'ESRI Shapefile', overwrite_layer = T)

goodco <- joinedgrid[!is.na(joinedgrid$GridID),] #select records that were in a grid cell 
nrow(goodco) # 1696 sets were in the Gulf Study Area

write.csv(goodco, "Data/GulfGriddedData.csv")

## set up data to view how many observations are in each grid cell
NumGrid <- data.frame(table(goodco$GridID))
colnames(NumGrid) <- c("Grid","Frequency")
NumGrid$Grid <- as.numeric(as.character(NumGrid$Grid))
subgrid@data$freq <- NumGrid[match(subgrid@data$id,NumGrid$Grid),"Frequency"]
GulfGrid@data$freq <- NumGrid[match(GulfGrid@data$id, NumGrid$Grid), "Frequency"]
GulfGrid@data$freq[is.na(GulfGrid@data$freq)] = 0

####Examine allocation and distribution of sets among grids

#Histogram
hist(NumGrid$Frequency, main = 'Allocation of sets among grid cells', xlab = 'Sets per grid', col = 'grey90')
abline(v = mean(NumGrid$Frequency), col = 'red', lwd = 2)
abline(v = median(NumGrid$Frequency), col = 'blue', lwd = 2, lty = 2)

mean(NumGrid$Frequency) # 1.28
median(NumGrid$Frequency) # 1
max(NumGrid$Frequency) # 5
sd(NumGrid$Frequency) # 0.59

#fortify for ggplot mapping
FortData <- fortify(subgrid, region = "id")
FortData.df <- join(FortData, subgrid@data, by = "id")
FortGulf <- fortify(Gulf_SA)

p1 <- ggplot() + 
  geom_polygon(data=FortGulf,
               aes(long, lat, group = group),fill=NA,col="black")+
  geom_polygon(data = FortData.df, 
               aes(long, lat, group = group, fill = freq)) + 
  coord_equal()+
  scale_fill_gradient(low="yellow", high="red")+
  theme_bw()+
  labs(x="UTM X (Zone 21)",y="UTM Y (Zone 21)",fill="Number of sets");p1

ggsave("Output/SetsByGrid_Gulf_4km_Analysis.png",p1)

#Spatial autocorrelation of survey sets

library("spdep", lib.loc="~/R/win-library/3.5")
#nb.GulfGrid <- poly2nb(GulfGrid, row.names = GulfGrid$id)
nb.GulfGrid <- knn2nb(knearneigh(coordinates(GulfGrid), k = 8,RANN = FALSE))
nb.GulfGrid.wts <- nb2listw(nb.GulfGrid, style = "B")
moran(GulfGrid$freq, nb.GulfGrid.wts, n = length(nb.GulfGrid.wts$neighbours), S0 = Szero(nb.GulfGrid.wts))
moran.mc(GulfGrid$freq, nb.GulfGrid.wts, nsim = 99) # p = 0.01, I = 0.03

#############################################
# 3. melt and cast to get site x species matrix -- We saved time and did this up front. The joined data
# up at this point is already in wide format. We do however, need to group the observations by GridID and 
#ID the cells as Presence - Absence

IDvars <-  c("ID","year","month","day","set.number","stratum","vertical.position","GridID")
Species <- setdiff(names(goodco),IDvars)

## Trim out barren sites -- those sites with only one species

SiteSummary <- goodco%>%select_(.,.dots=c("GridID",Species))%>%group_by(GridID)%>%
  summarise_all(funs(sum(.)))%>%ungroup()%>%data.frame()

Indval <- SiteSummary$GridID
SpeciesSiteSummary <- SiteSummary[,Species]
SpeciesSiteSummary[SpeciesSiteSummary>0]=1
unique(unlist(SpeciesSiteSummary)) ## just 0 and 1 so we can now add across rows to figure which Grids are barren (1 or 0 species)
speccount <- rowSums(SpeciesSiteSummary)
range(speccount) ## 7 - 62
table(speccount) ## 0 barren sites

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
length(IncludedSpecies) #135 unique taxa

### Trim the data ----------

Output <- SiteSummary[SiteSummary$GridID%in%IncludedSites,IncludedSpecies] %>% 
  dplyr::select(-Clupea.harengus) # remove herring (134 unique taxa)
rownames(Output) <- SiteSummary[SiteSummary$GridID%in%IncludedSites,"GridID"] 
Output[Output>0] <- 1 # convert finally to precence (1) and absence (0)

##Save the data
write.csv(Output, "Data/ClusterData4km_Gulf.csv", row.names=T)

#Move onto cluster analysis
