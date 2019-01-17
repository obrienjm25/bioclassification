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
#fourregions_merged<-read.csv("Data/fourregions_merged_130318.csv",stringsAsFactors = F) #962256 observations
#length(unique(fourregions_merged$common_name)) #183 unique taxa
#Nfld <- fourregions_merged[fourregions_merged$region=="NEWFOUNDLAND",] #191784 observations
#length(unique(Nfld$common_name)) #92 unique taxa
#rm(fourregions_merged);gc()
#save.image("Data/NfldData.RData")


##
#Import data 
#This data file contains all records of all taxa for all habitats and years from DFO's groundfish 
#biological surveys in Newfoundlan Region
#the first step of this file is paring it down to the dataset we want
load("Data/NfldData.RData")

#Examine distribution of entries among years and seasons
table(Nfld$year, Nfld$season)

#Examine dataframe for duplicates
Nfld$ID <- paste(Nfld$year, Nfld$set, Nfld$month, Nfld$day, Nfld$common_name, Nfld$longitude, Nfld$latitude, sep = "_") #unique ID for each set/species
length(unique(Nfld$ID)) #191748
length(Nfld$ID) #191784 ## 36 Duplicate entries exist

length(unique(Nfld$common_name)) # 92 ## equal to number of unique species
length(unique(Nfld$spec_code)) #91
length(unique(Nfld$spec)) #2
#unique spec_code and unique taxa don't match up. Perhaps incorrect reclassification?

#Examine which unique ID's are duplicated
Freq_UniqueID <- data.frame(table(Nfld$ID))
names(Freq_UniqueID) <- c("ID", "Freq")
Freq_UniqueID$ID <- as.character(Freq_UniqueID$ID)
Freq_UniqueID <- Freq_UniqueID[order(Freq_UniqueID$Freq, decreasing = T),]
Duplicates <- Freq_UniqueID[which(Freq_UniqueID$Freq > 1),]
DuplicateNames <- Duplicates %>% separate(., ID, c("Year", "Set", "Month", "Day", "common_name", "long", "lat"), sep="_") %>%
  dplyr::select(.,common_name)
DuplicateNames <- sort(unique(DuplicateNames$common_name)) #21 Duplicated species
DuplicateNames
#duplicates appear to be due to lumping redfish into one category or repeating invalid sets
#Duplicate entries will collapse (sum of weights) when the dataframe is reshaped

#correct redundant categories in horizontal.position and vertical.position

Nfld$horizontal.position <- gsub("brackish;marine", "brackish; marine", Nfld$horizontal.position)
Nfld$vertical.position <- gsub("pelagic-oceanic", "pelagic,oceanic", Nfld$vertical.position)

###### 
# 2. Filter dataframe to include observations 2007 (maybe 2003 since we only have data to 2013) 
#and after, groundfish and inverts (exclude pelagic fish)
#separate fall and spring survey sets and focus on fall for now
#####

BottomHabitats <- c('bathydemersal', 'benthic', 'benthopelagic', 'demersal', 'sessile', 'reef-associated', 'bathypelagic')

NfldFall <- Nfld %>% 
  filter(Nfld$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) %>%
  filter(season == 'Fall') #39653 observations

NfldSpring <- Nfld %>% 
  filter(Nfld$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) %>%
  filter(season == 'Spring') #21820

length(unique(NfldFall$common_name)) #85 unique taxa
length(unique(NfldSpring$common_name)) #82 unique taxa

#Unique ID that can be used to go from long to wide format, where each row is a sample sation.
NfldFall$ID <- paste(NfldFall$year,NfldFall$month,NfldFall$day, NfldFall$longitude, NfldFall$latitude, sep="_")
length(unique(NfldFall$ID)) #3407 unique sets

#cast to wide format

WideNfldFall_wgt <- dcast(NfldFall[,c("ID","species","totwgt")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(NfldFall[!duplicated(NfldFall$ID),],ID,latitude,longitude,year,month,day,
                     set,strat,season,vertical.position),by="ID")%>%
  dplyr::select(.,ID,longitude,latitude,year,month,day,set,strat,season,vertical.position,unique(NfldFall$species))%>%
  data.frame()


###biotic data to point shapefile - convert to coordinate system that matches the grid
coordinates(WideNfldFall_wgt)<-~longitude+latitude #tell R what the coordates are
proj4string(WideNfldFall_wgt)<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system
Nflddata <- spTransform(WideNfldFall_wgt,CRS("+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
#writeOGR(Nflddata, dsn = "Data/Shapefiles", layer = 'WideNL_surveydata', driver = 'ESRI Shapefile')

#Nfld survey strata boundaries
#RVstrataNL <- readOGR('Data/Shapefiles/NF_SamplingStrata_20140514.shp')
#RVstrataNL <- spTransform(RVstrataNL, CRS('+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
#NAFO divisions
#NAFO <- readOGR("Data/Shapefiles/NAFO_Divisions.shp")
#proj4string(NAFO) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
#NAFO <- spTransform(NAFO, CRS('+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))

#Examine which NAFO divisions contain survey data
#plot(Nflddata, cex = 0.1, col = 'red')
#plot(NAFO, add = T) #Little to no sampling in divisions 2G, 3M, 3Pn, 3Ps
#rm_list <- factor(c('2G','3M','3Pn','3Ps','0B','4Vn','4Vs','4R'))
#NAFO_rm <- NAFO[NAFO@data$ZONE %in% rm_list,]
#NL_SA <- aggregate(RVstrataNL) # aggegrate strata polygons to form boundaries of study area
#NL_SA <- erase(NL_SA,NAFO_rm)#mask study area polygon with 2G, 3M, 3Pn, 3Ps NAFO divisions
#writeOGR(as(NL_SA, 'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'NL_RVsurveyAgg', driver = 'ESRI Shapefile')

#load polygon for Nfld study area
NL_SA <- readOGR('Data/Shapefiles/NL_RVsurveyAgg.shp')

#Calculate density of survey pts and avg distance to nearest neighbour in survey area

Survey.owin <- as.owin(NL_SA) #Need to change class of survey polygon & pts to use spatstat functions
Survey.pts <- ppp(coordinates(Nflddata)[,1], coordinates(Nflddata)[,2], window = Survey.owin)
intensity(Survey.pts) # ~1 pts per 167 sq. km
NearNeigh <- nndist(Survey.pts, k = 1)
stats::quantile(NearNeigh, probs = c(0.1,0.25,0.5,0.75,0.8,0.90,0.95, na.rm =TRUE)) 
#95% of surveys with NN within 13 km, median of 5.6 km
mean(NearNeigh) # ~6.1 km on average to nearest neighbour
#standard tow is 0.8 nautical miles (1.48 km)
#Grid size should be no smaller than 2 X 2 km and no larger than 13 X 13 km

#join wide format survey dataframe with grid
load("Data/NLGridsWholeCell.RData")
NLGrid <- Griddata[[20]] # 4 km resolution 
GridID<- over(Nflddata,sp::geometry(NLGrid))
joinedgrid<-cbind(Nflddata@data,GridID ) #Bind the data 
head(joinedgrid)

plot(NL_SA,axes=T)
points(Nflddata@coords[!is.na(GridID),1],Nflddata@coords[!is.na(GridID),2],col="red",pch=19,cex=0.5)
points(Nflddata@coords[is.na(GridID),1],Nflddata@coords[is.na(GridID),2],col="black",pch=19,cex=0.5)

#each record got assigned to a gridcell - each gridcell has an ID ("GridID"). 
#if gridID=NA, then that point was outside of our study area
populated<-(unique(joinedgrid$GridID))[complete.cases(unique(joinedgrid$GridID))] #get list of populated grid cells

NLGrid@data$id <- 1:nrow(NLGrid@data)
length(NLGrid) #35390 grid cells initiallly 
subgrid<-NLGrid[populated,] #select only those grid cells that our data populates

#View it and see it makes sense
plot(NL_SA,axes=T)
plot(subgrid,add=T)
points(Nflddata@coords[!is.na(GridID),1],Nflddata@coords[!is.na(GridID),2],col="red",cex=0.1)
points(Nflddata@coords[is.na(GridID),1],Nflddata@coords[is.na(GridID),2],col="black",pch = 16,cex=0.5)
length(populated) #   3135 grid cells are populated by all benthic megainverts + fish species (~ 9% of cells in regional grid)

#Write this shapefile for later use
#writeOGR(subgrid, dsn = "Data/NfldGrids", layer = "NL_4km_Grid", driver = 'ESRI Shapefile')

goodco <- joinedgrid[!is.na(joinedgrid$GridID),] #select records that were in a grid cell 
nrow(goodco) # 3390 sets were in the Nfld Study Area

#write.csv(goodco, "Data/NfldGriddedData.csv")

## set up data to view how many observations are in each grid cell
NumGrid <- data.frame(table(goodco$GridID))
colnames(NumGrid) <- c("Grid","Frequency")
NumGrid$Grid <- as.numeric(as.character(NumGrid$Grid))
subgrid@data$freq <- NumGrid[match(subgrid@data$id,NumGrid$Grid),"Frequency"]
NLGrid@data$freq <- NumGrid[match(NLGrid@data$id, NumGrid$Grid), "Frequency"]
NLGrid@data$freq[is.na(NLGrid@data$freq)] = 0

####Examine allocation and distribution of sets among grids

#Histogram
hist(NumGrid$Frequency, main = 'Allocation of sets among grid cells', xlab = 'Sets per grid', col = 'grey90')
abline(v = mean(NumGrid$Frequency), col = 'red', lwd = 2)
abline(v = median(NumGrid$Frequency), col = 'blue', lwd = 2, lty = 2)

mean(NumGrid$Frequency) # 1.08
median(NumGrid$Frequency) # 1
max(NumGrid$Frequency) # 4
sd(NumGrid$Frequency) # 0.29

#fortify for ggplot mapping
FortData <- fortify(subgrid, region = "id")
FortData.df <- join(FortData, subgrid@data, by = "id")
FortNfld <- fortify(NL_SA)

p1 <- ggplot() + 
  geom_polygon(data=FortNfld,
               aes(long, lat, group = group),fill=NA,col="black")+
  geom_polygon(data = FortData.df, 
               aes(long, lat, group = group, fill = freq)) + 
  coord_equal()+
  scale_fill_gradient(low="yellow", high="red")+
  theme_bw()+
  labs(x="UTM X (Zone 21)",y="UTM Y (Zone 21)",fill="Number of sets");p1

#ggsave("Output/SetsByGrid_Nfld_4km_Analysis.png",p1)

#Spatial autocorrelation of survey sets

library("spdep", lib.loc="~/R/win-library/3.5")
#nb.NLGrid <- poly2nb(MaritimeGrid, row.names = MaritimeGrid$id)
nb.NLGrid <- knn2nb(knearneigh(coordinates(NLGrid), k = 8,RANN = FALSE))
nb.NLGrid.wts <- nb2listw(nb.NLGrid, style = "B")
moran(NLGrid$freq, nb.NLGrid.wts, n = length(nb.NLGrid.wts$neighbours), S0 = Szero(nb.NLGrid.wts))
moran.mc(NLGrid$freq, nb.NLGrid.wts, nsim = 99)

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
range(speccount) ## 2 - 39
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
removedSpecies #Is is weird that pollock and witch flounder among these>
IncludedSpecies <- setdiff(Species,removedSpecies)
length(IncludedSpecies) #69 unique species

### Trim the data ----------

Output <- SiteSummary[SiteSummary$GridID%in%IncludedSites,IncludedSpecies]
rownames(Output) <- SiteSummary[SiteSummary$GridID%in%IncludedSites,"GridID"] 
Output[Output>0] <- 1 # convert finally to precence (1) and absence (0)

##Save the data
write.csv(Output, "Data/ClusterData4km_Nfld.csv", row.names=T)

#Move onto cluster analysis
