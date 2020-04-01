####Climate Forecasting & Hindcasting####

#Load required packages
pkgs <- list("maptools", "plyr", "rgdal", "sp", "reshape", "reshape2", "tidyr", "dplyr", "raster", "rgeos","purrr","data.table")
invisible(lapply(pkgs, library, character.only = T))

####Data Prep####

#load Maritimes Data

load('Data/MaritimesData.RData')

#Filter dataframe to include observations groundfish and inverts

BottomHabitats <- c('bathydemersal', 'benthic', 'benthopelagic', 'demersal', 'sessile', 'reef-associated', 'bathypelagic')

#split data into 2007-2011 & >2012

MaritimesPre2012 <- Maritimes %>% 
  filter(Maritimes$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007 & year < 2012) #25812 observations

MaritimesPost2012 <- Maritimes %>% 
  filter(Maritimes$vertical.position %in% BottomHabitats) %>%
  filter(year >= 2012) #25497 observations

#make list object with pre- and post-2012 dataframes
MaritimesDFs <- list(MaritimesPre2012, MaritimesPost2012)   

map(MaritimesDFs, ~ length(unique(.$common_name))) #156 pre 2012 and 146 unique taxa post 2012

#Unique ID that can be used to go from long to wide format (rows are sample stations)
MaritimesDFs <- MaritimesDFs %>% 
  map(~ data.frame(.,ID = paste(.$year,.$month,.$day, .$set, sep="_")))

map(MaritimesDFs, ~ length(unique(.$ID))) #1734 & 1678 unique sets respectively

#Convert each to Spatial Points Dataframe

WideMaritimes_wgt <- MaritimesDFs %>% 
  map(~ inner_join(dplyr::select(.[!duplicated(.$ID),],ID,longitude,latitude,latitude,year,month,day, set,strat,season,vertical.position),
                   dcast(.[,c("ID","species","totwgt")],ID~species, fun.aggregate = sum), by = "ID")) %>% #cast to wide format
  map(.,data.frame) %>% #convert to dataframe
  map(~ `coordinates<-`(.,c('longitude','latitude'))) %>% #convert to SPDF
  map(~ `proj4string<-`(.,CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))) %>% #set crs
  map(~ spTransform(.,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))) #transform to cartesian projection

#join wide format survey dataframe with grid
load("Data/MaritimeGridsWholeCell.RData") #load Maritimes grid data
MaritimeGrid <- Griddata[[20]] # 4 km resolution grid

GridID <- map(WideMaritimes_wgt,~ over(.,sp::geometry(MaritimeGrid))) #get grid ID that 
joinedgrid <- map2(WideMaritimes_wgt, GridID, ~ cbind(.x@data, .y)) %>%  #Bind the data
  map(., ~ dplyr::rename(., GridID = .y)) 

#each record got assigned to a gridcell - each gridcell has an ID ("GridID"). 
#if gridID=NA, then that point was outside of our study area
#remove unpopulated grid cells and records outside study area from Spatial DF and DF respectively

populated <- map(joinedgrid, ~ (unique(.$GridID))[complete.cases(unique(.$GridID))]) #get list of populated grid cells

MaritimeGrid@data$id <- 1:nrow(MaritimeGrid@data)
length(MaritimeGrid) #14780 grid cells initiallly 
subgrid <- populated %>% map(~ MaritimeGrid[.,]) #select only those grid cells that our data populates
map(subgrid, nrow) #1409 & 1358 populated grid cells pre- and post- 2012, respectively

goodco <- joinedgrid %>% map(~ .[!is.na(.$GridID),]) #select records that were in a grid cell 
map(goodco, nrow) # 1529 & 1489 sets were in the Maritime Study Area pre- and post-2012 respectively

#Get site x species matrix -- Group the observations by GridID and 
#ID the cells as Presence - Absence

IDvars <-  c("ID","year","month","day","set","strat","season","vertical.position","GridID")
Species <- map(goodco, ~ setdiff(names(.),IDvars))


#Summarize by gridID
SiteSummary <- map2(goodco, Species, ~ dplyr::select(.x, c("GridID",.y))) %>%
  map(.,~ group_by(.,GridID)) %>%
  map(., ~ summarise_all(.,funs(sum(.)))) %>% 
  map(.,ungroup)%>% map(.,data.frame)

#Determine which grid cells had trawl sets before and after 2012
overlap <- intersect(SiteSummary[[1]]$GridID, SiteSummary[[2]]$GridID) #238 overlapping sites
saveRDS(overlap, file = 'Data/Maritimes_overlappingSite.rds')
#Create 2 versions of SiteXSpecies dataframe (1 using pre 2012 data for overlapping sites, 1 using post 2012 data for overlapping sites)

SiteSummary2 <- vector('list',2)
SiteSummary2[[1]] <- bind_rows(SiteSummary[[1]], SiteSummary[[2]][-which(SiteSummary[[2]]$GridID %in% overlap),])  
SiteSummary2[[1]][is.na(SiteSummary2[[1]])] <- 0
SiteSummary2[[2]] <- bind_rows(SiteSummary[[1]][-which(SiteSummary[[1]]$GridID %in% overlap),], SiteSummary[[2]]) 
SiteSummary2[[2]][is.na(SiteSummary2[[2]])] <- 0

sp.order <- sort(names(SiteSummary2[[1]])[-1]) #ordered species names encompassing both time periods

SiteSummary2  <-  map(SiteSummary2, ~dplyr::select(.,GridID,sp.order)) #reorder columns by species names

#Convert to presence/absence
Indval <- map(SiteSummary2, ~.$GridID)
SpeciesSiteSummary <- map(SiteSummary2, ~dplyr::select(., -GridID)) %>% 
  map(., ~transmute_all(.,as.logical)) %>% 
  map(., ~transmute_all(.,as.integer))

## Trim out barren sites -- those sites with only one species

map(SpeciesSiteSummary, ~unique(unlist(.))) ## just 0 and 1 so we can now add across rows to figure which Grids are barren (1 or 0 species)
speccount <- map(SpeciesSiteSummary, ~ rowSums(.))
map(speccount, range) ## 1 - 39 per site
map(speccount, table) ## 2 barren sites to be removed

BarrenSites <- map2(Indval, speccount, ~ .x[.y<2])
IncludedSites <- map2(Indval, BarrenSites, ~setdiff(.x,.y))

##############################
# drop "Rare species"
############################

## Trim out species that are only observed in fewer than 1% of stations

dropThreshold <- map(SpeciesSiteSummary, ~(ceiling(nrow(.)*0.01)))

map2(dropThreshold, SpeciesSiteSummary, ~print(paste("Dropthreshold is",.x,"out of",nrow(.y),"sites, or ~1%")))

SpeciesCountSummary <- map(SpeciesSiteSummary, colSums)

removedSpecies <- map2(SpeciesCountSummary,dropThreshold, ~ names(.x)[which(.x<.y)])
removedSpecies #51 removed species pre 2012, 48 removed post 2012
IncludedSpecies <- map2(SpeciesSiteSummary,removedSpecies,~ setdiff(names(.x),.y))
map(IncludedSpecies,length) #109 unique species pre 2012, 112 post 2012

### Trim the data ----------

rnames <- map2(SiteSummary2, IncludedSites, ~ .x[.x$GridID %in% .y, 'GridID'])

Output <- pmap(list(SiteSummary2,IncludedSites,IncludedSpecies), ~ ..1[..1$GridID %in% ..2, names(..1) %in% ..3]) %>% 
  map(., ~transmute_all(.,as.logical)) %>% #convert to presence (1) and absence (0)
  map(., ~transmute_all(.,as.integer)) %>% 
  map(.,~ dplyr::select(.,-Clupea.harengus)) %>%  #remove Atlantic herring
  map2(.,rnames,  ~ `rownames<-`(x = .x,value = .y)) #give gridID to rownames

##Save the data
write.csv(Output[[1]], "Data/ClusterData4kmPre2012.csv", row.names=T)
write.csv(Output[[2]], "Data/ClusterData4kmPost2012.csv", row.names=T)

####Examine multivariate relationship between pre and post 2012 data at sites sampled in both time periods####

#combine overlapping sites into one dataframe of presence/absence data
OverlappingSites <- bind_rows(SiteSummary[[1]][which(SiteSummary[[1]]$GridID %in% overlap),], SiteSummary[[2]][which(SiteSummary[[2]]$GridID %in% overlap),])  
OverlappingSites[is.na(OverlappingSites)] <- 0
rownames(OverlappingSites) <- paste(OverlappingSites$GridID, c(rep('pre2012',238),rep('post2012',238)),sep = "-")
OverlappingSites <- OverlappingSites[-1] %>% 
  transmute_all(.,as.logical) %>% 
  transmute_all(.,as.integer) %>% 
  dplyr::select(.,-Clupea.harengus)

#compute simpson dissimilarity matrix

library(simba)
library(vegan)

Overlap_diss <- sim(OverlappingSites, method = 'simpson')
OverlapMDS <- nmds(Overlap_diss, k = 2)
mdspts <- (OverlapMDS$points)
plot(mdspts[1:238,1],mdspts[1:238,2],col = 'red')
points(mdspts[239:476,1],mdspts[239:476,2],col = 'blue', add = T)

####Hierarchical clustering####

#load required packages

pkgs <- list('purrr','vegan', 'simba', 'maptools', 'pvclust', 'dendroextras', 'dendextend', 'reshape', 'reshape2', 'dplyr', 'tidyr', 'NbClust', 'rgdal')
invisible(lapply(pkgs, library, character.only = T))

#read in data if needed

SxSpfiles <- list.files(path = 'Data/', pattern = 'Cluster.+2012', full.names = T)
SiteXSpecies<- map(SxSpfiles, ~ read.csv( ., stringsAsFactors = F, row.names = 1))

#Compute dissimilarity matrices and cluster data using hierarchical clustering with average linkages
#Compare multiple measures of dissimilarity

diss.measure <-c('simpson','soerensen','jaccard','gower')

dist.matrix <- cross2(diss.measure, SiteXSpecies) %>% #create list of all combinations of pre/post 2012 data and dissimilarity measures
  map(., ~sim(.[[2]], method = .[[1]])) #compute dissimilarity matrix for each combination

dendro.obj <- map(dist.matrix, ~hclust(., method = 'average')) #create dendrogram object from each diss. matrix


#Calculate cophenetic correlation value for each tree
#strongest correlation using simpson coefficient. Use simpson.
copCor <- map(dendro.obj,cophenetic) %>% 
  map2(dist.matrix,., ~ paste(attributes(.x)$method, cor(.x,.y), sep = ' ')) %>% flatten_chr(.) %>% 
  print()

#compare fit of multiple linkage methods

lnkmtd <- c('ward.D', 'ward.D2','sin','com', 'ave','mcq','med','cen') # 8 common cluster algorithms

copCor.lnkmtd <- map(SiteXSpecies, ~ sim(., method = 'simpson')) %>% #create dissimilarity matrices for pre/post 2012 subsets
  cross2(lnkmtd, .) %>% #create list of all combinations input diss matrix and linkage method
  map(., ~ list(dm = .[[2]], do = hclust(.[[2]], method = .[[1]]))) %>% #for each combination, make list of [[1]] diss matrix, [[2]] dendrogram object using that linkage method   
  map_chr(., ~paste(.$do$method,cor(.$dm, cophenetic(.$do)), sep = ' ')) %>% print() #calculate and print cophenetic correlation for each combination
# average linkage method is best 

#Create dendrogram object using Simpson's distance and average linkage method for both pre & post 2012 subsets

benthtree <- map(SiteXSpecies, ~ sim(., method = 'simpson')) %>% 
  map(., ~hclust(.,method = 'average'))

#save trees
saveRDS(benthtree[[1]], file="Data/benthtree_Post2012.rds")
saveRDS(benthtree[[2]], file='Data/benthtree_Pre2012.rds')

####Choose height cut-off for cluster membership

pkgs <- list('vegan', 'simba', 'maptools', 'pvclust', 'dendroextras', 'dendextend', 'reshape', 'reshape2', 'dplyr', 'tidyr', 'NbClust','ggplot2')
invisible(lapply(pkgs, library, character.only = T))

benthtreePre2012<-readRDS("Data/benthtree_Pre2012.rds")
SiteXSpeciesPre2012<-read.csv("Data/ClusterData4kmPre2012.csv", stringsAsFactors = F, row.names = 1)
dissPre2012<-sim(SiteXSpeciesPre2012,  method='simpson')

benthtreePost2012<-readRDS("Data/benthtree_Post2012.rds")
SiteXSpeciesPost2012<-read.csv("Data/ClusterData4kmPost2012.csv", stringsAsFactors = F, row.names = 1)
dissPost2012<-sim(SiteXSpeciesPost2012,  method='simpson')

#Examine internal cluster validity index as a function of cluster size
#calculate various CVI's for partititions with k=2-20 

cvi <- c('ch', 'cindex', 'ptbiserial','db','silhouette')
nclust <- c(2:20)
vals <- list()

for (i in cvi) {
  
  vals[[i]] <- NbClust(data = SiteXSpeciesPost2012, diss = dissPost2012, distance = NULL, 
                       method = 'average',index = i, min.nc = 2, max.nc = 20)$All.index
  
}

cvi.df <- data.frame(NC = rep(nclust, length(cvi)), Value = unlist(vals))
cvi.df$Index <- gsub('[.][0-9]+','',row.names(cvi.df))
p.cvi <- ggplot(cvi.df) +
  geom_point(aes(x = NC, y = Value)) +
  labs(x = 'Number of clusters', y = 'Index value') +
  facet_wrap(~ Index, scales = 'free_y')+
  theme_bw(); p.cvi

#2, 4, 7, 11 clusters all possible Pre 2012
#3,5-6,#10-11, 13 all possible Post 2012

#Examine top clusters at various splits
ncells1<-list()
ncells2<-list()
ncells3<-list()
ncells4<-list()
ncells5<-list()
ncells6<-list()
ncells7<-list()
ncells8<-list()
ncells9<-list()
ncells10<-list()
nClusters <- list()

s<-seq(30,100, by=0.1)
for (i in 1:length(s)){
  clus<-dendroextras::slice(benthtreePost2012, h=(s[i]/100))
  cluscount<-as.data.frame(table(clus))
  cluscount2<-cluscount[order(-cluscount$Freq),]
  nClusters[i] <- (length(cluscount2$Freq))
  ncells1[i]<-sum(cluscount2$Freq[1:1])
  ncells2[i]<-sum(cluscount2$Freq[2:2])
  ncells3[i]<-sum(cluscount2$Freq[3:3])
  ncells4[i]<-sum(cluscount2$Freq[4:4])
  ncells5[i]<-sum(cluscount2$Freq[5:5])
  ncells6[i]<-sum(cluscount2$Freq[6:6])
  ncells7[i]<-sum(cluscount2$Freq[7:7])
  ncells8[i]<-sum(cluscount2$Freq[8:8])
  ncells9[i]<-sum(cluscount2$Freq[9:9])
  ncells10[i]<-sum(cluscount2$Freq[10:10])
}
table<-cbind.data.frame(s/100, unlist(ncells1),unlist(ncells2),unlist(ncells3),unlist(ncells4),unlist(ncells5),unlist(ncells6),unlist(ncells7),unlist(ncells8),unlist(ncells9),unlist(ncells10),unlist(nClusters))
colnames(table)<-c("distance_Bsim", "ncells_top_1","ncells_top_2","ncells_top_3","ncells_top_4","ncells_top_5","ncells_top_6","ncells_top_7","ncells_top_8","ncells_top_9","ncells_top_10","nClusters")

#Choose Bsim 0.607, 11 clusters for pre 2012
#Choose Bsim 0.607), 13 clusters for post 2012

##slice trees

dendro_files <- list.files('Data/', pattern = 'benth.+2012', full.names = T)
benthtree <- map(dendro_files, readRDS)
seth<-0.607 #use my choice - more informative biologically
cl<-map(benthtree, ~ dendroextras::slice(., h=seth)) #Cut tree

###
colorcount<- map(cl, ~ as.data.frame(table(.))) #get Table of cluster memberships
colorcount #Number of sites in each cluster # 11 & 13 clusters with 0.607 cut-off for Pre & Post 2012 data, respectively
colorcount <-map(colorcount, ~ .[order(-.$Freq),]) %>% #order by frequency
  map(., ~ dplyr::rename(., cl = '.'))
colorcount
setfreq<-map(colorcount, ~ min(.$Freq[1:6], na.rm=T)) #Select the first 6 clusters to be color-coded

#Count the number of sites in the top clusters
map(colorcount, ~ sum(.$Freq[1:6])) # 2505 & 2504 sites in top 6 clusters at 0.607 cut-off for Pre & Post 2012 data, respectively

# #Set Colors for top 6 clusters (SubBiomes)
colorscheme<- list(c("#a6cee3","#1f78b4","#ff7f00", "#33a02c", "#6a3d9a", "#e6ab02"), c("#1f78b4","#a6cee3","#ff7f00", "#33a02c", "#6a3d9a", "#e6ab02")) 

#Assign colours to clusters

for (i in 1:length(colorcount)){
  colorcount[[i]]$assigned[colorcount[[i]]$Freq<setfreq[[i]]]<-"#e0e0e0" #set all clusters that occur at less than the assigned cut off (6 clusters) to grey
  colorcount[[i]]$assigned[1:6] <- colorscheme[[i]]
  colorcount[[i]] <- colorcount[[i]][order(colorcount[[i]]$cl),]
}

colorcount #Each cluster is assigned a color
map(colorcount, ~ table(.$assigned))

#Create color-coded dendrogram
colortree<- map2(benthtree, colorcount, ~ colour_branches(.x, h=seth,col=as.character(.y$assigned)))

#Plot color-coded dendrogram
plot(colortree[[2]],xlab="", ylab="Simpson Distance", leaflab = "none") #Pre 2012 dendrogram
plot(colortree[[1]],xlab="", ylab="Simpson Distance", leaflab = "none") #Post 2012 dendrogram

###plot this cluster analysis on map
#import the subsetted grid (grid cells for which we had data) that we wrote earlier
grid<-readShapeSpatial("Data/MaritimeGrids/Maritimes_4km_Grid.shp") 
head(grid@data)

#Join colors and clusters to map
colorcountmap<-colorcount
id <-map(cl, ~ as.integer(names(.)))
colors<-map2(id, cl, ~as.data.frame(cbind(.x, .y))) %>% 
  map(., ~ dplyr::rename(., id = .x, cl = .y))
map(colors, head) #each grid cell is assigned its cluster
colorsmerge<-map(colors, ~ merge(grid, ., by="id"))
map(colorsmerge, ~head(.@data)) #join this to the shapefile
colorsmerge2<- map2(colorsmerge, colorcountmap, ~ merge(.x, .y, by="cl"))
map(colorsmerge2, ~head(.@data))#each grid cell is assigned its color to match the dendrogram
colorsgrid<- map(colorsmerge2, ~subset(.,!is.na(.$cl))) %>% #drop grid cells for which we had data, but which got dropped at the rare species/barren sites stage
  map(., ~`proj4string<-`(., CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))) #set coordinate reference system
map(colorsgrid, ~head(.@data)) 
colorGridData<- map(colorsgrid, ~.@data)

#Plot color-coded grid
par(mfrow = c(1,2), usr = c(76051.06, 1065725.49, 4502270.20, 5300000))
Land <- readOGR('Data/Shapefiles/Maritimes_prov_borders.shp')
strata <- readOGR('Data/Shapefiles/MaritimesRegionStrataAgg.shp')
plot(strata) #2007-2011
plot(crop(Land,extent(76051.06, 1065725.49, 4502270.20, 5300000)), add = T, col = 'grey')
plot(colorsgrid[[2]], col = as.character(colorsgrid[[2]]$assigned), border = NA, add = T)
plot(extent(76051.06, 1065725.49, 4502270.20, 5300000), add = T)
plot(strata) #2012-2016
plot(crop(Land,extent(76051.06, 1065725.49, 4502270.20, 5300000)), add = T, col = 'grey')
plot(colorsgrid[[1]], col = as.character(colorsgrid[[1]]$assigned), border = NA, add = T)
plot(extent(76051.06, 1065725.49, 4502270.20, 5300000), add = T)

#Examine frequency of clusters for those grid cells overlapping both time periods


colors[[2]]$cl <- recode(colors[[2]]$cl, `6`='WSS: Banks/Inner BoF',`7` = 'WSS/Outer BoF',
                         `8` = 'ESS: Banks',`11`='ESS',`10`= 'LC/Shelf Break', `1`='Slope', .default ='Other' )

colors[[1]]$cl <- recode(colors[[1]]$cl, `11`='WSS: Banks/Inner BoF',`9` = 'WSS/Outer BoF',
                         `2` = 'ESS: Banks',`4`='ESS',`3`= 'LC/Shelf Break', `1`='Slope', .default ='Other' )

freqCl <- bind_rows(colors[[2]][colors[[2]]$id %in% overlap, ],
                    colors[[1]][colors[[1]]$id %in% overlap, ],.id = 'Period') #create df w/ grid Id and cluster membership
freqCl$Period <- recode(freqCl$Period, '1' = '2007-2011', '2'='2012-2016') #recode groups
saveRDS(freqCl, 'Output/Maritimes_OverlappingSites_Assignment.rds')#save cluster assignments for overlapping grid cells
freqCl <- freqCl %>% group_by(Period, cl) %>% 
  dplyr::count(., sort = F) %>% bind_rows(.,data.frame(Period='2012-2016',cl='Slope',n=0)) %>% 
  mutate(cl = factor(cl, levels = c('Slope','LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF','WSS: Banks/Inner BoF')))

#Plot
my.colors <- MAR.palette$assigned[c(3,4,5,6,1,2)] #color scheme
names(my.colors) <- unique(freqCl$cl)

p1 <- ggplot(freqCl) + 
  geom_bar(aes(x = cl, y = n, group = Period, fill = cl, alpha = Period, color = Period), stat = 'identity', position = position_dodge2(padding=0.1)) +
  labs(x = NULL, y = 'Number of sites') +
  theme_classic() +
  scale_alpha_manual(values = c(1,0.75)) + 
  scale_fill_manual(name = 'Assemblage', values = my.colors) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        text = element_text(size = 14),
        panel.border = element_rect(fill = NA, colour = 'black')) +
  scale_color_manual(values = c('black','black'))

ggsave(plot = p1, 'Output/Maritimes_splitComparison.tiff', width = 6, height = 4.5, units = 'in', dpi = 300, compression = 'lzw')

####Split Environmental Layers####

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
library(sf)

#functions to be used to calculate annual summaries
rowMin <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  min = min(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12))
  return(min)}

rowMax <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  max = max(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12))
  return(max)}

#Objects and functions required to convert pt shapefiles of environmental
#variable summaries to raster format (4km resolution)

source("Code/VoronoiResample.R")

Spth <- 'U:/MPA_group/Data/Shapefiles/'
Rpth <- 'U:/MPA_group/Data/Rasters/'
extentMar <- extent(-68,-55,41,48)
MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")

#Create list of files with bottom temperature and salinity data
#Data provided as monthly predictions for pts on irregular grid extending from SW NS halfway up Labrador (2007-2015) 

bnam.files_period1 <- list.files(path = "Data/Zeliang_OceanographicData/", full.names = T,pattern = "(20[01][78901]).+Bottom_UVTS")
names(bnam.files_period1) <- "Pre2012"

bnam.files_period2 <- list.files(path = "Data/Zeliang_OceanographicData/", full.names = T,pattern = "(201[2-5]).+Bottom_UVTS")
names(bnam.files_period2) <- "Post2012"

#Bottom Temperature

filename <- "BT_NWA_pts_pre2012"
varname <- "BT_pre2012"

df.long.pre2012 <- bnam.files_period1 %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20[01][78901])", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T))) %>% 
  data.table::rbindlist(.) %>% 
  data.frame() %>% 
  select(-c('V3','V4','V6')) %>% #7813260 observations of 5 variables 
  `names<-`(.,c("x","y",varname,"year","month")) 
  df.long.pre2012$BT_pre2012[which(df.long.pre2012$BT_pre2012==0)] <- NA
  
df.wide.pre2012 <- unite(df.long.pre2012, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(max_ann = pmap_dbl(.,rowMax),
         min_ann = pmap_dbl(.,rowMin)) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame() %>% 
  `coordinates<-`(.,c('x','y')) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(df.wide.pre2012, dsn = 'U:/MPA_group/Data/Shapefiles', layer = filename, driver = "ESRI Shapefile")

VoronoiResample(Spth,Rpth,paste0(filename,'.shp'),extentMar,MaritimesSA,"Mar",4000,varname) #Pt shapefile to raster with 4 km res

filename <- "BT_NWA_pts_post2012"
varname <- "BT_post2012"

df.long.post2012 <- bnam.files_period2 %>% 
  map(~ cbind(fread(.),strapplyc(.,"(201[2-5])", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T))) %>% 
  data.table::rbindlist(.) %>% 
  data.frame() %>% 
  select(-c('V3','V4','V6')) %>% #6250608 observations of 5 variables 
  `names<-`(.,c("x","y",varname,"year","month")) 
  df.long.post2012$BT_post2012[which(df.long.post2012$BT_post2012==0)] <- NA

df.wide.post2012 <- unite(df.long.post2012, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(max_ann = pmap_dbl(.,rowMax),
         min_ann = pmap_dbl(.,rowMin)) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame() %>% 
  `coordinates<-`(.,c('x','y')) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(df.wide.post2012, dsn = 'U:/MPA_group/Data/Shapefiles', layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

VoronoiResample(Spth,Rpth,paste0(filename,'.shp'),extentMar,MaritimesSA,"Mar",4000,varname) #Pt shapefile to raster with 4 km res

#Max Annual Bottom Salinity

filename <- "BSal_NWA_pts_pre2012"
varname <- "BSal_pre2012"

df.long.pre2012 <- bnam.files_period1 %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20[01][78901])", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T))) %>% 
  data.table::rbindlist(.) %>% 
  data.frame() %>% 
  select(-c('V3','V4','V5')) %>% #7813260 observations of 5 variables 
  `names<-`(.,c("x","y",varname,"year","month")) 
  df.long.pre2012$BSal_pre2012[which(df.long.pre2012$BSal_pre2012==0)] <- NA

df.wide.pre2012 <- unite(df.long.pre2012, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(max_ann = pmap_dbl(.,rowMax)) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame() %>% 
  `coordinates<-`(.,c('x','y')) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(df.wide.pre2012, dsn = 'U:/MPA_group/Data/Shapefiles', layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

VoronoiResample(Spth,Rpth,paste0(filename,'.shp'),extentMar,MaritimesSA,"Mar",4000,varname) #Pt shapefile to raster with 4 km res

filename <- "BSal_NWA_pts_post2012"
varname <- "BSal_post2012"

df.long.post2012 <- bnam.files_period2 %>% 
  map(~ cbind(fread(.),strapplyc(.,"(201[2-5])", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T))) %>% 
  data.table::rbindlist(.) %>% 
  data.frame() %>% 
  select(-c('V3','V4','V5')) %>% #6250608 observations of 5 variables 
  `names<-`(.,c("x","y",varname,"year","month")) 
  df.long.post2012$BSal_post2012[which(df.long.post2012$BSal_post2012==0)] <- NA

df.wide.post2012 <- unite(df.long.post2012, "lonlat",c("x","y"),sep = "/") %>% 
  spread(month,varname,drop = T) %>%  
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate(max_ann = pmap_dbl(.,rowMax)) %>% 
  select(-(year:M_9)) %>% 
  group_by(x,y) %>% 
  summarise_all(funs(mean)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.double) %>% 
  data.frame() %>% 
  `coordinates<-`(.,c('x','y')) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(df.wide.post2012, dsn = 'U:/MPA_group/Data/Shapefiles', layer = filename, driver = "ESRI Shapefile", overwrite_layer = T)

VoronoiResample(Spth,Rpth,paste0(filename,'.shp'),extentMar,MaritimesSA,"Mar",4000,varname) #Pt shapefile to raster with 4 km res

####Random Forest Predictions pre and post 2012####

library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(dplyr)
library(tidyr)
library(purrr)
library(randomForest)
library(ggplot2)
library(gridExtra)
library(gtable)
library(caret)
library(pROC)

#Read in required files (model object, shapefiles, rasters)

RF.mod <- readRDS('Output/Maritimes_RF_model.rds') #random forest model object
time_periods <- c('pre2012','post2012')
Mar_SA <- readOGR("Data/Shapefiles/MaritimesRegionStrataAgg.shp") # Mar_SA (RV survey boundaries)
LandBuffer <- readOGR('Data/Shapefiles/Maritimes_LandBuffer_5km.shp') # 5km buffer around land points

raster.list <- map(time_periods, ~ list.files(path = paste0('U:/MPA_group/Data/Rasters/',.), 
                                              pattern = 'Mar', full.names = T)) #get list of raster files
env_pred <- map(raster.list, ~stack(.)) %>%  #stack rasters in list
  map(., ~ `names<-`(.,gsub('Mar_','',names(.)))) %>% #remove 'Mar_' from layer names
  map(., ~ `names<-`(.,gsub('_p.+2012','',names(.)))) %>% #remove time period suffix from layer names
  map(., ~mask(., Mar_SA)) %>%  #mask raster cells of predictors outside RV survey boundaries
  map(., ~mask(., LandBuffer, inverse = T)) #also mask cells overlapped by 5km land buffer

#Predict cluster membership for both time periods over entire Scotian Shelf
predict.map <- map(env_pred, ~raster::predict(., RF.mod)) # predict cluster membership over raster surfaces
my.colors = c('#e6ab02','#6a3d9a','#33a02c','#ff7f00', '#a6cee3','#1f78b4') # color scheme for rasters
par(mfrow = c(1,2))
map(predict.map, ~plot(., col = my.colors))
writeRaster(predict.map[[1]], filename = 'Output/Maritimes_Hindcast2007_2011_map.tif', overwrite = T) #write predictions to file  
writeRaster(predict.map[[2]], filename = 'Output/Maritimes_Hindcast2012_2015_map.tif', overwrite = T) #write predictions to file 

#Predict cluster membership for cells overlapping both time periods only

load("Data/MaritimeGridsWholeCell.RData") #load Maritimes grid data
MaritimeGrid <- Griddata[[20]] # 4 km resolution grid
MaritimeGrid@data$id <- 1:nrow(MaritimeGrid@data) # assign grid cells id
overlapID <- readRDS('Data/Maritimes_overlappingSite.rds')# GridID of cells overlapping both periods
overlapCells <- MaritimeGrid[MaritimeGrid@data$id %in% overlapID,] #spatial polygons for overlapping grid cells
overlap.data <- map(env_pred, ~raster::extract(.,overlapCells, factors = T, nl = nlayers(.), df = T)) %>% #extract data for overlapping cells
  map(., ~ bind_cols(id = overlapCells@data$id, dplyr::select(.,-ID))) %>% #add grid IDs to dataframe
  map(., data.frame)  
overlap.predictions <- map(overlap.data, ~predict(object = RF.mod, newdata = ., type = 'response')) %>% #predicted cluster assignment
  map(., ~recode_factor(., '1'='Slope', '2'= 'LC/Shelf Break','3'='ESS','4' = 'ESS: Banks',
                        '5' = 'WSS/Outer BoF','6'='WSS: Banks/Inner BoF')) %>%  #recode clusters to match observations
  map(.,data.frame)
overlap.prob <- map(overlap.data, ~predict(object = RF.mod, newdata = ., type = 'prob')) %>% #model predictions for overlapping cells as probabilities
  map(.,data.frame) %>% #convert to dataframe
  map(., ~ `colnames<-`(.,as.character(c(1:6)))) #rename columns
  
#Compare predictions with cluster assignments from hierarchical clustering

assignments <- readRDS('Output/Maritimes_OverlappingSites_Assignment.rds') %>% #read in cluster assignments
  arrange(.,id) %>% #order rows by grid ID ascending
  arrange(.,Period) %>% #group by time period
  mutate_at(.,'cl',as.factor) %>% #convert cluster to factor
  bind_cols(.,data.table::rbindlist(overlap.predictions)) %>% #bind observed cluster assignments with model predictions
  dplyr::rename(.,observed = cl, predicted = .x..i..) %>% #rename columns
  filter_all(.,all_vars(!is.na(.))) #remove incomplete cases (grid cells without complete suite of environmental predictors)

#reorder factor levels of predicted assignment to match order of observed factor levels

assignments$predicted <- factor(assignments$predicted, levels = levels(assignments$observed))

confusionHindcast <- confusionMatrix(assignments$predicted, assignments$observed) #compute confusion matrix for all data
#accuracy = 74.8%, more accurately predicts true positives for slope, and banks clusters
saveRDS(confusionHindcast, 'Output/ConfusionMatrixAllTimes.rds') #save matrix
#confusion matrix for 2007-2011 only (71.2% accuracy)
confusionHindcastpre2012 <- confusionMatrix(data = assignments$predicted[assignments$Period == '2007-2011'], 
                                            reference = assignments$observed[assignments$Period == '2007-2011'])
saveRDS(confusionHindcastpre2012 , 'Output/ConfusionMatrixPre2012.rds') #save matrix
#confusion matrix for 2012-2016 only (78.4 % accuracy)
confusionHindcastpost2012 <- confusionMatrix(data = assignments$predicted[assignments$Period == '2012-2016'], 
                                            reference = assignments$observed[assignments$Period == '2012-2016'])
saveRDS(confusionHindcastpost2012, 'Output/ConfusionMatrixPost2012.rds') #save matrix

#Evaluate accuracy of predictions with AUC

overlap.prob <- map(overlap.data, ~predict(object = RF.mod, newdata = ., type = 'prob')) %>% #model predictions for overlapping cells as probabilities
  map(.,data.frame) %>% #convert to dataframe
  map(., ~ `colnames<-`(.,c('Slope', 'LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF','WSS: Banks/Inner BoF'))) %>% #rename columns to match output of RF model
  data.table::rbindlist(., idcol = 'time_period') %>% #row bind probabilities for pre 2012 and post 2012 data 
  filter_all(., all_vars(!is.na(.))) #remove rows with NA's
  
mROC <- multiclass.roc(assignments$observed, dplyr::select(overlap.prob,-time_period)) #calculate multiclass ROC (and AUC)
mROC$auc #multiclass AUC = 0.938. Therefore RF model distinguishes between clusters well in hindcast

####Future environmental layers####

#load required packages
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(raster)
library(data.table)
library(rgdal)
library(rgeos)
library(sp)
library(R.matlab)
library(gsubfn)

#functions to be used to calculate annual summaries
rowMin <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  min = min(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12))
  return(min)}

rowMax <- function(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12,...){
  max = max(c(M_1,M_2,M_3,M_4,M_5,M_6,M_7,M_8,M_9,M_10,M_11,M_12))
  return(max)}

#read in anomolies from present climatology
BNAM_MAT <- readMat('Data/dTbtm_dSbtm_for_RS.mat') #list with lon, lat, monthly btm temp/sal in separate matrices or arrays

col.names <- paste0('M_', 1:12) #monthly col names
lon <- BNAM_MAT[[1]] %>% `dim<-`(.,NULL)  #flatten longitude array 
lat <- BNAM_MAT[[2]] %>% `dim<-`(.,NULL)  #flatten latitude array 

#Isolate btm temp and sal data, reshape, and convert to spatial df

fut_layers <- list(BNAM_MAT[[8]], BNAM_MAT[[9]]) %>% 
  map(., ~aperm(.,c(2,3,1))) %>%  #reshape array to place month index in 3rd dimension
  map(., ~`dim<-`(.,c(321201,12))) %>% #convert to 2d array, 12 cols (months) per location  
  map(., ~`dimnames<-`(.,list(NULL,col.names))) %>% #rename columns
  map(., ~data.frame(x = lon, y = lat, .)) %>%  #cbind coordinates for each location
  map(., ~`coordinates<-`(., c('x','y'))) %>% #convert to spatial df
  map(., ~`proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))) #set CRS

fut_bt <- fut_layers[[1]] #btm temperature monthly anomalies
fut_bsal <- fut_layers[[2]] #btm salinity monthly anomalies
fut_bt$id <- 1:nrow(fut_bt) #create id variable for each location
fut_bsal$id <- 1:nrow(fut_bsal) #create id variable for each location

#Since geometries of present and future climate are slightly different in precision, spatial join
#or attribute join based on lon/lat will not work
#Instead, create voronoi polygons around future pts and use spatial overlay to join present pts with future data

source("Code/Voronoipolygons2.R") #function to create voronoi polygons around points

vor_bt <- voronoipolygons(fut_bt) #voronoi polygons for btm temp
proj4string(vor_bt) <- proj4string(fut_bt) #set CRS

vor_bsal <- voronoipolygons(fut_bsal) #voronoi polygons for btm sal
proj4string(vor_bsal) <- proj4string(fut_bsal) #set CRS

#Import data for present climatology

#Create list of files with bottom temperature and salinity contemporary data
#Data provided as monthly predictions for pts on irregular grid extending from SW NS halfway up Labrador (2007-2015) 

bnam.files <- list.files(path = "Data/Zeliang_OceanographicData/", full.names = T,pattern = "Bottom_UVTS")

#Objects and functions required to convert pt shapefiles of environmental
#variable summaries to raster format (4km resolution)

source("Code/VoronoiResample.R")

Spth <- 'U:/MPA_group/Data/Shapefiles/'
Rpth <- 'U:/MPA_group/Data/Rasters/FutureClimate/'
extentMar <- extent(-68,-55,41,48)
MaritimesSA <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")

#Bottom Salinity

filename <- "Bsal_NWA_pts_future"
varname <- "BSal"

#create long format data frame of salinity data from text files
df.long <- bnam.files %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20..)", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T))) %>% 
  data.table::rbindlist(.) %>% 
  data.frame() %>% 
  dplyr::select(-c('V3','V4','V5')) %>% #14063868 observations of 5 variables 
  `names<-`(.,c("x","y",varname,"year","month")) 
df.long$BSal[which(df.long$BSal==0)] <- NA

#Reshape to wide format (months in columns) take monthly means across years and convert to spatial df
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>%  
  spread(month,varname,drop = T) %>%
  group_by(lonlat) %>% 
  summarise_at(.,vars(M_1:M_9), mean, na.rm = T) %>% 
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate_at(., c("x","y"), as.numeric) %>% 
  `coordinates<-`(.,c("x","y")) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) 

#Extract anomalies from voronoi polygons that present climatology pts overlay

anomalies <- over(df.wide, vor_bsal)

#rbind future anomalies with present climatology df and add anomalies to monthly values 

df.wide$id <- anomalies$id
df.wide2 <- bind_rows(dplyr::select(df.wide@data, colnames(anomalies)),anomalies) %>% 
  bind_cols(x = rep(df.wide@coords[,1],2), y = rep(df.wide@coords[,2],2),.) %>% 
  group_by(id,x,y) %>% 
  summarise_at(., vars(contains('M_')),sum, na.rm = T) %>% 
  ungroup() %>% 
  mutate(max_ann = pmap_dbl(.,rowMax)) %>% 
  dplyr::select(x,y,max_ann) %>%  
  data.frame() %>% 
  `coordinates<-`(., c("x","y")) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(df.wide2, dsn = 'U:/MPA_group/Data/Shapefiles', layer = filename, driver = "ESRI Shapefile")

VoronoiResample(Spth,Rpth,paste0(filename,'.shp'),extentMar,MaritimesSA,"Mar",4000,varname) #Pt shapefile to raster with 4 km res

#Bottom Temperature

filename <- "BT_NWA_pts_future"
varname <- "BT"

#create long format data frame of temperature data from text files
df.long <- bnam.files %>% 
  map(~ cbind(fread(.),strapplyc(.,"(20..)", simplify = T),strapplyc(.,"(M_.[0-9]?)",simplify = T))) %>% 
  data.table::rbindlist(.) %>% 
  data.frame() %>% 
  dplyr::select(-c('V3','V4','V6')) %>% #14063868 observations of 5 variables 
  `names<-`(.,c("x","y",varname,"year","month")) 
df.long$BT[which(df.long$BT==0)] <- NA

#Reshape to wide format (months in columns) take monthly means across years and convert to spatial df
df.wide <- unite(df.long, "lonlat",c("x","y"),sep = "/") %>%  
  spread(month,varname,drop = T) %>%
  group_by(lonlat) %>% 
  summarise_at(.,vars(M_1:M_9), mean, na.rm = T) %>% 
  separate(.,lonlat,c("x","y"), sep = "/") %>% 
  mutate_at(., c("x","y"), as.numeric) %>% 
  `coordinates<-`(.,c("x","y")) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) 

#Extract anomalies from voronoi polygons that present climatology pts overlay

anomalies <- over(df.wide, vor_bt)

#rbind future anomalies with present climatology df and add anomalies to monthly values 

df.wide$id <- anomalies$id
df.wide2 <- bind_rows(dplyr::select(df.wide@data, colnames(anomalies)),anomalies) %>% 
  bind_cols(x = rep(df.wide@coords[,1],2), y = rep(df.wide@coords[,2],2),.) %>% 
  group_by(id,x,y) %>% 
  summarise_at(., vars(contains('M_')),sum, na.rm = T) %>% 
  ungroup() %>% 
  mutate(max_ann = pmap_dbl(.,rowMax),
         min_ann = pmap_dbl(.,rowMin)) %>% 
  dplyr::select(x,y,max_ann, min_ann) %>%  
  data.frame() %>% 
  `coordinates<-`(., c("x","y")) %>% 
  `proj4string<-`(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(df.wide2, dsn = 'U:/MPA_group/Data/Shapefiles', layer = filename, driver = "ESRI Shapefile")

VoronoiResample(Spth,Rpth,paste0(filename,'.shp'),extentMar,MaritimesSA,"Mar",4000,varname) #Pt shapefile to raster with 4 km res

####Random forest predictions for 2075 RCP85 scenario####

library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(dplyr)
library(tidyr)
library(purrr)
library(randomForest)
library(ggplot2)
library(gridExtra)
library(gtable)

#Read in required files (model object, shapefiles, rasters)

RF.mod <- readRDS('Output/Maritimes_RF_model.rds') #random forest model object
Mar_SA <- readOGR("Data/Shapefiles/MaritimesRegionStrataAgg.shp") # Mar_SA (RV survey boundaries)
LandBuffer <- readOGR('Data/Shapefiles/Maritimes_LandBuffer_5km.shp') # 5km buffer around land points

raster.list <- list.files(path = 'U:/MPA_group/Data/Rasters/FutureClimate',
                          pattern = 'Mar', full.names = T) #get list of raster files
env_pred <- stack(raster.list) %>%  #stack rasters in list
  `names<-`(.,gsub('Mar_','',names(.))) %>% #remove 'Mar_' from layer names
  mask(., Mar_SA) %>%  #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) #also mask cells overlapped by 5km land buffer

#Predict cluster membership for 2075 under RCP85 scenario over entire Scotian Shelf
predict.map <- raster::predict(env_pred, RF.mod) # predict cluster membership over raster surfaces
my.colors = c('#e6ab02','#6a3d9a','#33a02c','#ff7f00', '#a6cee3','#1f78b4') # color scheme for rasters
plot(predict.map, col = my.colors)
writeRaster(predict.map, filename = 'Output/Maritimes_Forecast2075_map.tif', overwrite = T) #write predictions to file  

#Identify anomalies from present classification

present.class <- raster('Output/Maritimes_PredClust_map.tif')
grid_anomalies <- predict.map != present.class #raster layer w/ binary values 1 = classification changed, 0 = no change
writeRaster(grid_anomalies, filename = 'Output/Maritimes_GridAnomalies_RCP85_2075.tif', overwrite = T)
ClimateSensitive <- rasterToPolygons(grid_anomalies, fun = function(x){x>0}, na.rm = T) #isolate Climate Sensitive areas
ClimateSensitive <- aggregate(ClimateSensitive)
plot(ClimateSensitive) #verify correct grid cells isolated
writeOGR(as(ClimateSensitive,'SpatialPolygonsDataFrame'), dsn = 'Output', layer = 'Mar_ClimateSensitive', driver = "ESRI Shapefile")

#Determine change in area of each group

source('Code/SPERA_colour_palettes.R')

present_cl <- raster('Output/Maritimes_PredClust_map.tif') #raster of present predicted cluster membership
future_cl <- raster('Output/Maritimes_Forecast2075_map.tif') #raster of future predicted cluster membership under RCP 8.5

Cl_area <- data.frame(Cluster = c('Slope','Laurentian Channel/Shelf Break', 'ESS', 'ESS: Banks', 'WSS/Outer BoF','WSS: Banks/Inner BoF'),
                      area_present = as.vector(table(raster::values(present_cl))),
                      area_future = as.vector(table(raster::values(future_cl)))) %>% #total grid cells in each cluster
  mutate_if(.,is.integer, function(x){x*16}) %>% #calculate area of each cluster
  mutate(., area_change = area_future - area_present, #calculate change in area of each cluster
         perc_change = ((area_future - area_present)/area_present)*100, #percent change in area
         relative_area = (area_future/area_present)*100) #relative size of each cluster to present conditions
Cl_area$Cluster <- factor(Cl_area$Cluster, levels = Cl_area$Cluster[order(Cl_area$perc_change)])

my.colors <- MAR.palette$assigned[c(6,5,3,4,1,2)] #color scheme
my.colors2 <- MAR.palette$assigned[c(5,4,3,6,1,2)] #color scheme

tiff('Output/Maritimes_ChangeClusterArea.tiff', width = 8, height = 8, units = 'in', res = 300, compression = 'lzw')

p1 <- ggplot(Cl_area, aes(x = Cluster, y = perc_change, fill = Cluster, color = Cluster)) +
  geom_bar(show.legend = F, stat = 'identity') +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  labs(x = NULL, y = 'Change in area (%)') +
  scale_x_discrete(labels = c('Slope','LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF', 'WSS: Banks/Inner BoF')) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1), text = element_text(size = 16), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,12,'pt'), hjust = 0.5),
        axis.ticks.length = unit(2, 'mm')) +
  scale_color_manual(values = rep('black',6)); p1
dev.off()

# alternative format
p2 <- ggplot(Cl_area, aes(x = Cluster, y = perc_change, fill = Cluster, color = Cluster)) +
  geom_bar(show.legend = F, stat = 'identity') +
  theme_classic() +
  scale_fill_manual(values = my.colors2) +
  labs(x = NULL, y = 'Change in area (%)') +
  scale_x_discrete(labels = c('LC/Shelf Break','ESS Banks','ESS','Slope','WSS/Outer BoF', 'WSS Banks/Inner BoF')) +
  coord_flip() +
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(margin = ggplot2::margin(12,0,0,0,'pt'), hjust = 0.5),
        axis.ticks.length = unit(2, 'mm')) +
  scale_color_manual(values = rep('black',6)); p2

saveRDS(p2, 'Output/Maritimes_ChangeClusterArea_flipped.rds')
ggsave(p2, filename = 'Output/Maritimes_ChangeClusterArea_flipped.tiff', width = 8, height = 8, units = 'in', dpi = 300, compression = 'lzw')

####Predict classification under RCP 8.5 in 2075 w/out depth####

RF.mod_NoDepth <- readRDS('Output/Maritimes_RF_model_NoDepth.rds') #random forest model object
env_pred <- dropLayer(env_pred,'bathy')
predict.map <- raster::predict(env_pred, RF.mod_NoDepth) # predict cluster membership over raster surfaces
my.colors = c('#e6ab02','#6a3d9a','#33a02c','#ff7f00', '#a6cee3','#1f78b4') # color scheme for rasters
plot(predict.map, col = my.colors)
writeRaster(predict.map, filename = 'Output/Maritimes_Forecast2075_map_NoDepth.tif', overwrite = T) #write predictions to file  

#Identify anomalies from present classification
present.class <- raster('Output/Maritimes_PredClust_map_NoDepth.tif')
grid_anomalies <- predict.map != present.class #raster layer w/ binary values 1 = classification changed, 0 = no change
writeRaster(grid_anomalies, filename = 'Output/Maritimes_GridAnomalies_RCP85_2075_NoDepth.tif', overwrite = T)
ClimateSensitive <- rasterToPolygons(grid_anomalies, fun = function(x){x>0}, na.rm = T) #isolate Climate Sensitive areas
ClimateSensitive <- aggregate(ClimateSensitive)
plot(ClimateSensitive) #verify correct grid cells isolated
writeOGR(as(ClimateSensitive,'SpatialPolygonsDataFrame'), dsn = 'Output', layer = 'Mar_ClimateSensitive_NoDepth', driver = "ESRI Shapefile")
