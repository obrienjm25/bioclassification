#Import data 

load("Data/NfldData.RData")

#correct redundant categories in horizontal.position and vertical.position

Nfld$horizontal.position <- gsub("brackish;marine", "brackish; marine", Nfld$horizontal.position)
Nfld$vertical.position <- gsub("pelagic-oceanic", "pelagic,oceanic", Nfld$vertical.position)

#Filter dataframe to include observations 2007 (maybe 2003 since we only have data to 2013) 
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

#Unique ID that can be used to go from long to wide format, where each row is a sample sation.
NfldFall$ID <- paste(NfldFall$year,NfldFall$month,NfldFall$day, NfldFall$longitude, NfldFall$latitude, sep="_")
length(unique(NfldFall$ID)) #3407 unique sets
NfldSpring$ID <- paste(NfldSpring$year,NfldSpring$month,NfldSpring$day, NfldSpring$longitude, NfldSpring$latitude, sep="_")
length(unique(NfldSpring$ID)) #2228 unique sets

#cast to wide format

WideNfldFall_wgt <- dcast(NfldFall[,c("ID","species","totwgt")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(NfldFall[!duplicated(NfldFall$ID),],ID,latitude,longitude,year,month,day,
                            set,strat,season,vertical.position),by="ID")%>%
  dplyr::select(.,ID,longitude,latitude,year,month,day,set,strat,season,vertical.position,unique(NfldFall$species))%>%
  data.frame()

WideNfldSpring_wgt <- dcast(NfldSpring[,c("ID","species","totwgt")],ID~species, fun.aggregate = sum)%>%
  left_join(.,dplyr::select(NfldSpring[!duplicated(NfldSpring$ID),],ID,latitude,longitude,year,month,day,
                            set,strat,season,vertical.position),by="ID")%>%
  dplyr::select(.,ID,longitude,latitude,year,month,day,set,strat,season,vertical.position,unique(NfldSpring$species))%>%
  data.frame()

#biotic data to point shapefile - convert to coordinate system that matches the grid
coordinates(WideNfldFall_wgt)<-~longitude+latitude #tell R what the coordates are
proj4string(WideNfldFall_wgt)<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system

coordinates(WideNfldSpring_wgt)<-~longitude+latitude #tell R what the coordates are
proj4string(WideNfldSpring_wgt)<-CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system

#Examine whether multivariate differences between spring and fall surveys

#reproject to UTM zone 21
Fall <- spTransform(WideNfldFall_wgt,CRS("+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
Spring <- spTransform(WideNfldSpring_wgt,CRS("+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
NAFO <- readOGR("Data/Shapefiles/NAFO_Divisions.shp")
proj4string(NAFO) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
NAFO <- spTransform(NAFO, CRS('+proj=utm +zone=21 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
NAFO_3LON <- NAFO[NAFO@data$ZONE %in% c('3L','3O','3N'),] #NAFO boundaries for divisions containing both spring and fall surveys

Fall <- Fall[NAFO_3LON,] #clip fall and spring survey spatial dataframes to within
Spring <- Spring[NAFO_3LON,] #NAFO divisions containing surveys from both seasons
combined_seasons <- bind(Fall, Spring) # combine fall and spring surveys into one dataframe

combined_seasons@data <- replace_na(data = combined_seasons@data, #replace NA values (species not caught in both seasons) with 0 
                                    list(Cryptopsaras.couesii = 0, 
                                         Petromyzon.marinus = 0, 
                                         Ulvaria.subbifurcata = 0))

#Test for multivariate differences using PERMANOVA and compare with ANOSIM

library('vegan')
library('simba')
perm.season <- adonis(combined_seasons@data[,9:93]~season, data = combined_seasons@data)
perm.season$aov.tab #significant effect of season (p = 0.001)
saveRDS(perm.season,'Output/permanova_seasons_nfld.RDS') #save PERMANOVA results
#compare PERMANOVA based on Bray-Curtis similarity with simpson dissimilarity
sitexspecies <- combined_seasons@data[,9:93] #extract site x species matrix
sitexspecies[sitexspecies>0] <- 1 #convert weights to P-A
Nfld_simp_diss <- sim(sitexspecies, method = 'simpson')
perm.season.simp <- adonis(Nfld_simp_diss ~ season, data = combined_seasons@data)
perm.season.simp$aov.tab #significant effect of season (p = 0.001) but small effect R2 = 0.013
Nfld_bc_diss <- vegdist(combined_seasons@data[,9:93], method = 'bray') #compute Bray-Curtis dissimilarity matrix for nlfd data
anosim_groups <- factor(combined_seasons@data$season) #create factor vector for seasons
Nfld.ano <- anosim(Nfld_bc_diss, anosim_groups) #ANOSIM comparing seasons
summary(Nfld.ano) #significant (p = 0.001) of season
hist(Nfld.ano$perm) #distribution of permuted R values is shift to left. Could indicate difference in dispersion rather than location
mds_seasons <- labdsv::nmds(Nfld_bc_diss, k = 2) #visualize with nmds
saveRDS(mds_seasons2, 'Output/Nfld_mds_k3.rds')
mds_seasons2b <- metaMDS(Nfld_bc_diss, k = 3, trymax = 100, sratmax = 0.999999) #visualize with nmds
mds_seasons2b <- metaMDS(Nfld_bc_diss, k = 3, trymax = 100, previous.best = mds_seasons2, sratmax = 0.999999)
points(mds_seasons[['points']][1:1848,], col = 'red') #seems to be a lot of overlap
points(mds_seasons[['points']][1848:3883,], col = 'blue') #small difference in location, but large samples?

#try PERMANOVA with smaller sample (n = 50) for each season

Spring.id <- sample(Spring@data$ID,50) #ID of spring subset
Fall.id <- sample(Fall@data$ID,50) #ID of fall subset
Spring.sub <- Spring[Spring@data$ID %in% Spring.id,] # subset of spring survey sets
Fall.sub <- Fall[Fall@data$ID %in% Fall.id,] #subset of fall survey sets
combined_seasons.sub <- bind(Fall.sub, Spring.sub) # combine fall and spring surveys into one dataframe
combined_seasons.sub@data <- replace_na(data = combined_seasons.sub@data, #replace NA values (species not caught in both seasons) with 0 
                                        list(Cryptopsaras.couesii = 0, 
                                             Petromyzon.marinus = 0, 
                                             Ulvaria.subbifurcata = 0))
sitexspecies <- combined_seasons.sub@data[,9:93] #extract site x species matrix
sitexspecies[sitexspecies>0] <- 1 #convert weights to P-A
Nfld_simp_diss <- sim(sitexspecies, method = 'simpson')
perm.season.simp <- adonis(Nfld_simp_diss ~ season, data = combined_seasons.sub@data)
perm.season.simp$aov.tab #(p = 0.435), differences in centroid of spring and fall groups is likely small