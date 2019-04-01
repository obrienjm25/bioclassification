############################################################################################
#############################SPERA Data Cleaning and Prep###################################
############################################################################################

#load required packages

library(lubridate)
library(dplyr)
library(rfishbase)
library(taxize)
library(tidyr)
library(reshape2)
library(reshape)
library(purrr)
library(data.table)
library(worrms)

####GULF####

gulf<-read.csv("Data/sGSL-Obrien-February-2019(RV data).csv",stringsAsFactors = F) #read in survey data
spec.code <- read.csv('Data/sGULF_RV_TAXA_List_2007to2017.csv', stringsAsFactors = F)

min(gulf$year) #2007
max(gulf$year) #2017
table(gulf$duration) #all tows < 40min and > 20min (valid sets)

#reshape dataframe to long format & join with species codes
gulf_wgt <- dplyr::select(gulf, -contains('.number.caught')) %>% #make separate df for biomass
  `colnames<-`(., gsub('.weight.caught','',colnames(.))) #standardize colnames

gulf_tot <- dplyr::select(gulf, -contains('.weight.caught')) %>% #make separate df for abundance
  `colnames<-`(., gsub('.number.caught','',colnames(.))) #standardize colnames

gulf2 <- bind_rows(gulf_wgt, gulf_tot,.id = 'Metric') %>% #rowbind biomass and abundance dataframes
  dplyr::select(.,-X) #remove indexing column
gulf2$Metric <- recode(gulf2$Metric, '1' = 'biomass', '2' = 'catch') #recode factor indicating metric

names(gulf2)[23:367] <- spec.code$english #rename columns to match species code index

gulf2 <- gather(gulf2, 'common.name','value',23:367) %>%  #reshape to long format (entry for each species per set)
  spread(.,Metric,value) %>%   #spread multiple trawl measures (biomass, abundance) across cols
  inner_join(.,spec.code, by = c('common.name' = 'english')) %>%  #join species code table with survey data
  dplyr::select(., -latin,-common.name) #remove old taxa and common names

gulf2$region<-"GULF" #add region identifier
head(gulf2)
names(gulf2)
length(unique(gulf2$species)) #222 Unique taxa

#change species and common names to lower case
gulf2$species <- paste(substring(gulf2$species,1,1),tolower(substring(gulf2$species,2)),sep='')
gulf2$common_name <- paste(substring(gulf2$common_name,1,1),tolower(substring(gulf2$common_name,2)),sep='')

#remove rows that had NA or 0 listed for abundance
gulf2$abundance <- rowSums(gulf2[22:23], na.rm = T) #590640 observations in original dataframe
gulf2 <- gulf2[is.na(gulf2$abundance) == FALSE,] #590640 (No NAs)
gulf2 <- gulf2[gulf2$abundance>0,] #52022 observations

#Remove taxa that aren't fish or benthic invertebrates

rm.spec <- as.integer(c(9000,9001,9200,9300,9302,9303,9304,9360,9400,9500,9620,9630,9991,9992,9993,9994,9995,9996,
             9997,9998,9999,1099,1100,1200,1221,1224,1297,1511,1600,3199,3999,4348,6000,7000,8000,8100,8113,8321,
             8200))
gulf2 <- filter(gulf2, !(spec_code %in% rm.spec)) #45348 observations

##create new column for species that are going to be taxized
gulf2$taxa <- tolower(gulf2$species)

#remove taxonomic resolution abbreviations and correct problematic taxa
abb <- c(" p\\."," s\\.p\\."," c\\."," s\\.c\\.", " o\\.", " f\\.", " s\\.f\\.", " sp\\.", " sp$", " spp\\.")
gulf2$taxa <- gsub(paste(abb,collapse = '|'),'',gulf2$taxa)
gulf2$taxa <- gsub("protobranchia, heterodonta","protobranchia",gulf2$taxa)

#Use Classify function to expand species names into taxonomic groups 
source("Code/ClassifyFunction.R")
SpeciesList <- unique(gulf2$taxa)
SpeciesList[131] <- 'mya truncata'
SpeciesList[144] <- 'bathypolypus'
outlist <- lapply(SpeciesList,FUN=Classify)
taxInfo <- do.call("rbind", outlist)
taxInfo[131,'newsciname'] <- 'mya'
taxInfo[144,'newsciname'] <- 'bathypolypus bairdii'
write.csv(taxInfo, 'Data/Gulf_taxonomy')
#merging column data from taxinfo to gulf
gulf2 <- inner_join(gulf2, dplyr::select(taxInfo,-species), by =c('taxa' = 'newsciname'))
write.csv(gulf2,'Data/GulfData_FishInverts.csv')

#merge data with functional traits
functional<-read.csv("Data/functionaldatabase_141217.csv",stringsAsFactors = F)
dupl_var <- intersect(names(functional),names(gulf2))[-2]
gulf_merged <- left_join(gulf2,dplyr::select(functional, -dupl_var),by="species")

#get functional data on verticial and horizontal position for groups not in database
functional.nc <- unique(gulf_merged$taxa[is.na(gulf_merged$vertical.position)]) #taxa with missing functional data
functional.nc2 <- paste(toupper(substring(functional.nc,1,1)),substring(functional.nc,2),sep='') #same vector as above, but capitalized first letter
wm.details <- wm_records_names(name = functional.nc, fuzzy  = F) %>% #get data on environment (marine, brackish, fresh)
  map(~ list(aphID = as.character(.$AphiaID), scientificname = as.character(.$scientificname), phylum = as.character(.$phylum),
             marine = as.integer(.$isMarine), brackish = as.integer(.$isBrackish), fw = as.integer(.$isFreshwater))) %>% 
  data.table::rbindlist(.) %>% 
  filter(., !(aphID %in% c('771447','602300')))
wm.details$taxa <- functional.nc
wm.details$horizontal.position <- factor(rowSums(dplyr::select(wm.details,marine:fw), na.rm = T)) %>% 
  recode(., `1` = 'marine', `2` = 'brackish; marine' , `3` = 'freshwater; brackish; marine')
wm.details[75,8] <- 'freshwater; brackish; marine'

wm.attr <- wm_attr_data_(id = as.integer(wm.details$aphID), include_inherited = T) %>% #get data on vertical position from WoRMS
  filter(., measurementTypeID == 4) %>%
  group_by(id) %>% 
  summarise_at(.,'measurementValue', first) %>% 
  dplyr::rename(vertical.position = measurementValue)

fish.attr <- species(validate_names(functional.nc2), fields = c('Species',species_fields$habitat)) %>% #get data on vertical position from fishbase #get data on vertical position from fishbase
  dplyr::rename(vertical.position = DemersPelag) 

#join functional data from different sources

for (i in 1:nrow(fish.attr)){
  fish.attr$Species[i] <- wm.details$taxa[agrepl(fish.attr$Species[i],wm.details$taxa,ignore.case = T, max.distance = 0.2)]
}

functional_update <- left_join(wm.details,wm.attr, by = c('aphID' = 'id')) %>% 
  left_join(., dplyr::select(fish.attr,c('Species','vertical.position')), by = c('taxa' = 'Species')) %>% 
  mutate(vertical.position = coalesce(vertical.position.x,vertical.position.y)) %>% 
  dplyr::select(-vertical.position.x, - vertical.position.y) %>% 
  mutate_if(is.factor,as.character)

#fill in NAs by looking up similar taxa in functional database, fix errors, & align spellings with functional database

functional_update$vertical.position[which(is.na(functional_update$vertical.position))] <- c('bathypelagic','pelagic,oceanic',
                                                                                            'demersal','demersal','demersal','benthic',
                                                                                            'benthic','benthic','benthic','benthopelagic',
                                                                                            'benthic','benthic','benthic','benthic',
                                                                                            'benthic','benthic','benthic')
functional_update$vertical.position <- gsub('os','ic', functional_update$vertical.position)
functional_update$vertical.position <- gsub('-o',',o', functional_update$vertical.position)
functional_update$vertical.position <- gsub('nekton','benthic', functional_update$vertical.position)
functional_update$vertical.position <- gsub('plankton > zooplankton','pelagic,oceanic', functional_update$vertical.position)

functional_update[80,9] <- 'pelagic,oceanic'

#join updated functional traits with gulf database
gulf_merged <- left_join(gulf_merged, dplyr::select(functional_update,taxa:vertical.position), by = 'taxa') %>% 
  mutate(vertical.position = coalesce(vertical.position.x,vertical.position.y),
         horizontal.position = coalesce(horizontal.position.x,horizontal.position.y)) %>% 
  dplyr::select(-vertical.position.x, - vertical.position.y,-horizontal.position.x,-horizontal.position.y)

#write updated gulf dataframe to RDS object
saveRDS(gulf_merged, 'Data/Gulf_invertsAdded.rds')

####QUEBEC####

#5424 sets from 1995 to 2018
QC.set <- read.csv("Data/Quebec_RVdata/NGSL_RV_SetMetadata_eng.csv",stringsAsFactors = F) #read in survey set metadata
QC.catch <- read.csv("Data/Quebec_RVdata/NGSL_RV_Catch_eng.csv",stringsAsFactors = F) #read in catch data
spec.code <- read.csv('Data/Quebec_RVdata/NGSL_RV_taxalist.csv', stringsAsFactors = F)
spec.code$spec_code <- as.integer(spec.code$spec_code)
spec.code$AphID <- as.integer(spec.code$AphID)

#reformat lon/lat from DD MM. to decimal degrees

QC.set$DD_lat <- as.numeric(substr(QC.set$lat,1,2)) + (as.numeric(substr(QC.set$lat, 3,7))/60)
QC.set$DD_lon <- as.numeric(substr(QC.set$lon,1,2)) + (as.numeric(substr(QC.set$lon, 3,7))/60)

#Join set, catch, and species tables
QC <- inner_join(QC.catch,QC.set) #joining by vessel, survey, fishing number, and station
QC <- inner_join(QC, spec.code) #join by species code
#129027 observations

#Clean dataset

QC <- filter(QC, year >= 2007 & year <=2017) %>% #restrict year range
  filter(., gear_type == 6) %>% #restrict to sets with Campelan trawl
  filter(., tow_type == 1) %>% #restrict to fishing sets
  filter(., trawl_operation %in% c(1,2)) #restrict to valid sets/normal trawl operation
#74369 observations
length(unique(QC$spec_code)) #518 unique taxa
length(unique(paste0(QC$vessel,QC$survey,QC$nbpc,QC$stn,QC$year,QC$month,QC$day))) #1923 sets
QC_taxa2007_2017 <- unique(QC$spec_code) #create list of unique taxa within our dataset

#Create taxalist particular to our dataset for further cleaning/aggregating
# QC_taxalist_2007_2017 <- spec.code[spec.code$spec_code %in% QC_taxa2007_2017,]
# write.csv(QC_taxalist_2007_2017,'Data/Quebec_RVdata/NGSL_RV_taxalist_2007_2017.csv')

min(QC$year) #2007
max(QC$year) #2017
range(QC$dist_spd) #Towed distances between 0.5-1 nautical miles (within valid range)

#Rejoin filtered dataframe with updataed taxa list for 2007-2017

QC_taxalist_2007_2017 <- read.csv('Data/Quebec_RVdata/NGSL_RV_taxalist_2007_2017.csv', stringsAsFactors = F)

QC2 <-dplyr::select(QC, -c(AphID, species,common_name)) %>%  #select columns except placeholders for revised names
  inner_join(.,QC_taxalist_2007_2017) %>% #join species code table with survey data
  dplyr::select(., -latin,-english) %>% #remove old taxa and common names
  filter(., exclude != 'rm')#remove pelagic, non-metazoan, and uncertain taxa
  
QC2$region<-"QUEBEC" #add region identifier
head(QC2)
names(QC2)
length(unique(QC2$species)) #249 Unique taxa

#remove rows that had NA or 0 listed for abundance
QC2$abundance <- rowSums(cbind(QC2$biomass_corr,QC2$catch_corr), na.rm = T) #68836 observations in original dataframe
QC2 <- QC2[is.na(QC2$abundance) == FALSE,] #68849 (No NAs)
QC2 <- QC2[QC2$abundance>0,] #68834 observations

##create new column for species that are going to be taxized
QC2$taxa <- tolower(QC2$species)

#remove taxonomic resolution abbreviations
abb <- c(" p\\."," s\\.p\\."," c\\."," s\\.c\\.", " o\\.", " f\\.", " s\\.f\\.", " sp\\.", " sp$", " spp\\.")
QC2$taxa <- gsub(paste(abb,collapse = '|'),'',QC2$taxa)

#Use Classify function to expand species names into taxonomic groups 
source("Code/ClassifyFunction.R")
SpeciesList <- unique(QC2$taxa)
SpeciesList[46] <- 'bathypolypus'
SpeciesList[70] <- 'paraliparis copei'
SpeciesList[113] <- 'arctinula greenlandica'
SpeciesList[144] <- 'stichaeus punctatus'
SpeciesList[207] <- 'notoscopelus'
SpeciesList[223] <- 'epizoanthus'
SpeciesList[234] <- 'liponema'
SpeciesList[243] <- 'scabrotrophon'
outlist <- lapply(SpeciesList,FUN=Classify)
taxInfo <- do.call("rbind", outlist)
taxInfo[46,'newsciname'] <- 'bathypolypus bairdii'
taxInfo[70,'newsciname'] <- 'paraliparis copei copei'
taxInfo[113,'newsciname'] <- 'similipecten greenlandicus'
taxInfo[144,'newsciname'] <- 'stichaeus punctatus punctatus'
taxInfo[207,'newsciname'] <- 'notoscopelus kroyeri'
taxInfo[223,'newsciname'] <- 'epizoanthus erdmanni'
taxInfo[234,'newsciname'] <- 'liponema multicorne'
taxInfo[243,'newsciname'] <- 'scabrotrophon fabricii'
write.csv(taxInfo, 'Data/QC_taxonomy')
#merging column data from taxinfo to QC
QC2 <- inner_join(QC2, dplyr::select(taxInfo,-species), by =c('taxa' = 'newsciname'))
write.csv(QC2,'Data/QCData_FishInverts.csv')

#merge data with functional traits
functional<-read.csv("Data/functionaldatabase_141217.csv",stringsAsFactors = F)
dupl_var <- intersect(names(functional),names(QC2))[-2]
QC_merged <- left_join(QC2,dplyr::select(functional, -dupl_var),by="species")

#get functional data on verticial and horizontal position for groups not in database
functional.nc <- unique(QC_merged$taxa[is.na(QC_merged$vertical.position)]) #141 taxa with missing functional data
functional.nc2 <- paste(toupper(substring(functional.nc,1,1)),substring(functional.nc,2),sep='') #same vector as above, but capitalized first letter
functional.nc[[44]] <- "Flabellum (Ulocyathus) alabastrum"
wm.details <- wm_records_names(name = functional.nc, fuzzy  = F) %>% #get data on environment (marine, brackish, fresh)
  map(~ list(aphID = as.character(.$AphiaID), scientificname = as.character(.$scientificname), phylum = as.character(.$phylum),
             marine = as.integer(.$isMarine), brackish = as.integer(.$isBrackish), fw = as.integer(.$isFreshwater))) %>%     
  data.table::rbindlist(.) %>% 
  filter(., !(aphID %in% c('771447','602300','123506','146121','990530')))
wm.details$taxa <- unique(QC_merged$taxa[is.na(QC_merged$vertical.position)])
wm.details$horizontal.position <- factor(rowSums(dplyr::select(wm.details,marine:fw), na.rm = T)) %>% 
  recode(., `1` = 'marine', `2` = 'brackish; marine' , `3` = 'freshwater; brackish; marine')
wm.details[97,8] <- 'freshwater; brackish; marine'
wm.details[120,8] <- 'marine'

wm.attr <- wm_attr_data_(id = as.integer(wm.details$aphID), include_inherited = T) %>% #get data on vertical position from WoRMS
  filter(., measurementTypeID == 4) %>%
  group_by(id) %>% 
  summarise_at(.,'measurementValue', first) %>% 
  dplyr::rename(vertical.position = measurementValue)

fish.attr <- species(validate_names(functional.nc2), fields = c('Species',species_fields$habitat)) %>% #get data on vertical position from fishbase #get data on vertical position from fishbase
  dplyr::rename(vertical.position = DemersPelag) 

#join functional data from different sources

for (i in 1:nrow(fish.attr)){
  fish.attr$Species[i] <- wm.details$taxa[grep(fish.attr$Species[i],wm.details$taxa,ignore.case = T)]
}

functional_update <- left_join(wm.details,wm.attr, by = c('aphID' = 'id')) %>% 
  left_join(., dplyr::select(fish.attr,c('Species','vertical.position')), by = c('taxa' = 'Species')) %>% 
  mutate(vertical.position = coalesce(vertical.position.x,vertical.position.y)) %>% 
  dplyr::select(-vertical.position.x, - vertical.position.y) %>% 
  mutate_if(is.factor,as.character)

#fill in NAs by looking up similar taxa in functional database, fix errors, & align spellings with functional database

functional_update$vertical.position[which(is.na(functional_update$vertical.position))] <- c('benthic','benthic','bathypelagic',
                                                                                            'demersal','pelagic,oceanic','benthic',
                                                                                            'demersal','dermersal','benthic','demersal',
                                                                                            'benthic','benthic','benthic','demersal',
                                                                                            'benthic','benthic','demersal','benthic',
                                                                                            'benthopelagic','bathydemersal','benthic',
                                                                                            'pelagic,oceanic','demersal')

functional_update$vertical.position <- gsub('benthos > macrobenthos','benthic', functional_update$vertical.position)
functional_update$vertical.position <- gsub('os','ic', functional_update$vertical.position)
functional_update$vertical.position <- gsub('-o',',o', functional_update$vertical.position)
functional_update$vertical.position <- gsub('nekton','pelagic,neritic', functional_update$vertical.position)
functional_update$vertical.position <- gsub('plankton > zooplankton','pelagic,neritic', functional_update$vertical.position)

functional_update[functional_update$aphID == '107567',9] <- 'benthic'
functional_update[functional_update$aphID == '158357',9] <- 'benthic'
functional_update[functional_update$aphID == '107566',9] <- 'demersal'
functional_update[functional_update$aphID == '158359',9] <- 'demersal'
functional_update[functional_update$aphID == '107532',9] <- 'benthic'
functional_update[functional_update$aphID == '106986',9] <- 'benthic'
functional_update[functional_update$aphID == '158362',9] <- 'benthic'

#join updated functional traits with QC database
QC_merged <- left_join(QC_merged, dplyr::select(functional_update,taxa:vertical.position), by = 'taxa') %>% 
  mutate(vertical.position = coalesce(vertical.position.x,vertical.position.y),
         horizontal.position = coalesce(horizontal.position.x,horizontal.position.y)) %>% 
  dplyr::select(-vertical.position.x, - vertical.position.y,-horizontal.position.x,-horizontal.position.y)

#write updated QC dataframe to RDS object
saveRDS(QC_merged, 'Data/QC_invertsAdded.rds')
