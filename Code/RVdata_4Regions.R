#### Load packages

library(purrr)
library(data.table)
library(fread)
library(dplyr)
library(tidyr)
library(rlang)

#### Load regional RV datasets
load('Data/MaritimesData.RData')
load('Data/NfldData.RData')
QC <- readRDS('Data/QC_invertsAdded.rds')
Gulf <- readRDS('Data/Gulf_invertsAdded.rds')

#### Tidy datasets

# select columns with: year, month, day, set, depth, latitude, 
# longitude, species, common name, tow distance, total weight, total number, region

RV_vars <- c('year','month','day','set','depth','latitude','longitude','species',
             'common_name', 'phylum', 'subphylum','distance','totwgt','totno',
             'region','vertical.position')
# Maritimes
Maritimes$species[grep(pattern = 'Boltenia', x = Maritimes$species)] <- 'Boltenia ovifera' 
Maritimes <- dplyr::select(Maritimes, one_of(RV_vars)) %>% 
  mutate(., vertical.position = gsub('pelagic-oceanic','pelagic,oceanic', vertical.position)) %>% 
  mutate(., vertical.position = gsub('pelagic-neritic','pelagic,neritic', vertical.position)) 


# Newfoundland & Labrador
Nfld <- dplyr::select(Nfld, one_of(RV_vars)) %>% 
  mutate(., vertical.position = gsub('pelagic-oceanic','pelagic,oceanic', vertical.position)) %>% 
  mutate(., vertical.position = gsub('pelagic-neritic','pelagic,neritic', vertical.position)) 

# Quebec
QC <- QC %>%
  rename(., set = stn, 
         latitude = DD_lat, 
         longitude = DD_lon,
         distance = dist_spd,
         totno = catch_corr,
         totwgt = biomass_corr) %>% 
  mutate(., longitude = -(longitude)) %>% 
  dplyr::select(., one_of(RV_vars)) %>% 
  mutate(., vertical.position = gsub('pelagic-oceanic','pelagic,oceanic', vertical.position)) %>% 
  mutate(., vertical.position = gsub('pelagic-neritic','pelagic,neritic', vertical.position)) %>% 
  mutate(., vertical.position = gsub('dermersal','demersal',vertical.position))

# Gulf
Gulf <- Gulf %>% 
  rename(., set = set.number,
         totno = catch,
         totwgt = biomass) %>% 
  dplyr::select(., one_of(RV_vars)) %>% 
  mutate(., vertical.position = gsub('pelagic-oceanic','pelagic,oceanic', vertical.position)) %>% 
  mutate(., vertical.position = gsub('pelagic-neritic','pelagic,neritic', vertical.position)) 

# Join RV datasets

RV.data <- bind_rows(Maritimes, Nfld, QC, Gulf) %>% 
  arrange(., region, species, year, month, day, set)

# Export file

saveRDS(RV.data, 'R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVdata_4regions.rds')
