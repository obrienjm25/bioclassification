# Load packages

library(purrr)
library(data.table)
library(fread)
library(dplyr)
library(tidyr)
library(rlang)

# List of csv files with final SiteXSpecies matrices

filelist <- list.files(path = 'Data/', pattern = 'ClusterData4km', full.names = T)[1:4]

# Get species list for each region from SiteXSpecies matrices
species.ls <- map(filelist, ~fread(.)) %>% 
  map(., ~dplyr::select(., -V1)) %>%  
  map(., colnames) %>% 
  set_names(., c('MAR','SGSL','NL','NGSL')) %>% 
  map(., sort)

# lists of regex expressions and partially filled functions to edit taxa names by region
patterns <- c('\\.',rep('\\.(?=.+)',3))
sub.fn <- map2(patterns, species.ls, ~partial(.f = gsub, 
                                              pattern = !!.x, 
                                              perl = T, 
                                              x = !!.y, 
                                              replacement = ' '))

# Edit taxa names in each region
species.ls2 <- invoke_map(sub.fn) %>%  
  map(., ~gsub('Boltenia ovifera ', 'Boltenia ovifera', .x)) %>% 
  map(., ~gsub('Cirripedia s c.', 'Cirripedia s.c.', .x)) %>% 
  map(., ~gsub('Astarte sp .', 'Astarte sp. ', .x)) %>% 
  map(., ~gsub('Protobranchia  heterodonta', 'Protobranchia, heterodonta', .x))
  
# Load regional data frames with complete taxonomic classification and common names
load('Data/MaritimesData.RData')
Maritimes$species[grep(pattern = 'Boltenia', x = Maritimes$species)] <- 'Boltenia ovifera' 
load('Data/NfldData.RData')
QC <- readRDS('Data/QC_invertsAdded.rds')
Gulf <- readRDS('Data/Gulf_invertsAdded.rds')

# Filter regional data frames to include 1 record only for each unique taxon in final list
# And join into a common data frame
taxa.info.4region <- list(MAR = Maritimes, SGSL = Gulf, NL = Nfld, NGSL = QC) %>% 
  map2(., species.ls2, ~filter(.x, species %in% .y)) %>%
  map(., ~dplyr::select(., phylum, species, common_name, region)) %>% 
  map(., distinct) %>%
  rbind_list() %>% 
  mutate(., region = case_when(
    region == 'MARITIME' ~ 'MAR',
    region == 'GULF' ~ 'SGSL',
    region == 'NEWFOUNDLAND' ~ 'NL',
    region == 'QUEBEC' ~ 'NGSL'
  )) %>% 
  arrange(.,species, region)

# Categorize each taxa with new factor corresponding to list of all regions in which it's found
region.fct <- mutate(taxa.info.4region, 
                   NGSL = as.numeric(region == 'NGSL'),
                   SGSL = as.numeric(region == 'SGSL'),
                   NL = as.numeric(region == 'NL'),
                   MAR = as.numeric(region == 'MAR')) %>% 
  dplyr::select(-common_name) %>%
  group_by(species) %>% 
  summarise_at(vars(NGSL, SGSL, NL, MAR), sum) %>% 
  mutate(.,
         NGSL = if_else(NGSL == 1, 'NGSL',''),
         SGSL = if_else(SGSL == 1, 'SGSL',''),
         NL = if_else(NL == 1, 'NL', ''),
         MAR = if_else(MAR == 1, 'MAR','')) %>% 
  unite(., Region, NGSL:MAR, sep = ' ') %>%
  mutate(., Region = gsub(' +$', '' ,Region),
         Region = gsub('^ +', '' ,Region),
         Region = gsub(' +', ', ' ,Region))

# Join new factor with full data frame 
# Remove duplicate rows where one taxa is called by multiple names depending on region
Taxa.ls.full <- left_join(taxa.info.4region, region.fct, by = c('species')) %>% 
  rename(Taxa = species) %>% 
  dplyr::select(-region) %>%
  mutate(., common_name = toupper(common_name)) %>% 
  distinct(., Taxa, .keep_all = T) %>%  
  arrange(phylum, Taxa) 

# Export full taxa list

write.csv(Taxa.ls.full, 'Data/FinalTaxaList.csv', row.names = F)
  