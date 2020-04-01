####################################################
############Functional trait analysis###############
####################################################

####Outline####

# 1) Create regional functional trait databases
# 2) Calculate various functional diversity indices (FDisp, Frich) observed in contemporary trawl sets (2007-present)
# 3) Examine how functional diversity varies between major assemblages and as a function of species richness
# 4) Predict functional diversity over study area (interpolation or predictive modelling)
# 5) Identify hotspots with quantile breaks
# 6) Examine temporal trends in functional diversity for subset of species consistently sampled since 1970s
# 7) Examine functional responses to environmental and anthropogenic gradients (hierarchical bayesian models)

##### 1. Create regional functional traits databases ####

#### Load packages ####

library(tidyr)
library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(rfishbase)
library(ggfortify)
library(factoextra)
library(FD)
library(adiv)
library(sf)
library(tmap)
library(RColorBrewer)
library(raster)
library(ranger)

#### Load existing functional database ####

traits <- fread('Data/functionaldatabase_141217.csv') %>% #functional database
  mutate(., subphylum = if_else(species == 'Hippoglossina oblonga', 
                                'Vertebrata', 
                                subphylum)) %>% # add subphylum for 4-spot flounder 
  filter(., subphylum == 'Vertebrata') %>% #fish only
  # filter(., !(grepl('^pelagic', vertical.position))) %>% # remove pelagics
  mutate(., fb.names = unique(validate_names(species))) %>% #add column for validated fishbase names
  dplyr::select(., common_name:genus, resilience:fertilization, #select variables of interest
                vertical.position:fecundity.max,no.of.food3:fb.names) %>% 
  rename(., repro.mode = reproduction_Mode, #rename some vars
         shallow.m = shallow) %>% 
  `names<-`(., gsub('_','.',names(.))) #format var names

#### Download other fishbase tables ####

# fecund <- fecundity(traits$fb.names) %>% 
#   dplyr::select(Species,StockCode,Locality,FecundityMin,FecundityMax)

# Diet, Depth, and longevity
ModEstimates <- estimate(traits$fb.names) %>%
  dplyr::select(Species,DepthMin:DepthMax,AgeMax) %>% 
  # dplyr::select(Species:MaxLengthTL,Troph:seTroph,K,DepthMin:DepthMax,AgeMax)
  rename(shallow.m = DepthMin,
         deep.m = DepthMax,
         age = AgeMax)

# Larval characteristics
larvae <- larvae(traits$fb.names) %>% 
  dplyr::select(Species,StockCode,PlaceofDevelopment) %>% 
  # dplyr::select(Species,StockCode,LarvalArea,PlaceofDevelopment,LarvalDurationMod) %>% 
  distinct(., Species, .keep_all = T)
  

# Morphology
morph <- morphology(traits$fb.names) %>% 
  # dplyr::select(Species,StockCode,Females,Males,BodyShapeI:BodyShapeII,PosofMouth) %>% 
  dplyr::select(Species,PosofMouth) %>%
  distinct(., Species, .keep_all = T)

# Reproductive strategies
repro <- reproduction(traits$fb.names) %>% 
  dplyr::select(Species,StockCode,ReproMode,Fertilization) %>% 
  rename(repro.mode = ReproMode, 
         fertilization = Fertilization)

# Other species characterstics: Body Shape, Migration, Depth range, Longevity
species <- species(traits$fb.names) %>% 
  dplyr::select(Species,DepthRangeShallow:DepthRangeDeep,LongevityWild) %>% 
  # dplyr::select(Species,BodyShapeI,AnaCat,DepthRangeShallow:DepthRangeDeep,LongevityWild,
  #               Length,LTypeMaxM,Weight)
  rename(shallow.m = DepthRangeShallow,
         deep.m = DepthRangeDeep,
         age = LongevityWild)

# stocks <- stocks(traits$fb.names)

# Swimming characteristics: Broad adult swim type (e.g. movement of body or caudal fin) 
# and specific swim mode (e.g. anguilliform)
swim <- swimming(traits$fb.names) %>% 
  dplyr::select(Species,AdultType,AdultMode)

#### Merge fb tables and existing trait table ####

traits2 <- left_join(traits, ModEstimates, by = c('fb.names' = 'Species')) %>% #left-join by fishbase name
  left_join(., larvae, by = c('fb.names' = 'Species')) %>% #                      "           "
  left_join(., morph, by = c('fb.names' = 'Species')) %>% #                       "           "
  left_join(., repro, by = c('fb.names' = 'Species')) %>% #                       "           "
  left_join(., species, by = c('fb.names' = 'Species')) %>% #                     "           "
  left_join(., swim, by = c('fb.names' = 'Species')) %>% #                        "           "
  mutate_at(., vars(shallow.m.x, deep.m.x, age.x), as.double) %>% #coalesce duplicate variables to fill some NAs
  mutate(., shallow.m = coalesce(shallow.m.x,shallow.m.y,shallow.m),
         deep.m = coalesce(deep.m.x, deep.m.y, deep.m),
         age = coalesce(age, age.x, age.y),
         repro.mode = coalesce(repro.mode.x, repro.mode.y),
         fertilization = coalesce(fertilization.x, fertilization.y)) %>% 
  dplyr::select(-ends_with('.x')) %>% # remove duplicate variables
  dplyr::select(-ends_with('.y')) %>% # remove duplicate variables
  dplyr::select(-c(family,genus, resilience.category,salinity.tolerance:depth.category, # select subset of variables
                   trophic.level))

#### PCA of diet items ####

# Get diet data from functional trait table
diet.data <- dplyr::select(traits2, starts_with('f3')) %>% # extract cols with proportions of prey items in diet
  dplyr::select(-f3.pinnipeds) %>% # remove column for pinnepeds
  filter(., complete.cases(.)) %>% 
  `rownames<-`(., traits2$fb.names[!is.na(traits2$f3.amphipods)]) # species names as row names

# Run PCA on column centred data  
diet.pca <- prcomp(diet.data, center  = T, scale = F)

# Examine eigenvalues (first 9 principal dimensions explain ~ 80% of variance)

get_eigenvalue(diet.pca) # eigenvalues
fviz_eig(diet.pca, ncp = 20) #scree plot

# Examine quality and contribution of variables on PCA1-PCA9(i.e. diet items)
# bony fish, other benthic crustaceans, planktonic copepods, polychaetes, amphipods, other benthic invertebrates,
# other planktonic inverts, other planktonic crustaceans, euphausids, and bivalves are well represented and contribute
# most to variance of first 9 PCA axes

fviz_cos2(diet.pca, choice = 'var', axes = 1:9) # quality
fviz_contrib(diet.pca, choice = 'var', axes = 1:9) # contributions
map(1:9, ~fviz_contrib(diet.pca, choice = 'var', axes = .)) # look at diet item contributions axis by axis

# PCA1 bonyfish, other benthic crustaceans, polychaetes
# PCA2 other benthic crustaceans, planktonic copepods, amphipods, bonyfish
# PCA3 polychaetes, planktonic copepods, amphipods, cladocerans, other benthic crustaceans
# PCA4 other benthic invertebrates, planktonic copepods, polychaetes, bonyfish
# PCA5 polychaetes, amphipods, euphasiids, planktonic copepods, other benthic invertebrates, shrmimps/prawns, cladocerans
# PCA6 amphipods, other planktonic crustaceans, euphasiids, planktonic copepods, shrimps/prawns, other benthic invertebrates,
# polychaetes, other benthic crustaceans, bonyfish, bivalves
# PCA7 other planktonic invertebrates, euphausiids, other planktonic copepods
# PCA8 other planktonic crustaceans, euphausiids, bivalves, amphipods
# PCA9 bivalves, euphausiids, other planktonic invertebrates, polychaetes

# Biplot of first 2 PCA axes

guilds = filter(traits2, !is.na(f3.bonyfish)) %>% #categorical feeding guilds
  transmute(feeding.guild = case_when(
    feeding.guild.category == 1 ~ 'piscivore',
    feeding.guild.category == 2 ~ 'benthivore',
    feeding.guild.category == 3 ~ 'planktivore',
    feeding.guild.category == 4 ~ 'zoopiscivore'
  ))

autoplot(diet.pca, data = guilds, col = 'feeding.guild', #biplot with coloured by feeding guild
         loadings = T, loadings.label = T, loadings.colour = 'blue',  loadings.label.size = 4)
  
# Get species PCA scores on first 9 axes and join with functional traits

scores <- as.data.frame(diet.pca$x) %>% # extract scores as dataframe
  mutate(., fb.names = row.names(.)) %>% # add column for fishbase names
  dplyr::select(., fb.names, num_range('PC', 1:9)) # select fb names and scores on PCA1:9

#### Join PCA scores with traits table and write to csv file ####

left_join(traits2, scores, by = 'fb.names') %>% # join by fishbase names
  dplyr::select(-starts_with('f3')) %>% # remove raw diet items
  write.csv(., file = 'Data/functionaldatabase_180220.csv', row.names = F) # write to csv


#### Join our trait table with trait table from Centre for Ocean Life (Beukhof et al. 2019) ####

# Our trait table
traits.dfo <- fread('Data/functionaldatabase_180220.csv') %>% #129 unique taxa
  dplyr::select(-c(length.cm, PosofMouth,age:AdultMode)) %>% # remove unwanted vars
  rename(., RV.name = species) # species are names used in RV

# Centre for Ocean Life table
traits.col <- fread('Data/TraitCollectionFishNAtlanticNEPacificContShelf.csv') %>% # 1395 unique taxa
  # rename synonymous taxa
  mutate(., taxon = gsub('Scomberesox saurus saurus', 'Scomberesox saurus', taxon), 
         taxon = gsub('Osmerus mordax mordax', 'Osmerus mordax', taxon))

# Join tables
traits.join <- left_join(traits.col, traits.dfo, by = c('taxon' = 'fb.names')) %>% # left-join
  filter(., !is.na(common.name)) %>% # exclude taxa not caught in Atlantic RV trawls
  dplyr::select(., common.name:order, everything()) %>% # reorder vars
  arrange(., class, order, family, genus, species) %>% # order by phylogeny
  filter(., LME %in% c(8,9)) # restrict to rows for Large Marine Ecosystems in Canadian Atlantic

# Write to file

write.csv(traits.join, 'Data/TraitCollectionFishNWAtlanticContShelf.csv', row.names = F)

#### RV data prep ####

abb.new <- c('SGSL', 'MAR', 'NL', 'NGSL') # new region abbreviations

# Load RV Data
RV <- readRDS('Data/RVdata_4regions.rds') %>% 
  mutate(., subphylum = if_else(species == 'Hippoglossina oblonga', 
                                'Vertebrata', 
                                subphylum)) %>% # add subphylum for 4-spot flounder 
  filter(., subphylum == 'Vertebrata') %>% #fish only
  # filter(., !(grepl('^pelagic', vertical.position))) # remove pelagics
  # rename regions
  mutate(., region = case_when(
    region == 'GULF' ~ 'SGSL',
    region == 'QUEBEC' ~ 'NGSL',
    region == 'NEWFOUNDLAND' ~ 'NL',
    region == 'MARITIME' ~ 'MAR')) 

RV$totno[RV$species == 'Rouleina maderensis'] <- 2 # imputed abundance for Madeiran smooth-head

RV.taxa <- unique(RV$species) #167 taxa of fish

missing <- setdiff(RV.taxa, traits.join$RV.name) # 38 missing from functional trait tbl

# # Create species list for each region
# sp.list <- split(RV, RV$region) %>% # split RV df by region
#   map(., ~ pull(., species)) %>% # pull species column
#   map(., unique) %>% # only unique species
#   map(., sort) # sort alphabetically

# Create sf points objects with RV survey data (abundance) for each region
# and join with 4-km grid

epsg.codes <- c(26920, 26920, 26921, 26920) # projections for each region

RV.abun <- RV %>% # RV data
  filter(year >= 2007) %>% #filter obs. to post-2007
  # filter(month %in% c(6,7,8)) %>% # filter to summer survey
  mutate(ID = paste(year, month, day, set, sep = '_'), # unique id for each set
         # standardize Maritimes abundance by tow distance
         totabun = case_when(
           region == 'MAR' ~ totno*(1.75/distance), 
           region == 'NL'~ totwgt,
           TRUE ~ totno)) %>% 
  split(., .$region) %>% # split data frame by region
  map(., ~pivot_wider(., 
              id_cols = c(ID, longitude, latitude), # reshape to wide format
              names_from = species, 
              values_from = totabun, 
              values_fn = list(totabun = sum))) %>% 
  map(., ~dplyr::select(., -ID)) %>% 
  map(., ~st_as_sf(., coords = c('longitude','latitude'), # convert to sf object 
           crs = 4269)) %>% # NAD83
  map2(., epsg.codes, ~st_transform(.x, crs = .y)) %>% # transform to UTM projections 
  map(., ~mutate_if(., is.double, list(~ replace_na(., 0)))) # replace NAs with 0's
  # st_join(., grid) %>% 
  # filter(!is.na(layer)) %>% 
  # group_by(layer) %>%
  # summarise_if(is.double, mean)

# Load grids

gridfiles <- list.files('Data/', 
                        pattern = 'GridsWholeCell', 
                        full.names = T) # list of grid files
grid <- list() # empty list

for (i in 1:length(gridfiles)) {
  load(gridfiles[i])
  grid = append(grid, Griddata[[20]])
  rm(Griddata)
}

grid <- map(grid, st_as_sf) %>% # convert grid to sf
  set_names(., abb.new) %>% # set names
  .[sort(names(.))] # sort list by region name alphabetically


# spatial intersection returning logical indicating which cells contain RV data
sel <- map2(grid, RV.abun, ~st_intersects(.x, .y)) %>% 
  map(., ~lengths(.) > 0) 

# Join RV data with grid and summarize data by grid cell
popgrid <- grid %>% # 4-km grid
  map2(., sel, ~filter(.x, .y)) %>% # select cells populated with RV data
  map2(., RV.abun, ~st_join(.x, .y)) %>% # join with RV data
  map(., ~group_by(., layer)) %>% # group by cell
  map(., ~summarise_if(.,is.double, mean)) %>% # average abundance if 2+ obs. in a cell
  map(., ungroup) %>% # ungroup
  map(., ~`row.names<-`(., .$layer)) %>% # row names =  grid cell id
  map(., ~dplyr::select(., -layer)) # remove grid id col

# list of species
sp.list <- popgrid %>% 
  map(., st_drop_geometry) %>% 
  map(., colnames) 

#### Calculate various functional diversity indices from trait data ####

# functional traits to include in analysis
vars.keep <- c('RV.name', 'taxon', 'LME', 'tl', 'body.shape', 'fin.shape', 'AR',
               'offspring.size', 'spawning.type', 'age.maturity', 'fecundity',
               'repro.mode','fertilization', 'egg.location', 'PlaceofDevelopment',
               'growth.coefficient','length.max', 'age.max', 'resilience',
               'vertical.position','shallow.m','deep.m', 'depth.range', 'horizontal.position',
               'feeding.guild.category', 'no.of.food3')

# Import trait table and perform data prep

taxa.out <- map(popgrid, ~ colSums(st_drop_geometry(.)) == 0) %>% 
  map2(., popgrid, ~names(st_drop_geometry(.y))[.x]) # taxa with records outside grid only

trait <- fread('Data/TraitCollectionFishNWAtlanticContShelf_032520.csv') %>%  
  # Add new variables and edit existing ones
  mutate(LME = if_else(LME == 8, 'MAR-SGSL', 'NL-NGSL'),
         fecundity = as.double(fecundity),
         feeding.guild.category = as.factor(feeding.guild.category),
         depth.range = deep.m - shallow.m,
         fertilization = as.numeric(grepl('internal', fertilization)),
         resilience = factor(resilience, 
                            levels = c('Very low','Low','Medium','High'),
                            ordered = T),
         PlaceofDevelopment = case_when(
           grepl('^p.*', PlaceofDevelopment) == T ~ 'planktonic',
           grepl('^in .+', PlaceofDevelopment) == T ~ 'demersal',
           grepl('demersal', PlaceofDevelopment) == T ~ 'demersal',
           grepl('^internal$', PlaceofDevelopment) == T ~ 'internal',
           grepl('^internal.+', PlaceofDevelopment) == T ~ 'internal/planktonic',
           TRUE ~ NA_character_),
         vertical.position = gsub(',' ,'-', vertical.position),
         horizontal.position = gsub(' ', '', horizontal.position)) %>%  
  dplyr::select(one_of(vars.keep)) %>% # select only traits of interest 
  list(.) %>% # create list with trait table repeated 4X
  rep(., 4) %>% 
  # filter obs. to correct large marine ecosystem
  map2(names(popgrid),., ~filter(.y, grepl(pattern = .x, LME))) %>% 
  map2(., sp.list, ~filter(.x, RV.name %in% .y)) %>% # filter only taxa in RV data
  map2(., taxa.out, ~filter(.x, !(RV.name %in% .y))) %>% # filter out taxa found only outside study grid
  map(., ~arrange(., RV.name)) %>% # arrange by taxa
  map(., ~`row.names<-`(., .$RV.name)) %>% # taxa names as row names
  map(., ~dplyr::select(., -c(RV.name, taxon, LME))) # keep only columns taxa columns

# Calculate proportion of missing values in trait tables
map(trait, ~sapply(., function(x) sum(is.na(x))/length(x))*100) # Caudal fin aspect ratio and No. food items with most missing values

# taxa in RV survey missing from trait table
missing<- map2(sp.list, trait, ~setdiff(.x, row.names(.y))) 
map(missing, length) # MAR = 0, NGSL = 8, NL = 0, SGSL = 3, mostly rare taxa identified to family level or large pelagics

# Create matrix objects needed for calculating functional indices

abund <- map(popgrid, st_drop_geometry) %>%
  map2(., trait, ~ dplyr::select(.x, one_of(row.names(.y)))) %>% # to deal with missing taxa for now 
  map(., as.matrix) # abundance matrix

gow <- trait %>% 
  map(., ~gowdis(., ord = 'podani')) # trait matrix based on gower dissimilarity

# Calculate various indices:
# Functional dispersion, richness, eveness, functional group richness, redundancy

fdisper <- map2(gow, abund, ~fdisp(d = .x, a = .y))
frich <- map2(gow, abund, ~dbFD(x = .x, a = .y, calc.FRic = T, corr = 'cailliez'))
# fgr <- map2(gow, abund, ~dbFD(x = .x, a = .y, calc.FRic = F, corr = 'cailliez', 
                              # calc.FGR = T, clust.type = 'ward'))
red <- map2(abund, gow, ~uniqueness(comm = .x, dis = .y))

# Join indices with gridded data

popgrid <- pmap(list(popgrid, fdisper, red, frich), 
                ~bind_cols(..1, 
                           FD = ..2$FDis, 
                           FRed = ..3$red$R, 
                           FRic = ..4$FRic,
                           FEve = ..4$FEve))

# Map distribution of Functional dispersion

# bathymetric layers
bathy <- list.files('Data/Rasters/', 'bathy\\.tif$', full.names = T) %>% 
  map(., raster) %>% 
  set_names(., abb.new) %>% 
  .[sort(names(.))]

tmap_mode('view') # switch to interactive viewer 

map2(bathy, popgrid,
     ~tm_shape(.x)  + 
       tm_raster(palette = rev(brewer.pal(7, 'Blues')),
            breaks = c(-6000,-1000,-500, -200,-100,-50,0,1000)) + 
       tm_shape(.y) + 
       tm_fill(col = 'FD', 
               style = 'quantile',
               n = 5))

# Examine variation in functional diversity between major assemblages

source('Code/SPERA_colour_palettes.R') # load color palettes
palettes <- list(MAR.palette, QC.palette, NL.palette, GULF.palette)
# read in shapefiles with grid cells classified to assemblage types
assemblages <- list.files('Data/Shapefiles/','GridClusterAssignment.*.shp', full.names = T) %>% 
  map(., st_read) %>% 
  map2(., epsg.codes, ~st_transform(.x, crs = .y)) %>% 
  set_names(., abb.new) %>%
  .[sort(names(.))]

# Spatial join grid cells with functional diversity indices 
# with assemblage classification and calculate taxa richness per cell
popgrid2 <- map2(popgrid, assemblages, ~st_join(.x, .y, left = F, join = st_nearest_feature)) %>%
  # convert taxa abundances to p(1)/a(0)
  map2(., sp.list, ~mutate_at(.x, 1:length(.y), list(~ as.numeric(. > 0)))) %>%
  # calculate richness from p/a 
  map2(., sp.list, ~mutate(.x, SpRich = rowSums(st_drop_geometry(.x)[,.y]))) %>% 
  map(., ~mutate(., SpRich = as.integer(SpRich))) %>% 
  map(., ~bind_cols(., long = st_coordinates(st_centroid(.))[,1], lat = st_coordinates(st_centroid(.))[,2]))

# save joined grids to shapefiles
map2(popgrid2, names(popgrid2), 
     ~st_write(.x, dsn = 'Output', 
               layer = paste(.y,'FDiv','sf', sep = '_'),
               driver = 'ESRI Shapefile',
               update = T))

# calculate mean and max functional dispersion in each assemblage type
FD.assemblage <- popgrid2 %>% 
  map(., ~group_by(., name)) %>% 
  map(., ~summarise(., FD_mean = mean(FD, na.rm = T),
            red_mean = mean(FRed, na.rm = T),
            MaxRic = round(max(SpRich, na.rm = T)),
            MeanRic = round(mean(SpRich, na.rm = T)))) %>% 
  map(., ungroup) %>% 
  map2(., palettes, ~left_join(.x, .y))

# Plot functional dispersion as function of assemblage type
map(FD.assemblage, 
    ~ggplot(., aes(x = name, y = FD_mean)) +
      geom_bar(aes(fill = assigned), stat= 'identity') + 
      scale_fill_identity())

# Plot functional richness as a function of Taxa Richness
map(popgrid2, 
    ~ggplot(., aes(x = SpRich, y = FRic)) + 
      geom_point(aes(col = assigned)) +
      scale_colour_identity())

# Plot functional redundancy as a function of taxa richness
map(FD.assemblage, 
    ~ggplot(., aes(x = MaxRic, y = red_mean)) + 
      geom_point(aes(col = assigned)) +
      scale_colour_identity())

#### Predict Functional Dispersion over study domain ####

#import predictor variables retained for modelling

#list of retained variables
retained_vars <- list.files('Output/', 'retain_list\\.txt', full.names = T) %>% 
  map(., ~readLines(file(.))) %>% 
  set_names(., abb.new) %>% 
  .[sort(names(.))]

# get list of raster files
patterns <- c('Mar_','QC_','NL_','Gulf_')
raster.list <- map(patterns,
                   ~list.files(path = 'Data/Rasters/', 
                               pattern = paste0(., '.+\\.tif$'), 
                               full.names = T))
# get study area boundary .shp
RV.bound <- list.files('Data/Shapefiles/', 'Agg\\.shp$', full.names = T) %>% 
  map(., st_read) %>% 
  map2(., epsg.codes, ~st_transform(.x, crs = .y)) %>% 
  set_names(., abb.new) %>% 
  .[sort(names(.))]

# get 5-km land buffer .shp
LandBuffer <- list.files('Data/Shapefiles/', 'LandBuffer_5km\\.shp$', full.names = T) %>%  
  map(., st_read) %>% 
  map2(., epsg.codes, ~st_transform(.x, crs = .y)) %>% 
  set_names(., abb.new) %>% 
  .[sort(names(.))]

env_pred <- map(raster.list, stack) %>%  #stack rasters in list
  map2(., patterns, ~`names<-`(.x, gsub(.y,'',names(.x)))) %>% #remove region abb. from layer names
  map2(., retained_vars, ~.x[[.y]]) %>%  #subset raster stack to include only retained variables
  map2(., RV.bound, ~mask(.x, .y)) %>%  #mask raster cells of predictors outside RV survey boundaries
  map2(., LandBuffer, ~mask(.x, .y, inverse = T)) #also mask cells overlapped by 5km land buffer

#Extract raster values at each populated grid cell

RF.data <- map2(env_pred, popgrid2,
                ~raster::extract(.x, .y, 
                                 factors = T, 
                                 nl = nlayers(.x), 
                                 df = T)) %>% #extract values into df
                  map(., ~dplyr::select(., -ID)) %>% 
  map2(., popgrid2,
       ~bind_cols(.x, dplyr::select(st_drop_geometry(.y), FD, SpRich, name, long, lat))) %>% 
  map(., ~filter(., complete.cases(.))) 

# Build models

# Specify RF model formula
formula <-  map(retained_vars, ~as.formula(paste('FD','~',paste(., collapse = '+')))) 

# Fit model and predict over study domain (R^2 = 0.21)
RF.mod <- map2(formula, RF.data,
               ~ranger(.x, data = .y, num.trees = 10000, importance = 'impurity'))
pred <- map2(env_pred, RF.mod, 
             ~raster::predict(.x, .y, type = 'response', 
                        fun = function(model, ...) predict(model, ...)$predictions)) %>% 
  set_names(c('MAR','NGSL','NL','SGSL'))

# write predictions to tiff files

map2(pred, names(pred), ~writeRaster(.x, filename = paste0('Output/',.y,'_pred_FDiv.tiff')))

# Plot predictions with quantile breaks

# Palette with quantile breaks
map(pred, ~tm_shape(.) + 
      tm_raster(style = 'quantile', 
                palette = '-RdYlBu', 
                n = 5,
                contrast = c(0,1)))

# With continuous color palette
map(pred, ~tm_shape(.) + 
      tm_raster(style = 'cont',
                palette = '-RdYlBu', 
                contrast = c(0,1)))

#### Figure prep: spatial predictions of functional dispersion ####

#### Maritimes

pal <- rev(brewer.pal(n = 5, 'RdYlBu'))

# Shapefiles
StudyArea <- st_read('Data/Shapefiles/MaritimesRasterBoundary.shp') #study area boundaries
landfiles <- list.files(pattern = 'gadm', full.names = T) %>% 
  discard(~ grepl('FRA', .x)) %>% 
  discard(~ grepl('_0_', .x))
admin.borders <- map(landfiles, readRDS) %>% 
  map(., st_as_sf) %>% 
  rbindlist() %>% 
  st_as_sf() %>% 
  st_crop(., c(xmin = -71, ymin = 39, xmax = -45, ymax = 59)) %>% 
  st_transform(., crs = 26920)

# tmap version

map <-  tm_shape(pred$MAR, bbox = st_bbox(c(xmin =20000, 
                                           ymin = 4562134, 
                                           xmax = 1029435, 
                                           ymax = 5300000))) +
  tm_raster(legend.show = F,
            style = 'quantile', 
            palette = pal, 
            n = 5, 
            contrast = c(0,1)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +    
  tm_shape(admin.borders) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_scale_bar(breaks = c(0,50,100,200),
               color.dark = 'grey40',
               color.light = 'grey10',
               text.size = 4,
               position = c(0.01,0.90)) +
  tm_graticules(x = c(-67,-63,-59),
                y = c(47,45,43,41),
                alpha = 0.2,
                labels.size = 1.5) +
  tm_layout(title = 'MAR',
            title.position = c('right','top'),
            title.fontface = 'bold',
            title.size = 2,
            bg.color = '#e0e0e0')

tmap_save(map, filename = 'Output/MAR_FD_RFpred.tiff', compression = 'lzw')

#### Northern Gulf of St. Lawrence 

# Shapefiles
StudyArea <- st_read('Data/Shapefiles/QCRasterBoundary.shp') #study area boundaries

# tmap version
map <-  tm_shape(pred$NGSL, bbox = st_bbox(StudyArea), ext = 1.04) + 
  tm_raster(legend.show = F,
            style = 'quantile', 
            palette = pal, 
            n = 5, 
            contrast = c(0,1)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +  
  tm_shape(admin.borders) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_scale_bar(breaks = c(0,50,100,200),
               color.dark = 'grey40',
               color.light = 'grey10',
               text.size = 4,
               position = c(0.01,0.85)) +
  tm_graticules(x = c(-70,-66,-62,-58),
                y = c(48,50,52),
                alpha = 0.2,
                labels.size = 1.5) +
  tm_layout(title = 'NGSL',
            title.position = c(0.77,0.9),
            title.fontface = 'bold',
            title.size = 2,
            bg.color = '#e0e0e0')

tmap_save(map, filename = 'Output/NGSL_FD_RFpred.tiff', compression = 'lzw')

# Newfoundland and Labrador

# Shapefiles
StudyArea <- st_read('Data/Shapefiles/NLRasterBoundary.shp') #study area boundaries
Land <- st_read('Data/Shapefiles/NL_landborders.shp') # Land polygons

# tmap version

map <-  tm_shape(pred$NL, 
                 bbox = st_bbox(c(xmin =270169.8,
                                  ymin = 4742510.1,
                                  xmax = 1314169.8,
                                  ymax = 6398510.1)), 
                 ext = 1.04) + 
  tm_raster(legend.show = F,
            style = 'quantile', 
            palette = pal, 
            n = 5, 
            contrast = c(0,1)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +
  tm_shape(Land) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_scale_bar(breaks = c(0,100,200,400),
               color.dark = 'grey40',
               color.light = 'grey10',
               text.size = 4,
               position = c(0.01,0.005)) +
  tm_graticules(x = c(-48,-54,-60),
                y = c(44,48,52,56),
                alpha = 0.2,
                labels.size = 1.5) +
  tm_layout(title = 'NL',
            title.position = c('right','top'),
            title.fontface = 'bold',
            title.size = 2,
            bg.color = '#e0e0e0')

tmap_save(map, filename = 'Output/NL_FD_RFpred.tiff', compression = 'lzw')

# Southern Gulf of St. Lawrence

# Shapefiles

StudyArea <- st_read('Data/Shapefiles/GulfRasterBoundary.shp') #study area boundaries
Land <- st_read('Data/Shapefiles/Gulf_landborders.shp') # Land polygons

# tmap version

map <-  tm_shape(pred$SGSL, 
                 bbox = st_bbox(StudyArea), 
                 ext = 1.04) + 
  tm_raster(legend.show = F,
            style = 'quantile', 
            palette = pal, 
            n = 5, 
            contrast = c(0,1)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +
  tm_shape(Land) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_scale_bar(breaks = c(0,25,50,100),
               color.dark = 'grey40',
               color.light = 'grey10',
               text.size = 4,
               position = c(0.01,0.05)) +
  tm_graticules(x = c(-65,-63,-61),
                y = c(49,48,47,46),
                alpha = 0.2,
                labels.size = 1.5) +
  tm_layout(title = 'SGSL',
            title.position = c('right','top'),
            title.fontface = 'bold',
            title.size = 2,
            bg.color = '#e0e0e0')

tmap_save(map, filename = 'Output/SGSL_FD_RFpred.tiff', compression = 'lzw')

legend <- tm_shape(pred$SGSL) + 
  tm_raster(legend.show = F) +
  tm_add_legend(type = 'fill',
                col = pal,
                labels = c('1-20th percentile',
                           '21-40th percentile',
                           '41-60th percentile',
                           '61-80th percentile',
                           '81-100th percentile')) +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = T,
            legend.width = 10,
            legend.height = 1,
            legend.title.fontface = 'bold'); legend

tmap_save(legend, filename = 'Output/legend_FD_RFpred.tiff', compression = 'lzw')


# load('Data/MaritimesData.RData')
# 
# year.breaks <- c(1970,1975,1980,1985,1990,1995,2000,2005,2010,2015,2020)
# 
# RV <- Maritimes %>% 
#   filter(., !(vertical.position %in% c('pelagic', 'pelagic,neritic','pelagic,oceanic'))) %>%
#   mutate(., subphylum = if_else(species == 'Hippoglossina oblonga', 'Vertebrata', subphylum)) %>% 
#   filter(., subphylum == 'Vertebrata') %>% 
#   mutate(., year.bin = cut(year, 
#                            year.breaks,
#                            right = F, 
#                            labels = c('1970-1974','1975-1979','1980-1984','1985-1989',
#                                       '1990-1994','1995-1999','2000-2004','2005-2009',
#                                       '2010-2014','2015-2019')))
# 
# SpXYr <- table(RV$species, RV$year.bin) %>% 
#   data.frame() %>% 
#   rename(species = Var1, year.bin = Var2) %>% 
#   mutate(., Freq = if_else(Freq == 0, NA_integer_, Freq))
# 
# p1 <- ggplot(SpXYr, aes(x = year.bin, y = reorder(species, desc(species)), fill = Freq)) +
#   geom_tile(width = 1, height = 1, colour = 'grey50', size = 0.5) +
#   scale_fill_viridis_c(guide = 'colorbar') +
#   scale_y_discrete(expand = c(0,0)) +
#   scale_x_discrete(expand = c(0,0)) +
#   labs(y = 'Taxon') +
#   theme(panel.background = element_rect(fill = 'grey85', colour = 'grey50'),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.y = element_text(face = 'italic'))
# 
# ggsave(plot = p1, 'C:/Users/obrienjoh/Desktop/MAR_matrix_RVcatch.tiff', 
#        width = 10,height = 20, unit = 'in',compression = 'lzw')
