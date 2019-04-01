# Model cluster assignment based on environmental correlates using Random Forest
# John O'Brien (John.O'brien@dfo-mpo.gc.ca) _ Created 2019
#
## Code for Biological classification analysis - Eastern Canadian groundfish - Stanley, Heaslip, Jeffery, O'Brien

#Outline of code workflow

#1) Resample environmental rasters to resolution of grids used in biological classification (See 4Environmental_layers_SPERA.R)
#2) Decide which variables to include using criteria of biological relevance &
#3) collinearity (examine correlation coefficients among variables w/ 'corrplot')
#4)	Choose size of training dataset and out-of-bag samples
#5) Fit a classification RF model to data 
# - use SpatialPolygonsDataframe with cluster assignment as attribute
# - extract predictor values from rasters for each populated grid cell polygon
# - make sure rasters are aligned (same extent, origin and resolution)
# - make sure grid cells align with raster grid
#6)	Cross-validation of best model
#7)	Examine variable importance plots
#8) Predict cluster membership in non-survey grids
#9) Include maximum cluster probability or proportion of votes as measure of
#uncertainty
#10) Project cluster membership under future climate scenarios

#compare with species archtype model? (see Murillo et al. 2018 for example)

####load required packages####
library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(dplyr)
library(tidyr)
library(purrr)
library(corrplot)
library(tcltk2)
library(randomForest)
library(pROC)
library(ggplot2)
library(gridExtra)
library(gtable)

####Remove correlated variables####

LandBuffer <- readOGR('Data/Shapefiles/Maritimes_LandBuffer_5km.shp')

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Mar', full.names = T)
env_pred <- stack(raster.list)
env_pred <- mask(env_pred, LandBuffer, inverse = T)
pearson_corr <- layerStats(env_pred, 'pearson', na.rm = T)
corr_matrix <- pearson_corr[[1]]
corr_plot <- corrplot(corr_matrix)

#layerStats fx can calculate pearson correlation b/w layers, but Spearman may be more informative
#will need to extract values from each raster and create correlation matrix to do this

vals <-  c()#create empty list for values
for (i in 1:nlayers(env_pred)){ #for each raster layer
  v <-  getValues(raster(env_pred, i)) #get 1D array of values
  vals <-  cbind(vals, v) #add 1D array to col of vals
}
colnames(vals) <- gsub('Mar_','',names(env_pred))  #rename columns to variable names

#Compare Pearson vs. Spearman correlation
corr.Pearson <- cor(vals, use='pairwise.complete.obs', method = 'pearson')
corr.Spearman <- cor(vals, use='pairwise.complete.obs', method = 'spearman')
plot(corr.Spearman ~ corr.Pearson) #similar outcome (Rsq = 0.88)

diag(corr.Spearman) <- NA #set same-variable correlations to NA
corr.Spearman <- abs(corr.Spearman) # coerce correlations to absolute value

#Run GUI to select variable pair with highest correlation, subjectively remove one,
#continue until no variables remain with correlation higher than threshold (0.5)
threshold <- 0.65
#Set up filenames
tmpfilename <- 'temporary.txt'
removefilename <- 'eliminate_list.txt'
retainfilename <- 'retain_list.txt'

eliminated <- c('') #empty list for eliminated variables
maxcor <- max(corr.Spearman, na.rm = T) # find highest correlation b/w predictors
maxrow <- which(corr.Spearman == maxcor, arr.ind = T)[1,1] # row index of maxcor
maxcol <- which(corr.Spearman == maxcor, arr.ind = T)[1,2] # col index of maxcor

#write functions for 2 buttons (variable names of most highly correlated variables)

write1 <- function(){#if button 1 is pushed
  tmpfile <- file(tmpfilename) #create temporary file
  writeLines(colnames(corr.Spearman)[maxrow], tmpfile) #write name of maxrow (variable to be eliminated) to tmpfile
  close(tmpfile)
}
write2 <- function(){#if button 2 is pushed
  tmpfile <- file(tmpfilename)
  writeLines(colnames(corr.Spearman)[maxcol], tmpfile) #write name of maxcol (variable to be elimnated) to tmpfile
  close(tmpfile)
}

#create GUI widget
root <- tktoplevel()
btn1 <- tk2button(root, text = paste(colnames(corr.Spearman)[maxrow]), command = write1)
tkpack(btn1)
btn2 <- tk2button(root, text = paste(colnames(corr.Spearman)[maxcol]), command = write2)
tkpack(btn2)

#create 'while' loop that stops running when maxcor falls below 0.5
while(maxcor > threshold){
  tkconfigure(btn1, text = paste(colnames(corr.Spearman)[maxrow]), command = write1)
  tkconfigure(btn2, text = paste(colnames(corr.Spearman)[maxcol]), command = write2)
  while(!file.exists(tmpfilename)){a = 1} #wait until button is pushed/tmpfile is written
  tmpfile <- file(tmpfilename, 'rw') #open tmpfile
  eliminate <- readLines(tmpfile,1) #read one line
  close(tmpfile) #close tmpfile
  file.remove(tmpfilename) #delete tmpfile
  eliminated <- c(eliminated, eliminate) #eliminate variable from correlation matrix
  corr.Spearman <- corr.Spearman[-which(colnames(corr.Spearman)==eliminate),-which(colnames(corr.Spearman)==eliminate)]
  if(dim(corr.Spearman)[1]>1){ #if at least 2 variables remain, recalculate maxcor/row/col
    maxcor = max(corr.Spearman, na.rm = T)
    maxrow <- which(corr.Spearman == maxcor, arr.ind = T)[1,1] 
    maxcol <- which(corr.Spearman == maxcor, arr.ind = T)[1,2]
  } else maxcor = 0 #otherwise end the loop
}

tkdestroy(root) #remove widget
cornames <- colnames(corr.Spearman) #write results to files
write.table(cornames, file = retainfilename, row.names = F)
rmlist <- file(removefilename)
writeLines(eliminated, rmlist)
close(rmlist)

####Random Forest Model fitting####

#import predictor variables retained for modelling

retained_vars <- readLines(file('Output/Mar_retain_list.txt'))#list of retained variables
raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Mar', full.names = T) #get list of raster files
RV_strata <- readOGR('Data/Shapefiles/MaritimesRegionStrataAgg.shp')
LandBuffer <- readOGR('Data/Shapefiles/Maritimes_LandBuffer_5km.shp') # 5km buffer around land points

env_pred <- stack(raster.list) %>%  #stack rasters in list
  `names<-`(.,gsub('Mar_','',names(.))) %>% #remove 'Mar_' from layer names
   .[[retained_vars]] %>%  #subset raster stack to include only retained variables
  mask(., RV_strata) %>%  #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) #also mask cells overlapped by 5km land buffer

#import spatial polygons dataframe of populated 4km grid cells with cluster assignment as 
#attribute variable

MarGrid_populated <- readOGR('Data/Shapefiles/GridClusterAssignment4km_Maritimes.shp')
minor_cl <- factor(which(summary(MarGrid_populated$cl) < 20)) #minor clusters
#exclude minor clusters
MarGrid_populated <- MarGrid_populated[-which(MarGrid_populated@data$cl %in% minor_cl),] 

#Extract raster values at each populated grid cell
p.data <- raster::extract(env_pred,MarGrid_populated, factors = T, nl = nlayers(env_pred), df = T) #extract values into df
RF.data <- data.frame(coordinates(MarGrid_populated), MarGrid_populated@data, p.data[,-1]) #join environmental data with grid df
RF.data <- RF.data %>% rename(., x = X1, y = X2)
missing.data <- missingdata <- RF.data[!complete.cases(RF.data),]#incomplete cases
table(missing.data$cl) #cluster 8 mostly effected by incomplete cases (23) because location in inner BoF
RF.data <- RF.data[complete.cases(RF.data),] #remove incomplete cases
RF.data <- droplevels(RF.data) #drop unused factor levels from dataframe & recode retained factors
RF.data$cl <- recode_factor(RF.data$cl, `1` = '1', `5` = '2', `4` = '3', `2` = '4', `9` = '5', `8` = '6') #Slope, LC,ESS, ESS:Banks, WSS/Inner BoF, WSS:Banks/Inner BoF

#build model

# specify model formula
formula <-  as.formula(paste('cl','~',paste(retained_vars, collapse = '+'))) 
#fit model, 20.02% OOB error rate (higher for ESS 30.94%)
RF.mod1 <- randomForest(formula, data = RF.data, ntree = 10000, importance = T) 
saveRDS(RF.mod1, 'Output/Maritimes_RF_model.rds')

####Evaluate model accuracy w/ 10-fold cross-validation####

nfold <-  10 #number of splits of the data
n <- nrow(RF.data) #number of observations in data
groups <- sample(rep(1:nfold, length = n), n) #assign observation to 1 of 10 splits of the data (groups)
CalVal.list <- 1:nfold %>% map(~ list(cal = which(groups != .), val = which(groups == .))) #indices for 10 calibration and validation datasets

cv10 <- CalVal.list %>% 
  map(~ list(model = randomForest(formula, data = RF.data[.$cal,], ntree = 10000), val = RF.data[.$val,])) %>% 
  map(~ list(obs = as.numeric(.$val$cl), predictions = predict(.$model, newdata = .$val, type = 'prob'))) %>%  
  map(~ multiclass.roc(.$obs, .$predictions)$auc) #compute multi-class AUC using each calibration/validation pair 
AUC <- mean(unlist(cv10)) #average AUC for each 90:10 split of data
#High AUC 0.962 suggests high predictive power of model

####Predict classification for grids without biological data and plot####
predict.map <- raster::predict(env_pred, RF.mod1) # predict cluster membership over raster surfaces
my.colors = c('#e6ab02','#6a3d9a','#33a02c','#ff7f00', '#a6cee3','#1f78b4') # color scheme for rasters
plot(predict.map, col = my.colors)
writeRaster(predict.map, filename = 'Output/Maritimes_PredClust_map.tif', overwrite = T) #write prediction to file
# raster_boundary <- aggregate(rasterToPolygons(predict.map))
# writeOGR(as(raster_boundary, 'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'MaritimesRasterBoundary', driver = 'ESRI Shapefile')

####Variable importance plots####

varImp <- data.frame(importance(RF.mod1, scale = F))
varImp$predictor <- c('Aspect','Avg max PP (spring/summer)','Depth','DO','BPI (0.5 km)', 'BPI (10 km)',
                      'Avg max salinity','Avg max temperature','Avg min temperature',
                      'Avg mean bottom stress','Slope')
varImp <- varImp[order(-varImp$MeanDecreaseAccuracy),]
varImp$predictor <- factor(varImp$predictor, 
                           levels = varImp$predictor[order(varImp$MeanDecreaseAccuracy)])
colnames(varImp)[1:7] <- c('Slope','Laurentian Channel/\nShelf Break', 'ESS',
                           'ESS: Banks','WSS/Outer BoF', 'WSS: Banks/\nInner BoF','Whole Model')
varImp <- gather(varImp, 'class', 'MeanDecreaseAccuracy', 1:7)
varImp$class <- factor(varImp$class, levels = c('Whole Model','Slope','Laurentian Channel/\nShelf Break', 'ESS',
                                                'ESS: Banks','WSS/Outer BoF', 'WSS: Banks/\nInner BoF'))

tiff('Output/Maritimes_VarImpPlot.tiff', width = 9, height = 2.5, units = 'in', res = 300)                                                                                                
p1 <- ggplot(varImp, aes(x = predictor, y = MeanDecreaseAccuracy)) +
  geom_point(stat = 'identity') +
  facet_grid(~class) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = 'Mean Decrease in Accuracy') + 
  theme(axis.text.x = element_text(angle = 90), 
        text = element_text(size = 10)); p1
dev.off()

####Highlighting areas of lower model certainty####

#extract environmental raster cell values into dataframe
env_pred_df <- data.frame(getValues(env_pred)) 
#get RF predictions for environmental layers cell values as probability
RF.predictions <- predict(object = RF.mod1, newdata = env_pred_df, type = 'prob')

#extract proportion of vote counts for cluster predicted by model
#for each raster cell (i.e. extract maximum probability)

VoteCounts <- data.frame(RF.predictions) %>%
  rowwise() %>%
  mutate(maxVC = max(c(X1,X2,X3,X4,X5,X6))) %>% 
  dplyr::select(maxVC) %>%
  data.frame()

# create raster layer with max proportion of vote counts
VoteCountsRaster <- raster(env_pred[[1]]) 
VoteCountsRaster <- setValues(VoteCountsRaster, as.vector(VoteCounts$maxVC))
#convert raster layer to a spatial polygons data frame (subset cells < 0.7 only, and <0.5 only)
VoteCountsShape0.5 <- rasterToPolygons(VoteCountsRaster, fun = function(x){x<0.5}, na.rm = T)
VoteCountsShape0.7 <- rasterToPolygons(VoteCountsRaster, fun = function(x){x<0.7}, na.rm = T)
#write SPDFs to new shapefiles
writeOGR(VoteCountsShape0.5, dsn = 'Output', layer = 'Maritimes_RF_uncertainty_0.5', driver = 'ESRI Shapefile', overwrite_layer = T) 
writeOGR(VoteCountsShape0.7, dsn = 'Output', layer = 'Maritimes_RF_uncertainty_0.7', driver = 'ESRI Shapefile', overwrite_layer = T) 

####Evaluate uncertainty threshold  by examining relationship between assignment probability and success of prediction####

uncertainty <- data.frame(RF.mod1$votes) %>% 
  rowwise() %>% 
  mutate(maxVC = max(c(X1,X2,X3,X4,X5,X6))) %>% 
  dplyr::select(maxVC) %>%
  data.frame(., predicted = RF.mod1$predicted, actual = RF.mod1$y) %>% 
  mutate(correct = as.integer(predicted == actual)) %>% 
  arrange(., predicted)

roc(uncertainty$correct,uncertainty$maxVC, plot = T) #examine ROC. Threshold should favour highlighting potential
                                                     #false positives at the expense ofthe true positive rate

#examine the true negative/false positive rates at the 0.5 & 0.7 thresholds proposed above

uncertainty %>%#
  split(.$correct) %>% #split data by whether classification was correct
  map(~ stats::quantile(.$maxVC, probs = seq(0.1,0.8,0.05))) #calculate range of percentiles
#Thresholds of 0.5 & 0.7, result in true negative rate of ~25% and ~70% respectively (i.e. FPR = 75% & 30%)
uncertainty %>% filter(maxVC < 0.7) %>% count(.,correct) #False negative rate of ~30%

plot(uncertainty$correct ~ uncertainty$maxVC) #examine relationship between assignment probability & assignment success

####Distribution of important variables across cluster####

my.colors = c('#e6ab02','#6a3d9a','#33a02c','#ff7f00', '#a6cee3','#1f78b4') # color scheme for clusters

#boxplot for depth
box_depth <- ggplot(RF.data, aes(x = cl,y = bathy, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  labs(x = NULL, y = 'Depth (m)') +
  scale_x_discrete(labels = c('Slope','LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF', 'WSS: Banks/Inner BoF')) +
  theme(axis.text.x = element_text(colour = 'white'), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,0,'pt')),
        axis.ticks.length = unit(2, 'mm')); box_depth

#boxplot for min ann bottom temp  
box_minT <- ggplot(RF.data, aes(x = cl,y = min_ann_BT, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  labs(x = NULL, y = expression(paste('   Average min\ntemperature (\u00B0C)'))) +
  scale_x_discrete(labels = c('Slope','LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF', 'WSS: Banks/Inner BoF')) +
  theme(axis.text.x = element_text(colour = 'white'), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,0,'pt')),
        axis.ticks.length = unit(2, 'mm')); box_minT

#boxplot for max ann bottom temp
box_maxT <- ggplot(RF.data, aes(x = cl,y = max_ann_BT, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  scale_y_continuous(limits = c(0,16)) +
  labs(x = NULL, y = expression(paste('   Average max\ntemperature (\u00B0C)'))) +
  scale_x_discrete(labels = c('Slope','LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF', 'WSS: Banks/Inner BoF')) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,0,'pt'), hjust = 0.5),
        axis.ticks.length = unit(2, 'mm')); box_maxT

#boxplot for max ann bottom salinity
box_maxSal <- ggplot(RF.data, aes(x = cl,y = max_ann_BSal, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  labs(x = NULL, y = expression(paste('Average max salinity (\u2030)'))) +
  scale_x_discrete(labels = c('Slope','LC/Shelf Break','ESS','ESS: Banks','WSS/Outer BoF', 'WSS: Banks/Inner BoF')) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,0,'pt')),
        axis.ticks.length = unit(2, 'mm')); box_maxSal

#coerce ggplot objects to graphical objects
g1 <- ggplotGrob(box_depth)
g2 <- ggplotGrob(box_maxSal)
g3 <- ggplotGrob(box_minT)
g4 <- ggplotGrob(box_maxT)

#Arrange grobs and plot in 2X2 array
r1 <- cbind(g1, g3, size = 'first') #bind/align plot elements of row 1
r2 <- cbind(g2, g4, size = 'first') #bind/align plot elements of row 2
g <- rbind(r1, r2, size = 'first') #bind/align plot elements of both rows

tiff('Output/Maritimes_EnvVariation.tiff', width = 8, height = 8, units = 'in', res = 300)
plot(g)
dev.off()

####Build Model without depth####

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Mar', full.names = T) #get list of raster files
env_pred <- stack(raster.list) #stack rasters in list
names(env_pred) <- gsub('Mar_','',names(env_pred)) #remove 'Mar_' from layer names
retained_vars <- readLines(file('Output/Mar_retain_list.txt'))#list of retained variables
retained_vars <- retained_vars[-3] #remove depth from predictors
env_pred <- env_pred[[retained_vars]] #subset raster stack to include only retained variables
RV_strata <- readOGR('Data/Shapefiles/MaritimesRegionStrataAgg.shp')
LandBuffer <- readOGR('Data/Shapefiles/Maritimes_LandBuffer_5km.shp') # 5km buffer around land points
env_pred <- mask(env_pred, RV_strata) #mask raster cells of predictors outside RV survey boundaries
env_pred <- mask(env_pred, LandBuffer, inverse = T) #also mask cells overlapped by 5km land buffer

#import spatial polygons dataframe of populated 4km grid cells with cluster assignment as 
#attribute variable

MarGrid_populated <- readOGR('Data/Shapefiles/GridClusterAssignment4km_Maritimes.shp')
minor_cl <- factor(names(table(MarGrid_populated$cl))[which(table(MarGrid_populated$cl) < 20)]) #minor clusters
#exclude minor clusters
MarGrid_populated <- MarGrid_populated[-which(MarGrid_populated@data$cl %in% minor_cl),] 

#Extract raster values at each populated grid cell
p.data <- raster::extract(env_pred,MarGrid_populated, factors = T, nl = nlayers(env_pred), df = T) #extract values into df
RF.data <- data.frame(coordinates(MarGrid_populated), MarGrid_populated@data, p.data[,-1]) #join environmental data with grid df
RF.data <- RF.data %>% rename(., x = X1, y = X2)
missing.data <- missingdata <- RF.data[!complete.cases(RF.data),]#incomplete cases
table(missing.data$cl) #cluster 8 mostly effected by incomplete cases (23) because location in inner BoF
RF.data <- RF.data[complete.cases(RF.data),] #remove incomplete cases
RF.data <- droplevels(RF.data) #drop unused factor levels from dataframe & recode retained factors
RF.data$cl <- recode_factor(RF.data$cl, `1` = '1', `5` = '2', `4` = '3', `2` = '4', `9` = '5', `8` = '6') #Slope, LC,ESS, ESS:Banks, WSS/Inner BoF, WSS:Banks/Inner BoF

#build model

# specify model formula
formula <-  as.formula(paste('cl','~',paste(retained_vars, collapse = '+'))) 
#fit model, 22.75% OOB error rate (higher for ESS 44.84%)
RF.mod1 <- randomForest(formula, data = RF.data, ntree = 10000, importance = T) 
saveRDS(RF.mod1, 'Output/Maritimes_RF_model_NoDepth.rds')

#Evaluate model accuracy w/ 10-fold cross-validation

nfold <-  10 #number of splits of the data
n <- nrow(RF.data) #number of observations in data
groups <- sample(rep(1:nfold, length = n), n) #assign observation to 1 of 10 splits of the data (groups)
CalVal.list <- 1:nfold %>% map(~ list(cal = which(groups != .), val = which(groups == .))) #indices for 10 calibration and validation datasets

cv10 <- CalVal.list %>% 
  map(~ list(model = randomForest(formula, data = RF.data[.$cal,], ntree = 10000), val = RF.data[.$val,])) %>% 
  map(~ list(obs = as.numeric(.$val$cl), predictions = as.numeric(predict(.$model, newdata = .$val)))) %>% 
  map(~ multiclass.roc(.$obs, .$predictions)$auc) #compute multi-class AUC using each calibration/validation pair 
AUC <- mean(unlist(cv10)) #average AUC for each 90:10 split of data
#High AUC 0.902 suggests high predictive power of model

#Predict classification for grids without biological data and plot#
predict.map <- raster::predict(env_pred, RF.mod1) # predict cluster membership over raster surfaces
my.colors = c('#e6ab02','#6a3d9a','#33a02c','#ff7f00', '#a6cee3','#1f78b4') # color scheme for rasters
plot(predict.map, col = my.colors)
writeRaster(predict.map, 'Output/Maritimes_PredClust_map_NoDepth.tif', overwrite = T) #save raster