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

LandBuffer <- readOGR('Data/Shapefiles/Gulf_LandBuffer_5km.shp')

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Gulf', full.names = T)
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
colnames(vals) <- gsub('Gulf_','',names(env_pred))  #rename columns to variable names

#Compare Pearson vs. Spearman correlation
corr.Pearson <- cor(vals, use='pairwise.complete.obs', method = 'pearson')
corr.Spearman <- cor(vals, use='pairwise.complete.obs', method = 'spearman')
plot(corr.Spearman ~ corr.Pearson) #similar outcome (Rsq = 0.9259)
summary(lm(as.vector(corr.Spearman)~as.vector(corr.Pearson)))

diag(corr.Spearman) <- NA #set same-variable correlations to NA
corr.Spearman <- abs(corr.Spearman) # coerce correlations to absolute value
diag(corr.Pearson) <- NA
corr.Pearson <- abs(corr.Pearson)

#Run GUI to select variable pair with highest correlation, subjectively remove one,
#continue until no variables remain with correlation higher than threshold (0.7)
threshold <- 0.7
#Set up filenames
tmpfilename <- 'temporary.txt'
removefilename <- 'eliminate_list.txt'
retainfilename <- 'retain_list.txt'

eliminated <- c('') #empty list for eliminated variables
maxcor <- max(corr.Pearson, na.rm = T) # find highest correlation b/w predictors
maxrow <- which(corr.Pearson == maxcor, arr.ind = T)[1,1] # row index of maxcor
maxcol <- which(corr.Pearson == maxcor, arr.ind = T)[1,2] # col index of maxcor

#write functions for 2 buttons (variable names of most highly correlated variables)

write1 <- function(){#if button 1 is pushed
  tmpfile <- file(tmpfilename) #create temporary file
  writeLines(colnames(corr.Pearson)[maxrow], tmpfile) #write name of maxrow (variable to be eliminated) to tmpfile
  close(tmpfile)
}
write2 <- function(){#if button 2 is pushed
  tmpfile <- file(tmpfilename)
  writeLines(colnames(corr.Pearson)[maxcol], tmpfile) #write name of maxcol (variable to be elimnated) to tmpfile
  close(tmpfile)
}

#create GUI widget
root <- tktoplevel()
btn1 <- tk2button(root, text = paste(colnames(corr.Pearson)[maxrow]), command = write1)
tkpack(btn1)
btn2 <- tk2button(root, text = paste(colnames(corr.Pearson)[maxcol]), command = write2)
tkpack(btn2)

#create 'while' loop that stops running when maxcor falls below 0.7
while(maxcor > threshold){
  tkconfigure(btn1, text = paste(colnames(corr.Pearson)[maxrow]), command = write1)
  tkconfigure(btn2, text = paste(colnames(corr.Pearson)[maxcol]), command = write2)
  while(!file.exists(tmpfilename)){a = 1} #wait until button is pushed/tmpfile is written
  tmpfile <- file(tmpfilename, 'rw') #open tmpfile
  eliminate <- readLines(tmpfile,1) #read one line
  close(tmpfile) #close tmpfile
  file.remove(tmpfilename) #delete tmpfile
  eliminated <- c(eliminated, eliminate) #eliminate variable from correlation matrix
  corr.Pearson <- corr.Pearson[-which(colnames(corr.Pearson)==eliminate),-which(colnames(corr.Pearson)==eliminate)]
  if(dim(corr.Pearson)[1]>1){ #if at least 2 variables remain, recalculate maxcor/row/col
    maxcor = max(corr.Pearson, na.rm = T)
    maxrow <- which(corr.Pearson == maxcor, arr.ind = T)[1,1] 
    maxcol <- which(corr.Pearson == maxcor, arr.ind = T)[1,2]
  } else maxcor = 0 #otherwise end the loop
}

tkdestroy(root) #remove widget
cornames <- colnames(corr.Pearson) #write results to files
write.table(cornames, file = retainfilename, row.names = F)
rmlist <- file(removefilename)
writeLines(eliminated, rmlist)
close(rmlist)

####Random Forest Model fitting####

#import predictor variables retained for modelling

raster.list <- list.files(path = 'Data/Rasters/', pattern = 'Gulf', full.names = T) #get list of raster files
Gulf_SA <- readOGR('Data/Shapefiles/Gulf_RVsurveyAgg.shp') # Gulf_SA (RV survey boundaries)
LandBuffer <- readOGR('Data/Shapefiles/Gulf_LandBuffer_5km.shp') # 5km buffer around land points
retained_vars <- readLines(file('Output/Gulf_retain_list.txt'))#list of retained variables
env_pred <- stack(raster.list) %>% #stack rasters in list
  `names<-`(.,gsub('Gulf_','',names(.))) %>% #remove 'Gulf_' from layer names
  .[[retained_vars]] %>% #subset raster stack to include only retained variables
  mask(., Gulf_SA) %>% #mask raster cells of predictors outside RV survey boundaries
  mask(., LandBuffer, inverse = T) #also mask cells overlapped by 5km land buffer

#import spatial polygons dataframe of populated 4km grid cells with cluster assignment as 
#attribute variable

GulfGrid_populated <- readOGR('Data/Shapefiles/GridClusterAssignment4km_Gulf.shp')
minor_cl <- factor(names(table(GulfGrid_populated$cl))[which(table(GulfGrid_populated$cl) < 20)]) #minor clusters
#exclude minor clusters
GulfGrid_populated <- GulfGrid_populated[-which(GulfGrid_populated@data$cl %in% minor_cl),] 

#Extract raster values at each populated grid cell
p.data <- raster::extract(env_pred,GulfGrid_populated, factors = T, nl = nlayers(env_pred), df = T) #extract values into df
RF.data <- data.frame(coordinates(GulfGrid_populated), GulfGrid_populated@data, p.data[,-1]) #join environmental data with grid df
RF.data <- RF.data %>% rename(., x = X1, y = X2)
missing.data <- missingdata <- RF.data[!complete.cases(RF.data),]#incomplete cases
table(missing.data$cl) #Inshore/MI and NS/SGB mostly effected by incomplete cases (11 & 8)
RF.data <- RF.data[complete.cases(RF.data),] #remove incomplete cases
RF.data <- droplevels(RF.data) #drop unused factor levels from dataframe & recode retained factors
RF.data$cl <- recode_factor(RF.data$cl, `7` = '1', `6` = '2', `9` = '3', `3` = '4')

#build model

# specify model formula
formula <-  as.formula(paste('cl','~',paste(retained_vars, collapse = '+'))) 
#fit model, 9.14% OOB error rate although both nearshore clusters have higher OOB (22/20.8%)
RF.mod1 <- randomForest(formula, data = RF.data, ntree = 10000, importance = T) 
RF.output <- saveRDS(RF.mod1, 'Output/Gulf_RF_model.rds')

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
#High AUC 0.979 suggests high predictive power of model

####Predict classification for grids without biological data and plot####
predict.map <- raster::predict(env_pred, RF.mod1) # predict cluster membership over raster surfaces
my.colors <- c("#a6cee3", "#ff7f00", "#6a3d9a", "#1f78b4") # color scheme for rasters
plot(predict.map, col = my.colors)
writeRaster(predict.map, filename = 'Output/Gulf_PredClust_map.tif', overwrite = T) #write prediction to file
raster_boundary <- aggregate(rasterToPolygons(predict.map))
writeOGR(as(raster_boundary, 'SpatialPolygonsDataFrame'), dsn = 'Data/Shapefiles', layer = 'GulfRasterBoundary', driver = 'ESRI Shapefile', overwrite_layer = T)

####Variable importance plots####

varImp <- data.frame(importance(RF.mod1, scale = F))
varImp$predictor <- c('Aspect','Avg max PP (spring/summer)', 'BPI (1 km)', 'BPI (20 km)',
                      'Avg max temperature', 'Avg max current velocity (NS)',
                      'Avg min temperature','Avg min SST','Avg mean temperature',
                      'Avg mean bottom stress', 'Average mean MLD', 'Avg mean current velocity (EW)', 
                      'Avg mean current velocity (NS)', 'Average mean MLD (summer)',
                      'Avg mean SST (winter)', 'Slope')
varImp <- varImp[order(-varImp$MeanDecreaseAccuracy),]
varImp$predictor <- factor(varImp$predictor, 
                           levels = varImp$predictor[order(varImp$MeanDecreaseAccuracy)])
colnames(varImp)[1:5] <- c("Magdalen Shallows/\nChaleur Bay","Inshore/\nMagdalen Is.","Laurentian Channel", 
                           "Northumberland Strait/\nSt. George's Bay","Whole Model")
varImp <- gather(varImp, 'class', 'MeanDecreaseAccuracy', 1:5)
varImp$class <- factor(varImp$class, levels = c('Whole Model',"Magdalen Shallows/\nChaleur Bay","Inshore/\nMagdalen Is.",
                                                "Laurentian Channel", "Northumberland Strait/\nSt. George's Bay"))

tiff('Output/Gulf_VarImpPlot.tiff', width = 9, height = 3.5, units = 'in', res = 300)                                                                                                
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

#extract environemtnal raster cell values into dataframe
env_pred_df <- data.frame(getValues(env_pred)) 
#get RF predictions for environmental layers cell values as probability
RF.predictions <- predict(object = RF.mod1, newdata = env_pred_df, type = 'prob')

#extract proportion of vote counts for cluster predicted by model
#for each raster cell (i.e. extract maximum probability)

VoteCounts <- data.frame(RF.predictions) %>%
  dplyr::rowwise() %>% 
  mutate(maxVC = max(c(X1,X2,X3,X4))) %>% 
  dplyr::select(maxVC) %>%
  data.frame()
  
# create raster layer with max proportion of vote counts
VoteCountsRaster <- raster(env_pred[[1]]) 
VoteCountsRaster <- setValues(VoteCountsRaster, as.vector(VoteCounts$maxVC))
#convert raster layer to a spatial polygons data frame (subset cells < 0.7 only, and <0.5 only)
VoteCountsShape0.5 <- rasterToPolygons(VoteCountsRaster, fun = function(x){x<0.5}, na.rm = T)
VoteCountsShape0.7 <- rasterToPolygons(VoteCountsRaster, fun = function(x){x<0.7}, na.rm = T)
#write SPDFs to new shapefiles
writeOGR(VoteCountsShape0.5, dsn = 'Output', layer = 'Gulf_RF_uncertainty_0.5', driver = 'ESRI Shapefile', overwrite_layer = T) 
writeOGR(VoteCountsShape0.7, dsn = 'Output', layer = 'Gulf_RF_uncertainty_0.7', driver = 'ESRI Shapefile', overwrite_layer = T) 

####Evaluate uncertainty threshold  by examining relationship between assignment probability and success of prediction####

uncertainty <- data.frame(RF.mod1$votes) %>% 
  rowwise() %>% 
  mutate(maxVC = max(c(X1,X2,X3,X4))) %>% 
  dplyr::select(maxVC) %>%
  data.frame(., predicted = RF.mod1$predicted, actual = RF.mod1$y) %>% 
  mutate(correct = as.integer(predicted == actual)) %>% 
  arrange(., predicted)

roc(uncertainty$correct,uncertainty$maxVC, plot = T) #examine ROC. Threshold should favour highlighting potential
#false positives at the expense ofthe true positive rate

#examine the true negative/false positive rates at the 0.5 & 0.7 thresholds proposed above

uncertainty %>%
  split(.$correct) %>% #split data by whether classification was correct
  map(~ stats::quantile(.$maxVC, probs = seq(0.1,0.7,0.05))) #calculate range of percentiles
#0.5 & 0.7 threshold results in true negative rate of ~10% and ~55% respectively (i.e. FPR = 90% & 45%)
uncertainty %>% filter(maxVC < 0.7) %>% count(.,correct) #False negative rate of ~9%

plot(uncertainty$correct ~ uncertainty$maxVC) #examine relationship between assignment probability & assignment success

####Distribution of important variables across cluster####

my.colors <- c("#a6cee3", "#ff7f00", "#6a3d9a", "#1f78b4") # color scheme for clusters

#boxplot for max ann bottom temp
box_maxT <- ggplot(RF.data, aes(x = cl,y = max_ann_BT, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  scale_y_continuous(limits = c(0,25)) +
  labs(x = NULL, y = expression(paste('   Avg max\ntemperature (\u00B0C)'))) +
  theme(axis.text.x = element_blank(), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,12,'pt'), hjust = 0.5),
        axis.ticks.length = unit(2, 'mm')); box_maxT

#boxplot for mean ann bottom temp
box_meanT <- ggplot(RF.data, aes(x = cl,y = mn_ann_BT, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  scale_y_continuous(limits = c(0,10)) +
  labs(x = NULL, y = expression(paste('   Avg mean\ntemperature (\u00B0C)'))) +
  theme(axis.text.x = element_blank(), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,12,'pt'), hjust = 0.5),
        axis.ticks.length = unit(2, 'mm')); box_meanT

#boxplot for min ann bottom temp  
box_minT <- ggplot(RF.data, aes(x = cl,y = min_ann_BT, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  labs(x = NULL, y = expression(paste('   Avg min\ntemperature (\u00B0C)'))) +
  scale_x_discrete(labels = c("Magdalen Shallows\n& Chaleur Bay", "Inshore/Magdalen Is.",
                              "Laurentian Channel","N. Strait &\n St. George's Bay")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,12,'pt')),
        axis.ticks.length = unit(2, 'mm')); box_minT

#boxplot for Avg max PP (spring/summer)  
box_PP <- ggplot(RF.data, aes(x = cl,y = avg_max_sprsum_PP, fill = cl)) +
  geom_boxplot(show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = my.colors) +
  labs(x = NULL, y = expression(Avg~max~PP~'('*mg~C~m^{-2}~day^{-1}*')')) +
  scale_x_discrete(labels = c("Magdalen Shallows\n& Chaleur Bay", "Inshore/Magdalen Is.",
                              "Laurentian Channel","N. Strait &\n St. George's Bay")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), text = element_text(size = 12), 
        axis.title.y = element_text(margin = ggplot2::margin(0,12,0,12,'pt')),
        axis.ticks.length = unit(2, 'mm')); box_PP

#coerce ggplot objects to graphical objects
g1 <- ggplotGrob(box_maxT)
g2 <- ggplotGrob(box_meanT)
g3 <- ggplotGrob(box_minT)
g4 <- ggplotGrob(box_PP)

#Arrange grobs and plot in 2 X 2 array

c1 <- rbind(g1, g3, size = 'first') #bind/align plot elements of column 1
c2 <- rbind(g2, g4, size = 'first') #bind/align plot elements of column 2
g <- cbind(c1, c2, size = 'first') #bind/align plot elements of both rows

tiff('Output/Gulf_EnvVariation.tiff', width = 14, height = 8, units = 'in', res = 300)
plot(g)
dev.off()

