########################################
# Create dendrogram and obtain species groupings for benthic assemblages 
# Katie Gale (katie.gale@dfo-mpo.gc.ca) _ Created 2015
# Formatted and uploaded 2017 
#
# As published in 
# Rubidge, E., Gale, K.S.P., Curtis, J.M.R., McClelland, E., Feyrer, L., Bodtker, K., and Robb, C.
# 2016. Methodology of the Pacific Marine Ecological Classification System and its Application
# to the Northern and Southern Shelf Bioregions. DFO Can. Sci. Advis. Sec. Res. Doc.
# 2016/035. xi + 124 p.
#
#
## Adapted code for Biological classification analysis - Eastern Canadian groundfish - Stanley, Heaslip, Jeffery, O'Brien
#########################################

pkgs <- list('vegan', 'simba', 'maptools', 'pvclust', 'dendroextras', 'dendextend', 'raster', 'reshape', 'reshape2', 'dplyr', 'tidyr', 'NbClust', 'rgdal','purrr')
invisible(lapply(pkgs, library, character.only = T))

##############################
# CLUSTER ANALYSIS STARTS
##############################
# start comparisons

#read in data if needed
SiteXSpecies<-read.csv( "Data/ClusterData4km_QC.csv", stringsAsFactors = F, row.names = 1)

#Compute dissimilarity matrices and cluster data using hierarchical clustering with average linkages
#Compare multiple measuers of dissimilarity

diss.measure <-c('simpson','soerensen','jaccard','gower')

dist.matrix <- map(diss.measure, ~sim(SiteXSpecies, method = .)) #compute dissimilarity matrix for each combination

dendro.obj <- map(dist.matrix, ~hclust(., method = 'average')) #create dendrogram object from each diss. matrix

#Calculate cophenetic correlation value for this tree
#strongest correlation using simpson coefficient (0.770). Use simpson.

copCor <- map(dendro.obj,cophenetic) %>% 
  map2(dist.matrix,., ~ paste(attributes(.x)$method, cor(.x,.y), sep = ' ')) %>% flatten_chr(.) %>% 
  print()

#compare fit of multiple linkage methods

lnkmtd <- c('ward.D', 'ward.D2','sin','com', 'ave','mcq','med','cen') # 8 common cluster algorithms

copCor.lnkmtd <- map(lnkmtd, ~ list(dm = dist.matrix[[1]], 
                                    do = hclust(dist.matrix[[1]], method = .))) %>% #for each combination, make list of [[1]] diss matrix, [[2]] dendrogram object using that linkage method   
  map_chr(., ~paste(.$do$method,cor(.$dm, cophenetic(.$do)), sep = ' ')) %>% print() #calculate and print cophenetic correlation for each combination
# average linkage method is best 

saveRDS(dendro.obj[[1]], file="Data/benthtree4km_QC.rds")

### Need to choose a cutoff - central to the whole analysis. 
## plots to help choose cutoff in 2aCutoffPlotsJOB_NFLD.R
## Once you've picked, move on

#
benthtree <- readRDS('Data/benthtree4km_QC.rds')
plot(benthtree, hang=-1) #dendrogram
rect.hclust(benthtree, h=0.499, border="red") #Cutoff as calculated by me
seth<-0.499 #use my choice - more informative biologically
cl<-dendroextras::slice(benthtree, h=seth) #Cut tree

###
colorcount<-as.data.frame(table(cl)) #get Table of cluster memberships
colorcount #Number of sites in each cluster. 3 clusters with 0.524 cut-off
colorcount <-colorcount[order(-colorcount$Freq),] #order by frequency
colorcount
setfreq<-min(colorcount$Freq[1:3], na.rm=T) #Select the first 3 clusters to be color-coded

#Count the number of sites in the top clusters
sum(colorcount$Freq[1:3]) # 1575 sites in top 4 clusters at 0.462 cut-off

# #Set Colors for top 3 clusters (SubBiomes)
colorscheme<-colorscheme<-c("#6a3d9a", #purple
                            "#1f78b4", #dblue
                            "#ff7f00") #orange
     
#Assign colours to clusters

colorcount$assigned[colorcount$Freq<setfreq]<-"#e0e0e0" #set all clusters that occur at less than the assigned cut off (3 clusters) to grey
colorcount$assigned[1:3] <- colorscheme
colorcount <- colorcount[order(colorcount$cl),]

colorcount #Each cluster is assigned a color
table(colorcount$assigned)

#Create color-coded dendrogram

colortree<-colour_branches(benthtree, h=seth,col=as.character(colorcount$assigned))

#create legend
colorcount2<-subset(colorcount, Freq>=setfreq) #Select top 3 clusters
colorcount2$cl<-as.character(colorcount2$cl)
colorcount2$Freq<-as.numeric(colorcount2$Freq)
colorcount2<-colorcount2[order(-colorcount2$Freq),] #order by frequency
colorcount2$Freq<-as.character(colorcount2$Freq) 
new<-c("others", paste("<", setfreq), "#e0e0e0") #add a column for the other clusters (less than top 3 - color them grey)
legendcluster<-rbind(colorcount2, new)
legendcluster$cl_name <- c("Deep Channels", "Shallow Banks & Straits","Channel Heads & Slopes", "Other") 
legendcluster #legend

#Print dendrogram
# #small - pub 
# tiff("Output/Dendrogram_QC4km_large.tiff", res=350, height=8, width=14, units="cm", pointsize=7)
#large - csas
 tiff("Output/Dendrogram_QC4km_large.tiff", res=300, height=2000, width=3800)

plot(colortree,xlab="", ylab="Simpson Distance", leaflab = "none")

 dev.off()


###plot this cluster analysis on map
#import the subsetted grid (grid cells for which we had data) that we wrote earlier
grid<-readOGR("Data/QCGrids/QC_4km_Grid.shp") 
head(grid@data)

#Join colors and clusters to map
colorcountmap<-colorcount
id <-names(cl)
colors<-as.data.frame(cbind(id, cl))
head(colors) #each grid cell is assigned its cluster
colorsmerge<-merge(grid, colors, by="id")
head(colorsmerge@data) #join this to the shapefile
colorsmerge2<-merge(colorsmerge, colorcountmap, by="cl")
head(colorsmerge2@data)#each grid cell is assigned its color to match the dendrogram
colorsgrid<-subset(colorsmerge2,!is.na(colorsmerge2$cl)) #drop grid cells for which we had data, but which got dropped at the rare species/barren sites stage
head(colorsgrid@data)
colorGridData<-colorsgrid@data

#write shapefile with cluster and color assignment attributes for 4 km subgrid
writeOGR(colorsgrid, dsn = 'Data/Shapefiles', layer = 'GridClusterAssignment4km_QC', driver = "ESRI Shapefile", overwrite_layer = T)

#Mapping shapefiles
QC_SA <- readOGR('Data/Shapefiles/QC_StudyArea.shp')
Land <- readOGR('Data/Shapefiles/QC_landborders.shp')

#adjust legend for map
legendmap<-legendcluster

#Plot map in projected coordinate system
tiff("Output/ColourGrid_QCkm.tiff", res=300, height=3000, width=2400)
par(mar = c(6,3,2,1), usr = c(2171.388,1047275,5043630,5840451), xpd = NA)
plot(QC_SA)
plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
plot(Land, col = 'gray97', border = 'gray10', add = T)
axis(1,at = c(250000,500000,750000), pos = 5043630, tick = T)
axis(2,at = c(5200000,5400000,5600000,5800000), pos = 2171.388, tick = T)
plot(extent(Land), add = T)
dev.off()


#write csv file for SiteXSpecies matrix with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(dplyr::select(SiteXSpecies_ClustAssigned, -layer,-Freq), "Data/Coloredgrid4km_cluster_assignment_QC.csv")


###Move  on to indicator species analysis


