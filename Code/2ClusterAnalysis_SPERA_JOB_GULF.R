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
SiteXSpecies<-read.csv( "Data/ClusterData4km_Gulf.csv", stringsAsFactors = F, row.names = 1)

#Compute dissimilarity matrices and cluster data using hierarchical clustering with average linkages
#Compare multiple measuers of dissimilarity

diss.measure <-c('simpson','soerensen','jaccard','gower')

dist.matrix <- map(diss.measure, ~sim(SiteXSpecies, method = .)) #compute dissimilarity matrix for each combination

dendro.obj <- map(dist.matrix, ~hclust(., method = 'average')) #create dendrogram object from each diss. matrix

#Calculate cophenetic correlation value for this tree
#strongest correlation using simpson coefficient (0.698). Use simpson.

copCor <- map(dendro.obj,cophenetic) %>% 
  map2(dist.matrix,., ~ paste(attributes(.x)$method, cor(.x,.y), sep = ' ')) %>% flatten_chr(.) %>% 
  print()

#compare fit of multiple linkage methods

lnkmtd <- c('ward.D', 'ward.D2','sin','com', 'ave','mcq','med','cen') # 8 common cluster algorithms

copCor.lnkmtd <- map(lnkmtd, ~ list(dm = dist.matrix[[1]], 
                                    do = hclust(dist.matrix[[1]], method = .))) %>% #for each combination, make list of [[1]] diss matrix, [[2]] dendrogram object using that linkage method   
  map_chr(., ~paste(.$do$method,cor(.$dm, cophenetic(.$do)), sep = ' ')) %>% print() #calculate and print cophenetic correlation for each combination
# average linkage method is best 

saveRDS(dendro.obj[[1]], file="Data/benthtree4km_Gulf.rds")

### Need to choose a cutoff - central to the whole analysis. 
## plots to help choose cutoff in 2aCutoffPlotsJOB_NFLD.R
## Once you've picked, move on

#
benthtree <- readRDS('Data/benthtree4km_Gulf.rds')
plot(benthtree, hang=-1) #dendrogram
rect.hclust(benthtree, h=0.498, border="red") #Cutoff as calculated by me
seth<-0.498 #use my choice - more informative biologically
cl<-dendroextras::slice(benthtree, h=seth) #Cut tree

###
colorcount<-as.data.frame(table(cl)) #get Table of cluster memberships
colorcount #Number of sites in each cluster. 9 clusters with 0.498 cut-off
colorcount <-colorcount[order(-colorcount$Freq),] #order by frequency
colorcount
setfreq<-min(colorcount$Freq[1:4], na.rm=T) #Select the top clusters to be color-coded

#Count the number of sites in the top clusters
sum(colorcount$Freq[1:4]) # 1316 sites in top 4 clusters at 0.498 cut-off

# #Set Colors for top 4 clusters (SubBiomes)

colorscheme<-c("#a6cee3", #lblue
                "#ff7f00", #orange
               "#6a3d9a", #purple
               "#1f78b4") #dblue

#a6cee3, lblue
#1f78b4, dblue
#e6ab02, yellow                

#Assign colours to clusters

colorcount$assigned[colorcount$Freq<setfreq]<-"#e0e0e0" #set all clusters that occur at less than the assigned cut off (4 clusters) to grey
colorcount$assigned[1:4] <- colorscheme
colorcount <- colorcount[order(colorcount$cl),]

colorcount #Each cluster is assigned a color
table(colorcount$assigned)

#Create color-coded dendrogram

colortree<-colour_branches(benthtree, h=seth,col=as.character(colorcount$assigned))

#create legend
colorcount2<-subset(colorcount, Freq>=setfreq) #Select top 4 clusters
colorcount2$cl<-as.character(colorcount2$cl)
colorcount2$Freq<-as.numeric(colorcount2$Freq)
colorcount2<-colorcount2[order(-colorcount2$Freq),] #order by frequency
colorcount2$Freq<-as.character(colorcount2$Freq) 
new<-c("others", paste("<", setfreq), "#e0e0e0") #add a column for the other clusters (less than top 4 - color them grey)
legendcluster<-rbind(colorcount2, new)
legendcluster$cl_name <- c("Magdalen Shallows" ,"Inshore/Magdalen Is.","Laurentian Channel","Northumberland Strait/\nSt. George's Bay","Other") 
legendcluster #legend

#Print dendrogram
# #small - pub 
# tiff("Output/Dendrogram_Gulf4km_large.tiff", res=350, height=8, width=14, units="cm", pointsize=7)
#large - csas
 tiff("Output/Dendrogram_Gulf4km_large.tiff", res=300, height=2000, width=3800)

plot(colortree,xlab="", ylab="Simpson Distance", leaflab = "none")

 dev.off()


###plot this cluster analysis on map
#import the subsetted grid (grid cells for which we had data) that we wrote earlier
grid<-readOGR("Data/GulfGrids/Gulf_4km_Grid.shp") 
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
writeOGR(colorsgrid, dsn = 'Data/Shapefiles', layer = 'GridClusterAssignment4km_Gulf', driver = "ESRI Shapefile", overwrite_layer = T)

#Mapping shapefiles
Gulf_SA <- readOGR('Data/Shapefiles/Gulf_RVsurveyAgg.shp')
Land <- readOGR('Data/Shapefiles/Gulf_landborders.shp')

#adjust legend for map
legendmap<-legendcluster

#Plot map in projected coordinate system
tiff("Output/ColourGrid_Gulf4km.tiff", res=300, height=3000, width=2400)
par(mar = c(5,5,8,1), usr = c(255798,746002.1,5033918,5472640), xpd = NA)
plot(Gulf_SA)
plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
plot(Land, col = 'gray97', border = 'gray10', add = T)
axis(1,at = c(300000,500000,700000), pos = 5033918, tick = T)
axis(2,at = c(5100000,5250000,5400000), pos = 255798, tick = T)
plot(extent(Land), add = T)
legend(x = 590000, y = 5455000 , legend = legendmap$cl_name,
       fill = as.character(legendmap$assigned),  bty = "n", cex = 0.85)
dev.off()


#write csv file for SiteXSpecies matrix with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(SiteXSpecies_ClustAssigned[,!(names(SiteXSpecies_ClustAssigned) %in% c('layer','Freq'))], "Data/Coloredgrid4km_cluster_assignment_Gulf.csv")


###Move  on to indicator species analysis


