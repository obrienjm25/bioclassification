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

pkgs <- list('vegan', 'simba', 'maptools', 'pvclust', 'dendroextras', 'dendextend', 'raster', 'reshape', 'reshape2', 'dplyr', 'tidyr', 'NbClust', 'rgdal')
invisible(lapply(pkgs, library, character.only = T))


##############################
# CLUSTER ANALYSIS STARTS
##############################
# start comparisons

#read in data if needed
SiteXSpecies<-read.csv( "Data/ClusterData4km_Nfld.csv", stringsAsFactors = F, row.names = 1)

#Compute dissimilarity matrices and cluster data using hierarchical clustering with average linkages
#Compare multiple measuers of dissimilarity
grid.pa.simp<-sim(SiteXSpecies,  method='simpson') #calculate B-sim (simpson) distance on the site x species matrix #if this is taking a long time - check that the matrix is good! if the grid names are in a column it won't run 
benthtree4km.simp<-hclust(grid.pa.simp, method='average') #cluster the simpson distance using average grouping

grid.pa.soer<-sim(SiteXSpecies, method = 'soerensen') #Soerensen coefficient
benthtree4km.soer <- hclust(grid.pa.soer, method = 'average')

grid.pa.jac<-sim(SiteXSpecies, method = 'jaccard') #Jaccard coefficient
benthtree4km.jac <- hclust(grid.pa.jac, method = 'average')

grid.pa.gow <- sim(SiteXSpecies, method = 'gower') #Gower coefficient
benthtree4km.gow <- hclust(grid.pa.gow, method = 'average')


#Calculate cophenetic correlation value for this tree
#strongest correlation using simpson coefficient. Use simpson.
copCor.simp<-cophenetic(benthtree4km.simp)
cor(grid.pa.simp,copCor.simp) #0.716

copCor.soer <- cophenetic(benthtree4km.soer)
cor(grid.pa.soer, copCor.soer) #0.308

copCor.jac <- cophenetic(benthtree4km.jac)
cor(grid.pa.jac, copCor.jac) #0.284

copCor.gow <-cophenetic(benthtree4km.gow) #0.461
cor(grid.pa.gow, copCor.gow)

#compare fit of multiple cluster algorithms

cl.alg <- c('ward.D', 'ward.D2','sin','com', 'ave','mcq','med','cen')

for (i in 1:6){
  cl.alg.temp <- cl.alg[i]
  benthtree.temp <- hclust(grid.pa.simp, method = cl.alg.temp)
  print(paste(cl.alg.temp, cor(grid.pa.simp, cophenetic(benthtree.temp)), sep = " "))
  rm(benthtree.temp)
  rm(cl.alg.temp)
} # average linkage method is best 

saveRDS(benthtree4km.simp, file="Data/benthtree4km_Nfld.rds")

### Need to choose a cutoff - central to the whole analysis. 
## plots to help choose cutoff in 2aCutoffPlotsJOB_NFLD.R
## Once you've picked, move on

#
benthtree <- readRDS('Data/benthtree4km.rds')
plot(benthtree, hang=-1) #dendrogram
rect.hclust(benthtree, h=0.494, border="red") #Cutoff as calculated by me
seth<-0.494 #use my choice - more informative biologically
cl<-dendroextras::slice(benthtree, h=seth) #Cut tree

###
colorcount<-as.data.frame(table(cl)) #get Table of cluster memberships
colorcount #Number of sites in each cluster # 12 clusters with 0.484 cut-off
colorcount <-colorcount[order(-colorcount$Freq),] #order by frequency
colorcount
setfreq<-min(colorcount$Freq[1:5], na.rm=T) #Select the first 5 clusters to be color-coded

#Count the number of sites in the top clusters
sum(colorcount$Freq[1:5]) # 3116 sites in top 5 clusters at 0.484 cut-off

# #Set Colors for top 5 clusters (SubBiomes)
colorscheme<-c("#a6cee3", #lblue
               "#1f78b4", #dblue
               "#ff7f00", #orange
               "#e6ab02", #yellow
               "#6a3d9a") #purple
                

#Assign colours to clusters

colorcount$assigned[colorcount$Freq<setfreq]<-"#e0e0e0" #set all clusters that occur at less than the assigned cut off (6 clusters) to grey
colorcount$assigned[1:5] <- colorscheme
colorcount <- colorcount[order(colorcount$cl),]

colorcount #Each cluster is assigned a color
table(colorcount$assigned)

#Create color-coded dendrogram
colortree<-colour_branches(benthtree, h=seth,col=as.character(colorcount$assigned))

#create legend
colorcount2<-subset(colorcount, Freq>=setfreq) #Select top 5 clusters
colorcount2$cl<-as.character(colorcount2$cl)
colorcount2$Freq<-as.numeric(colorcount2$Freq)
colorcount2<-colorcount2[order(-colorcount2$Freq),] #order by frequency
colorcount2$Freq<-as.character(colorcount2$Freq) 
new<-c("others", paste("<", setfreq), "#e0e0e0") #add a column for the other clusters (less than top 5 - color them grey)
legendcluster<-rbind(colorcount2, new)
legendcluster$cl_name <- c("Inner Shelf","Outer Shelf","Grand Banks","Slope", "Laurentian Channel/Shelf Break","Other") 
legendcluster #legend

#Print dendrogram
# #small - pub 
# tiff("Output/Dendrogram_Nfld4km_large.tiff", res=350, height=8, width=14, units="cm", pointsize=7)
#large - csas
 #tiff("Output/Dendrogram_Nfld4km_large.tiff", res=300, height=2000, width=3800)

plot(colortree,xlab="", ylab="Simpson Distance", leaflab = "none")

 #dev.off()


###plot this cluster analysis on map
#import the subsetted grid (grid cells for which we had data) that we wrote earlier
grid<-readOGR("Data/NfldGrids/NL_4km_Grid.shp") 
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
writeOGR(colorsgrid, dsn = 'Data/Shapefiles', layer = 'GridClusterAssignment4km_Nfld', driver = "ESRI Shapefile")

#Mapping shapefiles
NL_SA <- readOGR('Data/Shapefiles/NL_RVsurveyAgg.shp')
Land <- readOGR('Data/Shapefiles/NL_landborders.shp')

#adjust legend for map
legendmap<-legendcluster

#Plot map in projected coordinate system
tiff("Output/ColourGrid_Nfld4km.tiff", res=300, height=3000, width=2400)
par(mar = c(5,2,8,1), usr = c(67854.73,1500000,4700000,6800000), xpd = NA)
plot(NL_SA)
plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
plot(Land, col = 'gray97', border = 'gray10', add = T)
axis(1,at = c(250000,750000,1250000), pos = 4700000, tick = T)
axis(2,at = c(5000000,5500000,6000000,6500000), pos = 67854.73, tick = T)
plot(extent(67854.73,1500000,4700000,6800000), add = T)
legend(x = 750000, y = 6750000, legend = legendmap$cl_name,
       fill = as.character(legendmap$assigned),  bty = "n", cex = 1)
dev.off()


#write csv file for SiteXSpecies matrix with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(SiteXSpecies_ClustAssigned[,c(-74,-75)], "Data/Coloredgrid4km_cluster_assignment_Nfld.csv")


###Move  on to indicator species analysis


