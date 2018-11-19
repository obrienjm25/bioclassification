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

pkgs <- list('vegan', 'simba', 'maptools', 'pvclust', 'dendroextras', 'dendextend', 'reshape', 'reshape2', 'dplyr', 'tidyr', 'NbClust', 'rgdal')
invisible(lapply(pkgs, library, character.only = T))


##############################
# CLUSTER ANALYSIS STARTS
##############################
# start comparisons

#read in data if needed
SiteXSpecies<-read.csv( "Data/ClusterData4km.csv", stringsAsFactors = F, row.names = 1)

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
cor(grid.pa.simp,copCor.simp) #0.606

copCor.soer <- cophenetic(benthtree4km.soer)
cor(grid.pa.soer, copCor.soer) #0.349

copCor.jac <- cophenetic(benthtree4km.jac)
cor(grid.pa.jac, copCor.jac) #0.333

copCor.gow <-cophenetic(benthtree4km.gow) #0.549
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

saveRDS(benthtree4km.simp, file="Data/benthtree4km.rds")

### Need to choose a cutoff - central to the whole analysis. 
## plots to help choose cutoff in 2aCutoffPlotsJOB.R
## Once you've picked, move on

#
benthtree <- readRDS('Data/benthtree4km.rds')
plot(benthtree, hang=-1) #dendrogram
rect.hclust(benthtree, h=0.585, border="red") #Cutoff as calculated by me
seth<-0.585 #use my choice - more informative biologically
cl<-dendroextras::slice(benthtree, h=seth) #Cut tree

###
colorcount<-as.data.frame(table(cl)) #get Table of cluster memberships
colorcount #Number of sites in each cluster # 17 clusters with 0.58 cut-off
colorcount <-colorcount[order(-colorcount$Freq),] #order by frequency
colorcount
setfreq<-min(colorcount$Freq[1:6], na.rm=T) #Select the first 6 clusters to be color-coded

#Count the number of sites in the top clusters
sum(colorcount$Freq[1:6]) # 2494 sites in top 6 clusters at 0.585 cut-off

# #Set Colors for top 6 clusters (SubBiomes)
colorscheme<-c("#a6cee3", #lblue
               "#1f78b4",#dblue
               "#ff7f00", #orange
               "#33a02c", #green
               "#6a3d9a", #purple
               "#e6ab02") #yellow

#Assign colours to clusters

colorcount$assigned[colorcount$Freq<setfreq]<-"#e0e0e0" #set all clusters that occur at less than the assigned cut off (6 clusters) to grey
colorcount$assigned[1:6] <- colorscheme
colorcount <- colorcount[order(colorcount$cl),]

colorcount #Each cluster is assigned a color
table(colorcount$assigned)

#Create color-coded dendrogram
colortree<-colour_branches(benthtree, h=seth,col=as.character(colorcount$assigned))

#create legend
colorcount2<-subset(colorcount, Freq>=setfreq) #Select top 6 clusters
colorcount2$cl<-as.character(colorcount2$cl)
colorcount2$Freq<-as.numeric(colorcount2$Freq)
colorcount2<-colorcount2[order(-colorcount2$Freq),] #order by frequency
colorcount2$Freq<-as.character(colorcount2$Freq) 
new<-c("others", paste("<", setfreq), "grey") #add a column for the other clusters (less than top 6 - color them grey)

legendcluster<-rbind(colorcount2, new)
legendcluster #legend

#Print dendrogram
# #small - pub 
# tiff("Output/Dendrogram_Maritimes4km_large.tiff", res=350, height=8, width=14, units="cm", pointsize=7)
#large - csas
# tiff("Output/Dendrogram_Maritimes4km_large.tiff", res=300, height=2000, width=3800)

plot(colortree,xlab="", ylab="Simpson Distance", leaflab = "none")

# dev.off()


###plot this cluster analysis on map
#import the subsetted grid (grid cells for which we had data) that we wrote earlier
grid<-readShapeSpatial("Data/MaritimeGrids/Maritimes_4km_Grid.shp") 
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

#set coordinate system and project
proj4string(colorsgrid)<-CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system
MaritimeRegion <- readOGR("Data/MaritimesPlanningRegion/MaritimesPlanningArea.shp")
MaritimeProj <- spTransform(MaritimeRegion,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
RVSurvey <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")
Land <- readOGR("Data/Shapefiles/LandBordersMaritimes.shp")

#adjust legend for map
legendmap<-legendcluster

#Plot map
tiff("Output/ColourGrid_Maritimes4km.tiff", res=300, height=3000, width=2400)
plot(RVSurvey, col = 'gray97')
plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
plot(Land, col = 'beige', add = T)
dev.off()

#write
writeSpatialShape(colorsgrid, "Combined_2015.04.30/3.Cluster_04.30_GfBioTanCrabClam_RemCoast_NSBWCVI/0.55/Coloredgrid_0.55_NSB_WCVI_2015_05_07.shp")

#Get lat long for each Grid record
head(bdata1)		
bdata1$GridID<-rownames(bdata1)
bdata2<-merge(bdata1, colorGridData, by="GridID")		
head(bdata2)		
write.csv(bdata2,"Combined_2015.04.30/3.Cluster_04.30_GfBioTanCrabClam_RemCoast_NSBWCVI/pa_matrix_3615sitesx174sp_withLatLongx.csv" )		
write.csv(bdata2[names(bdata2) %in% c("cl","SP_ID","Id","MidLat"  ,"MidLong" ,"Freq" ,"assigned"), "Combined_2015.04.30/3.Cluster_04.30_GfBioTanCrabClam_RemCoast_NSBWCVI/0.65/Coloredgrid_cluster_assignment_2015.05.07.csv"


###Move  on to indicator species analysis


