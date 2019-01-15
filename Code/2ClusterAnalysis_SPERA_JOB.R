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
new<-c("others", paste("<", setfreq), "#e0e0e0") #add a column for the other clusters (less than top 6 - color them grey)
legendcluster<-rbind(colorcount2, new)
legendcluster$cl_name <- c("WSS/Outer BoF","WSS: Banks/Inner BoF","ESS","ESS: Banks","Laurentian Channel/Shelf Break","Slope","Other") 
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
#write shapefile with cluster and color assignment attributes for 4 km subgrid
writeOGR(colorsgrid, dsn = 'Data/Shapefiles', layer = 'GridClusterAssignment4km_Maritimes', driver = "ESRI Shapefile")

#Mapping shapefiles
MaritimeRegion <- readOGR("Data/MaritimesPlanningRegion/MaritimesPlanningArea.shp")
MaritimeProj <- spTransform(MaritimeRegion,CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
RVSurvey <- readOGR("Data/Shapefiles/MaritimesStudyArea.shp")
Land <- readOGR("Data/Shapefiles/LandBordersMaritimes.shp")

#adjust legend for map
legendmap<-legendcluster

#Plot map in projected coordinate system
tiff("Output/ColourGrid_Maritimes4km.tiff", res=300, height=3000, width=2400)
plot(RVSurvey, col = 'gray97')
plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
plot(Land, col = 'beige', add = T)
dev.off()

#plot map unprojected

library(raster)
colorsgrid.Proj <- spTransform(colorsgrid, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
RVsurvey.Proj <- spTransform(RVSurvey, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))                              
MaritimesProj2 <- spTransform(MaritimeRegion, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
Land.Proj <- spTransform(Land, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
Land.Proj <- crop(Land.Proj, extent(c(-68.25,-54.5,39.25,48.6)))
canada <- raster::getData("GADM", country = "CAN", level = 1)
QC.NFLD.PEI <- c("Newfoundland and Labrador", "Prince Edward Island", "Qu\u{e9}bec")
canada <- canada[canada$NAME_1 %in% QC.NFLD.PEI,]
canada <- crop(canada, extent(c(-68.25,-54.5,39.25,48.6)))

library(prettymapr)
par(mar = c(4,5,1,1), usr = c(-68.25,-54.5,39.25,48.6), xpd = NA)
plot(MaritimesProj2, border = NA)
plot(RVsurvey.Proj, add = T)
axis(1,at = c(-68,-64,-60,-56), pos = 39.25, tick = T,
     labels = c(expression(paste("68",degree," W",sep='')),expression(paste("64",degree," W",sep='')),expression(paste("60",degree," W",sep='')),expression(paste("56",degree," W",sep=''))))
axis(2,at = c(40,42,44,46,48), pos = -68.25, tick = T,
     labels = c(expression(paste("40",degree," N",sep='')),expression(paste("42",degree," N",sep='')),expression(paste("44",degree," N",sep='')),expression(paste("46",degree," N",sep='')),expression(paste("48",degree," N",sep=''))))
plot(colorsgrid.Proj, col = as.character(colorsgrid.Proj$assigned), border = NA,add=T)
plot(canada, col = 'gray97', border = 'gray10', add = T)
plot(Land.Proj, col = 'gray97', border = 'gray10', add = T)
plot(extent(-68.25,-54.5,39.25,48.6), add = T)
addnortharrow(pos = 'bottomleft', scale = 0.35)
addscalebar(pos = 'bottomright', htin =0.1,widthhint=0.3)
legend(x = -61.8, y = 42.6, legend = legendmap$cl_name,
       fill = as.character(legendmap$assigned),  bty = "n", cex = 0.85)

#write csv file for SiteXSpecies matriz with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(SiteXSpecies_ClustAssigned[,c(-117,-118)], "Data/Coloredgrid4km_cluster_assignment_Maritimes.csv")


###Move  on to indicator species analysis


