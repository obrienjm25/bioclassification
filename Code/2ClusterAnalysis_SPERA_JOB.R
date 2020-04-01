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

library(vegan)
library(simba)
library(maptools)
library(pvclust)
library(dendroextras)
library(dendextend)
library(dendroextras)
library(reshape)
library(reshape2)
library(dplyr)
library(tidyr)
library(rgdal)
library(NbClust)
library(ggplot2)
library(ggdendro)
library(sf)
library(tmap)
library(raster)
library(ggplot2)
library(ggsn)

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
cor(grid.pa.simp,copCor.simp) #0.585

copCor.soer <- cophenetic(benthtree4km.soer)
cor(grid.pa.soer, copCor.soer) #0.335

copCor.jac <- cophenetic(benthtree4km.jac)
cor(grid.pa.jac, copCor.jac) #0.316

copCor.gow <-cophenetic(benthtree4km.gow) #0.554
cor(grid.pa.gow, copCor.gow)

#compare fit of multiple cluster algorithms

cl.alg <- c('ward.D2','sin','com', 'ave','mcq','med','cen')

for (i in 1:7){
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
rect.hclust(benthtree, h=0.587, border="red") #Cutoff as calculated by me
seth<-0.587 #use my choice - more informative biologically
cl<-dendroextras::slice(benthtree, h=seth) #Cut tree

###
colorcount<-as.data.frame(table(cl)) %>%  #get Table of cluster memberships
  arrange(desc(Freq)) %>% #order by frequency
  print() #Number of sites in each cluster # 9 clusters with 0.587 cut-off

filter(colorcount, Freq >= 20) %>% nrow() # 6 major clusters (>= 20 sites)

#Count the number of sites in the top clusters

sum(colorcount$Freq[colorcount$Freq>= 20]) # 2520 sites in top 6 clusters at 0.587 cut-off

# Set Colors for top 6 clusters (SubBiomes)

source('Code/SPERA_colour_palettes.R')

#Assign colours to clusters

colorcount <- left_join(colorcount, MAR.palette) %>% 
  replace_na(list(assigned = "#e0e0e0", name = "Unclassified")) %>% #set all clusters that occur < cut-off to grey
  arrange(cl) %>% #order by cluster ID
  print() #Each cluster is assigned a color

#Create color-coded dendrogram
colortree<-colour_branches(benthtree, h=seth,col=colorcount$assigned) %>%   
  as.ggdend() %>% 
  segment() %>% 
  replace_na(list(col = '#000000')) %>% 
  left_join(., colorcount, by = c('col' = 'assigned')) %>% 
  mutate(lwd = if_else(col == '#000000', 0.75, 0.5))

gg.colortree <- ggplot(colortree) + 
  theme_classic() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, colour = col, size = lwd),
               lineend = 'square',
               linejoin = 'bevel') +
  scale_size_identity() +
  scale_color_identity('', 
                       guide = 'legend', 
                       breaks = MAR.palette$assigned, 
                       labels = MAR.palette$name) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  scale_y_continuous(name = 'Simpson Distance', limits = c(0,1), expand = c(0,0)) +
  labs(tag = 'MAR') +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(fill = 'white', colour = 'white'),
        legend.position = 'bottom',
        legend.justification = c(0,1),
        legend.direction = 'vertical',
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(colour = 'white'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.tag = element_text(face = 'bold', size = 16),
        plot.tag.position = c(0.8,0.92)); gg.colortree
 
saveRDS(gg.colortree, 'Output/ggDendrogram_MAR.rds') 

###plot this cluster analysis on map
#import the subsetted grid (grid cells for which we had data) that we wrote earlier
grid<-readShapeSpatial("Data/MaritimeGrids/Maritimes_4km_Grid.shp") 
head(grid@data)

#Join colors and clusters to map
colorcountmap<-colorcount
id <-names(cl)
colors<-as.data.frame(cbind(id, cl))
head(colors) #each grid cell is assigned its cluster
colorsmerge<-merge(grid, colors, by="id") #join this to the shapefile
head(colorsmerge@data) 
colorsmerge2<-merge(colorsmerge, colorcountmap, by="cl") # merge with assigned cluster colours and names
head(colorsmerge2@data)#each grid cell is assigned its color to match the dendrogram
colorsgrid<-subset(colorsmerge2,!is.na(colorsmerge2$cl)) #drop grid cells for which we had data, but which got dropped at the rare species/barren sites stage
head(colorsgrid@data)
colorGridData<-colorsgrid@data

#set coordinate system and project
crs(colorsgrid)<-CRS("+proj=utm +zone=20 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") #set coordinate reference system
#write shapefile with cluster and color assignment attributes for 4 km subgrid
writeOGR(colorsgrid, dsn = 'Data/Shapefiles', layer = 'GridClusterAssignment4km_Maritimes', driver = "ESRI Shapefile", overwrite_layer = T)

# Read in other mapping shapefiles

StudyArea <- st_read('Data/Shapefiles/MaritimesRasterBoundary.shp')
landfiles <- list.files(pattern = 'gadm', full.names = T) %>% 
  discard(~ grepl('FRA', .x)) %>% 
  discard(~ grepl('_0_', .x))
admin.borders <- map(landfiles, readRDS) %>% 
  map(., st_as_sf) %>% 
  rbindlist() %>% 
  st_as_sf() %>% 
  st_crop(., c(xmin = -71, ymin = 39, xmax = -45, ymax = 59)) %>% 
  st_transform(., crs = 26920)

#convert colorsgrid to simple features object
colorsgrid.sf <- st_as_sf(colorsgrid) %>% # convert to sf
  st_cast("MULTIPOLYGON") %>% #cast to multipolygon geometry
  st_transform(., crs = 26920) %>% 
  st_intersection(StudyArea, .) # remove grid cells outside study boundary

# Plot map of coloured grid cells showing distribution of cluster assignments

grid.map <- tm_shape(colorsgrid.sf, 
                     bbox = st_bbox(c(xmin =20000, 
                                      ymin = 4562134, 
                                      xmax = 1029435, 
                                      ymax = 5300000))) + 
  tm_fill(col = 'assigned') + 
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') + 
  tm_shape(admin.borders) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  # tm_compass(type = '4star',
  #            text.color = 'grey10',
  #            color.dark = 'grey10',
  #            color.light = 'grey20',
  #            position = c(0.08,0.8),
  #            size = 1,
  #            text.size = 1.2) +
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
            bg.color = '#e0e0e0',
            inner.margins = c(0,0,0,0)); grid.map
  # tm_add_legend(type = 'fill', 
  #               labels = MAR.palette$name, 
  #               col = MAR.palette$assigned) +
  # tm_layout(legend.text.size = 0.95,
  #           legend.position = c(0.01, 0.61),
  #           title = 'MAR',
  #           title.position = c('right','top'),
  #           title.fontface = 'bold',
  #           title.size = 2,
  #           bg.color = '#e0e0e0',
  #           legend.width = 3)

saveRDS(grid.map, 'Output/ColorGrid_MAR_tmap.rds')

tmap_save(grid.map, filename = 'Output/ColourGrid_Maritimes4km_legend.tiff', compression = 'lzw')

legend.map <- tm_shape(colorsgrid.sf) +
  tm_polygons() +
  tm_add_legend(type = 'fill',
                labels = MAR.palette$name,
                col = MAR.palette$assigned,
                title = 'MAR') +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = T,
            legend.width = 10,
            legend.height = 5,
            legend.title.fontface = 'bold'); legend.map

saveRDS(legend.map, 'Output/tmap_legend_MAR.rds')

tmap_save(legend.map, filename = 'Output/tmap_legend_MAR.tiff', width = 10, height = 5, unit = 'in', compression = 'lzw')

#ggplot version

ggmap <- ggplot() +
  geom_sf(data = colorsgrid.sf, 
          aes(fill = assigned),
          col = NA,
          show.legend = F,
          na.rm = T) + 
  scale_fill_identity() +
  geom_sf(data = StudyArea,
          fill = NA, 
          col = 'black') +
  geom_sf(data = admin.borders, 
          fill = 'grey30',
          col = 'grey20', 
          lwd = 0.25) +
  coord_sf(label_graticule = "SE",
           expand = F,
           xlim = c(0, 1040000),
           ylim = c(4530000, 5300000)) +
  scale_x_continuous(name = '', breaks = c(-67,-63,-59)) +
  scale_y_continuous(name = '',breaks = c(47,45,43,41)) +
  annotate('text', 
           label = 'bold(MAR)',
           parse = T,
           size = 6,
           x = 920000,
           y = 5220000) +
  annotation_scale(bar_cols = c('grey40','grey10'),
                   text_cex = 1,
                   line_col = 'grey10',
                   location = 'tl',
                   pad_y = unit(0.5, 'cm')) +
  theme(panel.background = element_rect(fill = '#e0e0e0'),
        panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        panel.border = element_rect(fill = NA, colour = 'black'),
        axis.text = element_text(size = 12.5),
        axis.text.y.right = element_text(hjust = 0.5, angle = 270))

saveRDS(ggmap, 'Output/ColorGrid_MAR_ggmap.rds')
ggsave(ggmap, 
       filename = 'Output/ColourGrid_Maritimes4km_legend.tiff',
       compression = 'lzw', dpi = 300,
       width = 8, height = 5.8, unit = 'in')

#write csv file for SiteXSpecies matrix with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(SiteXSpecies_ClustAssigned[,c(-116,-117,-118)], "Data/Coloredgrid4km_cluster_assignment_Maritimes.csv")


###Move  on to indicator species analysis


