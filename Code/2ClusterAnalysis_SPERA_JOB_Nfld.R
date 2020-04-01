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
library(purrr)

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

for (i in 1:8){
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
benthtree <- readRDS('Data/benthtree4km_Nfld.rds')
plot(benthtree, hang=-1) #dendrogram
rect.hclust(benthtree, h=0.512, border="red") #Cutoff as calculated by me
seth<-0.512 #use my choice - more informative biologically
cl<-dendroextras::slice(benthtree, h=seth) #Cut tree

###
colorcount<-as.data.frame(table(cl)) %>%  #get Table of cluster memberships
  arrange(desc(Freq)) %>% #order by frequency
  print() #Number of sites in each cluster # 8 clusters with 0.512

filter(colorcount, Freq >= 20) %>% nrow() # 5 major clusters (>= 20 sites)

#Count the number of sites in the top clusters

sum(colorcount$Freq[colorcount$Freq>= 20]) # 4897 sites in top 5 clusters at 0.494 cut-off

# Set Colors for top 5 clusters (SubBiomes)

source('Code/SPERA_colour_palettes.R')

#Assign colours to clusters

colorcount <- left_join(colorcount, NL.palette) %>% 
  replace_na(list(assigned = "#e0e0e0", name = "Unclassified")) %>% #set all clusters that occur < cut-off to grey
  arrange(as.numeric(cl)) %>% #order by cluster ID
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
                       breaks = NL.palette$assigned, 
                       labels = NL.palette$name) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  scale_y_continuous(name = 'Simpson Distance', limits = c(0,1), expand = c(0,0)) +
  labs(tag = 'NL') +
  theme(legend.text = element_text(size = 12),
        legend.position = 'top',
        legend.justification = c(0,0),
        legend.direction = 'vertical',
        legend.key = element_rect(fill = 'white', colour = 'white'),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(colour = 'white'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.tag = element_text(face = 'bold', size = 16),
        plot.tag.position = c(0.8,0.5)); gg.colortree

saveRDS(gg.colortree, 'Output/ggDendrogram_NL.rds') 

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
writeOGR(colorsgrid, dsn = 'Data/Shapefiles', layer = 'GridClusterAssignment4km_Nfld', driver = "ESRI Shapefile", overwrite_layer = T)

# Read in other mapping shapefiles

Land <- st_read('Data/Shapefiles/NL_landborders.shp')
LandBuffer <- st_read('Data/Shapefiles/NL_LandBuffer_5km.shp')
load("Data/NLGridsWholeCell.RData") #load grid files
grid.4km <- Griddata[[20]] %>% st_as_sf() %>% st_transform(., st_crs(LandBuffer)) # 4 km resolution 

sel.touches.buffer <- st_disjoint(grid.4km, LandBuffer, sparse = F)[,1] # identify grid cells that intersect land buffer

StudyArea <- dplyr::filter(grid.4km, sel.touches.buffer == T) %>% # filter out grid cells intersecting land buffer
  st_union() %>% # aggregate polygons
  st_transform(., crs = 26921)
  
#convert colorsgrid to simple features object
colorsgrid.sf <- st_as_sf(colorsgrid) %>% # convert to sf
  st_cast("MULTIPOLYGON") %>% #cast to multipolygon geometry
  st_transform(., crs= 26921)

sel_cells_SA <- st_within(colorsgrid.sf, StudyArea, sparse = F)[,1] # identify grid cells within study area

colorsgrid.sf2 <- filter(colorsgrid.sf, sel_cells_SA == T)  

# Plot map of coloured grid cells showing distribution of cluster assignments

grid.map <- tm_shape(colorsgrid.sf2,
                     bbox = st_bbox(StudyArea)) +
  tm_fill(col = 'assigned') +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') + 
  tm_shape(Land) +
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
  # tm_add_legend(type = 'fill',
  #               labels = NL.palette$name,
  #               col = NL.palette$assigned) +
  # tm_layout(legend.text.size = 0.85,
  #           legend.position = c(0.42,0.81),
  #           title = 'NL',
  #           title.position = c('right','top'),
  #           title.fontface = 'bold',
  #           title.size = 2,
  #           bg.color = '#e0e0e0',
  #           legend.width = 3)
  
saveRDS(grid.map, 'Output/ColorGrid_NL_tmap.rds')

tmap_save(grid.map, filename = 'Output/ColourGrid_Nfld_4km_legend.tiff', compression = 'lzw')

legend.map <- tm_shape(colorsgrid.sf2) +
  tm_polygons() +
  tm_add_legend(type = 'fill',
                labels = NL.palette$name,
                col = NL.palette$assigned,
                title = 'NL') +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = T,
            legend.width = 10,
            legend.height = 5,
            legend.title.fontface = 'bold'); legend.map

saveRDS(legend.map, 'Output/tmap_legend_Nfld.rds')

tmap_save(legend.map, filename = 'Output/tmap_legend_Nfld.tiff', width = 10, height = 5, unit = 'in', compression = 'lzw')

# ggplot version

ggmap <- ggplot() +
  geom_sf(data = colorsgrid.sf2, 
          aes(fill = assigned),
          col = NA,
          show.legend = F,
          na.rm = T) + 
  scale_fill_identity() +
  geom_sf(data = StudyArea, 
          fill = NA, 
          col = 'black') +
  geom_sf(data = Land, 
          fill = 'grey30', 
          col = 'grey20', 
          lwd = 0.25) +
  coord_sf(label_graticule = "NE",
           expand = T,
           xlim = c(270169.8, 1314169.8),
           ylim = c(4742510.1, 6398510.1)) +
  scale_x_continuous(name = '', 
                     breaks = c(-48,-54,-60),
                     expand = c(0.02,0)) +
  scale_y_continuous(name = '',
                     breaks = c(44,48,52,56),
                     expand = c(0.02,0)) +
  annotation_scale(bar_cols = c('grey40','grey10'),
                   text_cex = 1,
                   line_col = 'grey10',
                   location = 'bl') +
  annotate('text', 
           label = 'bold(NL)',
           parse = T,
           size = 6,
           x = 1200000,
           y = 6300000) +
  theme(panel.background = element_rect(fill = '#e0e0e0'),
        panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
        panel.border = element_rect(fill = NA, colour = 'black'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        axis.text = element_text(size = 12.5),
        axis.text.y.right = element_text(hjust = 0.5, angle = 270))

saveRDS(ggmap, 'Output/ColorGrid_NL_ggmap.rds')
ggsave(ggmap, 
       filename = 'Output/ColourGrid_Nfld_4km_legend.tiff', 
       compression = 'lzw', dpi = 300,
       width = 5.6, height = 8.8, unit = 'in')

# #Mapping shapefiles
# NL_SA <- readOGR('Data/Shapefiles/NL_RVsurveyAgg.shp')
# Land <- readOGR('Data/Shapefiles/NL_landborders.shp')
# 
# #adjust legend for map
# legendmap<-legendcluster
# 
# #Plot map in projected coordinate system
# tiff("Output/ColourGrid_Nfld4km.tiff", res=300, height=3000, width=2400)
# par(mar = c(5,2,8,1), usr = c(67854.73,1500000,4700000,6800000), xpd = NA)
# plot(NL_SA)
# plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
# plot(Land, col = 'gray97', border = 'gray10', add = T)
# axis(1,at = c(250000,750000,1250000), pos = 4700000, tick = T)
# axis(2,at = c(5000000,5500000,6000000,6500000), pos = 67854.73, tick = T)
# plot(extent(67854.73,1500000,4700000,6800000), add = T)
# legend(x = 750000, y = 6750000, legend = legendmap$cl_name,
#        fill = as.character(legendmap$assigned),  bty = "n", cex = 1)
# dev.off()


#write csv file for SiteXSpecies matrix with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(SiteXSpecies_ClustAssigned[,c(-74,-75)], "Data/Coloredgrid4km_cluster_assignment_Nfld.csv")


###Move  on to indicator species analysis


