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
library(nngeo)
library(purrr)

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
colorcount<-as.data.frame(table(cl)) %>%  #get Table of cluster memberships
  arrange(desc(Freq)) %>% #order by frequency
  print() #Number of sites in each cluster # 9 clusters with 0.498 cut-off

filter(colorcount, Freq >= 20) %>% nrow() # 4 major clusters (>= 20 sites)

#Count the number of sites in the top clusters

sum(colorcount$Freq[colorcount$Freq>= 20]) # 1316 sites in top 4 clusters at 0.498 cut-off

# Set Colors for top 4 clusters (SubBiomes)

source('Code/SPERA_colour_palettes.R')

#Assign colours to clusters

colorcount <- left_join(colorcount, GULF.palette) %>% 
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
                       breaks = GULF.palette$assigned, 
                       labels = GULF.palette$name) +
  guides(col = guide_legend(override.aes = list(size = 2))) +
  scale_y_continuous(name = 'Simpson Distance', limits = c(0,1), expand = c(0,0)) +
  labs(tag = 'SGSL') +
  theme(legend.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.justification = c(0,1),
        legend.direction = 'vertical',
        legend.key = element_rect(fill = 'white', colour = 'white'),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.line.x = element_line(colour = 'white'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.tag = element_text(face = 'bold', size = 16),
        plot.tag.position = c(0.8,0.92)); gg.colortree

saveRDS(gg.colortree, 'Output/ggDendrogram_Gulf.rds') 

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

# Read in other mapping shapefiles

Land <- st_read('Data/Shapefiles/Gulf_landborders.shp')
load("Data/GulfGridsWholeCell.RData") #load grid files
grid.4km <- Griddata[[20]] %>% st_as_sf() %>% st_transform(., st_crs(Land)) # 4 km resolution 

StudyArea <- grid.4km %>% st_union() %>% st_transform(crs = 26920)# aggregate grid cells

#convert colorsgrid to simple features object
colorsgrid.sf <- st_as_sf(colorsgrid) %>% # convert to sf
  st_transform(., crs = 26920) # match crs with study area

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
            bg.color = '#e0e0e0'); grid.map
  # tm_add_legend(type = 'fill', 
  #               labels = GULF.palette$name, 
  #               col = GULF.palette$assigned) +
  # tm_layout(legend.text.size = 0.78,
  #           legend.position = c(0.57,0.78),
  #           title = 'SGSL',
  #           title.position = c('right','top'),
  #           title.fontface = 'bold',
  #           title.size = 2,
  #           bg.color = '#e0e0e0',
  #           legend.width = 3)

saveRDS(grid.map, 'Output/ColorGrid_Gulf_tmap.rds')

tmap_save(grid.map, filename = 'Output/ColourGrid_Gulf_4km_legend.tiff', compression = 'lzw')

legend.map <- tm_shape(colorsgrid.sf2) +
  tm_polygons() +
  tm_add_legend(type = 'fill',
                labels = GULF.palette$name,
                col = GULF.palette$assigned,
                title = 'SGSL') +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = T,
            legend.width = 10,
            legend.height = 5,
            legend.title.fontface = 'bold'); legend.map

saveRDS(legend.map, 'Output/tmap_legend_Gulf.rds')

tmap_save(legend.map, filename = 'Output/tmap_legend_Gulf.tiff', width = 12, height = 5, unit = 'in', compression = 'lzw')

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
  coord_sf(label_graticule = "SW",
           expand = T,
           xlim = c(280798, 720798),
           ylim = c(5063640, 5447640)) +
  scale_x_continuous(name = '', 
                     breaks = c(-65,-63,-61),
                     expand = c(0.02,0)) +
  scale_y_continuous(name = '',
                     breaks = c(49,48,47,46),
                     expand = c(0.02,0)) +
  annotation_scale(bar_cols = c('grey40','grey10'),
                   text_cex = 1,
                   line_col = 'grey10',
                   location = 'bl',
                   width_hint = 0.2,
                   pad_y = unit(1.2, 'cm')) +
  annotate('text', 
           label = 'bold(SGSL)',
           parse = T,
           size = 6,
           x = 685000,
           y = 5420000) +
  theme(panel.background = element_rect(fill = '#e0e0e0'),
        panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
        panel.border = element_rect(fill = NA, colour = 'black'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        axis.text = element_text(size = 12.5),
        axis.text.y.left = element_text(hjust = 0.5, angle = 90))

saveRDS(ggmap, 'Output/ColorGrid_Gulf_ggmap.rds')
ggsave(ggmap,
       filename = 'Output/ColourGrid_Gulf_4km_legend.tiff', 
       compression = 'lzw',  dpi = 300,
       width = 6.7, height = 5.8, units = 'in')

# #Mapping shapefiles
# Gulf_SA <- readOGR('Data/Shapefiles/Gulf_RVsurveyAgg.shp')
# Land <- readOGR('Data/Shapefiles/Gulf_landborders.shp')
# 
# #adjust legend for map
# legendmap<-legendcluster
# 
# #Plot map in projected coordinate system
# tiff("Output/ColourGrid_Gulf4km.tiff", res=300, height=3000, width=2400)
# par(mar = c(5,5,8,1), usr = c(255798,746002.1,5033918,5472640), xpd = NA)
# plot(Gulf_SA)
# plot(colorsgrid, col = as.character(colorsgrid$assigned), border = NA,add=T)
# plot(Land, col = 'gray97', border = 'gray10', add = T)
# axis(1,at = c(300000,500000,700000), pos = 5033918, tick = T)
# axis(2,at = c(5100000,5250000,5400000), pos = 255798, tick = T)
# plot(extent(Land), add = T)
# legend(x = 590000, y = 5455000 , legend = legendmap$cl_name,
#        fill = as.character(legendmap$assigned),  bty = "n", cex = 0.85)
# dev.off()

#write csv file for SiteXSpecies matrix with cluster and color assignments
SiteXSpecies2 <- SiteXSpecies
SiteXSpecies2$id<-rownames(SiteXSpecies2)
SiteXSpecies_ClustAssigned<-merge(SiteXSpecies2, colorGridData, by="id")		
head(SiteXSpecies_ClustAssigned)		
write.csv(SiteXSpecies_ClustAssigned[,!(names(SiteXSpecies_ClustAssigned) %in% c('layer','Freq'))], "Data/Coloredgrid4km_cluster_assignment_Gulf.csv")


###Move  on to indicator species analysis


