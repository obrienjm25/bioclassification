####Figure Prep for SPERA Bioclassification####
#---------------------------------------------#

# load required packages

library(tmap)
library(grid)
library(gridExtra)
library(purrr)
library(ggplot2)
library(sf)
library(raster)
library(rgdal)
library(cowplot)
library(data.table)
library(dplyr)
library(tidyr)
library(wesanderson)
library(ggspatial)
library(ggsn)

# load colour palette

source('Code/SPERA_colour_palettes.R')

#### Dendrograms ####

# Read in dendrograms for 4 regions
dendro.ls <- list.files(path = 'Output/', pattern = 'ggDendro', full.names = T) %>% # file list
  map(readRDS) %>% # Read files in list
  map(., ~ . + theme(plot.tag = element_text(face = 'bold', size = 12), # format plot tag font
                     panel.border = element_rect(fill = NA, size = 1), # add panel borders
                     panel.background = element_rect(fill = 'grey85'), # add light background
                     axis.line.x = element_line(colour = 'black'))) %>% # add x-axis
  set_names(c('SGSL','MAR','NL','NGSL')) %>% # name elements in list
  .[c('NGSL','NL','SGSL','MAR')] # re-order list

# Remove ylabs for panels in 2nd column
dendro.ls$NL <- dendro.ls$NL + theme(axis.title.y = element_blank(),
                                     axis.text.y = element_blank())
dendro.ls$MAR <- dendro.ls$MAR + theme(axis.title.y = element_blank(),
                                       axis.text.y = element_blank())

# convert gg objects to grobs

dendro.ls <- map(dendro.ls, ggplotGrob) 

# Create first 2-panel row
g.r1 <- cbind(dendro.ls$NGSL, dendro.ls$NL, size = 'last') # bind panels
# g.c1$heights <- unit.pmax(dendro.ls$NGSL$heights, dendro.ls$SGSL$heights) # ensure consisent widths of plotting elements

# Create 2nd 2-panel row
g.r2 <- cbind(dendro.ls$SGSL, dendro.ls$MAR, size = 'last')
# g.c2$heights <- unit.pmax(dendro.ls$NL$heigths, dendro.ls$MAR$heigths)

# Bind columns
dendro.4panel <- rbind(g.r1, g.r2, size = 'first')

# Save 4 panel figure as tiff with lzw compression
ggsave(dendro.4panel, filename = 'Output/Dendrograms_4regions.tiff', width = 7, height = 10, compression = 'lzw')

#### Map with predicted distribution of major clusters ####

# Maritimes--------------------

# Load spatial layers

RF.pred <- raster('Output/Maritimes_PredClust_map.tif') # predicted distribution of 6 major clusters
RF.pred.sf <- as(RF.pred, 'SpatialPolygonsDataFrame') %>% # convert to sf object
  st_as_sf() %>% 
  dplyr::rename(cl = Maritimes_PredClust_map) %>% 
  mutate(., assigned = case_when(
    cl == 1 ~ '#5A283E',
    cl == 2 ~ '#F2300F',
    cl == 3 ~ '#649373',
    cl == 4 ~ '#1B5656',
    cl == 5 ~ '#9FA682',
    cl == 6 ~ '#E1BD6D'))
StudyArea <- st_read('Data/Shapefiles/MaritimesRasterBoundary.shp') #study area boundaries
landfiles <- list.files(pattern = 'gadm', full.names = T) %>% 
  discard(~ grepl('FRA', .x)) %>% 
  discard(~ grepl('_0_', .x))
admin.borders <- map(landfiles, readRDS) %>% 
  map(., st_as_sf) %>% 
  rbindlist() %>% 
  st_as_sf() %>% 
  st_crop(., c(xmin = -71, ymin = 39, xmax = -45, ymax = 59)) %>% 
  st_transform(., crs = 26920)
uncertainty <- st_read('Output/Maritimes_RF_uncertainty_0.7.shp') # grid cells with probability of assignment < 0.70

#ggplot version

ggmap <- ggplot() +
  geom_sf(data = RF.pred.sf, 
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
  geom_sf(data = uncertainty, 
          col = 'black', 
          fill = NA,
          show.legend = F) +
  coord_sf(label_graticule = "SE",
           expand = F,
           xlim = c(0, 1040000),
           ylim = c(4530000, 5300000)) +
  scale_x_continuous(name = '', breaks = c(-67,-63,-59)) +
  scale_y_continuous(name = '',breaks = c(47,45,43,41)) +
  annotation_scale(bar_cols = c('grey40','grey10'),
                   text_cex = 1,
                   line_col = 'grey10',
                   location = 'tl',
                   pad_y = unit(0.5, 'cm')) +
  annotate('text', 
           label = 'bold(MAR)',
           parse = T,
           size = 6,
           x = 920000,
           y = 5220000) +
  theme(panel.background = element_rect(fill = '#e0e0e0'),
        panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
        panel.border = element_rect(fill = NA, colour = 'black'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        axis.text = element_text(size = 12.5),
        axis.text.y.right = element_text(hjust = 0.5, angle = 270))

saveRDS(ggmap, 'Output/Maritimes_RF_predictions_ggmap_uncertainty_0.7.rds')
ggsave(ggmap, 
       filename = 'Output/Maritimes_RF_predictions_ggmap_uncertainty_0.7.tiff',
       compression = 'lzw', dpi = 300,
       width = 8, height = 5.8, unit = 'in')

# tmap version

map <-  tm_shape(RF.pred, bbox = st_bbox(c(xmin =20000, 
                                           ymin = 4562134, 
                                           xmax = 1029435, 
                                           ymax = 5300000))) +
  tm_raster(legend.show = F,
            style = 'cat',
            palette = c(NA, MAR.palette$assigned[c(6,5,3,4,1,2)])) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +    
  tm_shape(admin.borders) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_shape(uncertainty) + 
  tm_borders(col = 'black') +
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
            bg.color = '#e0e0e0'); map

saveRDS(map, 'Output/MaritimesRF_predictions_tmap_uncertainty_0.7.rds')

tmap_save(map, filename = 'Output/MaritimesRF_predictions_tmap_uncertainty_0.7.tiff', compression = 'lzw')

# Legend 

legend <- tm_shape(uncertainty) +
  tm_borders(col = 'black') +
  tm_add_legend(type = 'fill',
                col = 'white',
                labels = 'Prob. of Assignment < 0.70') +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = T,
            legend.width = 10,
            legend.height = 1,
            legend.title.fontface = 'bold'); legend

saveRDS(legend, 'Output/tmap_legend_RF_prediction.rds')

tmap_save(legend, filename = 'Output/tmap_legend_RF_prediction.tiff', width = 10, height = 1, unit = 'in', compression = 'lzw')

#### Map climate projections - Maritimes

# Load additional spatial layers

grid.anom <- raster('Output/Maritimes_GridAnomalies_RCP85_2075.tif') # 2075 anomalies (i.e. classification changed)
climate.sensitive <- st_read('Output/Mar_ClimateSensitive.shp') # Border of climate susceptible areas
RF.2075.RCP8.5 <- raster('Output/Maritimes_Forecast2075_map.tif') # Predictions for 2075 under RCP 8.5
RF.2075.RCP8.5.sf <- as(RF.2075.RCP8.5, 'SpatialPolygonsDataFrame') %>% # convert to sf object
  st_as_sf() %>% 
  dplyr::rename(cl = Maritimes_Forecast2075_map) %>% 
  mutate(., assigned = case_when(
    cl == 1 ~ '#5A283E',
    cl == 2 ~ '#F2300F',
    cl == 3 ~ '#649373',
    cl == 4 ~ '#1B5656',
    cl == 5 ~ '#9FA682',
    cl == 6 ~ '#E1BD6D'))

# Map forecasted distribution and highlight climate sensitive areas

#ggplot version

ggmap <- ggplot() +
  geom_sf(data = RF.2075.RCP8.5.sf, 
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
  geom_sf(data = climate.sensitive, 
          aes(col = 'blue'), 
          fill = 'black',
          alpha = 0.15) +
  scale_color_identity(guide = 'legend',
                       name = '',
                       labels = c('Climate susceptible (anomalies)')) +
  coord_sf(label_graticule = "SE",
           expand = F,
           xlim = c(0, 1040000),
           ylim = c(4530000, 5300000)) +
  scale_x_continuous(name = '', breaks = c(-67,-63,-59)) +
  scale_y_continuous(name = '',breaks = c(47,45,43,41)) +
  # scalebar(RF.2075.RCP8.5.sf, 
  #          location = 'topleft',
  #          dist = 100,
  #          dist_unit = 'km',
  #          transform = F) +
  annotation_scale(bar_cols = c('grey40','grey10'),
                   text_cex = 1,
                   line_col = 'grey10',
                   location = 'tl',
                   pad_y = unit(0.5, 'cm')) +
  theme(panel.background = element_rect(fill = '#e0e0e0'),
        panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
        panel.border = element_rect(fill = NA, colour = 'black'),
        axis.text = element_text(size = 12.5),
        axis.text.y.right = element_text(hjust = 0.5, angle = 270),
        legend.text = element_text(size = 12),
        legend.position = c(0.70,0.1),
        legend.background = element_blank())

saveRDS(ggmap, 'Output/MaritimesRF_predictions_ggmap_RCP85_2075.rds')

# tmap version

map.2075 <-  tm_shape(RF.2075.RCP8.5, bbox = st_bbox(c(xmin =20000, 
                                           ymin = 4562134, 
                                           xmax = 1029435, 
                                           ymax = 5300000))) + 
  tm_raster(legend.show = F, 
            style = 'cat', 
            palette = MAR.palette$assigned[c(6,5,3,4,1,2)]) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +    
  tm_shape(admin.borders) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_shape(grid.anom) +
  tm_raster(legend.show = F,
            style = 'cat',
            palette = c(NA, 'black'),
            alpha = 0.15) +
  tm_shape(climate.sensitive) + 
  tm_borders(col = 'blue') +
  tm_scale_bar(breaks = c(0,50,100,200),
               color.dark = 'grey40',
               color.light = 'grey10',
               text.size = 4,
               position = c(0.01,0.90)) +
  tm_graticules(x = c(-67,-63,-59),
                y = c(47,45,43,41),
                alpha = 0.2,
                labels.size = 1) +
  # tm_add_legend(type = 'fill',
  #               labels = MAR.palette$name,
  #               col = MAR.palette$assigned) +
  tm_add_legend(type = 'fill',
                labels = 'Climate susceptible (anomalies)',
                border.col = 'blue',
                col = 'white') +
  tm_layout(legend.outside = F,
            legend.position = c('right','bottom'),
            bg.color = '#e0e0e0',
            legend.text.size = 1); map.2075

saveRDS(map.2075, 'Output/MaritimesRF_predictions_tmap_RCP85_2075.rds')

tmap_save(map.2075, 
          filename = 'Output/MaritimesRF_predictions_tmap_RCP85_2075.tiff', 
          compression = 'lzw',
          dpi = 300,
          width = 11)

# Combine panels for change in area and projected distributions

panelA <- readRDS('Output/Maritimes_ChangeClusterArea_flipped.rds')
panelB <- readRDS('Output/MaritimesRF_predictions_ggmap_RCP85_2075.rds')

panel.fig <- plot_grid(panelA, panelB, 
          align = 'h', 
          axis = 'bt', 
          rel_widths = c(1.05,1.5),
          labels = 'AUTO')

ggsave(panel.fig, filename = 'Output/Maritimes_2075_2panel.tiff',
       width = 9.2, height = 4.5, dpi = 300, compression = 'lzw')

# Quebec--------------------

# load spatial layers

RF.pred <- raster('Output/QC_PredClust_map.tif') # predicted distribution of 3 major clusters
RF.pred.sf <- as(RF.pred, 'SpatialPolygonsDataFrame') %>% # convert to sf object
  st_as_sf() %>% 
  dplyr::rename(cl = QC_PredClust_map) %>% 
  mutate(., assigned = case_when(
    cl == 1 ~ '#C93312',
    cl == 2 ~ '#FAEFD1',
    cl == 3 ~ '#DC863B'))
StudyArea <- st_read('Data/Shapefiles/QCRasterBoundary.shp') #study area boundaries
Land <- st_read('Data/Shapefiles/QC_landborders.shp') # Land polygons
uncertainty <- st_read('Output/QC_RF_uncertainty_0.7.shp') # grid cells with probability of assignment < 0.70

# ggplot version

#ggplot alternative

ggmap <- ggplot() +
  geom_sf(data = RF.pred.sf, 
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
  geom_sf(data = uncertainty, 
          col = 'black', 
          fill = NA,
          show.legend = F) +
  coord_sf(label_graticule = "NW",
           expand = T,
           xlim = c(30174, 1026175),
           ylim = c(5224423, 5820424)) +
  scale_x_continuous(name = '', 
                     breaks = c(-70,-66,-62,-58), 
                     expand = c(0.02,0)) +
  scale_y_continuous(name = '',
                     breaks = c(48,50,52),
                     expand = c(0.02,0)) +
  annotation_scale(bar_cols = c('grey40','grey10'),
                   text_cex = 1,
                   line_col = 'grey10',
                   location = 'tl') +
  annotate('text', 
           label = 'bold(NGSL)',
           parse = T,
           size = 6,
           x = 900000,
           y = 5800000) +
  theme(panel.background = element_rect(fill = '#e0e0e0'),
        panel.grid.major = element_line(colour = gray(0.4, alpha = 0.2)),
        panel.border = element_rect(fill = NA, colour = 'black'),
        plot.margin = unit(c(0,0,0,0), 'lines'),
        axis.text = element_text(size = 12.5),
        axis.text.y.left = element_text(hjust = 0.5, angle = 90));ggmap

saveRDS(ggmap, 'Output/QC_RF_predictions_ggmap_uncertainty_0.7.rds')
ggsave(ggmap, 
       filename = 'Output/QC_RF_predictions_ggmap_uncertainty_0.7.tiff', 
       compression = 'lzw', dpi = 300,
       width = 9, height = 5.4, unit = 'in')

# tmap version

map <-  tm_shape(RF.pred, bbox = st_bbox(StudyArea), ext = 1.04) + 
  tm_raster(legend.show = F, palette = c(NA, QC.palette$assigned)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +  
  tm_shape(Land) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_shape(uncertainty) + 
  tm_borders(col = 'black') +
  tm_scale_bar(breaks = c(0,50,100,200),
               color.dark = 'grey40',
               color.light = 'grey10',
               text.size = 4,
               position = c(0.01,0.85)) +
  tm_graticules(x = c(-70,-66,-62,-58),
                y = c(48,50,52),
                alpha = 0.2,
                labels.size = 1.5) +
  tm_layout(title = 'NGSL',
            title.position = c(0.77,0.9),
            title.fontface = 'bold',
            title.size = 2,
            bg.color = '#e0e0e0');map

saveRDS(map, 'Output/QC_RF_predictions_tmap_uncertainty_0.7.rds')

tmap_save(map, filename = 'Output/QC_RF_predictions_tmap_uncertainty_0.7.tiff', compression = 'lzw')

# Newfoundland & Labrador--------------------

# Load spatial layers

RF.pred <- raster('Output/NL_PredClust_map.tif') # predicted distribution of 3 major clusters
RF.pred.sf <- as(RF.pred, 'SpatialPolygonsDataFrame') %>% # convert to sf object
  st_as_sf() %>% 
  dplyr::rename(cl = NL_PredClust_map) %>% 
  mutate(., assigned = case_when(
    cl == 1 ~ '#85B5C0',
    cl == 2 ~ '#3F7F92',
    cl == 3 ~ '#E1AF00',
    cl == 4 ~ '#3F0A25',
    cl == 5 ~ '#E23B35'))
StudyArea <- st_read('Data/Shapefiles/NLRasterBoundary.shp') #study area boundaries
Land <- st_read('Data/Shapefiles/NL_landborders.shp') # Land polygons
uncertainty <- st_read('Output/NL_RF_uncertainty_0.7.shp') # grid cells with probability of assignment < 0.70


# ggplot version

ggmap <- ggplot() +
  geom_sf(data = RF.pred.sf, 
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
  geom_sf(data = uncertainty, 
          col = 'black', 
          fill = NA,
          show.legend = F) +
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

saveRDS(ggmap, 'Output/NL_RF_predictions_ggmap_uncertainty_0.7.rds')
ggsave(ggmap, 
       filename = 'Output/NL_RF_predictions_ggmap_uncertainty_0.7.tiff', 
       compression = 'lzw', dpi = 300,
       width = 5.6, height = 8.8, unit = 'in')

# tmap version

map <-  tm_shape(RF.pred, 
                 bbox = st_bbox(c(xmin =270169.8,
                                  ymin = 4742510.1,
                                  xmax = 1314169.8,
                                  ymax = 6398510.1)), 
                 ext = 1.04) + 
  tm_raster(legend.show = F, palette = c(NA, NL.palette$assigned)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +
  tm_shape(Land) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_shape(uncertainty) + 
  tm_borders(col = 'black') +
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
            bg.color = '#e0e0e0'); map

saveRDS(map, 'Output/NL_RF_predictions_tmap_uncertainty_0.7.rds')

tmap_save(map, filename = 'Output/NL_RF_predictions_tmap_uncertainty_0.7.tiff', compression = 'lzw')

# Gulf--------------------

# Load spatial layers

RF.pred <- raster('Output/Gulf_PredClust_map.tif') # predicted distribution of 3 major clusters
RF.pred.sf <- RF.pred.sf <- as(RF.pred, 'SpatialPolygonsDataFrame') %>% # convert to sf object
  st_as_sf() %>% 
  dplyr::rename(cl = Gulf_PredClust_map) %>% 
  mutate(., assigned = case_when(
    cl == 1 ~ '#899DA4',
    cl == 2 ~ '#002F2F',
    cl == 3 ~ '#7F1100',
    cl == 4 ~ '#725900'))

StudyArea <- st_read('Data/Shapefiles/GulfRasterBoundary.shp') #study area boundaries
Land <- st_read('Data/Shapefiles/Gulf_landborders.shp') # Land polygons
uncertainty <- st_read('Output/Gulf_RF_uncertainty_0.7.shp') # grid cells with probability of assignment < 0.70

# ggplot version

ggmap <- ggplot() +
  geom_sf(data = RF.pred.sf, 
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
  geom_sf(data = uncertainty, 
          col = 'black', 
          fill = NA,
          show.legend = F) +
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

saveRDS(ggmap, 'Output/SGSL_RF_predictions_ggmap_uncertainty_0.7.rds')
ggsave(ggmap,
       filename = 'Output/SGSL_RF_predictions_ggmap_uncertainty_0.7.tiff', 
       compression = 'lzw',  dpi = 300,
       width = 6.7, height = 5.8, units = 'in')

# Map panel

map <-  tm_shape(RF.pred, 
                 bbox = st_bbox(StudyArea), 
                 ext = 1.04) + 
  tm_raster(legend.show = F, palette = c(NA, GULF.palette$assigned)) +
  tm_shape(StudyArea) + 
  tm_borders(lwd = 1, 
             col = 'black') +
  tm_shape(Land) +
  tm_polygons(col = 'grey30',
              border.col = 'grey20',
              lwd = 0.25) +
  tm_shape(uncertainty) + 
  tm_borders(col = 'black') +
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
            bg.color = '#e0e0e0'); map

saveRDS(map, 'Output/Gulf_RF_predictions_tmap_uncertainty_0.7.rds')

tmap_save(map, filename = 'Output/Gulf_RF_predictions_tmap_uncertainty_0.7.tiff', compression = 'lzw')


#### Variable Importance Plots ####

# Faceted by cluster

varImp.ls <- list.files('Output/', pattern = 'ImpPlot.rds$', full.names = T) %>% # list of files
  map(readRDS) %>% # read in ggplots
  map(., ~ . + theme(axis.title.x = element_blank(), # remove x-axis label
                     axis.text = element_text(size = 12), # increase axis text size
                     strip.text = element_text(size = 12)) + # increase strip text size
        scale_y_continuous(limits = c(-0.05,0.55))) %>% # apply consistent scale across regions
  set_names(c('SGSL', 'MAR', 'NL', 'NGSL')) # name list elements by region


# Plot variable importance plots in 4 rows; 1 per region
p <- ggdraw() +
  draw_plot(varImp.ls$NGSL, x = 0, y = 0.75, width = 0.625, height = 0.25) +
  draw_plot(varImp.ls$SGSL, x = 0, y = 0.5, width = 0.75, height = 0.25) +
  draw_plot(varImp.ls$NL, x = 0, y = 0.25, width = 0.875, height = 0.25) +
  draw_plot(varImp.ls$MAR, x = 0.01, y = 0, width = 0.99, height = 0.25) +
  draw_plot_label(label = c('NGSL','SGSL','NL','MAR'), size = 15, hjust = 1,
                  x = c(0.12,0.12,0.12,0.12), y = c(1,0.75,0.5,0.25))
# Save output
ggsave(plot = p, 'Output/VarImp_4regions_ByClass.tiff', height = 17, width = 16, compression = 'lzw')

# 4 panel figure - Variable Importance - whole model only

varImp.wm <- list.files('Output/', pattern = 'WM.rds$', full.names = T) %>% # list of files
  map(readRDS) %>% # read in ggplots
  map(., ~ . + theme(axis.title.x = element_blank(),# remove x-axis label
                     strip.text = element_text(size = 14, face = 'bold'), # adjust font of strip text
                     axis.text = element_text(size = 12)) + # increase axis text size 
        scale_y_continuous(limits = c(-0.05,0.25))) %>% # apply consistent scale across regions
  set_names(c('SGSL', 'MAR', 'NL', 'NGSL')) %>%  # name list elements by region
  .[c('NGSL','NL','SGSL','MAR')] # re-arrange list elements


# Plot on 2 X 2 grid

p <- plot_grid(plotlist = varImp.wm, 
               nrow = 2, ncol = 2, 
               align = 'hv', 
               axis = 'lb', 
               labels = 'Mean Decrease in Accuracy',
               label_x = c(0.7),
               label_y = c(-1),
               hjust = -0.5) + 
  theme(plot.margin = margin(0,0,40,0))

# Save output

ggsave(plot = p, 'Output/VarImp_4regions_WholeModel.tiff', height = 9, width = 11, unit = 'in', compression = 'lzw', dpi = 300)


#### Matrix plot for Indicator Species ####

# colour palettes

source('Code/SPERA_colour_palettes.R') 
palettes <- bind_rows(GULF.palette, MAR.palette, NL.palette, QC.palette, .id = 'Region') %>% 
  mutate(Region = case_when(
    Region == 1 ~ 'SGSL',
    Region == 2 ~ 'MAR',
    Region == 3 ~ 'NL',
    Region == 4 ~ 'NGSL'
  ))

# Read in tables with IndVal and Freq

indval.tbl <- list.files('Output/', pattern = 'Indi.+txt$', full.names = T) %>% # list of text files
  map(fread) %>% # read in tables 
  set_names(c('SGSL','MAR','NL','NGSL')) %>% # name list elements
  rbindlist(., idcol = 'Region') %>% # row bind & add id column for region
  select(-5) %>% # remove Frequency
  #label clusters
  mutate(Ecoregion2 = case_when(
    Region == 'SGSL' & Ecoregion == 1 ~ 'Magdalen Shallows',
    Region == 'SGSL' & Ecoregion == 2 ~ 'Inshore/Magdalen Is.',
    Region == 'SGSL' & Ecoregion == 3 ~ 'Laurentian Channel',
    Region == 'SGSL' & Ecoregion == 4 ~ "Northumberland Strait/St. George's Bay",
    Region == 'MAR' & Ecoregion == 1 ~ 'WSS/Outer BoF',
    Region == 'MAR' & Ecoregion == 2 ~ 'WSS: Banks/Inner BoF',
    Region == 'MAR' & Ecoregion == 3 ~ 'ESS: Banks',
    Region == 'MAR' & Ecoregion == 4 ~ 'ESS',
    Region == 'MAR' & Ecoregion == 5 ~ 'Laurentian Channel/Shelf Break',
    Region == 'MAR' & Ecoregion == 6 ~ 'Slope',
    Region == 'NL' & Ecoregion == 1 ~ 'Inner Shelf',
    Region == 'NL' & Ecoregion == 2 ~ 'Outer Shelf',
    Region == 'NL' & Ecoregion == 3 ~ 'Grand Banks',
    Region == 'NL' & Ecoregion == 4 ~ 'Slope',
    Region == 'NL' & Ecoregion == 5 ~ 'Laurentian Channel/Shelf Break',
    Region == 'NGSL' & Ecoregion == 1 ~ 'Deep Channels',
    Region == 'NGSL' & Ecoregion == 2 ~ 'Shallow Banks & Straits',
    Region == 'NGSL' & Ecoregion == 3 ~ 'Channel Heads & Slopes',
    TRUE ~ NA_character_
  )) %>% 
  inner_join(., palettes, by = c('Region','Ecoregion2' = 'name')) %>% # join colour palettes
  # fix inconsistent taxa spellings
  mutate(Species = gsub('^Gorgonocephalus sp $', 'Gorgonocephalus spp.', Species),
         Species = gsub('sp $', 'spp.', Species),
         Species = gsub('c $', '', Species),
         Species = gsub('o $', '', Species),
         Species = gsub('o\\.', '', Species),
         Species = gsub('f\\.', '', Species),
         Species = gsub('c\\.', '', Species),
         Species = gsub('sp\\.', 'spp.', Species)) %>%
  mutate(Dupl = duplicated(Species)) %>% # Identify indicator species that are shared across 2+ regions
  mutate(Region = factor(Region, levels = c('NGSL','SGSL','NL','MAR'))) %>% # re-order factor levels for Region
  arrange(Species) # arrange by species alphabetically # 121 taxa across all regions

dupl <- filter(indval.tbl, Dupl == T) %>% pull(Species) %>% unique() # 47 unique taxa shared by 2+ regions

# Matrix of full taxa list

mat.ind <- ggplot(indval.tbl, aes(x = Region, y = reorder(Species, desc(Species)), fill = assigned)) +
  geom_tile(width = 1, height = 1, colour = 'grey50', size = 0.5) +
  scale_fill_identity(guide = 'none') +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(y = 'Indicator taxa') +
  theme(panel.background = element_rect(fill = 'grey85', colour = 'grey50'),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(face = 'italic')); mat.ind

ggsave(plot = mat.ind, 'Output/IndVal_MatrixPlot_FullList.tiff', width = 4.25, height = 12, compression = 'lzw')

# Matrix with indicator taxa shared by at least 2 regions

mat.ind2 <- ggplot(filter(indval.tbl, Species %in% dupl),
                   aes(x = Region, 
                       y = reorder(Species, desc(Species)),
                       fill = assigned)) +
  geom_tile(width = 1, height = 1, colour = 'grey50', size = 0.5) +
  scale_fill_identity(guide = 'none') +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(y = 'Indicator taxa') +
  theme(panel.background = element_rect(fill = 'grey85', colour = 'grey50'),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = 'bold', size = 14),
        axis.text.y = element_text(face = 'italic')); mat.ind2

ggsave(plot = mat.ind2, 'Output/IndVal_MatrixPlot_SharedTaxa.tiff', width = 5.5, height = 12, compression = 'lzw')

#### Cluster validity Indices - 4 regions ####

# Read in files

p.cvi <- list.files('Output/', pattern = 'CVI.+rds', full.names = T) %>% # list of rds objects
  map(readRDS) %>% # read rds objects
  map(., ~. + theme(axis.title = element_blank())) %>% # remove axis labels
  set_names(c('SGSL', 'MAR', 'NL', 'NGSL')) %>%  # name list elements by region
  .[c('NGSL','SGSL','NL','MAR')] # re-arrange list elements

# Plot on 4 column grid

p <- plot_grid(plotlist = p.cvi, 
               nrow = 1, ncol = 4, 
               align = 'hv', 
               axis = 'lb', 
               labels = c('Number of Clusters'),
               label_x = c(0.5),
               label_y = c(0),
               hjust = -1) + 
  theme(plot.margin = margin(40,0,40,40)) +
  draw_label('Index value', x = 0, y = 0.5, angle = 90, fontface = 'bold', size = 14, vjust = -1) + 
  annotate(geom = 'text', 
           x = c(0.125,0.375,0.635,0.875), 
           y = c(1,1,1,1), label = c('NGSL','SGSL','NL','MAR'),
           vjust = -1,
           hjust = 0,
           fontface = 'bold')

# Save output

ggsave(plot = p, 'Output/CVI_4regions.tiff', height = 7.5, width = 7, unit = 'in', compression = 'lzw', dpi = 300)


#### Study Domain Map ####

# list shapefiles with land borders
landfiles <- list.files(pattern = 'gadm', full.names = T) %>% 
  discard(~ grepl('FRA', .x)) %>% 
  discard(~ grepl('_0_', .x))

# read in shapefiles, convert to sf df, combine, and crop
admin.borders <- map(landfiles, readRDS) %>% 
  map(., st_as_sf) %>% 
  rbindlist() %>% 
  st_as_sf() %>% 
  st_crop(., c(xmin = -69.5, ymin = 39, xmax = -45, ymax = 59)) %>% 
  st_simplify(., dTolerance = 0.025)

  
# read in shapefiles with study domains

names <- c('MAR','NGSL','NL','SGSL')
StudyArea <- list.files(pattern = '^StudyArea.+.shp', path = 'Data/Shapefiles', full.names = T) %>% 
  map(., st_read) %>%
  map(., ~ st_transform(., 4326)) %>%
  set_names(., names) %>% 
  map2(., names, ~ mutate(.x, region = .y)) %>% 
  map(., ~select(., region, geometry)) %>% 
  map(., ~st_cast(., 'MULTIPOLYGON')) %>% 
  rbindlist() %>% 
  st_as_sf() %>% 
  arrange((region))

# points for trawl sets

RVdata <- readRDS('R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVdata_4regions.rds') %>% 
  st_as_sf(., coords = c('longitude','latitude'), crs = 4326) %>% 
  filter(., year > 2007) %>%  
  distinct(., region, year, geometry, .keep_all = T) %>% 
  filter(., st_intersects(., st_union(StudyArea), sparse = F)) %>% 
  mutate(., region2 = case_when(
    region == 'GULF' ~ 'SGSL',
    region == 'MARITIME' ~ 'MAR',
    region == 'QUEBEC' ~ 'NGSL',
    region == 'NEWFOUNDLAND' ~ 'NL'
  ))

# map study area and trawl locations

pal <- wes_palette('Darjeeling1', n = 4, type = 'discrete') #colour palette
names(pal) = c('NGSL','NL','SGSL', 'MAR') # name palette elements

study.dom.gg <- ggplot() +
  geom_sf(data = admin.borders, fill = 'grey30', col = 'grey20',lwd = 0.25) +
  geom_sf(data = StudyArea, aes(fill = region), col = NA, alpha = 0.3) +
  geom_sf(data = RVdata, aes(col = region2), size = 0.5, shape = 16) + 
  scale_color_manual(values = pal, 
                     breaks = c('NGSL','SGSL','NL','MAR'),
                     aesthetics = c('fill','colour'), 
                     guide = guide_legend('Region', 
                                          override.aes = list(alpha = 0.3),
                                          title.theme = element_text(face = 'bold'))) +
  coord_sf(xlim = c(-69.5,-45.5), ylim = c(40.5,58),expand = F) +
  # annotation_scale(aes(location = 'tl',
  #                      line_col = 'grey40',
  #                      text_col = 'grey10',
  #                      width_hint = 1/5),
  #                  bar_cols = c('grey10','grey40')) +
  theme(panel.border = element_rect(colour = 'black', fill = NA),
        plot.margin = unit(c(0,0.5,0,0), 'lines'),
        panel.grid.major = element_blank()); study.dom.gg

ggsave(plot = study.dom.gg, 'Output/StudyDomain.tiff', height = 6, width = 6.5, unit = 'in', compression = 'lzw', dpi = 300)
  
#### Map with grid cells coloured by cluster assignment - 4 regions combined 

# Read in individual maps and legends

# map.pdf.ls <- list.files(path = 'Output/', pattern = 'Col.+pdf', full.names = T)
# map.rds.ls <- list.files(path = 'Output/', pattern = 'Col.+rds', full.names = T)
# leg.pdf.ls <- list.files(path = 'Output/', pattern = 'tmap.+pdf', full.names = T)
# leg.rds.ls <- list.files(path = 'Output/', pattern = 'tmap_leg.+rds', full.names = T)
# 
# maps <- map(map.rds.ls, readRDS) %>% 
#   set_names(c('SGSL','MAR','NL','NGSL')) %>% 
#   map(., ~ . + tm_layout(title = ''))
# 
# legends <- map(leg.rds.ls, readRDS) %>% 
#   set_names(c('SGSL','MAR','NL','NGSL'))
# 
# # Define layout with matrix
# 
# lay <- rbind(c(NA,NA,NA,NA,3,3,3),
#              c(4,4,4,4,3,3,3),
#              c(4,4,4,4,3,3,3),
#              c(4,4,4,4,3,3,3),
#              c(1,1,1,2,2,2,2),
#              c(1,1,1,2,2,2,2),
#              c(1,1,1,2,2,2,2))
# pdf(file = 'Output/ColorGrid_4regions.pdf')
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(3,2)))
# print(maps$NGSL, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(maps$SGSL, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
# print(maps$NL, vp = viewport(layout.pos.row = c(1:2), layout.pos.col = 2))
# print(maps$MAR, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
# print(legends$NGSL, vp = viewport(x = 0.1, y = 0.8, height = 0.25, width = 1, just = c(0,0)))
# print(legends$SGSL, vp = viewport(x = 0.1, y = 0.6, height = 0.25, width = 1, just = c(0,0)))
# print(legends$NL, vp = viewport(x = 0.4, y = 0.8, height = 0.25, width = 1, just = c(0,0)))
# print(legends$MAR, vp = viewport(x = 0.4, y = 0.6, height = 0.25, width = 1, just = c(0,0)))
# dev.off()
