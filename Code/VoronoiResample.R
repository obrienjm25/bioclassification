#function takes a large spatial points df without a regular grid, clips
#to a smaller study area, and uses voronoi (thiessen) polygons to map points 
#to a regular grid with resolution of your choosing before conversion to
#raster layer

VoronoiResample <- function(shpth, raspth, inpFile, 
                            ext, SApoly, SAname, resol, EnVar) {
  
  source("Code/Voronoipolygons2.R") # voronoi (thiessen) polygon function
  
  NWA_pts <- readOGR(paste(shpth,inpFile,sep='')) #read in pt shapefile
  summary_stats <- colnames(NWA_pts@data)
  SA_pts <- crop(NWA_pts, ext) #crop pts to extent
  SA_pts <- spTransform(SA_pts, proj4string(SApoly)) #reproject to cartesian projection of study area
  vor_SA <- voronoipolygons(SA_pts) # create voronoi polygons for each point
  proj4string(vor_SA) <- proj4string(SApoly)

  # create fishnet grid of 4km X 4km polygons spanning study area
  r.template <- raster(extent(SApoly), res = resol, crs = proj4string(SApoly))
  fishnet <- rasterToPolygons(r.template)
  fishnet$layer <- c(1:length(fishnet$layer))
  
  # convert voronoi and fishnet polygon df's to sf objects
  vor_SA <- st_as_sf(vor_SA)
  fishnet <- st_as_sf(fishnet)
  
  # spatial join of fishnet grid with voronoi polygons
  #group by fishnet polygon id and aggregate
  SPDF_4km <- st_join(fishnet, vor_SA) %>% 
    group_by(layer) %>% 
    summarise_at(summary_stats, mean, na.rm = T) %>% 
    as(.,'Spatial')
  
  #rasterize each attribute in spatial polygons df & save rasters (4km res)
  
  for (i in summary_stats){
    
    rasterize(SPDF_4km, r.template, field = SPDF_4km@data[i], filename = paste(raspth,paste(SAname,i,EnVar,sep="_"),".tif",sep = ''), overwrite = T)
    
  }
  
  
}

