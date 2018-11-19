## Function that creates grid within study area of given resolution and excludes grid cells that intersect with land

GridFilter2<-function(regionshape, resol = 1, landshape){
  
  require(rgdal)
  require(raster)
  require(rgeos)
  
  grid <- raster(extent(regionshape))
  res(grid) <- resol
  proj4string(grid)<-proj4string(regionshape)
  grid <- rasterize(regionshape, grid)
  gridpolygon <- rasterToPolygons(grid)
  gridpolygon$layer <- c(1:length(gridpolygon$layer))
  Grid.land <- unlist(over(landshape, gridpolygon, returnList = TRUE), use.names = FALSE) #select GridID for cells that intersect land
  gridpolygonNoLand <- gridpolygon[!gridpolygon@data$layer %in% Grid.land,]}