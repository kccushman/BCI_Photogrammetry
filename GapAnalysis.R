# Do preliminary gap analyses

# NOTE (October 20, 2020): There are some areas where clouds have apparently 
# *not* be automatically masked by current QAQC routines. I think this is a 
# small issue but needs to be fixed.

#### READ DATA ####

# Read in canopy height rasters:
  dsm09 <- raster::raster("CHM_2009.tif")
  dsm15 <- raster::raster("CHM_2015_corrected.tif")
  dsm17 <- raster::raster("CHM_2017_corrected.tif")
  dsm18 <- raster::raster("CHM_2018_corrected.tif")
  dsm19 <- raster::raster("CHM_2019_corrected.tif")
  dsm20 <- raster::raster("CHM_2020_corrected.tif")

# Read in soil polygon layers
  soils <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
  soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
# Read in BCI DEM
  dem <- raster::raster("D:/BCI_Spatial/LidarDEM_BCI.tif")
  dem <- raster::resample(dem,dsm15)
  
# Read in QAQC data
  qaqc <- read.csv("gridInfo_QAQC.csv")
 
#### PROCESS DATA ####
  
# Remove raster areas outside soil polygon
  dsm09 <- raster::mask(dsm09, soils)
  dsm15 <- raster::mask(dsm15, soils)  
  dsm17 <- raster::mask(dsm17, soils)  
  dsm18 <- raster::mask(dsm18, soils)  
  dsm19 <- raster::mask(dsm19, soils)  
  dsm20 <- raster::mask(dsm20, soils)
  
# Subtract ground elevation
  chm09 <- dsm09-dem
  chm15 <- dsm15-dem
  chm17 <- dsm17-dem
  chm18 <- dsm18-dem
  chm19 <- dsm19-dem
  chm20 <- dsm20-dem
  
# Set areas that fail QAQC to "NA"
  for(i in 1:dim(qaqc)[1]){
    
    # make extent object for current tile
      x1 <- qaqc[i, "xmin"] 
      x2 <- qaqc[i, "xmax"]
      y1 <- qaqc[i, "ymin"] 
      y2 <- qaqc[i, "ymax"] 
      extent_i <- raster::extent(c(x1,x2,y1,y2))
      extent_i <- as(extent_i, 'SpatialPolygons')
    
    # set values within failed tiles to NA  
    if(qaqc$Use15[i]==F){
      chm15 <- raster::mask(chm15, extent_i, inverse = T)
    }    
      
    if(qaqc$Use17[i]==F){
      chm17 <- raster::mask(chm17, extent_i, inverse = T)
    }
      
    if(qaqc$Use18[i]==F){
      chm18 <- raster::mask(chm18, extent_i, inverse = T)
    }
    
    if(qaqc$Use19[i]==F){
      chm19 <- raster::mask(chm19, extent_i, inverse = T)
    }
    
    if(qaqc$Use20[i]==F){
      chm20 <- raster::mask(chm20, extent_i, inverse = T)
    }
      print(i)
  }

# Make sure all years have the same extent
  chm09 <- raster::crop(chm09, raster::extent(soils))
  chm15 <- raster::crop(chm15, raster::extent(soils))
  chm17 <- raster::crop(chm17, raster::extent(soils))
  chm18 <- raster::crop(chm18, raster::extent(soils))
  chm19 <- raster::crop(chm19, raster::extent(soils))
  chm20 <- raster::crop(chm20, raster::extent(soils))
  
#### MAKE GRAPHS ####
  
  hist(dsm09, breaks=seq(-200,200,1),
       xlim=c(-5,75),
       main = "Canopy height distribution 2009")
  hist(dsm15, breaks=seq(-200,200,1),
       xlim=c(-5,75),
       main = "Canopy height distribution 2015")
  hist(dsm17, breaks=seq(-200,200,1),
       xlim=c(-5,75),
       main = "Canopy height distribution 2017")
  hist(dsm18, breaks=seq(-200,200,1),
       xlim=c(-5,75),
       main = "Canopy height distribution 2018")
  hist(dsm19, breaks=seq(-200,200,1),
       xlim=c(-5,75),
       main = "Canopy height distribution 2019")
  hist(dsm20, breaks=seq(-200,200,1),
       xlim=c(-5,75),
       main = "Canopy height distribution 2020")
  
  
#### BINARY LOW CANOPY ####
  
  low09 <- chm09
  low09@data@values[low09@data@values<=5] <- 1
  low09@data@values[low09@data@values>5] <- 0
  
  low15 <- chm15
  low15@data@values[low15@data@values<=5] <- 1
  low15@data@values[low15@data@values>5] <- 0
  
  low17 <- chm17
  low17@data@values[low17@data@values<=5] <- 1
  low17@data@values[low17@data@values>5] <- 0
  
  low18 <- chm18
  low18@data@values[low18@data@values<=5] <- 1
  low18@data@values[low18@data@values>5] <- 0
  
  low19 <- chm19
  low19@data@values[low19@data@values<=5] <- 1
  low19@data@values[low19@data@values>5] <- 0
  
  low20 <- chm20
  low20@data@values[low20@data@values<=5] <- 1
  low20@data@values[low20@data@values>5] <- 0
  
  
  
#### BINARY NEW GAPS ####  