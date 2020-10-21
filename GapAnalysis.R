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
  
# Remove raster areas outside soil polygon with 25 m buffer
  soilsAll <- sp::aggregate(soils,dissolve=T)
  soilsBuff <- raster::buffer(soils, width=-25, dissolve = F)
  
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

# Mask weird really high values in 2015 -- REVISIT THIS eventually
  chm15@data@values[chm15@data@values>200 & !is.na(chm15@data@values)] <- NA
    
# Save
  raster::writeRaster(chm09, "CHM_2009_QAQC.tif")
  raster::writeRaster(chm15, "CHM_2015_QAQC.tif")
  raster::writeRaster(chm17, "CHM_2017_QAQC.tif")
  raster::writeRaster(chm18, "CHM_2018_QAQC.tif")
  raster::writeRaster(chm19, "CHM_2019_QAQC.tif")
  raster::writeRaster(chm20, "CHM_2020_QAQC.tif")
  
#### MAKE GRAPHS ####
  
  #USE A BUFFER AROUND EDGES
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
  chm09 <- raster::mask(chm09, buffer)
  chm15 <- raster::mask(chm15, buffer)
  chm17 <- raster::mask(chm17, buffer)
  chm18 <- raster::mask(chm18, buffer)
  chm19 <- raster::mask(chm19, buffer)
  chm20 <- raster::mask(chm20, buffer)
  
  dens09 <- density(chm09@data@values, na.rm=T)
  dens15 <- density(chm15@data@values, na.rm=T)
  dens17 <- density(chm17@data@values, na.rm=T)
  dens18 <- density(chm18@data@values, na.rm=T)
  dens19 <- density(chm19@data@values, na.rm=T)
  dens20 <- density(chm20@data@values, na.rm=T)
  
  plot(dens09, col="black", lwd=2,
       main="Canopy height distribution",
       xlab="Canopy height (m)",
       ylim=c(0,0.01),
       xlim=c(-1,10))
  lines(dens15, col="#99d8c9", lwd=2)
  lines(dens17, col="#66c2a4", lwd=2)
  lines(dens18, col="#41ae76", lwd=2)
  lines(dens19, col="#238b45", lwd=2)
  lines(dens20, col="#005824", lwd=2)
  
  legend(x=35,y=0.05,
         c("2009","2015","2017","2018","2019","2020"),
         col=c("black","#99d8c9","#66c2a4","#41ae76","#238b45","#005824"),
         lwd=2,
         bty="n")
  

#### BINARY LOW CANOPY ####
  
  htThresh <-5
  
  low09 <- chm09
  low09@data@values[low09@data@values<=htThresh & !is.na(low09@data@values)] <- 1
  low09@data@values[low09@data@values>htThresh & !is.na(low09@data@values)] <- 0
  
  low15 <- chm15
  low15@data@values[low15@data@values<=htThresh & !is.na(low15@data@values)] <- 1
  low15@data@values[low15@data@values>htThresh & !is.na(low15@data@values)] <- 0
  
  low17 <- chm17
  low17@data@values[low17@data@values<=htThresh & !is.na(low17@data@values)] <- 1
  low17@data@values[low17@data@values>htThresh & !is.na(low17@data@values)] <- 0
  
  low18 <- chm18
  low18@data@values[low18@data@values<=htThresh & !is.na(low18@data@values)] <- 1
  low18@data@values[low18@data@values>htThresh & !is.na(low18@data@values)] <- 0
  
  low19 <- chm19
  low19@data@values[low19@data@values<=htThresh & !is.na(low19@data@values)] <- 1
  low19@data@values[low19@data@values>htThresh & !is.na(low19@data@values)] <- 0
  
  low20 <- chm20
  low20@data@values[low20@data@values<=htThresh & !is.na(low20@data@values)] <- 1
  low20@data@values[low20@data@values>htThresh & !is.na(low20@data@values)] <- 0
  
  # Estimate maximum age of low canopy areas in 2020
  low20age <- low20
  low20age@data@values[low20@data@values==1 & !is.na(low20@data@values)] <- 0
  low20age@data@values[low20@data@values==1 & low19@data@values==0 & !is.na(low20@data@values) & !is.na(low19@data@values)] <- 1
  low20age@data@values[low20@data@values==1 & low19@data@values==1 & 
                         low18@data@values==0 & !is.na(low20@data@values) & !is.na(low18@data@values)] <- 2
  low20age@data@values[low20@data@values==1 & low19@data@values==1 & low18@data@values==1 & 
                         low17@data@values==0 & !is.na(low20@data@values) & !is.na(low17@data@values)] <- 3
  low20age@data@values[low20@data@values==1 & low19@data@values==1 & low18@data@values==1 & low17@data@values==1 & 
                         low15@data@values==0 & !is.na(low20@data@values) & !is.na(low15@data@values)] <- 5
  low20age@data@values[low20@data@values==1 & low19@data@values==1 & low18@data@values==1 & low17@data@values==1 & low15@data@values==1 & 
                         low09@data@values==0 & !is.na(low20@data@values) & !is.na(low09@data@values)] <- 11
  
  low20NA <- low20
  low20NA@data@values[low20@data@values==1 & !is.na(low20@data@values)] <- 0
  low20NA@data@values[low20@data@values==1 & low20age@data@values == 0 & !is.na(low20@data@values)] <- 1
  
  100*length(low20age@data@values[low20age@data@values==1 & !is.na(low20age@data@values)])/length(low20age@data@values[low20@data@values==1 & !is.na(low20age@data@values)])
  100*length(low20age@data@values[low20age@data@values==2 & !is.na(low20age@data@values)])/length(low20age@data@values[low20@data@values==1 &!is.na(low20age@data@values)])
  100*length(low20age@data@values[low20age@data@values==3 & !is.na(low20age@data@values)])/length(low20age@data@values[low20@data@values==1 &!is.na(low20age@data@values)])
  100*length(low20age@data@values[low20age@data@values==5 & !is.na(low20age@data@values)])/length(low20age@data@values[low20@data@values==1 &!is.na(low20age@data@values)])
  100*length(low20age@data@values[low20age@data@values==11 & !is.na(low20age@data@values)])/length(low20age@data@values[low20@data@values==1 &!is.na(low20age@data@values)])
  
#### BINARY NEW GAPS ####  