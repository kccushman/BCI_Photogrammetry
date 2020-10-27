# Do preliminary gap analyses

# NOTE (October 20, 2020): There are some areas where clouds have apparently 
# *not* be automatically masked by current QAQC routines. I think this is a 
# small issue but needs to be fixed.

#### READ RAW DATA ####

# Read in canopy height rasters:
  dsm09 <- raster::raster("DSM_2009.tif")
  dsm15 <- raster::raster("DSM_2015_corrected.tif")
  dsm17 <- raster::raster("DSM_2017_corrected.tif")
  dsm18 <- raster::raster("DSM_2018_corrected.tif")
  dsm19 <- raster::raster("DSM_2019_corrected.tif")
  dsm20 <- raster::raster("DSM_2020_corrected.tif")

# Read in soil polygon layers
  soils <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
  soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
# Read in forest age
  age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  ageUse <- age[!(age$TYPE=="Clearings"),]

# Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    
# Read in BCI DEM
  dem <- raster::raster("D:/BCI_Spatial/LidarDEM_BCI.tif")
  dem <- raster::resample(dem,dsm15)
  
# Read in QAQC data
  qaqc <- read.csv("gridInfo_QAQC.csv")
 
#### PROCESS DATA ####
  
# Remove raster areas outside BCI perimeter (exclude within 25 m of lake)
  dsm09 <- raster::mask(dsm09, buffer)
  dsm15 <- raster::mask(dsm15, buffer)  
  dsm17 <- raster::mask(dsm17, buffer)  
  dsm18 <- raster::mask(dsm18, buffer)  
  dsm19 <- raster::mask(dsm19, buffer)  
  dsm20 <- raster::mask(dsm20, buffer)

# Remove raster areas in clearings
  dsm09 <- raster::mask(dsm09, ageUse)
  dsm15 <- raster::mask(dsm15, ageUse)  
  dsm17 <- raster::mask(dsm17, ageUse)  
  dsm18 <- raster::mask(dsm18, ageUse)  
  dsm19 <- raster::mask(dsm19, ageUse)  
  dsm20 <- raster::mask(dsm20, ageUse)

# Crop to ensure each raster has same extent
  dsm09 <- raster::crop(dsm09, extent(ageUse))
  dsm15 <- raster::crop(dsm15, extent(ageUse))  
  dsm17 <- raster::crop(dsm17, extent(ageUse))  
  dsm18 <- raster::crop(dsm18, extent(ageUse))  
  dsm19 <- raster::crop(dsm19, extent(ageUse))  
  dsm20 <- raster::crop(dsm20, extent(ageUse))    
  
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
  chm09 <- raster::crop(chm09, raster::extent(ageUse))
  chm15 <- raster::crop(chm15, raster::extent(chm09))
  chm17 <- raster::crop(chm17, raster::extent(chm09))
  chm18 <- raster::crop(chm18, raster::extent(chm09))
  chm19 <- raster::crop(chm19, raster::extent(chm09))
  chm20 <- raster::crop(chm20, raster::extent(chm09))
  
  
# Cloud masks for 2017-2020
  mask17 <- raster::raster("CloudMask_2017.tif")
  mask18 <- raster::raster("CloudMask_2018.tif")
  mask19 <- raster::raster("CloudMask_2019.tif")
  mask20 <- raster::raster("CloudMask_2020.tif")
  
  # Resample to resolution of CHMs
  mask17 <- raster::resample(mask17, chm17)
  mask18 <- raster::resample(mask18, chm18)
  mask19 <- raster::resample(mask19, chm19)
  mask20 <- raster::resample(mask20, chm20)
  
  # Remove cloud pixels
  chm17[!(mask17==1)] <- NA
  chm18[!(mask18==1)] <- NA
  chm19[!(mask19==1)] <- NA
  chm20[!(mask20==1)] <- NA
    
# # Save
#   raster::writeRaster(chm09, "CHM_2009_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm15, "CHM_2015_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm17, "CHM_2017_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm18, "CHM_2018_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm19, "CHM_2019_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm20, "CHM_2020_QAQC.tif", overwrite=T)
  

#### READ PROCESSED DATA ####  
  
  # Canopy height models for all years
    chm09 <- raster::raster("CHM_2009_QAQC.tif")
    chm15 <- raster::raster("CHM_2015_QAQC.tif")
    chm17 <- raster::raster("CHM_2017_QAQC.tif")
    chm18 <- raster::raster("CHM_2018_QAQC.tif")
    chm19 <- raster::raster("CHM_2019_QAQC.tif")
    chm20 <- raster::raster("CHM_2020_QAQC.tif")
    
  
#### CANOPY HEIGHT DISTRIBUTION ####

  dens09 <- density(raster::values(chm09), na.rm=T)
  dens15 <- density(raster::values(chm15), na.rm=T)
  dens17 <- density(raster::values(chm17), na.rm=T)
  dens18 <- density(raster::values(chm18), na.rm=T)
  dens19 <- density(raster::values(chm19), na.rm=T)
  dens20 <- density(raster::values(chm20), na.rm=T)
  
  plot(dens09, col="black", lwd=2,
       main="Canopy height distribution",
       xlab="Canopy height (m)",
       ylim=c(0,0.05),
       xlim=c(-1,50))
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
  
  plot(dens09, col="black", lwd=2,
       main="Canopy height distribution",
       xlab="Canopy height (m)",
       ylim=c(0,0.03),
       xlim=c(-1,15))
  lines(dens15, col="#99d8c9", lwd=2)
  lines(dens17, col="#66c2a4", lwd=2)
  lines(dens18, col="#41ae76", lwd=2)
  lines(dens19, col="#238b45", lwd=2)
  lines(dens20, col="#005824", lwd=2)
  
  

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
  
  d17to18 <- chm18-chm17
  d18to19 <- chm19-chm18
  d19to20 <- chm20-chm19
  
  colBrks2 <- c(-100,-5,-1,-0.5,0.5,1,5,100)
  colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                "white",
                                "aliceblue","cornflowerblue","darkblue"))
  
  raster::plot(d17to18,
               col = colPal2(length(colBrks2)-1),
               breaks = colBrks2,
               main = "2017 to 2018 height change")
  
  raster::plot(d18to19,
               col = colPal2(length(colBrks2)-1),
               breaks = colBrks2,
               main = "2018 to 2019 height change")
  
  raster::plot(d19to20,
               col = colPal2(length(colBrks2)-1),
               breaks = colBrks2,
               main = "2019 to 2020 height change")

  # Use ForestGapR package to delineate new gaps (kind of hack-y)
    gaps17to18 <- ForestGapR::getForestGaps(d17to18,
                                            threshold = -5, size=c(10,10^4))
    gaps17to18sp <- ForestGapR::GapSPDF(gaps17to18)
    gaps17to18sp@data$area <- NA
    gaps17to18sp@data$perimeter <- NA
    for(i in 1:length(gaps17to18sp)){
      gaps17to18sp[gaps17to18sp$gap_id==i,"area"] <- raster::area(gaps17to18sp[gaps17to18sp$gap_id==i,])
      gaps17to18sp[gaps17to18sp$gap_id==i,"perimeter"] <- spatialEco::polyPerimeter(gaps17to18sp[gaps17to18sp$gap_id==i,])
      
    }
    
    gaps17to18sp@data$ratio <- gaps17to18sp@data$area/gaps17to18sp@data$perimeter
    gaps17to18_circ <- gaps17to18sp[gaps17to18sp@data$ratio<5,]
    
    
    gaps18to19 <- ForestGapR::getForestGaps(d18to19,
                                            threshold = -5, size=c(10,10^4))
    gaps18to19sp <- ForestGapR::GapSPDF(gaps18to19)
    gaps18to19sp@data$area <- NA
    gaps18to19sp@data$perimeter <- NA
    for(i in 1:length(gaps18to19sp)){
      gaps18to19sp[gaps18to19sp$gap_id==i,"area"] <- raster::area(gaps18to19sp[gaps18to19sp$gap_id==i,])
      gaps18to19sp[gaps18to19sp$gap_id==i,"perimeter"] <- spatialEco::polyPerimeter(gaps18to19sp[gaps18to19sp$gap_id==i,])
      
    }
    gaps18to19sp@data$ratio <- gaps18to19sp@data$area/gaps18to19sp@data$perimeter
    gaps18to19_circ <- gaps18to19sp[gaps18to19sp@data$ratio<5,]
    
    
    gaps19to20 <- ForestGapR::getForestGaps(d19to20,
                                            threshold = -5, size=c(10,10^4))
    gaps19to20sp <- ForestGapR::GapSPDF(gaps19to20)
    gaps19to20sp@data$area <- NA
    gaps19to20sp@data$perimeter <- NA
    for(i in 1:length(gaps19to20sp)){
      gaps19to20sp[gaps19to20sp$gap_id==i,"area"] <- raster::area(gaps19to20sp[gaps19to20sp$gap_id==i,])
      gaps19to20sp[gaps19to20sp$gap_id==i,"perimeter"] <- spatialEco::polyPerimeter(gaps19to20sp[gaps19to20sp$gap_id==i,])
      
    }
    gaps19to20sp@data$ratio <- gaps19to20sp@data$area/gaps19to20sp@data$perimeter
    gaps19to20_circ <- gaps19to20sp[gaps19to20sp@data$ratio<5,]
    save(gaps17to18sp, gaps18to19sp, gaps19to20sp, gaps17to18, gaps18to19, gaps19to20, file="PrelimGapLayers.RData")
  
  