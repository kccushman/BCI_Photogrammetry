#### Read data ####

  # Read grid info
  gridInfo <- read.csv("gridInfo_QAQC.csv")
  
  # Canopy height models for all years
  chm17 <- raster::raster("CHM_2017_QAQC.tif")
  chm18 <- raster::raster("CHM_2018_QAQC.tif")
  chm19 <- raster::raster("CHM_2019_QAQC.tif")
  chm20 <- raster::raster("CHM_2020_QAQC.tif")
  
  # Make canopy height change rasters
  dchm17to18 <- chm18-chm17
  dchm18to19 <- chm19-chm18
  dchm19to20 <- chm20-chm19
  
  # Gap rasters
  gaps17to18 <- raster::raster("newGaps17to18.tif")
  gaps18to19 <- raster::raster("newGaps18to19.tif")
  gaps19to20 <- raster::raster("newGaps19to20.tif")

  # Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
  
  # Read topographic variables
  slope <- raster::raster("slopeBCI.tif")
    slope <- raster::crop(slope,gaps17to18)
    slope <- raster::mask(slope, buffer)

  aspect <- raster::raster("aspectBCI.tif")
    aspect <- raster::crop(aspect,gaps17to18)
    aspect <- raster::mask(aspect, buffer)
  curve <- raster::raster("pcurvatureBCI.tif")
    curve <- raster::crop(curve,gaps17to18)
    curve <- raster::mask(curve, buffer)
    
  dem <- raster::raster("LidarDEM_BCI.tif")  
    dem <- raster::crop(dem,gaps17to18)
    dem <- raster::mask(dem, buffer)
    
    # make easting and northing rasters
    east <- aspect
    east@data@values <- cos(aspect@data@values*(pi/180))
    # remove values with aspect = 0 (no aspect because totally flat)
    east@data@values[aspect@data@values==0] <- NA
    
    north <- aspect
    north@data@values <- sin(aspect@data@values*(pi/180))
    # remove values with aspect = 0 (no aspect because totally flat)
    north@data@values[aspect@data@values==0] <- NA
  
  # Match topographic variables to canopy height rasters
  # NOTE: very preliminary based on one spatial scale so far
  # aspect is calculated counterclockwise from east
  raster::crs(slope) <- raster::crs(chm17)
  raster::crs(aspect) <- raster::crs(chm17)
  raster::crs(curve) <- raster::crs(chm17)
  slope <- raster::resample(slope,chm17)
  aspect <- raster::resample(aspect,chm17)
  curve <- raster::resample(curve,chm17)

#### 2017 to 2018 gaps ####

  # Mask NA cells in slope, curvature, and aspect plots
  aspect_mask <- aspect; aspect_mask[is.na(chm17) | is.na(chm18)] <- NA
  curvature_mask <- curve; curvature_mask[is.na(chm17) | is.na(chm18)] <- NA
  slope_mask <- slope; slope_mask[is.na(chm17) | is.na(chm18)] <- NA

# Only keep gap pixels
  aspect_gap <- aspect_mask; aspect_gap[is.na(gaps17to18)] <- NA
  curvature_gap <- curvature_mask; curvature_gap[is.na(gaps17to18)] <- NA
  slope_gap <- slope_mask; slope_gap[is.na(gaps17to18)] <- NA

# Compare proportion of aspect, curvature, slope in gaps compared to whole island
  # aspect
  aspect_all <- raster::getValues(aspect_mask)
  aspect_gap <- raster::getValues(aspect_gap)
  aspect_all_density <- density(aspect_all, na.rm = T, bw=5)
  aspect_all_gap <- density(aspect_gap, na.rm = T, bw=5 )
  
  plot(aspect_all_density,
       ylim=c(0,0.004))
  lines(aspect_all_gap, col="red")
  abline(v=c(0,90,180,270,360))
  
  # curvature
  curvature_all <- raster::getValues(curvature_mask)
  curvature_gap <- raster::getValues(curvature_gap)
  curvature_all_density <- density(curvature_all, na.rm = T, bw=0.00002)
  curvature_all_gap <- density(curvature_gap, na.rm = T, bw=0.00002)

  plot(curvature_all_density,
       xlim=c(-0.5,0.5))
  lines(curvature_all_gap, col="red")
  
  # slope
  slope_all <- raster::getValues(slope_mask)
  slope_gap <- raster::getValues(slope_gap)
  slope_all_density <- density(slope_all, na.rm = T)
  slope_all_gap <- density(slope_gap, na.rm = T)

  plot(slope_all_density)
  lines(slope_all_gap, col="red")
  
  # Compare with soil types
  soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
  soil <- sp::spTransform(soil,raster::crs(chm17))

  # Make a vector of unique soil types
  soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                          PropGap = NA,
                          nGap = NA,
                          AreaSampled = NA)
  
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::crop(dchm17to18,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- raster::crop(gaps17to18,soil[soil$SOIL==soilTypes$Soil[i],])
    values_areai <- raster::getValues(areai)
    values_gapsi <- raster::getValues(gapsi)
    soilTypes$PropGap[i] <- length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$nGap[i] <- length(unique(gapsi@data@values[!is.na(gapsi@data@values)]))/length(values_areai[!is.na(values_areai)])
    soilTypes$AreaSampled[i] <-length(values_areai[!is.na(values_areai)])/(10000)
  }
    soilTypes17to18 <- soilTypes

  # Forest age
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    
    # Make new age class
    age@data$AgeClass <- "Other"
    age@data$AgeClass[age@data$Mascaro_Co == "> 400"] <- "OldGrowth"
    age@data$AgeClass[age@data$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    
    oldGrowth <- raster::aggregate(age[age@data$AgeClass=="OldGrowth",])
    secondary <- raster::aggregate(age[age@data$AgeClass=="Secondary",])
    
    areaOld <- raster::crop(dchm17to18,oldGrowth)
    gapsOld <- raster::crop(gaps17to18,oldGrowth)
    values_areaOld <- raster::getValues(areaOld)
    values_gapsOld <- raster::getValues(gapsOld)
    length(values_gapsOld[!is.na(values_gapsOld)])/length(values_areaOld[!is.na(values_areaOld)])
    
    areaSec <- raster::crop(dchm17to18,secondary)
    gapsSec <- raster::crop(gaps17to18,secondary)
    values_areaSec <- raster::getValues(areaSec)
    values_gapsSec <- raster::getValues(gapsSec)
    length(values_gapsSec[!is.na(values_gapsSec)])/length(values_areaSec[!is.na(values_areaSec)])
    
#### 2018 to 2019 gaps ####
    
    # Mask NA cells in slope, curvature, and aspect plots
    aspect_mask <- aspect; aspect_mask[is.na(chm18) | is.na(chm19)] <- NA
    curvature_mask <- curve; curvature_mask[is.na(chm18) | is.na(chm19)] <- NA
    slope_mask <- slope; slope_mask[is.na(chm18) | is.na(chm19)] <- NA
    
    # Only keep gap pixels
    aspect_gap <- aspect_mask; aspect_gap[is.na(gaps18to19)] <- NA
    curvature_gap <- curvature_mask; curvature_gap[is.na(gaps18to19)] <- NA
    slope_gap <- slope_mask; slope_gap[is.na(gaps18to19)] <- NA
    
    # Compare proportion of aspect, curvature, slope in gaps compared to whole island
    # aspect
    aspect_all <- raster::getValues(aspect_mask)
    aspect_gap <- raster::getValues(aspect_gap)
    aspect_all_density <- density(aspect_all, na.rm = T, bw=5)
    aspect_all_gap <- density(aspect_gap, na.rm = T, bw=5 )
    
    plot(aspect_all_density,
         ylim=c(0,0.004))
    lines(aspect_all_gap, col="red")
    
    # curvature
    curvature_all <- raster::getValues(curvature_mask)
    curvature_gap <- raster::getValues(curvature_gap)
    curvature_all_density <- density(curvature_all, na.rm = T, bw=0.002)
    curvature_all_gap <- density(curvature_gap, na.rm = T, bw=0.002)
    
    plot(curvature_all_density,
         xlim=c(-0.5,0.5))
    lines(curvature_all_gap, col="red")
    
    # slope
    slope_all <- raster::getValues(slope_mask)
    slope_gap <- raster::getValues(slope_gap)
    slope_all_density <- density(slope_all, na.rm = T)
    slope_all_gap <- density(slope_gap, na.rm = T)
    
    plot(slope_all_density)
    lines(slope_all_gap, col="red")
    
    # Compare with soil types
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(chm18))
    
    # Make a vector of unique soil types
    soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                            PropGap = NA,
                            nGap = NA,
                            AreaSampled = NA)
    
    for(i in 1:dim(soilTypes)[1]){
      areai <- raster::crop(dchm18to19,soil[soil$SOIL==soilTypes$Soil[i],])
      gapsi <- raster::crop(gaps18to19,soil[soil$SOIL==soilTypes$Soil[i],])
      values_areai <- raster::getValues(areai)
      values_gapsi <- raster::getValues(gapsi)
      soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
      soilTypes$nGap[i] <- length(unique(gapsi@data@values[!is.na(gapsi@data@values)]))/length(values_areai[!is.na(values_areai)])
      soilTypes$AreaSampled[i] <-length(values_areai[!is.na(values_areai)])/(10000)
    }
    soilTypes18to19 <- soilTypes
    
    # Forest age
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    
    # Make new age class
    age@data$AgeClass <- "Other"
    age@data$AgeClass[age@data$Mascaro_Co == "> 400"] <- "OldGrowth"
    age@data$AgeClass[age@data$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    
    oldGrowth <- raster::aggregate(age[age@data$AgeClass=="OldGrowth",])
    secondary <- raster::aggregate(age[age@data$AgeClass=="Secondary",])
    
    areaOld <- raster::crop(dchm18to19,oldGrowth)
    gapsOld <- raster::crop(gaps18to19,oldGrowth)
    values_areaOld <- raster::getValues(areaOld)
    values_gapsOld <- raster::getValues(gapsOld)
    length(values_gapsOld[!is.na(values_gapsOld)])/length(values_areaOld[!is.na(values_areaOld)])
    
    areaSec <- raster::crop(dchm18to19,secondary)
    gapsSec <- raster::crop(gaps18to19,secondary)
    values_areaSec <- raster::getValues(areaSec)
    values_gapsSec <- raster::getValues(gapsSec)
    length(values_gapsSec[!is.na(values_gapsSec)])/length(values_areaSec[!is.na(values_areaSec)]) 
    
#### 2019 to 2020 gaps ####
    
    # Mask NA cells in slope, curvature, and aspect plots
    aspect_mask <- aspect; aspect_mask[is.na(chm19) | is.na(chm20)] <- NA
    curvature_mask <- curve; curvature_mask[is.na(chm19) | is.na(chm20)] <- NA
    slope_mask <- slope; slope_mask[is.na(chm19) | is.na(chm20)] <- NA
    
    # Only keep gap pixels
    aspect_gap <- aspect_mask; aspect_gap[is.na(gaps19to20)] <- NA
    curvature_gap <- curvature_mask; curvature_gap[is.na(gaps19to20)] <- NA
    slope_gap <- slope_mask; slope_gap[is.na(gaps19to20)] <- NA
    
    # Compare proportion of aspect, curvature, slope in gaps compared to whole island
    # aspect
    aspect_all <- raster::getValues(aspect_mask)
    aspect_gap <- raster::getValues(aspect_gap)
    aspect_all_density <- density(aspect_all, na.rm = T, bw=5)
    aspect_all_gap <- density(aspect_gap, na.rm = T, bw=5 )
    
    plot(aspect_all_density,
         ylim=c(0,0.004))
    lines(aspect_all_gap, col="red")
    
    # curvature
    curvature_all <- raster::getValues(curvature_mask)
    curvature_gap <- raster::getValues(curvature_gap)
    curvature_all_density <- density(curvature_all, na.rm = T, bw=0.002)
    curvature_all_gap <- density(curvature_gap, na.rm = T, bw=0.002)
    
    plot(curvature_all_density,
         xlim=c(-0.5,0.5))
    lines(curvature_all_gap, col="red")
    
    # slope
    slope_all <- raster::getValues(slope_mask)
    slope_gap <- raster::getValues(slope_gap)
    slope_all_density <- density(slope_all, na.rm = T)
    slope_all_gap <- density(slope_gap, na.rm = T)
    
    plot(slope_all_density)
    lines(slope_all_gap, col="red")
    
    # Compare with soil types
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(chm19))
    
    # Make a vector of unique soil types
    soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                            PropGap = NA,
                            nGap = NA,
                            AreaSampled = NA)
    
    for(i in 1:dim(soilTypes)[1]){
      areai <- raster::crop(dchm19to20,soil[soil$SOIL==soilTypes$Soil[i],])
      gapsi <- raster::crop(gaps19to20,soil[soil$SOIL==soilTypes$Soil[i],])
      values_areai <- raster::getValues(areai)
      values_gapsi <- raster::getValues(gapsi)
      soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
      soilTypes$nGap[i] <- length(unique(gapsi@data@values[!is.na(gapsi@data@values)]))/length(values_areai[!is.na(values_areai)])
      soilTypes$AreaSampled[i] <-length(values_areai[!is.na(values_areai)])/(10000)
    }
    soilTypes19to20 <- soilTypes
    
    # Forest age
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    
    # Make new age class
    age@data$AgeClass <- "Other"
    age@data$AgeClass[age@data$Mascaro_Co == "> 400"] <- "OldGrowth"
    age@data$AgeClass[age@data$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    
    oldGrowth <- raster::aggregate(age[age@data$AgeClass=="OldGrowth",])
    secondary <- raster::aggregate(age[age@data$AgeClass=="Secondary",])
    
    areaOld <- raster::crop(dchm19to20,oldGrowth)
    gapsOld <- raster::crop(gaps19to20,oldGrowth)
    values_areaOld <- raster::getValues(areaOld)
    values_gapsOld <- raster::getValues(gapsOld)
    length(values_gapsOld[!is.na(values_gapsOld)])/length(values_areaOld[!is.na(values_areaOld)])
    
    areaSec <- raster::crop(dchm19to20,secondary)
    gapsSec <- raster::crop(gaps19to20,secondary)
    values_areaSec <- raster::getValues(areaSec)
    values_gapsSec <- raster::getValues(gapsSec)
    length(values_gapsSec[!is.na(values_gapsSec)])/length(values_areaSec[!is.na(values_areaSec)])
    
#### Compare soil rates across years ####
    names(soilTypes17to18) <- c("Soil","PropGap17to18","nGap17to18","Area17to18")
    names(soilTypes18to19) <- c("Soil","PropGap18to19","nGap18to19","Area18to19")
    names(soilTypes19to20) <- c("Soil","PropGap19to20","nGap19to20","Area19to20")
    
    soilTypes <- merge(soilTypes17to18,soilTypes18to19)
    soilTypes <- merge(soilTypes, soilTypes19to20)
    
    soilTypes$Good <- ifelse(soilTypes$Area17to18>50 & soilTypes$Area18to19>50 & soilTypes$Area19to20>50,
                             T,
                             F)
    
    soilTypes <- soilTypes[order(-soilTypes$Good, -soilTypes$PropGap19to20),]
    
    
    plot(nGap17to18~PropGap17to18, data=soilTypes[soilTypes$Good==T,],
         pch = 20,
         xlab = "Proportion of area in new gaps",
         ylab = "New gaps/area",
         main = "2017 to 2018")
    
    plot(nGap18to19~PropGap18to19, data=soilTypes[soilTypes$Good==T,],
         pch = 20,
         xlab = "Proportion of area in new gaps",
         ylab = "New gaps/area",
         main = "2018 to 2019")
    
    plot(nGap19to20~PropGap19to20, data=soilTypes[soilTypes$Good==T,],
         pch = 20,
         xlab = "Proportion of area in new gaps",
         ylab = "New gaps/area",
         main = "2019 to 2020")
    
    
    par(mfrow=c(1,2),
        mar=c(4,4,1,1))
    plot(PropGap18to19~PropGap17to18, data=soilTypes[soilTypes$Good==T,],
         pch=20,
         xlim=c(0.01,0.04),
         ylim=c(0.01,0.04),
         ylab="2018-2019",
         xlab="2017-2018")
    abline(a=0,b=1)
    summary(lm(PropGap18to19~PropGap17to18,data=soilTypes[soilTypes$Good==T,]))
    summary(lm(nGap18to19~nGap17to18,data=soilTypes[soilTypes$Good==T,]))
    
    
    plot(PropGap19to20~PropGap18to19,data=soilTypes[soilTypes$Good==T,],
         pch=20,
         xlim=c(0.01,0.04),
         ylim=c(0.01,0.04),
         ylab="2019-2020",
         xlab="2018-2019")
    abline(a=0,b=1)
    summary(lm(PropGap19to20~PropGap18to19,data=soilTypes[soilTypes$Good==T,]))
    summary(lm(nGap19to20~nGap18to19,data=soilTypes[soilTypes$Good==T,]))
    
    plot(PropGap19to20~PropGap17to18, data=soilTypes[soilTypes$Good==T,],
         pch=20,
         ylab="2019-2020",
         xlab="2017-2018")
    abline(a=0,b=1)
    summary(lm(PropGap19to20~PropGap17to18,data=soilTypes[soilTypes$Good==T,]))
    summary(lm(nGap19to20~nGap17to18,data=soilTypes[soilTypes$Good==T,]))
    
    
#### Average gap formation rates in space ####
    soilTypes$avgRate <- (1/3)*(soilTypes$PropGap17to18+soilTypes$PropGap18to19+soilTypes$PropGap19to20)
    soilTypes <- soilTypes[order(-soilTypes$Good,-soilTypes$avgRate),]
    
    manual.col <- colorRampPalette(c("#57471d","#b5b5b5"))
    soilTypes$col <- manual.col(length(unique(soilTypes$avgRate)))
    
    soil$plotCol <- NA
    for(i in 1:dim(soilTypes)[1]){
      if(soilTypes$Good[i]==T){
        soil[soil$SOIL==soilTypes$Soil[i],"plotCol"] <- soilTypes$col[i]
      }
    }
    
    par(mfrow=c(1,1),mar=c(1,1,1,1))
    raster::plot(soil,
                 col=soil$plotCol)
    
    
#### Correlations between average formation rates and topography ####
    
  # Slope
    
  soilTypes$avgSlope <- NA
    for(i in 1:dim(soilTypes)[1]){
      slopei <- raster::crop(slope,soil[soil$SOIL==soilTypes$Soil[i],])
      values_slopei <- raster::getValues(slopei)
      soilTypes$avgSlope[i] <- mean(values_slopei)
    }
    
    par(mfrow=c(1,1),mar=c(4,4,1,1))
    plot(avgRate~avgSlope, data = soilTypes[soilTypes$Good==T,],
         pch=20,
         xlab="Average slope", ylab="Average gap formation rate")
    summary(lm(avgRate~avgSlope, data = soilTypes[soilTypes$Good==T,]))
  
  # Aspect
  # East-west (value of 1 is facing east, value of -1 is facing due west)
  # North-south (value of 1 is facing north, value of -1 is facing due south)
  soilTypes$avgEastWest <- NA
  soilTypes$avgNorthSouth <- NA
  
    for(i in 1:dim(soilTypes)[1]){
      aspecti <- raster::crop(aspect,soil[soil$SOIL==soilTypes$Soil[i],])
      
      valuesi <- raster::getValues(aspecti)
      valuesi <- cos(valuesi*(pi/180))
      soilTypes$avgEastWest[i] <- mean(valuesi)
      
      valuesi <- raster::getValues(aspecti)
      valuesi <- sin(valuesi*(pi/180))
      soilTypes$avgNorthSouth[i] <- mean(valuesi)
      
    }
    
    par(mfrow=c(1,2),mar=c(4,4,1,1))
    plot(avgRate~avgEastWest, data = soilTypes[soilTypes$Good==T,],
         pch=20,
         xlab="Average East-West aspect", ylab="Average gap formation rate")
    summary(lm(avgRate~avgEastWest, data = soilTypes[soilTypes$Good==T,]))
    
    plot(avgRate~avgNorthSouth, data = soilTypes[soilTypes$Good==T,],
         pch=20,
         xlab="Average North-South aspect", ylab=NA)
    summary(lm(avgRate~avgNorthSouth, data = soilTypes[soilTypes$Good==T,]))


#### Calculate slope and easting for each gap ####
  IDs18 <- raster::unique(gaps17to18)
  IDs19 <- raster::unique(gaps18to19)
  IDs20 <- raster::unique(gaps19to20)
  
  gapTopo <- data.frame(gapID = c(IDs18, IDs19, IDs20),
                        year = c(rep(2018, length(IDs18)),
                                 rep(2019, length(IDs19)),
                                 rep(2020, length(IDs20))),
                        gapArea = NA,
                        slope = NA,
                        easting = NA,
                        northing = NA)
  
  gapVals17to18 <- raster::values(gaps17to18)
  gapVals18to19 <- raster::values(gaps18to19)
  gapVals19to20 <- raster::values(gaps19to20)
  slopeVals <- raster::values(slope)
  eastVals <- raster::values (east)
  northVals <- raster::values(north)
  
  
  for(i in 1:dim(gapTopo)[1]){
    if(gapTopo$year[i]==2018){
      # calculate gap area
      gapTopo$gapArea[i] <- length(gapVals17to18[gapVals17to18==gapTopo$gapID[i]
                                                 & !is.na(gapVals17to18)])
      
      # calculate average slope in gap
      gapTopo$slope[i] <- mean(slopeVals[gapVals17to18==gapTopo$gapID[i]
                                             & !is.na(gapVals17to18)])
      
      # calculate average easting in gap
      gapTopo$easting[i] <- mean(eastVals[gapVals17to18==gapTopo$gapID[i]
                                               & !is.na(gapVals17to18)])
      
      # calcuate average northing in gap
      gapTopo$northing[i] <- mean(northVals[gapVals17to18==gapTopo$gapID[i]
                                                & !is.na(gapVals17to18)])
    }
    
    if(gapTopo$year[i]==2019){
      # calculate gap area
      gapTopo$gapArea[i] <- length(gapVals18to19[gapVals18to19==gapTopo$gapID[i]
                                                 & !is.na(gapVals18to19)])
      
      # calculate average slope in gap
      gapTopo$slope[i] <- mean(slopeVals[gapVals18to19==gapTopo$gapID[i]
                                         & !is.na(gapVals18to19)])
      
      # calculate average easting in gap
      gapTopo$easting[i] <- mean(eastVals[gapVals18to19==gapTopo$gapID[i]
                                          & !is.na(gapVals18to19)])
      
      # calcuate average northing in gap
      gapTopo$northing[i] <- mean(northVals[gapVals18to19==gapTopo$gapID[i]
                                            & !is.na(gapVals18to19)])
    }
    
    if(gapTopo$year[i]==2020){
      # calculate gap area
      gapTopo$gapArea[i] <- length(gapVals19to20[gapVals19to20==gapTopo$gapID[i]
                                                 & !is.na(gapVals19to20)])
      
      # calculate average slope in gap
      gapTopo$slope[i] <- mean(slopeVals[gapVals19to20==gapTopo$gapID[i]
                                         & !is.na(gapVals19to20)])
      
      # calculate average easting in gap
      gapTopo$easting[i] <- mean(eastVals[gapVals19to20==gapTopo$gapID[i]
                                          & !is.na(gapVals19to20)])
      
      # calcuate average northing in gap
      gapTopo$northing[i] <- mean(northVals[gapVals19to20==gapTopo$gapID[i]
                                            & !is.na(gapVals19to20)])
    }
    print(i)
    
  }
    
  write.csv(gapTopo, file="gapTopoStats.csv", row.names = F)
  
  
    
#### Set up functions for MCMC routine ####
  # Load data
  gapData <- read.csv("gapTopoStats.csv")
  
  
  # library(Rlab)
  require(VGAM)
  require(poweRlaw)
  
  # define minimum gap size
  minGap <- 10

  calcul_lambda <- function(theta, var)
  {  
    pred =1+exp(theta[1]+theta[2]*var)  #Univariate model cf (eq.8) of paper
    return(pred)
  }
  
  fct <- function(x)
    return(dpldis(x[1],minGap,x[2],log=T))
  
  loglikelihood_pareto= function(data, xmin =minGap ,theta,var){
    lambda_pred = calcul_lambda(theta,var)
    d_sim = cbind(data,lambda_pred)
    temp = apply(d_sim,1,fct)
    LL= as.matrix(temp)
    ll_area_tot <- colSums(LL)
    return(ll_area_tot)
  }

#### Do MCMC univariate: slope ####
  # Make a vector of gap areas (must convert to integer if not already)
  d <- gapData[gapData$gapArea>=minGap & !is.na(gapData$slope),"gapArea"]
  
  # Transform slope data
  p <- sqrt(gapData[gapData$gapArea>=minGap & !is.na(gapData$slope),"slope"])
  var <- c(scale(p,center= TRUE,scale= TRUE))
  
  #________________________________
  #       START with model paramter (THETA)
  # At this step all parameter must be adjusted, 
  #_________________________________ 
  theta <- c(0.5,0.4)
  nb_param <- length(theta)
  step <- c(0.05,0.05)
  count= 0
  nb_it <- 1000 #number of iterations
  nb_sample <- 2    
  m_prior <- 0.5
  sd_prior <- 500
  theta <- matrix(theta,nrow=nb_param,ncol = 2,byrow=F)
  v_e1 <- var
  save_res  = "TRUE" # TRUE or FALSE  # If we want to save results

  #_________________________________ 
  #      COMPUTE LILKELIHOO_PARETO for the fist elements of theta and var     
  #_________________________________ 
  llp = loglikelihood_pareto(data = d  ,xmin = minGap ,theta,var)
  llp_res=c(llp, rep(0, nb_it))
  
  #________________________________
  #      Prepare OUTPUT             
  #________________________________
  accept <- rep(0,nb_param)
  res <- array(dim=c(nb_it+1,nb_param,nb_sample))
  res[1,,] <- theta
  
  #________________________________
  #   CHAINE DE Metropolis_Hastings      
  #________________________________
  for(compt in 1:nb_it)
  {
    
    ordre <- sample(1 : nb_param)
    
    for (i in ordre)
    {
      
      theta_new <- theta #candidat
      theta_new[i,] <- rnorm(2,theta[i,],step[i])
      
      llp_new = loglikelihood_pareto(d,xmin=minGap,theta=theta_new, var = var)
      
      logr = llp_new - llp +dnorm(theta[i,],theta_new[i,],step[i],log=T) - 
        dnorm(theta_new[i,],theta[i,],step[i],log=T) +
        dnorm(theta_new[i,], m_prior, sd_prior, log=T) - 
        dnorm(theta[i,], m_prior, sd_prior, log=T)
      r <- exp(logr)
      
      nb_ok <- which(!is.na(r))
      r_ok <- na.omit(r)
      
      accepte <- nb_ok[(runif(length(r_ok),min=0,max=1) < r_ok)]
      llp_res[compt+1]=llp_new
      theta[,accepte] <- theta_new[,accepte]
      llp <- llp_new
    }
    
    res[compt+1,,] <- theta
    print(compt)
  }
  
  #________________________________
  #      Plot results 
  #________________________________
  windows(10,10)
  par(mfrow=c(2,2),mar=c(3,2,2,2))
  
  theta_est = res[,,1]
  plot(theta_est[100:1000,1],col="blue",type = "l",main ="Theta 0" )
  plot(theta_est[100:1000,2],col="magenta",type = "l",main = "Theta 1")
  theta.1 = theta_est[100:1000,1]
  plot(density(theta.1))
  theta.2 = theta_est[100:1000,2]
  plot(density(theta.2),
       main = "Slope effect")
  abline(v=0,lty=2)
  
  #_____________________________
  # Save result
  #_____________________________
  if (save_res == TRUE) {
    write.csv(res,file="slope_MCMC.csv")
  }
  


#### Do MCMC univariate: easting ####
  # Make a vector of gap areas (must convert to integer if not already)
  d <- gapData[gapData$gapArea>=minGap & !is.na(gapData$easting),"gapArea"]
  
  # Transform slope data
  p <- gapData[gapData$gapArea>=minGap & !is.na(gapData$easting),"easting"]
  var <- c(scale(p,center= TRUE,scale= TRUE))
  
  #________________________________
  #       START with model paramter (THETA)
  # At this step all parameter must be adjusted, 
  #_________________________________ 
  theta <- c(0.5,0.4)
  nb_param <- length(theta)
  step <- c(0.05,0.05)
  count= 0
  nb_it <- 1000 #number of iterations
  nb_sample <- 2    
  m_prior <- 0.5
  sd_prior <- 500
  theta <- matrix(theta,nrow=nb_param,ncol = 2,byrow=F)
  v_e1 <- var
  save_res  = "TRUE" # TRUE or FALSE  # If we want to save results
  
  #_________________________________ 
  #      COMPUTE LILKELIHOO_PARETO for the fist elements of theta and var     
  #_________________________________ 
  llp = loglikelihood_pareto(data = d  ,xmin = minGap ,theta,var)
  llp_res=c(llp, rep(0, nb_it))
  
  #________________________________
  #      Prepare OUTPUT             
  #________________________________
  accept <- rep(0,nb_param)
  res <- array(dim=c(nb_it+1,nb_param,nb_sample))
  res[1,,] <- theta
  
  #________________________________
  #   CHAINE DE Metropolis_Hastings      
  #________________________________
  for(compt in 1:nb_it)
  {
    
    ordre <- sample(1 : nb_param)
    
    for (i in ordre)
    {
      
      theta_new <- theta #candidat
      theta_new[i,] <- rnorm(2,theta[i,],step[i])
      
      llp_new = loglikelihood_pareto(d,xmin=minGap,theta=theta_new, var = var)
      
      logr = llp_new - llp +dnorm(theta[i,],theta_new[i,],step[i],log=T) - 
        dnorm(theta_new[i,],theta[i,],step[i],log=T) +
        dnorm(theta_new[i,], m_prior, sd_prior, log=T) - 
        dnorm(theta[i,], m_prior, sd_prior, log=T)
      r <- exp(logr)
      
      nb_ok <- which(!is.na(r))
      r_ok <- na.omit(r)
      
      accepte <- nb_ok[(runif(length(r_ok),min=0,max=1) < r_ok)]
      llp_res[compt+1]=llp_new
      theta[,accepte] <- theta_new[,accepte]
      llp <- llp_new
    }
    
    res[compt+1,,] <- theta
    print(compt)
  }
  
  #________________________________
  #      Plot results 
  #________________________________
  windows(10,10)
  par(mfrow=c(2,2),mar=c(3,2,2,2))
  
  theta_est = res[,,1]
  plot(theta_est[,1],col="blue",type = "l",main ="Theta 0" )
  plot(theta_est[,2],col="magenta",type = "l",main = "Theta 1")
  theta.1 = theta_est[100:1000,1]
  plot(density(theta.1))
  theta.2 = theta_est[100:1000,2]
  plot(density(theta.2),
       main = "Easting effect")
  abline(v=0,lty=2)
  
  #_____________________________
  # Save result
  #_____________________________
  if (save_res == TRUE) {
    write.csv(res,file="easting_MCMC.csv")
  }
  
  
  
  
#### Do MCMC univariate: northing ####
  # Make a vector of gap areas (must convert to integer if not already)
  d <- gapData[gapData$gapArea>=minGap & !is.na(gapData$northing),"gapArea"]
  
  # Transform slope data
  p <- gapData[gapData$gapArea>=minGap & !is.na(gapData$northing),"northing"]
  var <- c(scale(p,center= TRUE,scale= TRUE))
  
  #________________________________
  #       START with model paramter (THETA)
  # At this step all parameter must be adjusted, 
  #_________________________________ 
  theta <- c(0.5,0.4)
  nb_param <- length(theta)
  step <- c(0.05,0.05)
  count= 0
  nb_it <- 1000 #number of iterations
  nb_sample <- 2    
  m_prior <- 0.5
  sd_prior <- 500
  theta <- matrix(theta,nrow=nb_param,ncol = 2,byrow=F)
  v_e1 <- var
  save_res  = "TRUE" # TRUE or FALSE  # If we want to save results
  
  #_________________________________ 
  #      COMPUTE LILKELIHOO_PARETO for the fist elements of theta and var     
  #_________________________________ 
  llp = loglikelihood_pareto(data = d  ,xmin = minGap ,theta,var)
  llp_res=c(llp, rep(0, nb_it))
  
  #________________________________
  #      Prepare OUTPUT             
  #________________________________
  accept <- rep(0,nb_param)
  res <- array(dim=c(nb_it+1,nb_param,nb_sample))
  res[1,,] <- theta
  
  #________________________________
  #   CHAINE DE Metropolis_Hastings      
  #________________________________
  for(compt in 1:nb_it)
  {
    
    ordre <- sample(1 : nb_param)
    
    for (i in ordre)
    {
      
      theta_new <- theta #candidat
      theta_new[i,] <- rnorm(2,theta[i,],step[i])
      
      llp_new = loglikelihood_pareto(d,xmin=minGap,theta=theta_new, var = var)
      
      logr = llp_new - llp +dnorm(theta[i,],theta_new[i,],step[i],log=T) - 
        dnorm(theta_new[i,],theta[i,],step[i],log=T) +
        dnorm(theta_new[i,], m_prior, sd_prior, log=T) - 
        dnorm(theta[i,], m_prior, sd_prior, log=T)
      r <- exp(logr)
      
      nb_ok <- which(!is.na(r))
      r_ok <- na.omit(r)
      
      accepte <- nb_ok[(runif(length(r_ok),min=0,max=1) < r_ok)]
      llp_res[compt+1]=llp_new
      theta[,accepte] <- theta_new[,accepte]
      llp <- llp_new
    }
    
    res[compt+1,,] <- theta
    print(compt)
  }
  
  #________________________________
  #      Plot results 
  #________________________________
  windows(10,10)
  par(mfrow=c(2,2),mar=c(3,2,2,2))
  
  theta_est = res[,,1]
  plot(theta_est[,1],col="blue",type = "l",main ="Theta 0" )
  plot(theta_est[,2],col="magenta",type = "l",main = "Theta 1")
  theta.1 = theta_est[100:1000,1]
  plot(density(theta.1))
  theta.2 = theta_est[100:1000,2]
  plot(density(theta.2),
       main = "Northing effect")
  abline(v=0,lty=2)
  
  #_____________________________
  # Save result
  #_____________________________
  if (save_res == TRUE) {
    write.csv(res,file="northing_MCMC.csv")
  }
  
#### GRAVEYARD ####    
# Get soil type grap frequency from lidar
  soilTypes$PropGapStanding <- NA
  dsm <- raster::raster("~/Desktop/Postdoc/BCI lidar 2009/BCI_DSM_2009.tif")
  dsm <- raster::crop(dsm, raster::extent(dCHM))
  dsm_mask <- dsm; dsm_mask[is.na(dCHM_binary)] <- NA
  dem <- raster::raster("~/Desktop/Postdoc/BCI lidar 2009/DEM_clipMCH.tif")
  dem <- raster::resample(dem, dsm)
  dem <- raster::crop(dem, raster::extent(dCHM))
  chm <- dsm-dem
  chm[is.na(dCHM_binary)] <- NA
  
  for(i in 1:dim(soilTypes)[1]){
    soili <- raster::crop(chm,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- ForestGapR::getForestGaps(soili, threshold = 2)
    
    valuesi <- raster::getValues(soili)
    valuesi <- valuesi[!is.na(valuesi)]
    
    gapsArea<- sum(ForestGapR::GapStats(gapsi,soili)$gap_area)
    
    soilTypes$PropGapStanding[i] <- gapsArea/length(valuesi)
  }

  plot(PropGap~PropGapStanding, data = soilTypes, pch=20,
       xlab="Proportion standing gaps in 2009",
       ylab = "Proportion new gaps 2018-2019")
  summary(lm(PropGap~PropGapStanding, data = soilTypes))


