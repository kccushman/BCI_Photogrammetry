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

  # Read topographic variables
  slope <- raster::raster("slopeBCI.tif")
  aspect <- raster::raster("aspectBCI.tif")
  curve <- raster::raster("pcurvatureBCI.tif")
  
  # Match topographic variables to canopy height rasters
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
  soil <- sp::spTransform(soil,raster::crs(chm17))

  # Make a vector of unique soil types
  soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                          PropGap = NA,
                          AreaSampled = NA)
  
  for(i in 1:dim(soilTypes)){
    areai <- raster::crop(dchm17to18,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- raster::crop(gaps17to18,soil[soil$SOIL==soilTypes$Soil[i],])
    values_areai <- raster::getValues(areai)
    values_gapsi <- raster::getValues(gapsi)
    soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
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
                            AreaSampled = NA)
    
    for(i in 1:dim(soilTypes)){
      areai <- raster::crop(dchm18to19,soil[soil$SOIL==soilTypes$Soil[i],])
      gapsi <- raster::crop(gaps18to19,soil[soil$SOIL==soilTypes$Soil[i],])
      values_areai <- raster::getValues(areai)
      values_gapsi <- raster::getValues(gapsi)
      soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
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
                            AreaSampled = NA)
    
    for(i in 1:dim(soilTypes)){
      areai <- raster::crop(dchm19to20,soil[soil$SOIL==soilTypes$Soil[i],])
      gapsi <- raster::crop(gaps19to20,soil[soil$SOIL==soilTypes$Soil[i],])
      values_areai <- raster::getValues(areai)
      values_gapsi <- raster::getValues(gapsi)
      soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
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
    names(soilTypes17to18) <- c("Soil","PropGap17to18","Area17to18")
    names(soilTypes18to19) <- c("Soil","PropGap18to19","Area18to19")
    names(soilTypes19to20) <- c("Soil","PropGap19to20","Area19to20")
    
    soilTypes <- merge(soilTypes17to18,soilTypes18to19)
    soilTypes <- merge(soilTypes, soilTypes19to20)
    
    soilTypes$Good <- ifelse(soilTypes$Area17to18>100 & soilTypes$Area18to19>100 & soilTypes$Area19to20>100,
                             T,
                             F)
    
    soilTypes <- soilTypes[order(-soilTypes$Good, -soilTypes$PropGap19to20),]
    
    plot(PropGap18to19~PropGap17to18,data=soilTypes[soilTypes$Good==T,],
         pch=20,
         ylab="2018-2019",
         xlab="2017-2018")
    abline(a=0,b=1)
    summary(lm(PropGap19to20~PropGap18to19,data=soilTypes[soilTypes$Good==T,]))
    
    
    plot(PropGap19to20~PropGap18to19,data=soilTypes[soilTypes$Good==T,],
         pch=20,
         ylab="2019-2020",
         xlab="2018-2019")
    abline(a=0,b=1)
    summary(lm(PropGap19to20~PropGap18to19,data=soilTypes[soilTypes$Good==T,]))
    
    
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


