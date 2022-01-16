#### Read data ####

  # Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs") 
  
#### Cloud mask 2015 ####

  # Read albedo (because clouds are bright)
  albedo15 <- raster::raster("Data_QAQC/CloudMasks/Min_Albedo_2015.tif")
  albedo15 <- raster::mask(albedo15,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  redblu15 <- raster::raster("Data_QAQC/CloudMasks/Max_RBdiff_2015.tif")
  redblu15 <- raster::mask(redblu15,buffer)
  
  
  thresh_albedo <- 475
  thresh_redblu <- 20
  mask15 <- albedo15
  mask15@data@values[albedo15@data@values >= thresh_albedo & redblu15@data@values <= thresh_redblu] <- 0
  mask15@data@values[albedo15@data@values < thresh_albedo | redblu15@data@values > thresh_redblu] <- 1
  raster::plot(mask15)
  
  raster::writeRaster(mask15, "Data_QAQC/CloudMasks/CloudMask_2015.tif")
    
#### Cloud mask 2018 ####
  
  # Read albedo (because clouds are bright)
  albedo18 <- raster::raster("Data_QAQC/CloudMasks/Min_Albedo_2018.tif")
  albedo18 <- raster::mask(albedo18,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  redblu18 <- raster::raster("Data_QAQC/CloudMasks/Max_RBdiff_2018.tif")
  redblu18 <- raster::mask(redblu18,buffer)
  
  
  thresh_albedo <- 375
  thresh_redblu <- 20
  mask18 <- albedo18
  mask18@data@values[albedo18@data@values >= thresh_albedo & redblu18@data@values <= thresh_redblu] <- 0
  mask18@data@values[albedo18@data@values < thresh_albedo | redblu18@data@values > thresh_redblu] <- 1
  raster::plot(mask18)
  
  raster::writeRaster(mask18, "Data_QAQC/CloudMasks/CloudMask_2018.tif")      
    
#### Cloud mask 2020 ####
  
  # Read albedo (because clouds are bright)
  albedo20 <- raster::raster("Data_QAQC/CloudMasks/Min_Albedo_2020.tif")
  albedo20 <- raster::mask(albedo20,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  redblu20 <- raster::raster("Data_QAQC/CloudMasks/Max_RBdiff_2020.tif")
  redblu20 <- raster::mask(redblu20,buffer)
  
  
  thresh_albedo <- 375
  thresh_redblu <- 15
  mask20 <- albedo20
  mask20@data@values[albedo20@data@values >= thresh_albedo & redblu20@data@values <= thresh_redblu] <- 0
  mask20@data@values[albedo20@data@values < thresh_albedo | redblu20@data@values > thresh_redblu] <- 1
  raster::plot(mask20)
  
  raster::writeRaster(mask20, "Data_QAQC/CloudMasks/CloudMask_2020.tif")      
  
#### Data QAQC mask 2015 ####
  
  # Aggregate criteria 1: range of heights within 0.2 m pixels
  qaqc15a <- raster::raster("Data_QAQC/htRangeRaster_2015.tif")
  qaqc15a_ag <- raster::aggregate(qaqc15a, fact = 50, fun=median)
  raster::plot(qaqc15a_ag,
               breaks=c(0,1.5,50),
               col=c("grey","red"))
  
  # Aggregate criteria 2: number of non-empty cells
  qaqc15b <- raster::raster("Data_QAQC/DSM_noGapFill_2015.tif")
  
  # Make binary (1 if value, 0 if NA)
  qaqc15b_bin <- raster::values(qaqc15b)
  qaqc15b_bin[!is.na(qaqc15b_bin)] <- 1
  qaqc15b_bin[is.na(qaqc15b_bin)] <- 0
  raster::values(qaqc15b) <- qaqc15b_bin
  
  # Aggregate by adding non-NA values
  qaqc15b_ag <- raster::aggregate(qaqc15b, fact = 10, fun=sum)
  raster::plot(qaqc15b_ag,
               breaks=c(0,90,100),
               col=c("red","grey"))
  
  # Make binary map
  thresh_1 <- 1.5
  thresh_2 <- 90
  
  mask15 <- qaqc15a_ag
  # Set mask values
  mask_vals <- raster::values(mask15)
  vals_a <- raster::values(qaqc15a_ag)
  vals_b <- raster::values(qaqc15b_ag)
  
  # Bad pixels == 0
  mask_vals[vals_a >= thresh_1 | vals_b <= thresh_2] <- 0
  # Good pixels == 1
  mask_vals[vals_a < thresh_1 & vals_b > thresh_2] <- 1
  raster::values(mask15) <- mask_vals
  
  raster::plot(mask15)
  raster::writeRaster(mask15, "Data_QAQC/QAQCMask_2015.tif")
    
#### Data QAQC mask 2018 ####
    
    # Aggregate criteria 1: range of heights within 0.2 m pixels
    qaqc18a <- raster::raster("Data_QAQC/htRangeRaster_2018.tif")
    qaqc18a_ag <- raster::aggregate(qaqc18a, fact = 50, fun=median)
    raster::plot(qaqc18a_ag,
                 breaks=c(0,1.5,50),
                 col=c("grey","red"))
    
    # Aggregate criteria 2: number of non-empty cells
    qaqc18b <- raster::raster("Data_QAQC/DSM_noGapFill_2018.tif")
    
    # Make binary (1 if value, 0 if NA)
    qaqc18b_bin <- raster::values(qaqc18b)
    qaqc18b_bin[!is.na(qaqc18b_bin)] <- 1
    qaqc18b_bin[is.na(qaqc18b_bin)] <- 0
    raster::values(qaqc18b) <- qaqc18b_bin
    
    # Aggregate by adding non-NA values
    qaqc18b_ag <- raster::aggregate(qaqc18b, fact = 10, fun=sum)
    raster::plot(qaqc18b_ag,
                 breaks=c(0,90,100),
                 col=c("red","grey"))
    
    # Make binary map
    thresh_1 <- 1.5
    thresh_2 <- 90
    
    mask18 <- qaqc18a_ag
    # Set mask values
    mask_vals <- raster::values(mask18)
    vals_a <- raster::values(qaqc18a_ag)
    vals_b <- raster::values(qaqc18b_ag)
    
    # Bad pixels == 0
    mask_vals[vals_a >= thresh_1 | vals_b <= thresh_2] <- 0
    # Good pixels == 1
    mask_vals[vals_a < thresh_1 & vals_b > thresh_2] <- 1
    raster::values(mask18) <- mask_vals
    
    raster::plot(mask18)
    raster::writeRaster(mask18, "Data_QAQC/QAQCMask_2018.tif")
    
#### Data QAQC mask 2020 ####
    
    # Aggregate criteria 1: range of heights within 0.2 m pixels
    qaqc20a <- raster::raster("Data_QAQC/htRangeRaster_2020.tif")
    qaqc20a_ag <- raster::aggregate(qaqc20a, fact = 50, fun=median)
    raster::plot(qaqc20a_ag,
                 breaks=c(0,1.5,50),
                 col=c("grey","red"))
    
    # Aggregate criteria 2: number of non-empty cells
    qaqc20b <- raster::raster("Data_QAQC/DSM_noGapFill_2020.tif")
    
    # Make binary (1 if value, 0 if NA)
    qaqc20b_bin <- raster::values(qaqc20b)
    qaqc20b_bin[!is.na(qaqc20b_bin)] <- 1
    qaqc20b_bin[is.na(qaqc20b_bin)] <- 0
    raster::values(qaqc20b) <- qaqc20b_bin
    
    # Aggregate by adding non-NA values
    qaqc20b_ag <- raster::aggregate(qaqc20b, fact = 10, fun=sum)
    raster::plot(qaqc20b_ag,
                 breaks=c(0,90,100),
                 col=c("red","grey"))
    
    # Make binary map
    thresh_1 <- 1.5
    thresh_2 <- 90
    
    mask20 <- qaqc20a_ag
    # Set mask values
    mask_vals <- raster::values(mask20)
    vals_a <- raster::values(qaqc20a_ag)
    vals_b <- raster::values(qaqc20b_ag)
    
    # Bad pixels == 0
    mask_vals[vals_a >= thresh_1 | vals_b <= thresh_2] <- 0
    # Good pixels == 1
    mask_vals[vals_a < thresh_1 & vals_b > thresh_2] <- 1
    raster::values(mask20) <- mask_vals
    
    raster::plot(mask20)
    raster::writeRaster(mask20, "Data_QAQC/QAQCMask_2020.tif")
    
    