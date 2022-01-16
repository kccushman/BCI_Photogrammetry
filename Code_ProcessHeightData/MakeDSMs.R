#### Read data ####
  
  # Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
  
  # read lidar (2009) and drone photogrammetry (2015, 2018, 2020) point cloud catalog objects (aligned and trimmed tiles)
  cat09 <- lidR::catalog("PointClouds/Processed/Lidar/BCI09Tiles_trim/")
  cat15at <- lidR::catalog("PointClouds/Processed/DronePhotogrammetry/BCI15Tiles_alignedTrim/")
  cat18at <- lidR::catalog("PointClouds/Processed/DronePhotogrammetry/BCI18Tiles_alignedto15Trim/")
  cat20at <- lidR::catalog("PointClouds/Processed/DronePhotogrammetry/BCI20Tiles_alignedto18Trim/")

#### Make digital surface models (DSMs) for each year ####


  dsm09 <- lidR::grid_canopy(cat09,
                             res = 1,
                             algorithm = lidR::p2r(subcircle=0.01))
  
  dsm15 <- lidR::grid_canopy(cat15at,
                             res = 1,
                             algorithm = lidR::p2r(subcircle=0.01,
                                                   na.fill = lidR::tin()))
  dsm18 <- lidR::grid_canopy(cat18at,
                             res = 1,
                             algorithm = lidR::p2r(subcircle=0.01,
                                                   na.fill = lidR::tin()))
  dsm20 <- lidR::grid_canopy(cat20at,
                             res = 1,
                             algorithm = lidR::p2r(subcircle=0.01,
                                                   na.fill = lidR::tin()))

# Crop lidar data to extent of photogrammetry data
  dsm09 <- raster::crop(dsm09, raster::extent(dsm15))

#### Write data to .tif rasters ####  
  raster::writeRaster(dsm09, file = "Data_HeightRasters/DSM_2009.tif")
  raster::writeRaster(dsm15, file = "Data_HeightRasters/DSM_2015_corrected.tif")
  raster::writeRaster(dsm18, file = "Data_HeightRasters/DSM_2018_corrected.tif")
  raster::writeRaster(dsm20, file = "Data_HeightRasters/DSM_2020_corrected.tif")
  
#### Make uncorrected DSMs for photogrammetry data ####
  
  # Define catalog objects
    cat15_raw <- lidR::catalog("PointClouds/Processed/DronePhotogrammetry/BCI15Tiles/")
    cat18_raw <- lidR::catalog('PointClouds/Processed/DronePhotogrammetry/BCI18Tiles/')
    cat20_raw <- lidR::catalog('PointClouds/Processed/DronePhotogrammetry/BCI20Tiles/')
  
  # Make DSMs 
    dsm15_raw <- lidR::grid_canopy(cat15_raw,
                                   res = 1,
                                   algorithm = lidR::p2r(subcircle=0.01,
                                                         na.fill = lidR::tin()))
   
    dsm18_raw <- lidR::grid_canopy(cat18_raw,
                                   res = 1,
                                   algorithm = lidR::p2r(subcircle=0.01,
                                                         na.fill = lidR::tin()))
  
    dsm20_raw <- lidR::grid_canopy(cat20_raw,
                                   res = 1,
                                   algorithm = lidR::p2r(subcircle=0.01,
                                                         na.fill = lidR::tin()))
  # Write raster files
    raster::writeRaster(dsm15_raw, file = "Data_HeightRasters/DSM_2015_raw.tif")
    raster::writeRaster(dsm18_raw, file = "Data_HeightRasters/DSM_2018_raw.tif")
    raster::writeRaster(dsm20_raw, file = "Data_HeightRasters/DSM_2020_raw.tif")

  
#### Make DSMs for QAQC procedure ####
    
  # 1 m pixels with no points (no gap fill)
    qaqc15a <- lidR::grid_canopy(cat15_at,
                                 res = 1,
                                 algorithm = lidR::p2r(subcircle=0.01))
    qaqc18a <- lidR::grid_canopy(cat18_at,
                                 res = 1,
                                 algorithm = lidR::p2r(subcircle=0.01))
    qaqc20a <- lidR::grid_canopy(cat20_at,
                                 res = 1,
                                 algorithm = lidR::p2r(subcircle=0.01))
    
    raster::writeRaster(qaqc15a, "Data_QAQC/DSM_noGapFill_2015.tif")
    raster::writeRaster(qaqc18a, "Data_QAQC/DSM_noGapFill_2018.tif")
    raster::writeRaster(qaqc20a, "Data_QAQC/DSM_noGapFill_2020.tif")
    
    
  # 0.2 m pixels with height range of points
    qaqc15b <- lidR::grid_metrics(cat15_at,
                                  res = 0.2,
                                  func = ~max(Z)-min(Z))
    qaqc18b <- lidR::grid_metrics(cat18_at,
                                  res = 0.2,
                                  func = ~max(Z)-min(Z))
    qaqc20b <- lidR::grid_metrics(cat20_at,
                                  res = 0.2,
                                  func = ~max(Z)-min(Z))
    
    raster::writeRaster(qaqc15b, "Data_QAQC/htRangeRaster_2015.tif")
    raster::writeRaster(qaqc15b, "Data_QAQC/htRangeRaster_2018.tif")
    raster::writeRaster(qaqc15b, "Data_QAQC/htRangeRaster_2020.tif")