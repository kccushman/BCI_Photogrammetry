##### Read DEM at 1 m resolution ####

dem_1m <- raster::raster("Data_HeightRasters/LidarDEM_BCI.tif")  


#### Smooth over various scales ####

  # sigma = 1
  w1 <- raster::focalWeight(dem_1m, 1, "Gauss")
  dem_1m_1x <- raster::focal(dem_1m, w1)
  raster::writeRaster(dem_1m_1x, "Data_TopographyRasters/DEM_smooth_1.tif")
  
  # sigma = 2
  w2 <- raster::focalWeight(dem_1m, 2, "Gauss")
  dem_1m_2x <- raster::focal(dem_1m, w2)
  raster::writeRaster(dem_1m_2x, "Data_TopographyRasters/DEM_smooth_2.tif")
  
  # sigma = 3 
  w3 <- raster::focalWeight(dem_1m, 3, "Gauss")
  dem_1m_3x <- raster::focal(dem_1m, w3)
  raster::writeRaster(dem_1m_3x, "Data_TopographyRasters/DEM_smooth_3.tif")
  
  # sigma = 4
  w4 <- raster::focalWeight(dem_1m, 4, "Gauss")
  dem_1m_4x <- raster::focal(dem_1m, w4)
  raster::writeRaster(dem_1m_4x, "Data_TopographyRasters/DEM_smooth_4.tif")
  
  # sigma = 6
  w6 <- raster::focalWeight(dem_1m, 6, "Gauss")
  dem_1m_6x <- raster::focal(dem_1m, w6)
  raster::writeRaster(dem_1m_6x, "Data_TopographyRasters/DEM_smooth_6.tif",overwrite=T,
                      options=c("COMPRESS=NONE", "TFW=YES"))
  
  # sigma = 8
  w8 <- raster::focalWeight(dem_1m, 8, "Gauss")
  dem_1m_8x <- raster::focal(dem_1m, w8)
  raster::writeRaster(dem_1m_8x, "Data_TopographyRasters/DEM_smooth_8.tif")
  
  # sigma = 12
  w12 <- raster::focalWeight(dem_1m, 12, "Gauss")
  dem_1m_12x <- raster::focal(dem_1m, w12)
  raster::writeRaster(dem_1m_12x, "Data_TopographyRasters/DEM_smooth_12.tif",overwrite=T,
                      options=c("COMPRESS=NONE", "TFW=YES"))
  
  # sigma = 16
  w16 <- raster::focalWeight(dem_1m, 16, "Gauss")
  dem_1m_16x <- raster::focal(dem_1m, w16)
  raster::writeRaster(dem_1m_16x, "Data_TopographyRasters/DEM_smooth_16.tif",overwrite=T,
                      options=c("COMPRESS=NONE", "TFW=YES"))
  
  # sigma = 24
  w24 <- raster::focalWeight(dem_1m, 24, "Gauss")
  dem_1m_24x <- raster::focal(dem_1m, w24)
  raster::writeRaster(dem_1m_24x, "Data_TopographyRasters/DEM_smooth_24.tif")
  
  # sigma = 32
  w32 <- raster::focalWeight(dem_1m, 32, "Gauss")
  dem_1m_32x <- raster::focal(dem_1m, w32)
  raster::writeRaster(dem_1m_32x, "Data_TopographyRasters/DEM_smooth_32.tif")
  
  # sigma = 48
  w48 <- raster::focalWeight(dem_1m, 48, "Gauss")
  dem_1m_48x <- raster::focal(dem_1m, w48)
  raster::writeRaster(dem_1m_48x, "Data_TopographyRasters/DEM_smooth_48.tif")
  
  # sigma = 64
  w64 <- raster::focalWeight(dem_1m, 64, "Gauss")
  dem_1m_64x <- raster::focal(dem_1m, w64)
  raster::writeRaster(dem_1m_64x, "Data_TopographyRasters/DEM_smooth_64.tif")
  
  
  # NOTE: calculation of slope, curvature, and HAND was subsequently calculated in ArcGIS Pro from each smoothed DEM