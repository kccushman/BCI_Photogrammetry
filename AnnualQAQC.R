#### Read grid info and soils data ####
gridInfo <- read.csv("gridInfo.csv")

path1 <- "D:/BCI_Spatial/Lidar_Data/"
path2 <- "D:/BCI_Spatial/UAV_Data/TiledPointClouds/"

soils <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")


#### 2015 ####
  
  # Make variables to store two QAQC metrics
    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Median height range within single pixels (0.1 x 0.1 m resolution)

    gridInfo$QAQC15_A <- NA
    gridInfo$QAQC15_B <- NA
  
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::readLAS(paste0(path2, "BCI15Tiles_alignedTrim/BCI15at_",gridInfo$ID[i],".laz"))
      
      chmA <- lidR::grid_canopy(data,
                                 res = 1,
                                 algorithm = lidR::p2r(subcircle=0.01))
      
      gridInfo$QAQC15_A[i] <-  length(chmA@data@values[is.na(chmA@data@values)])/length(chmA@data@values)
      
      chmB <- lidR::grid_metrics(data,
                                res = 0.2,
                                func = ~max(Z)-min(Z))
      
      gridInfo$QAQC15_B[i] <- mean(chmB@data@values,na.rm=T)
      
    }
    
    plot(QAQC_15B~QAQC15_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=19)

#### 2017 ####
  cat17at <- lidR::catalog(paste0(path2,"BCI17Tiles_alignedTrim/"))
  
  gridInfo$QAQC17 <- NA
  
  chm17 <- lidR::grid_canopy(cat17at,
                             res = 0.5,
                             algorithm = lidR::p2r(subcircle=0.01))
  
  for(i in 1:dim(gridInfo)[1]){
    x1 <- gridInfo[i, "xmin"] 
    x2 <- gridInfo[i, "xmax"]
    y1 <- gridInfo[i, "ymin"] 
    y2 <- gridInfo[i, "ymax"] 
    
    cropi <- raster::crop(chm17,
                          raster::extent(c(x1,x2,y1,y2)))
    
    gridInfo$chm17[i] <-  length(cropi@data@values[is.na(cropi@data@values)])/length(cropi@data@values)
  }