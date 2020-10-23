#### Read grid info and soils data ####
gridInfo <- read.csv("gridInfo_QAQC.csv")

path1 <- "D:/BCI_Spatial/Lidar_Data/"
path2 <- "D:/BCI_Spatial/UAV_Data/TiledPointClouds/"

soils <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

# Read in forest age
age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
ageUse <- age[!(age$TYPE=="Clearings"),]

# Read polygon buffer 25 m inland from lake
buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  

#### 2015 ####
  
  # Make variables to store two QAQC metrics

    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Mean height range within single pixels (0.1 x 0.1 m resolution)

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
    
    plot(QAQC15_B~QAQC15_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=19,
         col = adjustcolor("black",alpha.f = 0.5))
    
    gridInfo$Use15 <- ifelse(gridInfo$QAQC15_A < 0.3 & gridInfo$QAQC15_B < 0.5,
                             T,
                             F)
    
  # Plot masked tiles on canopy height change raster
    
    chm09 <- raster::raster("CHM_2009.tif")
    chm15 <- raster::raster("CHM_2015_corrected.tif")
    
    dHeight09to15 <- chm15-chm09
    
    colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
    colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue","darkblue"))
    
    raster::plot(dHeight09to15,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2009 to 2015 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use15[i]==F){
        
        x1 <- gridInfo[i, "xmin"] 
        x2 <- gridInfo[i, "xmax"]
        y1 <- gridInfo[i, "ymin"] 
        y2 <- gridInfo[i, "ymax"] 
        
        raster::plot(raster::extent(c(x1,x2,y1,y2)), add=T, col="red", lwd=2)
      }
    }
    
    
  # Make reference folder based on QAQC metrics
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use15[i]==T){
        file.copy(from = paste0(path2, "BCI15Tiles_aligned/BCI15a_",gridInfo$ID[i],".las"),
                  to = paste0(path2, "BCI15Tiles_ref"))
        file.rename(from = paste0(path2, "BCI15Tiles_ref/BCI15a_",gridInfo$ID[i],".las"),
                    to = paste0(path2, "BCI15Tiles_ref/BCI15r_",gridInfo$ID[i],".las"))
      }
      
      if(gridInfo$Use15[i]==F){
         data <- lidR::readLAS(paste0(path1,"BCI09Tiles/BCI09_",gridInfo$ID[i],".laz"))
         lidR::writeLAS(data, paste0(path2, "BCI15Tiles_ref/BCI15r_",gridInfo$ID[i],".las"))
      }
    }
    
#### 2017 ####
    
    # Make variables to store two QAQC metrics
    
    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Mean height range within single pixels (0.1 x 0.1 m resolution)
    
    gridInfo$QAQC17_A <- NA
    gridInfo$QAQC17_B <- NA
    
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::readLAS(paste0(path2, "BCI17Tiles_alignedTrim/BCI17at_",gridInfo$ID[i],".laz"))
      
      chmA <- lidR::grid_canopy(data,
                                res = 1,
                                algorithm = lidR::p2r(subcircle=0.01))
      
      gridInfo$QAQC17_A[i] <-  length(chmA@data@values[is.na(chmA@data@values)])/length(chmA@data@values)
      
      chmB <- lidR::grid_metrics(data,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
      
      gridInfo$QAQC17_B[i] <- mean(chmB@data@values,na.rm=T)
      
    }
    
    plot(QAQC17_B~QAQC17_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=19,
         col = adjustcolor("black",alpha.f = 0.5))
    
    gridInfo$Use17 <- ifelse(gridInfo$QAQC17_A < 0.035 & gridInfo$QAQC17_B < 1.5,
                             T,
                             F)
    
    write.csv(gridInfo, "gridInfo_QAQC.csv", row.names = F)
    
  # Plot masked tiles on canopy height change raster
    
    chm15 <- raster::raster("CHM_2015_corrected.tif")
    chm17 <- raster::raster("CHM_2017_corrected.tif")
    
    dHeight15to17 <- chm17-chm15
    
    colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
    colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue","darkblue"))
    
    raster::plot(dHeight15to17,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2015 to 2017 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use17[i]==F){
        
        x1 <- gridInfo[i, "xmin"] 
        x2 <- gridInfo[i, "xmax"]
        y1 <- gridInfo[i, "ymin"] 
        y2 <- gridInfo[i, "ymax"] 
        
        raster::plot(raster::extent(c(x1,x2,y1,y2)), add=T, col="red", lwd=2)
      }
    }
    
    # Make reference folder based on QAQC metrics
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use17[i]==T){
        data <- lidR::readLAS(paste0(path2, "BCI17Tiles_aligned/BCI17a_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI17Tiles_ref/BCI17r_",gridInfo$ID[i],".laz"))
      }
      
      if(gridInfo$Use17[i]==F){
          data <- lidR::readLAS(paste0(path2,"BCI15Tiles_ref/BCI15r_",gridInfo$ID[i],".las"))
          lidR::writeLAS(data, paste0(path2, "BCI17Tiles_ref/BCI17r_",gridInfo$ID[i],".laz"))
      }
    }
 
#### 2018 ####
    
    # Make variables to store two QAQC metrics
    
    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Mean height range within single pixels (0.1 x 0.1 m resolution)
    
    gridInfo$QAQC18_A <- NA
    gridInfo$QAQC18_B <- NA
    
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::readLAS(paste0(path2, "BCI18Tiles_alignedTrim/BCI18at_",gridInfo$ID[i],".laz"))
      
      chmA <- lidR::grid_canopy(data,
                                res = 1,
                                algorithm = lidR::p2r(subcircle=0.01))
      
      gridInfo$QAQC18_A[i] <-  length(chmA@data@values[is.na(chmA@data@values)])/length(chmA@data@values)
      
      chmB <- lidR::grid_metrics(data,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
      
      gridInfo$QAQC18_B[i] <- mean(chmB@data@values,na.rm=T)
      
    }
    
    plot(QAQC18_B~QAQC18_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=19,
         col = adjustcolor("black",alpha.f = 0.5))
    
    gridInfo$Use18 <- ifelse(gridInfo$QAQC18_A < 0.035 & gridInfo$QAQC18_B < 1.5,
                             T,
                             F)
    
    write.csv(gridInfo, "gridInfo_QAQC.csv", row.names = F)
    
    # Plot masked tiles on canopy height change raster
    
    chm17 <- raster::raster("CHM_2017_corrected.tif")
    chm18 <- raster::raster("CHM_2018_corrected.tif")
    
    dHeight17to18 <- chm18-chm17
    
    colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
    colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue","darkblue"))
    
    raster::plot(dHeight17to18,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2017 to 2018 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use18[i]==F | gridInfo$Use17[i]==F){
        
        x1 <- gridInfo[i, "xmin"] 
        x2 <- gridInfo[i, "xmax"]
        y1 <- gridInfo[i, "ymin"] 
        y2 <- gridInfo[i, "ymax"] 
        
        raster::plot(raster::extent(c(x1,x2,y1,y2)), add=T, col="red", lwd=2)
      }
    }
    
    # Make reference folder based on QAQC metrics
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use18[i]==T){
        data <- lidR::readLAS(paste0(path2, "BCI18Tiles_aligned/BCI18a_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI18Tiles_ref/BCI18r_",gridInfo$ID[i],".laz"))
      }
      
      if(gridInfo$Use18[i]==F){
        data <- lidR::readLAS(paste0(path2,"BCI17Tiles_ref/BCI17r_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI18Tiles_ref/BCI18r_",gridInfo$ID[i],".laz"))
      }
    }
    
#### 2019 ####
    
    # Make variables to store two QAQC metrics
    
    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Mean height range within single pixels (0.1 x 0.1 m resolution)
    
    gridInfo$QAQC19_A <- NA
    gridInfo$QAQC19_B <- NA
    
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::readLAS(paste0(path2, "BCI19Tiles_alignedTrim/BCI19at_",gridInfo$ID[i],".laz"))
      
      chmA <- lidR::grid_canopy(data,
                                res = 1,
                                algorithm = lidR::p2r(subcircle=0.01))
      
      gridInfo$QAQC19_A[i] <-  length(chmA@data@values[is.na(chmA@data@values)])/length(chmA@data@values)
      
      chmB <- lidR::grid_metrics(data,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
      
      gridInfo$QAQC19_B[i] <- mean(chmB@data@values,na.rm=T)
      
    }
    
    plot(QAQC19_B~QAQC19_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=19,
         col = adjustcolor("black",alpha.f = 0.5))
    
    gridInfo$Use19 <- ifelse(gridInfo$QAQC19_A < 0.035 & gridInfo$QAQC19_B < 1.5,
                             T,
                             F)
    
    write.csv(gridInfo, "gridInfo_QAQC.csv", row.names = F)
    
    # Plot masked tiles on canopy height change raster
    
    chm18 <- raster::raster("CHM_2018_corrected.tif")
    chm19 <- raster::raster("CHM_2019_corrected.tif")
    
    dHeight18to19 <- chm19-chm18
    
    colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
    colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue","darkblue"))
    
    raster::plot(dHeight18to19,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2018 to 2019 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use19[i]==F | gridInfo$Use18[i]==F){
        
        x1 <- gridInfo[i, "xmin"] 
        x2 <- gridInfo[i, "xmax"]
        y1 <- gridInfo[i, "ymin"] 
        y2 <- gridInfo[i, "ymax"] 
        
        raster::plot(raster::extent(c(x1,x2,y1,y2)), add=T, col="red", lwd=2)
      }
    }
    
    # Make reference folder based on QAQC metrics
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use19[i]==T){
        data <- lidR::readLAS(paste0(path2, "BCI19Tiles_aligned/BCI19a_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI19Tiles_ref/BCI19r_",gridInfo$ID[i],".laz"))
      }
      
      if(gridInfo$Use19[i]==F){
        data <- lidR::readLAS(paste0(path2,"BCI18Tiles_ref/BCI18r_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI19Tiles_ref/BCI19r_",gridInfo$ID[i],".laz"))
      }
    }
    
#### 2020 ####
    
    # Make variables to store two QAQC metrics
    
    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Mean height range within single pixels (0.1 x 0.1 m resolution)
    
    gridInfo$QAQC20_A <- NA
    gridInfo$QAQC20_B <- NA
    
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::readLAS(paste0(path2, "BCI20Tiles_alignedTrim/BCI20at_",gridInfo$ID[i],".laz"))
      
      chmA <- lidR::grid_canopy(data,
                                res = 1,
                                algorithm = lidR::p2r(subcircle=0.01))
      
      gridInfo$QAQC20_A[i] <-  length(chmA@data@values[is.na(chmA@data@values)])/length(chmA@data@values)
      
      chmB <- lidR::grid_metrics(data,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
      
      gridInfo$QAQC20_B[i] <- mean(chmB@data@values,na.rm=T)
      
    }
    
    plot(QAQC20_B~QAQC20_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=20,
         col = adjustcolor("black",alpha.f = 0.5))
    
    gridInfo$Use20 <- ifelse(gridInfo$QAQC20_A < 0.035 & gridInfo$QAQC20_B < 1.5,
                             T,
                             F)
    
    write.csv(gridInfo, "gridInfo_QAQC.csv", row.names = F)
    
    # Plot masked tiles on canopy height change raster
    
    chm19 <- raster::raster("CHM_2019_corrected.tif")
    chm20 <- raster::raster("CHM_2020_corrected.tif")
    
    dHeight19to20 <- chm20-chm19
    
    colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
    colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue","darkblue"))
    
    raster::plot(dHeight19to20,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2019 to 2020 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use20[i]==F | gridInfo$Use18[i]==F){
        
        x1 <- gridInfo[i, "xmin"] 
        x2 <- gridInfo[i, "xmax"]
        y1 <- gridInfo[i, "ymin"] 
        y2 <- gridInfo[i, "ymax"] 
        
        raster::plot(raster::extent(c(x1,x2,y1,y2)), add=T, col="red", lwd=2)
      }
    }
    
    # Make reference folder based on QAQC metrics
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use20[i]==T){
        data <- lidR::readLAS(paste0(path2, "BCI20Tiles_aligned/BCI20a_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI20Tiles_ref/BCI20r_",gridInfo$ID[i],".laz"))
      }
      
      if(gridInfo$Use20[i]==F){
        data <- lidR::readLAS(paste0(path2,"BCI18Tiles_ref/BCI18r_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI20Tiles_ref/BCI20r_",gridInfo$ID[i],".laz"))
      }
    }
    
#### Find tile ID from coordinates ####
    
    coords <- locator(1)
    
    gridInfo[gridInfo$xmin<coords$x & gridInfo$xmax>coords$x & gridInfo$ymin<coords$y & gridInfo$ymax>coords$y,"ID"]
    
#### Cloud mask 2017 ####
    
  # Read albedo (because clouds are bright)
  albedo17 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Min_Albedo_2017.tif")
  albedo17 <- raster::mask(albedo17,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  
  redblu17 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Max_RBdiff_2017.tif")
  redblu17 <- raster::mask(redblu17,buffer)
  
  
  thresh_albedo <- 375
  thresh_redblu <- 20
  mask17 <- albedo17
  mask17@data@values[albedo17@data@values >= thresh_albedo & redblu17@data@values <= thresh_redblu] <- 0
  mask17@data@values[albedo17@data@values < thresh_albedo | redblu17@data@values > thresh_redblu] <- 1
  raster::plot(mask17)
  
  raster::writeRaster(mask17, "CloudMask_2017.tif")

#### Cloud mask 2019 ####
  
  # Read albedo (because clouds are bright)
  albedo17 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Min_Albedo_2017.tif")
  albedo17 <- raster::mask(albedo17,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  
  redblu17 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Max_RBdiff_2017.tif")
  redblu17 <- raster::mask(redblu17,buffer)
  
  
  thresh_albedo <- 375
  thresh_redblu <- 20
  mask17 <- albedo17
  mask17@data@values[albedo17@data@values >= thresh_albedo & redblu17@data@values <= thresh_redblu] <- 0
  mask17@data@values[albedo17@data@values < thresh_albedo | redblu17@data@values > thresh_redblu] <- 1
  raster::plot(mask17)
  
  raster::writeRaster(mask17, "CloudMask_2017.tif")      
    

  
  
  
  
      