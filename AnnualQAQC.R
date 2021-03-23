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

#### 2009 ####

dsm09 <- raster::raster("DSM_2009.tif")
dem09 <- raster::raster("D:/BCI_Spatial/BCI_Topo/LidarDEM_BCI.tif")
dem09 <- raster::crop(dem09,dsm09)
raster::crs(dem09) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
raster::crs(dsm09) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
dem09 <- raster::resample(dem09,dsm09)

chm09 <- dsm09-dem09
chm09 <- raster::mask(chm09, buffer)

threshLo <- 10


binLo <- chm09
binLo@data@values[!(is.na(binLo@data@values)) & binLo@data@values < threshLo] <- 0
binLo@data@values[!(is.na(binLo@data@values)) & binLo@data@values >= threshLo] <- 1


binHi <- chm09
binHi@data@values[binHi@data@values < threshHi] <- NA

raster::writeRaster(binLo, "binaryLoCanopy.tif", overwrite=T)
raster::writeRaster(binHi, "binaryHiCanopy.tif", overwrite=T)

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
    
    gridInfo$Use15 <- ifelse(gridInfo$QAQC15_A < 0.05 & gridInfo$QAQC15_B < 1.5,
                             T,
                             F)
    
  # Plot masked tiles on canopy height change raster
    
    chm09 <- raster::raster("DSM_2009.tif")
    chm15 <- raster::raster("DSM_2015_corrected_tin.tif")
    
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
                  to = paste0(path2, "BCI15Tiles_ref2"))
        file.rename(from = paste0(path2, "BCI15Tiles_ref2/BCI15a_",gridInfo$ID[i],".las"),
                    to = paste0(path2, "BCI15Tiles_ref2/BCI15r_",gridInfo$ID[i],".las"))
      }
      
      if(gridInfo$Use15[i]==F){
         data <- lidR::readLAS(paste0(path1,"BCI09Tiles/BCI09_",gridInfo$ID[i],".laz"))
         lidR::writeLAS(data, paste0(path2, "BCI15Tiles_ref2/BCI15r_",gridInfo$ID[i],".las"))
      }
    }
    
  # Make full raster of height range
    
    cat15_at <- lidR::catalog("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI15Tiles_alignedTrim/")
    qaqc15 <- lidR::grid_metrics(cat15_at,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
    raster::writeRaster(qaqc15, "htRangeRaster_2015.tif")
    rm(qaqc15)
    
#### 2018 ####
    
    # Rename transition matrix file
    for(i in 1:dim(gridInfo)[1]){
      file.rename(from = list.files("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI18Tiles_dec/",
                                    full.names = T, pattern = paste0("BCI18d_",i,"_REG")),
                  to = paste0("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI18Tiles_dec/BCI18mat2_",i,".txt"))
    }
    
    # Trim tiles
      for(i in 1:dim(gridInfo)[1]){
        
        data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI18Tiles_alignedto15Full/BCI18af_",gridInfo$ID[i],".las")),
                                     xleft=gridInfo$xmin[i],
                                     xright=gridInfo$xmax[i],
                                     ybottom=gridInfo$ymin[i],
                                     ytop=gridInfo$ymax[i])
        if(length(data@data$X)>0){
          lidR::writeLAS(las = data,
                         file = paste0(path2, "BCI18Tiles_alignedto15Trim/BCI18at_",gridInfo$ID[i],".laz"))
        }
      }
    
    # Make variables to store two QAQC metrics
    
    # A. Proportion of pixels (1 x 1 m resolution) with no data
    # B. Mean height range within single pixels (0.1 x 0.1 m resolution)
    
    gridInfo$QAQC18_A <- NA
    gridInfo$QAQC18_B <- NA
    
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::readLAS(paste0(path2, "BCI18Tiles_alignedto15Trim/BCI18at_",gridInfo$ID[i],".laz"))
      
      chmA <- lidR::grid_canopy(data,
                                res = 1,
                                algorithm = lidR::p2r(subcircle=0.01))
      
      gridInfo$QAQC18_A[i] <-  length(chmA@data@values[is.na(chmA@data@values)])/length(chmA@data@values)
      
      chmB <- lidR::grid_metrics(data,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
      
      gridInfo$QAQC18_B[i] <- mean(chmB@data@values,na.rm=T)
      print(i)
    }
    
    plot(QAQC18_B~QAQC18_A, data = gridInfo,
         xlab="Proportion of empty 1 m pixels",
         ylab="Mean height range in 0.2 m pixels",
         pch=19,
         col = adjustcolor("black",alpha.f = 0.5))
    
    gridInfo$Use18 <- ifelse(gridInfo$QAQC18_A < 0.05 & gridInfo$QAQC18_B < 1.5,
                             T,
                             F)
    
    write.csv(gridInfo, "gridInfo_QAQC.csv", row.names = F)
    
    # Plot masked tiles on canopy height change raster
    
    chm17 <- raster::raster("DSM_2017_corrected_tin.tif")
    dsm18 <- raster::raster("DSM_2018_corrected_tin.tif")
    
    dHeight15to18 <- chm18-chm15
    
    colBrks2 <- c(-100,-5,-0.5,0.5,5,100)
    colPal2 <- colorRampPalette(c("red","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue"))
    
    raster::plot(dHeight15to18,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2015 to 2018 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use18[i]==F | gridInfo$Use15[i]==F){
        
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
        data <- lidR::readLAS(paste0(path2, "BCI18Tiles_alignedto15/BCI18a_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI18Tiles_ref2/BCI18r_",gridInfo$ID[i],".laz"))
      }
      
      if(gridInfo$Use18[i]==F){
        data <- lidR::readLAS(paste0(path2,"BCI15Tiles_ref2/BCI15r_",gridInfo$ID[i],".las"))
        lidR::writeLAS(data, paste0(path2, "BCI18Tiles_ref2/BCI18r_",gridInfo$ID[i],".laz"))
      }
    }
    
    cat18_at <- lidR::catalog("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI18Tiles_alignedto15Trim/")
    qaqc18 <- lidR::grid_metrics(cat18_at,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
    raster::writeRaster(qaqc18, "htRangeRaster_2018.tif")
    rm(qaqc18)
    

#### 2020 ####
    
    # Rename transition matrix file
    for(i in 1:dim(gridInfo)[1]){
      file.rename(from = list.files("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI20Tiles_dec/",
                                    full.names = T, pattern = paste0("BCI20d_",i,"_REG")),
                  to = paste0("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI20Tiles_dec/BCI20mat3_",i,".txt"))
    }
    
    # Trim tiles
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI20Tiles_alignedto18Full/BCI20af_",gridInfo$ID[i],".las")),
                                   xleft=gridInfo$xmin[i],
                                   xright=gridInfo$xmax[i],
                                   ybottom=gridInfo$ymin[i],
                                   ytop=gridInfo$ymax[i])
      if(length(data@data$X)>0){
        lidR::writeLAS(las = data,
                       file = paste0(path2, "BCI20Tiles_alignedto18Trim/BCI20at_",gridInfo$ID[i],".laz"))
      }
      print(i)
    }
    
    
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
    
    gridInfo$Use20 <- ifelse(gridInfo$QAQC20_A < 0.05 & gridInfo$QAQC20_B < 1.5,
                             T,
                             F)
    
    write.csv(gridInfo, "gridInfo_QAQC.csv", row.names = F)
    
    # Plot masked tiles on canopy height change raster
    
    dsm19 <- raster::raster("DSM_2019_corrected_tin.tif")
    dsm20 <- raster::raster("DSM_2020_corrected_tin.tif")
    
    dHeight19to20 <- dsm20-dsm19
    
    colBrks2 <- c(-100,-20,-10,-5,-1,1,5,10,20,100)
    colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                                  "white",
                                  "aliceblue","cornflowerblue","darkblue"))
    
    raster::plot(dHeight19to20,
                 col = colPal2(length(colBrks2)-1),
                 breaks = colBrks2,
                 main = "Corrected 2019 to 2020 height change")
    
    for(i in 1:dim(gridInfo)[1]){
      if(gridInfo$Use20[i]==F | gridInfo$Use19[i]==F){
        
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
    
    cat20_at <- lidR::catalog("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI20Tiles_alignedto19Trim/")
    qaqc20 <- lidR::grid_metrics(cat20_at,
                                 res = 0.2,
                                 func = ~max(Z)-min(Z))
    raster::writeRaster(qaqc20, "htRangeRaster_2020.tif")
    rm(qaqc20)
    
    
#### Find tile ID from coordinates ####
    
    coords <- locator(1)
    ID <- gridInfo[gridInfo$xmin<coords$x & gridInfo$xmax>coords$x & gridInfo$ymin<coords$y & gridInfo$ymax>coords$y,"ID"]
    
    gridInfo[ID,c("QAQC18_A","QAQC18_B")]
    
#### Cloud mask 2015 ####

# Read albedo (because clouds are bright)
albedo15 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Min_Albedo_2015.tif")
albedo15 <- raster::mask(albedo15,buffer)

# Read red-blue difference (because clouds tend to be blue)
redblu15 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Max_RBdiff_2015.tif")
redblu15 <- raster::mask(redblu15,buffer)


thresh_albedo <- 475
thresh_redblu <- 20
mask15 <- albedo15
mask15@data@values[albedo15@data@values >= thresh_albedo & redblu15@data@values <= thresh_redblu] <- 0
mask15@data@values[albedo15@data@values < thresh_albedo | redblu15@data@values > thresh_redblu] <- 1
raster::plot(mask15)

raster::writeRaster(mask15, "CloudMask_2015.tif")
    
#### Cloud mask 2018 ####
  
  # Read albedo (because clouds are bright)
  albedo18 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Min_Albedo_2018.tif")
  albedo18 <- raster::mask(albedo18,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  redblu18 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Max_RBdiff_2018.tif")
  redblu18 <- raster::mask(redblu18,buffer)
  
  
  thresh_albedo <- 375
  thresh_redblu <- 20
  mask18 <- albedo18
  mask18@data@values[albedo18@data@values >= thresh_albedo & redblu18@data@values <= thresh_redblu] <- 0
  mask18@data@values[albedo18@data@values < thresh_albedo | redblu18@data@values > thresh_redblu] <- 1
  raster::plot(mask18)
  
  raster::writeRaster(mask18, "CloudMask_2018.tif")      
    

#### Cloud mask 2020 ####
  
  # Read albedo (because clouds are bright)
  albedo20 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Min_Albedo_2020.tif")
  albedo20 <- raster::mask(albedo20,buffer)
  
  # Read red-blue difference (because clouds tend to be blue)
  redblu20 <- raster::raster("D:/BCI_Spatial/UAV_Data/CloudMasks/Max_RBdiff_2020.tif")
  redblu20 <- raster::mask(redblu20,buffer)
  
  
  thresh_albedo <- 375
  thresh_redblu <- 15
  mask20 <- albedo20
  mask20@data@values[albedo20@data@values >= thresh_albedo & redblu20@data@values <= thresh_redblu] <- 0
  mask20@data@values[albedo20@data@values < thresh_albedo | redblu20@data@values > thresh_redblu] <- 1
  raster::plot(mask20)
  
  raster::writeRaster(mask20, "CloudMask_2020.tif")      
  
  
  
  
      