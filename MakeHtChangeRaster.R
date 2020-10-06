#### Read grid info data ####
gridInfo <- read.csv("gridInfo.csv")

#### 2009 to 2015 #### 


# Remove overlap from tiles

  #overlap <- 25
  
  path1 <- "D:/BCI_Spatial/Lidar_Data/"
  path2 <- "D:/BCI_Spatial/UAV_Data/TiledPointClouds/"
  
  # 2009
    for(i in 1:dim(gridInfo)[1]){
     
      data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path1,"BCI09Tiles/BCI09_",gridInfo$ID[i],".laz")),
                                   xleft=gridInfo$xmin[i],
                                   xright=gridInfo$xmax[i],
                                   ybottom=gridInfo$ymin[i],
                                   ytop=gridInfo$ymax[i])
      lidR::writeLAS(las = data,
                     file = paste0(path1,"BCI09Tiles_trim/BCI09at_",gridInfo$ID[i],".laz"))
    }
  
  # 2015
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI15Tiles_aligned/BCI15a_",gridInfo$ID[i],".las")),
                                   xleft=gridInfo$xmin[i],
                                   xright=gridInfo$xmax[i],
                                   ybottom=gridInfo$ymin[i],
                                   ytop=gridInfo$ymax[i])
      if(length(data@data$X)>0){
        lidR::writeLAS(las = data,
                       file = paste0(path2, "BCI15Tiles_alignedTrim/BCI15at_",gridInfo$ID[i],".laz"))
        }
    }


# Make canopy height rasters for each year

  cat09at <- lidR::catalog(paste0(path1,"BCI09Tiles_trim/"))
  cat15at <- lidR::catalog(paste0(path2,"BCI15Tiles_alignedTrim/"))
  
  chm09 <- lidR::grid_canopy(cat09at,
                             res = 2,
                             algorithm = lidR::p2r(subcircle=1))
  
  chm15 <- lidR::grid_canopy(cat15at,
                             res = 1,
                             algorithm = lidR::p2r(subcircle=0.5))
  
# Plot canopy height rasters and canopy height change
  colBrks1 <- seq(1,250)
  raster::plot(chm09,
               main="Canopy height 2009",
               col = terrain.colors(length(colBrks1), rev=T),
               breaks = colBrks1)
  raster::plot(chm15,
               main="Canopy height 2015",
               col = terrain.colors(length(colBrks1), rev=T),
               breaks = colBrks1)

dHeight17to18 <- chm18-chm17

colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))


raster::plot(dHeight17to18,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2017 to 2018 height change")
raster::plot(dHeight18to19,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2018 to 2019 height change")

raster::writeRaster(chm17, file = "CHM_2017_corrected.tif")
raster::writeRaster(chm18, file = "CHM_2018_corrected.tif")
raster::writeRaster(chm19, file = "CHM_2019_corrected.tif")

chm17 <- raster::raster("CHM_2017_corrected.tif")
chm18 <- raster::raster("CHM_2018_corrected.tif")
chm19 <- raster::raster("CHM_2019_corrected.tif")

# Compare to uncorrected rasters
cat17_old <- lidR::catalog("2017_06_22_AllBCI_EBEE.las")
cat18_old <- lidR::catalog("2018_Todos_BCI_Dipteryx_EBEE_PointCloud.las")
cat19_old <- lidR::catalog("Cloudpoint_2019_06_24_BCI_Dipteyx.las")

chm17_old <- lidR::grid_canopy(cat17_old,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.1))

chm18_old <- lidR::grid_canopy(cat18_old,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.1))

chm19_old <- lidR::grid_canopy(cat19_old,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.1))

chm17_old <- raster::crop(chm17_old[[1]],raster::extent(chm17))
chm18_old <- raster::crop(chm18_old[[1]],raster::extent(chm17))
chm19_old <- raster::crop(chm19_old[[1]],raster::extent(chm17))

dHeight17to18_old <- chm18_old-chm17_old
dHeight18to19_old <- chm19_old-chm18_old

colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))
raster::plot(dHeight17to18_old,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Uncorrected 2017 to 2018 height change")
raster::plot(dHeight18to19_old,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Uncorrected 2018 to 2019 height change")

## New alignment for 2018-2019
mask <- raster::raster("BCI_mask/BCI_mask.tif")
mask <- raster::crop(mask,raster::extent(chm19))


chm18_old <- raster::crop(chm18_old[[1]],raster::extent(chm19))


dHeight18to19 <- chm19-chm18_old
dHeight18to19[mask==0] <- NA


colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))

pdf(height=10, width=8, file="canopyHtchange.pdf")
raster::plot(dHeight18to19,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2018 to 2019 height change")
abline(v=gridInfo$xmin)
abline(h=gridInfo$ymin)
dev.off()

raster::writeRaster(chm19, file = "CHM_2019_corrected_2.tif")
raster::writeRaster(chm18_old, file = "CHM_2018_uncorrected.tif")


## try something out
gridInfo$QC <- NA

for(i in 1:dim(gridInfo)[1]){
x1 <- gridInfo[i,"xmin"]
x2 <- gridInfo[i,"xmax"]
y1 <- gridInfo[i,"ymin"]
y2 <- gridInfo[i,"ymax"]

cropi <- raster::crop(dHeight18to19,
                      raster::extent(c(x1,x2,y1,y2)))
gridInfo$QC[i] <- mean(abs(cropi@data@values),na.rm=T)
}

hist(gridInfo$QC)

dim(gridInfo[gridInfo$QC >4,])

