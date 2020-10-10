#### Read grid info data ####
gridInfo <- read.csv("gridInfo.csv")

#### 2009 to 2015 #### 

# Remove overlap from tiles

#overlap <- 30

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

cat09at <- lidR::catalog(paste0(path1,"lidar tiles 2009/lidar tiles 2009/"))
cat15at <- lidR::catalog(paste0(path2,"BCI15Tiles_alignedTrim/"))

chm09 <- lidR::grid_canopy(cat09at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.5))

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
chm09 <- raster::crop(chm09, raster::extent(chm15))

dHeight09to15 <- chm15-chm09

colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))

raster::plot(dHeight09to15,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2009 to 2015 height change")


raster::writeRaster(chm09, file = "CHM_2009.tif")
raster::writeRaster(chm15, file = "CHM_2015_corrected.tif")

#### 2015 to 2017 #### 

# Remove overlap from tiles

#overlap <- 30

path1 <- "D:/BCI_Spatial/Lidar_Data/"
path2 <- "D:/BCI_Spatial/UAV_Data/TiledPointClouds/"

# 2017
for(i in 1:dim(gridInfo)[1]){
  
  data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI17Tiles_aligned/BCI17a_",gridInfo$ID[i],".las")),
                               xleft=gridInfo$xmin[i],
                               xright=gridInfo$xmax[i],
                               ybottom=gridInfo$ymin[i],
                               ytop=gridInfo$ymax[i])
  if(length(data@data$X)>0){
    lidR::writeLAS(las = data,
                   file = paste0(path2, "BCI17Tiles_alignedTrim/BCI17at_",gridInfo$ID[i],".laz"))
  }
}

# Make canopy height raster for each year

cat17at <- lidR::catalog(paste0(path2,"BCI17Tiles_alignedTrim/"))

chm15 <- raster::raster("CHM_2015_corrected.tif")

chm17 <- lidR::grid_canopy(cat17at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.1))

# Plot canopy height rasters and canopy height change
colBrks1 <- seq(1,250)
raster::plot(chm15,
             main="Canopy height 2015",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)
raster::plot(chm17,
             main="Canopy height 2017",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)

dHeight15to17 <- chm17-chm15

colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))

raster::plot(dHeight15to17,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2015 to 2017 height change")

raster::writeRaster(chm17, file = "CHM_2017_corrected.tif")

#### QAQC #####

# Look for "holes" in each year
soils <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

# 2017
gridInfo$QAQC17 <- NA

chm17 <- lidR::grid_canopy(cat17at,
                           res = 0.5,
                           algorithm = lidR::p2r(subcircle=0.01))

for(i in 1:dim(gridInfo)[1]){
  x1 <- gridInfo[i, "xmin"] 
  x2 <- gridInfo[i, "xmax"]
  y1 <- gridInfo[i, "ymin"] 
  y2 <- gridInfo[i, "ymax"] 
  
  cropi <- raster::crop(dHeight15to17,
                        raster::extent(c(x1,x2,y1,y2)))
  
  gridInfo$QAQC1517[i] <-  length(cropi@data@values[is.na(cropi@data@values)])/length(cropi@data@values)
}



gridInfo$QAQC1517_mean <- NA
gridInfo$QAQC1517_median <- NA

for(i in 1:dim(gridInfo)[1]){
 x1 <- gridInfo[i, "xmin"] 
 x2 <- gridInfo[i, "xmax"]
 y1 <- gridInfo[i, "ymin"] 
 y2 <- gridInfo[i, "ymax"] 
 
 cropi <- raster::crop(dHeight15to17,
                       raster::extent(c(x1,x2,y1,y2)))
 
 gridInfo$QAQC1517_mean[i] <- mean(abs(cropi@data@values), na.rm = T)
 gridInfo$QAQC1517_median[i] <- median(abs(cropi@data@values), na.rm = T)
}

hist(gridInfo$QAQC1517_mean,
     breaks = seq(0,100,0.5),
     xlim=c(0,15),
     col="black", border="white")
hist(gridInfo$QAQC1517_median,
     breaks = seq(0,100,0.5),
     xlim=c(0,15),
     col="black", border="white")

head(gridInfo[gridInfo$QAQC1517_mean>4 & gridInfo$QAQC1517_mean<8,],100)

gridInfo[gridInfo$xmin<628039.3  & gridInfo$xmax>628039.3  & gridInfo$ymin<1011579  & gridInfo$ymax>1011579 ,]


