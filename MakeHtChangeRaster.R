#### Read grid info data ####
gridInfo <- read.csv("gridInfo.csv")

path1 <- "D:/BCI_Spatial/Lidar_Data/"
path2 <- "D:/BCI_Spatial/UAV_Data/TiledPointClouds/"

#### 2009 to 2015 #### 


# Make canopy height rasters for each year

cat09 <- lidR::catalog("D:/BCI_Spatial/Lidar_Data/lidar tiles 2009/lidar tiles 2009")

cat09at <- lidR::catalog(paste0(path1,"BCI09Tiles_trim/"))
cat15at <- lidR::catalog(paste0(path2,"BCI15Tiles_alignedTrim/"))

chm09 <- lidR::grid_canopy(cat09,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01))

chm09 <- raster::crop(chm09, soils)

chm15 <- lidR::grid_canopy(cat15at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01,
                                                 na.fill = lidR::tin()))

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


raster::writeRaster(chm09, file = "DSM_2009.tif")
raster::writeRaster(chm15, file = "DSM_2015_corrected_tin.tif")

#### 2015 to 2017 #### 

# Make canopy height raster for each year

cat17at <- lidR::catalog(paste0(path2,"BCI17Tiles_alignedTrim/"))

chm15 <- raster::raster("CHM_2015_corrected.tif")

chm17 <- lidR::grid_canopy(cat17at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01,
                                                 na.fill = lidR::tin()))


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

raster::writeRaster(chm17, file = "DSM_2017_corrected_tin.tif")


#### 2017 to 2018 #### 

# Make canopy height raster for each year

cat18at <- lidR::catalog(paste0(path2,"BCI18Tiles_alignedTrim/"))

chm17 <- raster::raster("DSM_2017_corrected.tif")

chm18 <- lidR::grid_canopy(cat18at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01,
                                                 na.fill = lidR::tin()))

# Plot canopy height rasters and canopy height change
colBrks1 <- seq(1,250)
raster::plot(chm17,
             main="Canopy height 2017",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)
raster::plot(chm18,
             main="Canopy height 2018",
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

raster::writeRaster(chm18, file = "DSM_2018_corrected_tin.tif")

#### 2018 to 2019 #### 

# Make canopy height raster for each year

cat19at <- lidR::catalog(paste0(path2,"BCI19Tiles_alignedTrim/"))

dsm18 <- raster::raster("DSM_2018_corrected.tif")

dsm19 <- lidR::grid_canopy(cat19at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01,
                                                 na.fill = lidR::tin()))

# Plot canopy height rasters and canopy height change
colBrks1 <- seq(1,250)
raster::plot(dsm18,
             main="Canopy height 2018",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)
raster::plot(dsm19,
             main="Canopy height 2019",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)

dHeight18to19 <- dsm19-chm18

colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))

raster::plot(dHeight18to19,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2018 to 2019 height change")

raster::writeRaster(dsm19, file = "DSM_2019_corrected_tin.tif")

#### 2019 to 2020 #### 

# Make canopy height raster for each year

cat20at <- lidR::catalog(paste0(path2,"BCI20Tiles_alignedTrim/"))

dsm19 <- raster::raster("DSM_2019_corrected.tif")

dsm20 <- lidR::grid_canopy(cat20at,
                           res = 1,
                           algorithm = lidR::p2r(subcircle=0.01,
                                                 na.fill = lidR::tin()))

# Plot canopy height rasters and canopy height change
colBrks1 <- seq(1,250)
raster::plot(dsm19,
             main="Canopy height 2019",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)
raster::plot(dsm20,
             main="Canopy height 2020",
             col = terrain.colors(length(colBrks1), rev=T),
             breaks = colBrks1)

dHeight19to20 <- dsm20-dsm19

colBrks2 <- c(-100,-20,-10,-5,-1,-0.5,0.5,1,5,10,20,100)
colPal2 <- colorRampPalette(c("red","darksalmon","yellow",
                              "white",
                              "aliceblue","cornflowerblue","darkblue"))

raster::plot(dHeight19to20,
             col = colPal2(length(colBrks2)-1),
             breaks = colBrks2,
             main = "Corrected 2019 to 2020 height change")

raster::writeRaster(dsm20, file = "DSM_2020_corrected_tin.tif")
