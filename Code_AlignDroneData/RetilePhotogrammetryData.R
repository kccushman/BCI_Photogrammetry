#### NOTE#### 

  # Initial and final photogrammetry point cloud data are archived
  # Intermediate, duplicated datasets (i.e. tiled with overlap, decimated point clouds) are not archived

  # Code is provided for reference

#### Packages ####

# Install the following packages, if needed

  # install.packages("lidR") # version 3.0.4 used
  # install.packages("raster") # version 3.3-13 used
  # install.packages("rgdal") # version 1.5-16 used
  # install.packages("sp") # version 1.4-2 used

#### Make lidR catalog entries for each dataset ####
  
  cat09 <- lidR::catalog("D:/BCI_Spatial/Lidar_Data/lidar tiles 2009/lidar tiles 2009/")
    sp::proj4string(cat09) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
  cat15 <- lidR::catalog("PointClouds/Raw/2015_07_30_All_BCI_PointCloud.las")
  cat18 <- lidR::catalog("PointClouds/Raw/2018_Todos_BCI_Dipteryx_EBEE_PointCloud.las")
  cat20 <- lidR::catalog("PointClouds/Raw/2020_07_31_BCI_Dipteryx_Ebee_DensePts.las")

#### Create grid for new .las tiles ####
  
  xmin <- 623400
  xmax <- 630100
  ymin <- 1009700
  ymax <- 1015200
  
  tileSz <- 150
  xmins <- seq(xmin,(xmax-tileSz),tileSz)
  xmaxs <- seq((xmin+tileSz),xmax,tileSz)
  ymins <- seq(ymin,(ymax-tileSz),tileSz)
  ymaxs <- seq((ymin+tileSz),ymax,tileSz)
  
  gridInfo <- data.frame(xmin=rep(xmins, length(ymins)),
                         xmax=rep(xmaxs, length(ymins)),
                         ymin=rep(ymins, each=length(xmins)),
                         ymax=rep(ymaxs, each=length(xmins)))
  
# Only keep grid cells that actually overlay BCI, using BCI soils shapefile
  soils <- rgdal::readOGR("Data_Ancillary/BCI_Soils/BCI_Soils.shp")
  soils <- sp::spTransform(soils,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
  gridInfo$Use <- NA
  
  for(i in 1:dim(gridInfo)[1]){
    poly_i <- as(raster::extent(as.numeric(gridInfo[i,1:4])),
                     'SpatialPolygons')
    sp::proj4string(poly_i) <- sp::proj4string(soils)
    
    test_i <- sp::over(x = poly_i,
                       y = soils)
    if(is.na(test_i$ID[1])){
      gridInfo$Use[i] <- F
    }
    
    if(!is.na(test_i$ID[1])){
      gridInfo$Use[i] <- T
    }
  }
  
  # Only keep tiles that overlap BCI soils shapefile
  gridInfo <- gridInfo[gridInfo$Use==T,]
  gridInfo$ID <- 1:dim(gridInfo)[1]
  
  write.csv(gridInfo, row.names = F, file = "Data_Ancillary/gridInfo.csv")
  
  overlap <- 30

#### Retile 2009 lidar with overlap, and only keep highest points ####
  
  # THESE FILES ARE NOT ARCHIVED
  
  for(i in 1:dim(gridInfo)[1]){
    data <- lidR::lasclipRectangle(cat09, 
                                   xleft = gridInfo$xmin[i] - overlap,
                                   ybottom = gridInfo$ymin[i] - overlap,
                                   xright = gridInfo$xmax[i] + overlap,
                                   ytop = gridInfo$ymax[i] + overlap)
    data <- lidR::decimate_points(data, algorithm = lidR::highest(res=2))
    
    if(length(data@data$X)>0){
      lidR::writeLAS(data, file=paste0("D:/BCI_Spatial/Lidar_Data/BCI09Tiles/BCI09_",i,".laz"))
    }
    
    print(i)
  }
  
#### Retile 2015 photogrammetry point cloud with overlap ####
  
  # THESE FILES ARE NOT ARCHIVED
  
  # full-resolution
  for(i in 1:dim(gridInfo)[1]){
    data <- lidR::lasclipRectangle(cat15, 
                                   xleft = gridInfo$xmin[i] - overlap,
                                   ybottom = gridInfo$ymin[i] - overlap,
                                   xright = gridInfo$xmax[i] + overlap,
                                   ytop = gridInfo$ymax[i] + overlap)
    if(length(data@data$X)>0){
      lidR::writeLAS(data, file=paste0("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI15Tiles/BCI15_",i,".laz"))
    }
    
    print(i)
  }

  #decimated
  for(i in 1:dim(gridInfo)[1]){
    data <- lidR::lasclipRectangle(cat15, 
                                   xleft = gridInfo$xmin[i] - overlap,
                                   ybottom = gridInfo$ymin[i] - overlap,
                                   xright = gridInfo$xmax[i] + overlap,
                                   ytop = gridInfo$ymax[i] + overlap)
    data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
    if(length(data@data$X)>0){
      lidR::writeLAS(data, file=paste0("C:/Users/cushmank/Desktop/BCI15Tiles_dec/BCI15d_",i,".laz"))
    }
    
    print(i)
  }



#### Retile 2018 photogrammetry point cloud with overlap ####
  
  # THESE FILES ARE NOT ARCHIVED
  
  for(i in 1:dim(gridInfo)[1]){
    data <- lidR::lasclipRectangle(cat18, 
                                   xleft = gridInfo$xmin[i] - overlap,
                                   ybottom = gridInfo$ymin[i] - overlap,
                                   xright = gridInfo$xmax[i] + overlap,
                                   ytop = gridInfo$ymax[i] + overlap)
    data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
    if(length(data@data$X)>0){
      lidR::writeLAS(data, file=paste0("C:/Users/cushmank/Desktop/BCI18Tiles_dec/BCI18d_",i,".laz"))
    }
    
    print(i)
  }

  for(i in 1:dim(gridInfo)[1]){
    data <- lidR::lasclipRectangle(cat18, 
                                   xleft = gridInfo$xmin[i] - overlap,
                                   ybottom = gridInfo$ymin[i] - overlap,
                                   xright = gridInfo$xmax[i] + overlap,
                                   ytop = gridInfo$ymax[i] + overlap)
    if(length(data@data$X)>0){
      lidR::writeLAS(data, file=paste0("C:/Users/cushmank/Desktop/BCI18Tiles/BCI18_",i,".laz"))
    }
    
    print(i)
  }

#### Retile 2020 photogrammetry point cloud with overlap ####
  
  # THESE FILES ARE NOT ARCHIVED
  
  
  for(i in 1:dim(gridInfo)[1]){
  data <- lidR::lasclipRectangle(cat20, 
                                 xleft = gridInfo$xmin[i] - overlap,
                                 ybottom = gridInfo$ymin[i] - overlap,
                                 xright = gridInfo$xmax[i] + overlap,
                                 ytop = gridInfo$ymax[i] + overlap)
  data <- lidR::decimate_points(data, algorithm = lidR::highest(res=0.5))
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("C:/Users/cushmank/Desktop/BCI20Tiles_dec/BCI20_",i,".laz"))
  }
  
  print(i)
}

for(i in 1:dim(gridInfo)[1]){
  data <- lidR::lasclipRectangle(cat20, 
                                 xleft = gridInfo$xmin[i] - overlap,
                                 ybottom = gridInfo$ymin[i] - overlap,
                                 xright = gridInfo$xmax[i] + overlap,
                                 ytop = gridInfo$ymax[i] + overlap)
  if(length(data@data$X)>0){
    lidR::writeLAS(data, file=paste0("C:/Users/cushmank/Desktop/BCI20Tiles/BCI20_",i,".laz"))
  }
  
  print(i)
}

#### Rename transition matrix files ####
  
  # renames output from running do2015to2009.bat before running apply2015transformation.bat
  
    # 2015
    
      for(i in 1:dim(gridInfo)[1]){
        file.rename(from = list.files("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI15Tiles_dec/",
                                      full.names = T, pattern = paste0("BCI15d_",i,"_REG")),
                    to = paste0("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI15Tiles_dec/BCI15mat_",i,".txt"))
      }

  # renames output from running do2018to2015.bat before running apply2018transformation.bat
  
    # 2018
    
    for(i in 1:dim(gridInfo)[1]){
      file.rename(from = list.files("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI18Tiles_dec/",
                                    full.names = T, pattern = paste0("BCI18d_",i,"_REG")),
                  to = paste0("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI18Tiles_dec/BCI18mat_",i,".txt"))
    }


  # renames output from running do2020to2018.bat before running apply2020transformation.bat
  
    # 2020
    
      for(i in 1:dim(gridInfo)[1]){
        file.rename(from = list.files("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI20Tiles_dec/",
                                      full.names = T, pattern = paste0("BCI20d_",i,"_REG")),
                    to = paste0("D:/BCI_Spatial/UAV_Data/TiledPointClouds/BCI20Tiles_dec/BCI20mat_",i,".txt"))
      }

#### Remove overlap from aligned point clouds ####

  # THESE DATA FOLDERS ARE ARCHIVED
  
  # 2015
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI15Tiles_alignedFull/BCI15af_",gridInfo$ID[i],".las")),
                                   xleft=gridInfo$xmin[i],
                                   xright=gridInfo$xmax[i],
                                   ybottom=gridInfo$ymin[i],
                                   ytop=gridInfo$ymax[i])
      if(length(data@data$X)>0){
        lidR::writeLAS(las = data,
                       file = paste0(path2, "PointsClouds/Processed/DronePhotogrammetry/BCI15Tiles_alignedTrim/BCI15at_",gridInfo$ID[i],".laz"))
      }
    }


  # 2018
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI18Tiles_alignedFull/BCI18af_",gridInfo$ID[i],".las")),
                                   xleft=gridInfo$xmin[i],
                                   xright=gridInfo$xmax[i],
                                   ybottom=gridInfo$ymin[i],
                                   ytop=gridInfo$ymax[i])
      if(length(data@data$X)>0){
        lidR::writeLAS(las = data,
                       file = paste0(path2, "PointsClouds/Processed/DronePhotogrammetry/BCI18Tiles_alignedto15Trim/BCI18at_",gridInfo$ID[i],".laz"))
      }
    }


  # 2020
    for(i in 1:dim(gridInfo)[1]){
      
      data <- lidR::clip_rectangle(las = lidR::readLAS(paste0(path2,"BCI20Tiles_alignedFull/BCI20af_",gridInfo$ID[i],".las")),
                                   xleft=gridInfo$xmin[i],
                                   xright=gridInfo$xmax[i],
                                   ybottom=gridInfo$ymin[i],
                                   ytop=gridInfo$ymax[i])
      if(length(data@data$X)>0){
        lidR::writeLAS(las = data,
                       file = paste0(path2, "PointsClouds/Processed/DronePhotogrammetry/BCI20Tiles_alignedto18Trim/BCI20at_",gridInfo$ID[i],".laz"))
      }
    }
