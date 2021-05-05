##### Read data ####

  gaps15to18 <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
  gaps18to20 <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")
  
  # Only keep gaps that pass area:perimeter threshold
  gaps15to18 <- gaps15to18[gaps15to18$use==T,]
  gaps18to20 <- gaps18to20[gaps18to20$use==T,]

# Make catalogs of aligned lidar tiles
  path2 <- "D:/BCI_Spatial/UAV_Data/TiledPointClouds/"
  
  cat15 <- lidR::catalog(paste0(path2,"BCI15Tiles_alignedTrim/"))
  cat18 <- lidR::catalog(paste0(path2,"BCI18Tiles_alignedto15Trim/"))
  cat20 <- lidR::catalog(paste0(path2,"BCI20Tiles_alignedto18Trim/"))
  
##### For 20 biggest gaps and 20 random gaps per interval, write las files of lidar data around gaps ####
  buff <- 20
  
  # 2015 - 2018
  
  big15to18 <- gaps15to18@data[order(-gaps15to18$area),"gap_id"][1:20]
  
  ran15to18 <- sample(gaps15to18$gap_id[!(gaps15to18$gap_id %in% big15to18)],
                      20,
                      replace = F)
  
  for(i in 1:length(big15to18)){
    gap_i <- gaps15to18[gaps15to18$gap_id == big15to18[i],]
    
    
    las15_i <- lidR::clip_rectangle(cat15,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    las18_i <- lidR::clip_rectangle(cat18,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    
    lidR::writeLAS(las15_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/15to18/Big_",i,"_15.las"))
    lidR::writeLAS(las18_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/15to18/Big_",i,"_18.las"))
  }
  
  for(i in 1:length(ran15to18)){
    gap_i <- gaps15to18[gaps15to18$gap_id == ran15to18[i],]
    
    
    las15_i <- lidR::clip_rectangle(cat15,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    las18_i <- lidR::clip_rectangle(cat18,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    
    lidR::writeLAS(las15_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/15to18/Random_",i,"_15.las"))
    lidR::writeLAS(las18_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/15to18/Random_",i,"_18.las"))
  }
  
  # 2018 - 2020
  
  big18to20 <- gaps18to20@data[order(-gaps18to20$area),"gap_id"][1:20]
  
  ran18to20 <- sample(gaps18to20$gap_id[!(gaps18to20$gap_id %in% big18to20)],
                      20,
                      replace = F)
  
  for(i in 1:length(big18to20)){
    gap_i <- gaps18to20[gaps18to20$gap_id == big18to20[i],]
    
    
    las18_i <- lidR::clip_rectangle(cat18,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    las20_i <- lidR::clip_rectangle(cat20,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    
    lidR::writeLAS(las18_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/18to20/Big_",i,"_18.las"))
    lidR::writeLAS(las20_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/18to20/Big_",i,"_20.las"))
  }
  
  for(i in 1:length(ran18to20)){
    gap_i <- gaps18to20[gaps18to20$gap_id == ran18to20[i],]
    
    
    las18_i <- lidR::clip_rectangle(cat18,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    las20_i <- lidR::clip_rectangle(cat20,
                                    xleft = raster::extent(gap_i)@xmin-buff,
                                    ybottom = raster::extent(gap_i)@ymin-buff,
                                    xright = raster::extent(gap_i)@xmax+buff,
                                    ytop = raster::extent(gap_i)@ymax+buff)
    
    lidR::writeLAS(las18_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/18to20/Random_",i,"_18.las"))
    lidR::writeLAS(las20_i, file = paste0("D:/BCI_Spatial/UAV_Data/ValidationSubsets/18to20/Random_",i,"_20.las"))
  }