#### Read data ####

    dchm17to18 <- raster::raster("dCHM17to18.tif")

  # Read in forest age
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    ageUse <- age[!(age$TYPE=="Clearings"),]
    
  # Read polygon buffer 25 m inland from lake
    buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
    buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
    buffer <- raster::intersect(buffer, ageUse)
  
#### Make polygons for 10% blocked LOO samples ####
  
  blocks <- raster::raster(nrows=12, ncol=12,
                           ext = raster::extent(buffer))
  blockPoly <- raster::rasterToPolygons(blocks)
  
  randomBlocks <- data.frame(block = 1:12^2,
                             xmin = NA,
                             xmax = NA,
                             ymin = NA,
                             ymax = NA,
                             overlap = NA,
                             area = NA)
  for(i in 1:dim(randomBlocks)[1]){
    randomBlocks$xmin[i] <- raster::extent(blockPoly[i,])@xmin
    randomBlocks$xmax[i] <- raster::extent(blockPoly[i,])@xmax
    randomBlocks$ymin[i] <- raster::extent(blockPoly[i,])@ymin
    randomBlocks$ymax[i] <- raster::extent(blockPoly[i,])@ymax
    
    overlap <- raster::crop(x=buffer,
                            y=raster::extent(blockPoly[i,]))
    
    if(is.null(overlap)){
      randomBlocks$overlap[i] <- F
    }
    
    if(!is.null(overlap)){
      
      for(j in 1:length(overlap)){
        if(j==1){
          area_i <- overlap@polygons[[j]]@area
        }
        
        if(j>1){
          area_i <- area_i + overlap@polygons[[j]]@area
        }
      }
      
      randomBlocks$overlap[i] <- T
      randomBlocks$area[i] <- area_i
    }
    
  }
  
  blockKeep <- which(randomBlocks$overlap==T)
  randomBlocks <- randomBlocks[blockKeep,]
  blockPoly <- blockPoly[blockKeep,]
  
  #plot to double check
  raster::plot(buffer)
  raster::plot(blockPoly, add=T)
  
  # make 10 blocks
  areaAll <- buffer@data$Shape_Area.1[1]
  
  sampleList <- list()
  set.seed(1)
  sample_blocks <- sample(1:length(blockPoly), replace=F)
  
  for(i in 1:10){
    areaSum <- rep(NA,length(blockPoly))
    for(j in 1:length(blockPoly)){
      areaSum[j] <- sum(randomBlocks[sample_blocks[1:j],"area"])
    }
    toOmit_start <- min(which(areaSum/areaAll > (0.1*(i-1))))
    toOmit_end <- max(which(areaSum/areaAll < (0.1*i)))
    sampleList[[i]] <- sample_blocks[toOmit_start:toOmit_end]
  }
  
  blank <- raster::raster("D:/BCI_Spatial/BCI_Topo/DEM_smooth_0.tif")
  blank <- raster::crop(blank, dchm17to18)

  for(i in 1:10){
    keepPoly <- blockPoly[-sampleList[[i]],]
    
    sampleRaster <- raster::mask(blank, keepPoly)
    
    sampleVals <- sampleRaster@data@values
    
    sampleVals[!is.na(sampleVals)] <- T
    sampleVals[is.na(sampleVals)] <- F
    
    write.csv(sampleVals,
              file=paste0("D:/BCI_Spatial/BlockSamples/Sample_",i,".csv"),
              row.names = F)
  }

  
  # Save .csv of random block info
  randomBlocks$block <- 1:dim(randomBlocks)[1]
  
  write.csv(randomBlocks, "bootstrapBlocks.csv", row.names = T)
  