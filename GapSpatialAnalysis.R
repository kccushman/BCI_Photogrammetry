#### Read data ####

  # Read grid info
  gridInfo <- read.csv("gridInfo_QAQC.csv")
  
  # Read in forest age
  age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  ageUse <- age[!(age$TYPE=="Clearings"),]
  
  # Canopy height models for all years
  chm17 <- raster::raster("CHM_2017_QAQC.tif")
  chm18 <- raster::raster("CHM_2018_QAQC.tif")
  chm19 <- raster::raster("CHM_2019_QAQC.tif")
  chm20 <- raster::raster("CHM_2020_QAQC.tif")
  
    # Mask all layers to exclude clearings
    chm17 <- raster::mask(chm17, ageUse)
    chm18 <- raster::mask(chm18, ageUse)
    chm19 <- raster::mask(chm19, ageUse)
    chm20 <- raster::mask(chm20, ageUse)
  
  # Make canopy height change rasters
  dchm17to18 <- chm18-chm17
  dchm18to19 <- chm19-chm18
  dchm19to20 <- chm20-chm19
  
  allLayerList <- list(dchm17to18,dchm18to19,dchm19to20)
  
  # # Gap rasters
  # gaps17to18 <- raster::raster("newGaps17to18.tif")
  # gaps18to19 <- raster::raster("newGaps18to19.tif")
  # gaps19to20 <- raster::raster("newGaps19to20.tif")
  
  # Gap rasters
  gaps17to18 <- raster::raster("newGaps17to18a.tif")
  gaps18to19 <- raster::raster("newGaps18to19a.tif")
  gaps19to20 <- raster::raster("newGaps19to20a.tif")
  
    # Mask gap layers to exclude clearings
    gaps17to18 <- raster::mask(gaps17to18, ageUse)
    gaps18to19 <- raster::mask(gaps18to19, ageUse)
    gaps19to20 <- raster::mask(gaps19to20, ageUse)
    
  gapLayerList <- list(gaps17to18, gaps18to19, gaps19to20)
  
  # Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  buffer <- raster::intersect(buffer, ageUse)
  
  # Low and high forest area from 2009
  
  lo09 <- raster::raster("binaryLoCanopy.tif")
  hi09 <- raster::raster("binaryHiCanopy.tif")
  
  lo09 <- raster::mask(lo09, buffer)
  hi09 <- raster::mask(hi09, buffer)
  
  # DEM
  dem09 <- raster::raster("D:/BCI_Spatial/BCI_Topo/LidarDEM_BCI.tif")
  dem09 <- raster::mask(dem09, buffer)
  dem09 <- raster::crop(dem09, chm17)
  
  # Read soil type polygons
  soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
  soil <- sp::spTransform(soil,raster::crs(chm17))
  
  # Basic stats for each year
  gapPct17to18 <- length(raster::values(gaps17to18)[!is.na(raster::values(gaps17to18))])/
    length(raster::values(dchm17to18)[!is.na(raster::values(dchm17to18))])
  
  gapPct18to19 <- length(raster::values(gaps18to19)[!is.na(raster::values(gaps18to19))])/
    length(raster::values(dchm18to19)[!is.na(raster::values(dchm18to19))])
  
  gapPct19to20 <- (12/13)*length(raster::values(gaps19to20)[!is.na(raster::values(gaps19to20))])/
    length(raster::values(dchm19to20)[!is.na(raster::values(dchm19to20))])
  
#### Make polygons for 10% blocked LOO samples ####
  
  blocks <- raster::raster(nrows=10, ncol=10,
                           ext = raster::extent(buffer))
  blockPoly <- raster::rasterToPolygons(blocks)
  
  randomBlocks <- data.frame(block = 1:100,
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
  
  # make a list of 1000 random samples leaving at least 10% out
  areaAll <- buffer@data$Shape_Area
  
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
  blank <- raster::crop(blank, chm17)

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

#### Define a function to test distribution of values ####
  
  topoTest <- function(topoLayer, gapLayer, allLayer, 
                       sampleDir,
                       mask, method="KS"){
    
    # Read topographic raster and convert to correct size
    topoRaster <- raster::raster(topoLayer)
    topoRaster <- raster::crop(topoRaster,gapLayer)
    topoRaster <- raster::mask(topoRaster,mask)
    
    # Make a vectors of all topo values, and gap values
    topoData <- raster::values(topoRaster)
    allData <- raster::values(allLayer)
    gapData <- raster::values(gapLayer)
    
    # # Make masked rasters of topographic values
    # allRaster <- topoRaster
    # allRaster@data@values[is.na(allData)] <- NA
    # 
    # gapRaster <- topoRaster
    # gapRaster@data@values[is.na(gapData)] <- NA
    
  
    # Do Kolmogorov-Smirnov or Watson-Williams test test
    nIter <- length(list.files(sampleDir))
    results <- data.frame(D=rep(NA,nIter),
                          p=rep(NA,nIter))

    for(i in 1:nIter){
      blockData <- read.csv(list.files(sampleDir,
                            pattern=paste0("_",i,".csv"),
                            full.names = T))
      
      allSample <- topoData
      allSample[is.na(allData)] <- NA
      allSample[blockData==0] <- NA
      allSample <- allSample[!is.na(allSample)]
      
      gapSample <- topoData
      gapSample[is.na(gapData)] <- NA
      gapSample[blockData==0] <- NA
      gapSample <- gapSample[!is.na(gapSample)]
      
      
      if(method=="KS"){
        KS <- ks.test(allSample, gapSample)
        results$D[i] <- KS$statistic
        results$p[i] <- KS$p.value
      }
      if(method=="WW"){
        WW <- circular::watson.wheeler.test(list(allSample, gapSample))
        results$D[i] <- WW$statistic
        results$p[i] <- WW$p.value
      }
      
      #print(i)
    }
    
    return(results)
  }
  
#### Define a function to calculate the ratio of distributions ####
  
  ratioCalc <- function(topoLayer, gapLayer, allLayer, 
                        sampleDir,
                        mask, 
                        nBins = 256){
    
    # Read topographic raster and convert to correct size
    topoRaster <- raster::raster(topoLayer)
    topoRaster <- raster::crop(topoRaster,gapLayer)
    topoRaster <- raster::mask(topoRaster,mask)
    
    
    # Make a vectors of all topo values, and gap values
    topoData <- raster::values(topoRaster)
    allData <- raster::values(allLayer)
    gapData <- raster::values(gapLayer)
    
    # # Make masked rasters of topographic values
    # allRaster <- topoRaster
    # allRaster@data@values[is.na(allData)] <- NA
    # 
    # gapRaster <- topoRaster
    # gapRaster@data@values[is.na(gapData)] <- NA
    
    
    # Do Kolmogorov-Smirnov or Watson-Williams test test
    nIter <- length(list.files(sampleDir))
    results <- matrix(nrow=nIter+1,
                      ncol=nBins,
                      data=NA)
    
    
    bins <- quantile(topoData,
                     probs = seq(0,1,length.out = nBins+1),
                     na.rm=T)
    
    results[1,] <- (bins[1:(nBins)]+bins[2:(nBins+1)])/2
    
    for(i in 1:nIter){
      blockData <- read.csv(list.files(sampleDir,
                                       pattern=paste0("_",i,".csv"),
                                       full.names = T))
      
      allSample <- topoData
      allSample[is.na(allData)] <- NA
      allSample[blockData==0] <- NA
      allSample <- allSample[!is.na(allSample)]
      
      gapSample <- topoData
      gapSample[is.na(gapData)] <- NA
      gapSample[blockData==0] <- NA
      gapSample <- gapSample[!is.na(gapSample)]
      
      allCount <- rep(NA,nBins)
      gapCount <- rep(NA,nBins)
      
      for(j in 1:(nBins)){
        allCount[j] <- length(allSample[allSample > bins[j] 
                                        & allSample <= bins[j+1]
                                        & !is.na(allSample)])
        
        gapCount[j] <- length(gapSample[gapSample > bins[j]
                                        & gapSample <= bins[j+1]
                                        & !is.na(gapSample)])
        if(gapCount[j]==0 | allCount[j]== 0){
          print("ERROR need fewer bins")
        }
        
      }
      
      #print(i)
      
      results[i+1,] <-gapCount/allCount

    }
    
    return(results)
  }
  
#### Test for variation with slope at different smoothing scales ####
  #smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
  smoothScales <- c(3,9,15,21,27,33,39,45,51)
  endYrs <- c("2018","2019","2020")
  
  for(i in 1:length(smoothScales)){
    
    for(j in 1:length(gapLayerList)){
      results_ij <- topoTest(topoLayer = paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",
                                                  smoothScales[i],".tif"),
                             gapLayer = gapLayerList[[j]],
                             allLayer = allLayerList[[j]],
                             sampleDir = "D:/BCI_Spatial/BlockSamples/",
                             mask = buffer)
      results_ij$Var <- "Slope"
      results_ij$Scale <- smoothScales[i]
      results_ij$Year <- endYrs[j]
      
      if(j==1 & i==1){
        resultsSlope <- results_ij
      }
      
      if(j>1 | i>1){
        resultsSlope <- rbind(resultsSlope, results_ij)
      }
      print(paste0(i,"_",j))
    }
    
  }
  
#### Test for variation with curvature at different smoothing scales ####
  #smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
  smoothScales <- c(3,9,15,21,27,33,39,45,51)
  
  for(i in 1:length(smoothScales)){
    
    for(j in 1:length(gapLayerList)){
      results_ij <- topoTest(topoLayer = paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",
                                                smoothScales[i],".tif"),
                             gapLayer = gapLayerList[[j]],
                             allLayer = allLayerList[[j]],
                             sampleDir = "D:/BCI_Spatial/BlockSamples/",
                             mask = buffer)
      results_ij$Var <- "Curvature"
      results_ij$Scale <- smoothScales[i]
      results_ij$Year <- endYrs[j]
      
      if(j==1 & i==1){
        resultsCurv <- results_ij
      }
      
      if(j>1|i>1){
        resultsCurv <- rbind(resultsCurv, results_ij)
      }
      print(paste0(i,"_",j))
    }
    
  }  
  

  
#### Test for variation with aspect at different smoothing scales ####
  #smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
  smoothScales <- c(3,15,27,39,51)
  endYrs <- c("2018","2019","2020")
  
  for(i in 1:length(smoothScales)){
    
    for(j in 1:length(gapLayerList)){
      results_ij <- topoTest(topoLayer = paste0("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_",
                                                smoothScales[i],".tif"),
                             gapLayer = gapLayerList[[j]],
                             allLayer = allLayerList[[j]],
                             sampleDir = "D:/BCI_Spatial/BlockSamples/",
                             mask = buffer,
                             method = "WW")
      results_ij$Var <- "Aspect"
      results_ij$Scale <- smoothScales[i]
      results_ij$Year <- endYrs[j]
      
      if(j==1 & i==1){
        resultsAspect <- results_ij
      }
      
      if(j>1|i>1){
        resultsAspect <- rbind(resultsAspect, results_ij)
      }
      print(paste0(i,"_",j))
    }
    
  }  
  
  
  
#### Test for variation with Laplacian convexity at different smoothing scales ####
  #smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
  smoothScales <- c(3,15,27,39,51)
  endYrs <- c("2018","2019","2020")
  
  for(i in 1:length(smoothScales)){
    
    for(j in 1:length(gapLayerList)){
      results_ij <- topoTest(topoLayer = paste0("D:/BCI_Spatial/BCI_Topo/Convex_smooth_",
                                                smoothScales[i],".tif"),
                             gapLayer = gapLayerList[[j]],
                             allLayer = allLayerList[[j]],
                             sampleDir = "D:/BCI_Spatial/BlockSamples/",
                             mask = buffer)
      results_ij$Var <- "Convexity"
      results_ij$Scale <- smoothScales[i]
      results_ij$Year <- endYrs[j]
      
      if(j==1 & i==1){
        resultsConvex <- results_ij
      }
      
      if(j>1|i>1){
        resultsConvex <- rbind(resultsConvex, results_ij)
      }
      print(paste0(i,"_",j))
    }
    
  }  
  
#### Summarize topographic results ####
  #smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
  smoothScales <- c(3,15,27,39,51)
  endYrs <- c("2018","2019","2020")
  
  summaryAspect <- data.frame(Scale = rep(smoothScales, length(endYrs)),
                              Year = rep(endYrs, each = length(smoothScales)),
                              minD = NA,
                              maxD = NA,
                              meanD = NA,
                              minP = NA,
                              maxP = NA,
                              meanP = NA)
  for(i in 1:dim(summaryAspect)[1]){
    summaryAspect$minD[i] <- quantile(resultsAspect[resultsAspect$Year==summaryAspect$Year[i] & 
                                                    resultsAspect$Scale==summaryAspect$Scale[i], "D"],
                                      probs = 0.025)
    summaryAspect$maxD[i] <- quantile(resultsAspect[resultsAspect$Year==summaryAspect$Year[i] & 
                                                      resultsAspect$Scale==summaryAspect$Scale[i], "D"],
                                      probs = 0.975)
    summaryAspect$meanD[i] <- mean(resultsAspect[resultsAspect$Year==summaryAspect$Year[i] & 
                                                      resultsAspect$Scale==summaryAspect$Scale[i], "D"])
    summaryAspect$minP[i] <- quantile(resultsAspect[resultsAspect$Year==summaryAspect$Year[i] & 
                                                      resultsAspect$Scale==summaryAspect$Scale[i], "p"],
                                      probs = 0.025)
    summaryAspect$maxP[i] <- quantile(resultsAspect[resultsAspect$Year==summaryAspect$Year[i] & 
                                                      resultsAspect$Scale==summaryAspect$Scale[i], "p"],
                                      probs = 0.975)
    summaryAspect$meanP[i] <- mean(resultsAspect[resultsAspect$Year==summaryAspect$Year[i] & 
                                                   resultsAspect$Scale==summaryAspect$Scale[i], "p"])
    
  }
  
  summarySlope <- data.frame(Scale = rep(smoothScales, length(endYrs)),
                              Year = rep(endYrs, each = length(smoothScales)),
                              minD = NA,
                              maxD = NA,
                              meanD = NA,
                              minP = NA,
                              maxP = NA,
                              meanP = NA)
  for(i in 1:dim(summarySlope)[1]){
    summarySlope$minD[i] <- quantile(resultsSlope[resultsSlope$Year==summarySlope$Year[i] & 
                                                      resultsSlope$Scale==summarySlope$Scale[i], "D"],
                                      probs = 0.025)
    summarySlope$maxD[i] <- quantile(resultsSlope[resultsSlope$Year==summarySlope$Year[i] & 
                                                      resultsSlope$Scale==summarySlope$Scale[i], "D"],
                                      probs = 0.975)
    summarySlope$meanD[i] <- mean(resultsSlope[resultsSlope$Year==summarySlope$Year[i] & 
                                                   resultsSlope$Scale==summarySlope$Scale[i], "D"])
    summarySlope$minP[i] <- quantile(resultsSlope[resultsSlope$Year==summarySlope$Year[i] & 
                                                      resultsSlope$Scale==summarySlope$Scale[i], "p"],
                                      probs = 0.025)
    summarySlope$maxP[i] <- quantile(resultsSlope[resultsSlope$Year==summarySlope$Year[i] & 
                                                      resultsSlope$Scale==summarySlope$Scale[i], "p"],
                                      probs = 0.975)
    summarySlope$meanP[i] <- mean(resultsSlope[resultsSlope$Year==summarySlope$Year[i] & 
                                                   resultsSlope$Scale==summarySlope$Scale[i], "p"])
    
  }
  
  summaryCurv <- data.frame(Scale = rep(smoothScales, length(endYrs)),
                              Year = rep(endYrs, each = length(smoothScales)),
                              minD = NA,
                              maxD = NA,
                              meanD = NA,
                              minP = NA,
                              maxP = NA,
                              meanP = NA)
  for(i in 1:dim(summaryCurv)[1]){
    summaryCurv$minD[i] <- quantile(resultsCurv[resultsCurv$Year==summaryCurv$Year[i] & 
                                                      resultsCurv$Scale==summaryCurv$Scale[i], "D"],
                                      probs = 0.025)
    summaryCurv$maxD[i] <- quantile(resultsCurv[resultsCurv$Year==summaryCurv$Year[i] & 
                                                      resultsCurv$Scale==summaryCurv$Scale[i], "D"],
                                      probs = 0.975)
    summaryCurv$meanD[i] <- mean(resultsCurv[resultsCurv$Year==summaryCurv$Year[i] & 
                                                   resultsCurv$Scale==summaryCurv$Scale[i], "D"])
    summaryCurv$minP[i] <- quantile(resultsCurv[resultsCurv$Year==summaryCurv$Year[i] & 
                                                      resultsCurv$Scale==summaryCurv$Scale[i], "p"],
                                      probs = 0.025)
    summaryCurv$maxP[i] <- quantile(resultsCurv[resultsCurv$Year==summaryCurv$Year[i] & 
                                                      resultsCurv$Scale==summaryCurv$Scale[i], "p"],
                                      probs = 0.975)
    summaryCurv$meanP[i] <- mean(resultsCurv[resultsCurv$Year==summaryCurv$Year[i] & 
                                                   resultsCurv$Scale==summaryCurv$Scale[i], "p"])
    
  }
  
  summaryConvex <- data.frame(Scale = rep(smoothScales, length(endYrs)),
                            Year = rep(endYrs, each = length(smoothScales)),
                            minD = NA,
                            maxD = NA,
                            meanD = NA,
                            minP = NA,
                            maxP = NA,
                            meanP = NA)
  for(i in 1:dim(summaryConvex)[1]){
    summaryConvex$minD[i] <- quantile(resultsConvex[resultsConvex$Year==summaryConvex$Year[i] & 
                                                  resultsConvex$Scale==summaryConvex$Scale[i], "D"],
                                    probs = 0.025)
    summaryConvex$maxD[i] <- quantile(resultsConvex[resultsConvex$Year==summaryConvex$Year[i] & 
                                                  resultsConvex$Scale==summaryConvex$Scale[i], "D"],
                                    probs = 0.975)
    summaryConvex$meanD[i] <- mean(resultsConvex[resultsConvex$Year==summaryConvex$Year[i] & 
                                               resultsConvex$Scale==summaryConvex$Scale[i], "D"])
    summaryConvex$minP[i] <- quantile(resultsConvex[resultsConvex$Year==summaryConvex$Year[i] & 
                                                  resultsConvex$Scale==summaryConvex$Scale[i], "p"],
                                    probs = 0.025)
    summaryConvex$maxP[i] <- quantile(resultsConvex[resultsConvex$Year==summaryConvex$Year[i] & 
                                                  resultsConvex$Scale==summaryConvex$Scale[i], "p"],
                                    probs = 0.975)
    summaryConvex$meanP[i] <- mean(resultsConvex[resultsConvex$Year==summaryConvex$Year[i] & 
                                               resultsConvex$Scale==summaryConvex$Scale[i], "p"])
    
  }
  
  

  par(mfrow=c(1,1),las=1, mar=c(3,5,1,1))
  # 2018 relative importance
  plot(meanD~Scale, data=summarySlope[summarySlope$Year=="2018",],
       type="l",
       ylim=c(0.02,0.14),
       xlab=NA,
       ylab=NA,
       cex.axis=1.5)
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summarySlope[summarySlope$Year=="2018","minD"],
              rev(summarySlope[summarySlope$Year=="2018","maxD"])),
          border=NA, col=adjustcolor("black",0.2))
  points(meanD~Scale, data=summarySlope[summarySlope$Year=="2018",],
         pch=19)

  
  lines(meanD~Scale, data=summaryCurv[summaryCurv$Year=="2018",],
        type="l", col="blue")
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryCurv[summaryCurv$Year=="2018","minD"],
              rev(summaryCurv[summaryCurv$Year=="2018","maxD"])),
          border=NA, col=adjustcolor("blue",0.2))
  points(meanD~Scale, data=summaryCurv[summaryCurv$Year=="2018",],
         pch=19, col="blue")
  
  lines(meanD~Scale, data=summaryConvex[summaryConvex$Year=="2018",],
        type="l", col="red")
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryConvex[summaryConvex$Year=="2018","minD"],
              rev(summaryConvex[summaryConvex$Year=="2018","maxD"])),
          border=NA, col=adjustcolor("red",0.2))
  points(meanD~Scale, data=summaryConvex[summaryConvex$Year=="2018",],
         pch=19, col="red")
  
  
  
  # 2019 relative importance
  plot(meanD~Scale, data=summarySlope[summarySlope$Year=="2019",],
       type="l",
       ylim=c(0.02,0.14),
       xlab=NA,
       ylab=NA,
       cex.axis=1.5)
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summarySlope[summarySlope$Year=="2019","minD"],
              rev(summarySlope[summarySlope$Year=="2019","maxD"])),
          border=NA, col=adjustcolor("black",0.2))
  points(meanD~Scale, data=summarySlope[summarySlope$Year=="2019",],
         pch=19)
  
  
  lines(meanD~Scale, data=summaryCurv[summaryCurv$Year=="2019",],
        type="l", col="blue")
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryCurv[summaryCurv$Year=="2019","minD"],
              rev(summaryCurv[summaryCurv$Year=="2019","maxD"])),
          border=NA, col=adjustcolor("blue",0.2))
  points(meanD~Scale, data=summaryCurv[summaryCurv$Year=="2019",],
         pch=19, col="blue")
  
  lines(meanD~Scale, data=summaryConvex[summaryConvex$Year=="2019",],
        type="l", col="red")
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryConvex[summaryConvex$Year=="2019","minD"],
              rev(summaryConvex[summaryConvex$Year=="2019","maxD"])),
          border=NA, col=adjustcolor("red",0.2))
  points(meanD~Scale, data=summaryConvex[summaryConvex$Year=="2019",],
         pch=19, col="red")
  
  # 2020 relative importance
  plot(meanD~Scale, data=summarySlope[summarySlope$Year=="2020",],
       type="l",
       ylim=c(0.02,0.14),
       xlab=NA,
       ylab=NA,
       cex.axis=1.5)
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summarySlope[summarySlope$Year=="2020","minD"],
              rev(summarySlope[summarySlope$Year=="2020","maxD"])),
          border=NA, col=adjustcolor("black",0.2))
  points(meanD~Scale, data=summarySlope[summarySlope$Year=="2020",],
         pch=19)
  
  
  lines(meanD~Scale, data=summaryCurv[summaryCurv$Year=="2020",],
        type="l", col="blue")
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryCurv[summaryCurv$Year=="2020","minD"],
              rev(summaryCurv[summaryCurv$Year=="2020","maxD"])),
          border=NA, col=adjustcolor("blue",0.2))
  points(meanD~Scale, data=summaryCurv[summaryCurv$Year=="2020",],
         pch=19, col="blue")
  
  lines(meanD~Scale, data=summaryConvex[summaryConvex$Year=="2020",],
        type="l", col="red")
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryConvex[summaryConvex$Year=="2020","minD"],
              rev(summaryConvex[summaryConvex$Year=="2020","maxD"])),
          border=NA, col=adjustcolor("red",0.2))
  points(meanD~Scale, data=summaryConvex[summaryConvex$Year=="2020",],
         pch=19, col="red")
  
  
  
  # Separate aspect plots
  plot(meanD~Scale, data=summaryAspect[summaryAspect$Year=="2018",],
       type="l",
       xlab=NA,
       ylab=NA,
       cex.axis=1.5)
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryAspect[summaryAspect$Year=="2018","minD"],
              rev(summaryAspect[summaryAspect$Year=="2018","maxD"])),
          border=NA, col=adjustcolor("black",0.2))
  points(meanD~Scale, data=summaryAspect[summaryAspect$Year=="2018",],
         pch=19)
  
  plot(meanD~Scale, data=summaryAspect[summaryAspect$Year=="2019",],
       type="l",
       xlab=NA,
       ylab=NA,
       cex.axis=1.5)
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryAspect[summaryAspect$Year=="2019","minD"],
              rev(summaryAspect[summaryAspect$Year=="2019","maxD"])),
          border=NA, col=adjustcolor("black",0.2))
  points(meanD~Scale, data=summaryAspect[summaryAspect$Year=="2019",],
         pch=19)
  
  plot(meanD~Scale, data=summaryAspect[summaryAspect$Year=="2020",],
       type="l",
       xlab=NA,
       ylab=NA,
       cex.axis=1.5)
  polygon(x=c(smoothScales,rev(smoothScales)),
          y=c(summaryAspect[summaryAspect$Year=="2020","minD"],
              rev(summaryAspect[summaryAspect$Year=="2020","maxD"])),
          border=NA, col=adjustcolor("black",0.2))
  points(meanD~Scale, data=summaryAspect[summaryAspect$Year=="2020",],
         pch=19)
  
#### Distribution in gaps vs non gaps for best smoothing scale: slope ####

  slopeComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  slopeComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[2]],
                           allLayer = allLayerList[[2]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  slopeComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[3]],
                           allLayer = allLayerList[[3]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  slopeCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryLoCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  slopeCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryHiCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  #load("prelimSlopeComp.RData")
  
  slope18_mean <- apply(slopeComp18[-1,], MARGIN = 2, mean)
  slope18_min <- apply(slopeComp18[-1,], MARGIN = 2, min)
  slope18_max <- apply(slopeComp18[-1,], MARGIN = 2, max)
  
  slope19_mean <- apply(slopeComp19[-1,], MARGIN = 2, mean)
  slope19_min <- apply(slopeComp19[-1,], MARGIN = 2, min)
  slope19_max <- apply(slopeComp19[-1,], MARGIN = 2, max)
  
  slope20_mean <- apply(slopeComp20[-1,], MARGIN = 2, mean)*(12/13)
  slope20_min <- apply(slopeComp20[-1,], MARGIN = 2, min)*(12/13)
  slope20_max <- apply(slopeComp20[-1,], MARGIN = 2, max)*(12/13)
  
  slopeLo_mean <- apply(slopeCompLo[-1,], MARGIN = 2, mean)
  slopeLo_min <- apply(slopeCompLo[-1,], MARGIN = 2, min)
  slopeLo_max <- apply(slopeCompLo[-1,], MARGIN = 2, max)
  
  slopeHi_mean <- apply(slopeCompHi[-1,], MARGIN = 2, mean)
  slopeHi_min <- apply(slopeCompHi[-1,], MARGIN = 2, min)
  slopeHi_max <- apply(slopeCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=slopeComp18[1,], y=slope18_mean,
       ylab=NA, cex.axis = 1.5,
       type="l", lwd=2,
       ylim=c(1e-2,0.04),
       log="y")
  polygon(x=c(slopeComp18[1,],rev(slopeComp18[1,])),
          y=c(slope18_min,
              rev(slope18_max)),
          border=NA, col=adjustcolor("black",0.2))

  lines(x=slopeComp19[1,], y=slope19_mean,
        lwd=2, col = "blue")
  polygon(x=c(slopeComp19[1,],rev(slopeComp19[1,])),
          y=c(slope19_min,
              rev(slope19_max)),
          border=NA, col=adjustcolor("blue",0.2))

  lines(x=slopeComp20[1,], y=slope20_mean,
       lwd=2,
       col="red")
  polygon(x=c(slopeComp20[1,],rev(slopeComp20[1,])),
          y=c(slope20_min,
              rev(slope20_max)),
          border=NA, col=adjustcolor("red",0.2))
  

  
  
  lines(x=slopeCompLo[1,], y=slopeLo_mean,
        lwd=2,
        col="orange")
  polygon(x=c(slopeCompLo[1,],rev(slopeCompLo[1,])),
          y=c(slopeLo_min,
              rev(slopeLo_max)),
          border=NA, col=adjustcolor("orange",0.2))

  lines(x=slopeCompHi[1,], y=slopeHi_mean,
        lwd=2,
        col="green")
  polygon(x=c(slopeCompHi[1,],rev(slopeCompHi[1,])),
          y=c(slopeHi_min,
              rev(slopeHi_max)),
          border=NA, col=adjustcolor("green",0.2))
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 2,
       xlim=range(slopeComp18[1,]),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5)
  lines(dens19, col = "blue", lwd = 2)
  lines(dens20, col="red", lwd = 2)
  
  
  
  #save(slopeComp18, slopeComp19, slopeComp20, slopeCompLo, slopeCompHi, file="prelimSlopeComp.RData")
  
#### Distribution in gaps vs non gaps for best smoothing scale: curvature ####
  endYrs <- c("2019","2019","2020")
  
  curvComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  curvComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                           gapLayer = gapLayerList[[2]],
                           allLayer = allLayerList[[2]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  curvComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                           gapLayer = gapLayerList[[3]],
                           allLayer = allLayerList[[3]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  curvCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                           gapLayer = raster::raster("binaryLoCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  curvCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                           gapLayer = raster::raster("binaryHiCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer)
  
  #load("prelimCurvComp.RData")
  
  curv18_mean <- apply(curvComp18[-1,], MARGIN = 2, FUN=mean)
  curv18_min <- apply(curvComp18[-1,], MARGIN = 2, min)
  curv18_max <- apply(curvComp18[-1,], MARGIN = 2, max)
  
  curv19_mean <- apply(curvComp19[-1,], MARGIN = 2, FUN=mean)
  curv19_min <- apply(curvComp19[-1,], MARGIN = 2, min)
  curv19_max <- apply(curvComp19[-1,], MARGIN = 2, max)
  
  curv20_mean <- apply(curvComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  curv20_min <- apply(curvComp20[-1,], MARGIN = 2, min)*(12/13)
  curv20_max <- apply(curvComp20[-1,], MARGIN = 2, max)*(12/13)
  
  curvLo_mean <- apply(curvCompLo[-1,], MARGIN = 2, FUN=mean)
  curvLo_min <- apply(curvCompLo[-1,], MARGIN = 2, min)
  curvLo_max <- apply(curvCompLo[-1,], MARGIN = 2, max)
  
  curvHi_mean <- apply(curvCompHi[-1,], MARGIN = 2, FUN=mean)
  curvHi_min <- apply(curvCompHi[-1,], MARGIN = 2, min)
  curvHi_max <- apply(curvCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=curvComp18[1,], y=curv18_mean,
       ylab=NA, cex.axis = 1.5,
       type="l", lwd=2,
       ylim=c(1e-2,0.04),
       log="y")
  
  polygon(x=c(curvComp18[1,],rev(curvComp18[1,])),
          y=c(curv18_min,
              rev(curv18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=curvComp19[1,], y=curv19_mean,
        lwd=2, col = "blue")
  polygon(x=c(curvComp19[1,],rev(curvComp19[1,])),
          y=c(curv19_min,
              rev(curv19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=curvComp20[1,], y=curv20_mean,
        lwd=2,
        col="red")
  polygon(x=c(curvComp20[1,],rev(curvComp20[1,])),
          y=c(curv20_min,
              rev(curv20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=1,lty=2)
  
  lines(x=curvCompLo[1,], y=curvLo_mean,
        lwd=2,
        col="orange")
  polygon(x=c(curvCompLo[1,],rev(curvCompLo[1,])),
          y=c(curvLo_min,
              rev(curvLo_max)),
          border=NA, col=adjustcolor("orange",0.2))

  lines(x=curvCompHi[1,], y=curvHi_mean,
        lwd=2,
        col="green")
  polygon(x=c(curvCompHi[1,],rev(curvCompHi[1,])),
          y=c(curvHi_min,
              rev(curvHi_max)),
          border=NA, col=adjustcolor("green",0.2))
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)], n=1024)
  dens19 <- density(topoVals[!is.na(all19)], n=1024)
  dens20 <- density(topoVals[!is.na(all20)], n=1024)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 2,
       xlim=range(curvComp18[1,]),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5)
  lines(dens19, col = "blue", lwd = 2)
  lines(dens20, col="red", lwd = 2)
  

  #save(curvComp18, curvComp19, curvComp20, curvCompLo, curvCompHi, file="prelimCurvComp.RData")

#### Distribution in gaps vs non gaps for best smoothing scale: aspect ####
  endYrs <- c("2019","2019","2020")
  
  aspectComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  aspectComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                          gapLayer = gapLayerList[[2]],
                          allLayer = allLayerList[[2]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  aspectComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                          gapLayer = gapLayerList[[3]],
                          allLayer = allLayerList[[3]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  aspectCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                          gapLayer = raster::raster("binaryLoCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  aspectCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                          gapLayer = raster::raster("binaryHiCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer)
  
  #load("prelimAspectComp.RData")
  
  aspect18_mean <- apply(aspectComp18[-1,], MARGIN = 2, FUN=mean)
  aspect18_min <- apply(aspectComp18[-1,], MARGIN = 2, min)
  aspect18_max <- apply(aspectComp18[-1,], MARGIN = 2, max)
  
  aspect19_mean <- apply(aspectComp19[-1,], MARGIN = 2, FUN=mean)
  aspect19_min <- apply(aspectComp19[-1,], MARGIN = 2, min)
  aspect19_max <- apply(aspectComp19[-1,], MARGIN = 2, max)
  
  aspect20_mean <- apply(aspectComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  aspect20_min <- apply(aspectComp20[-1,], MARGIN = 2, min)*(12/13)
  aspect20_max <- apply(aspectComp20[-1,], MARGIN = 2, max)*(12/13)
  
  aspectLo_mean <- apply(aspectCompLo[-1,], MARGIN = 2, FUN=mean)
  aspectLo_min <- apply(aspectCompLo[-1,], MARGIN = 2, min)
  aspectLo_max <- apply(aspectCompLo[-1,], MARGIN = 2, max)
  
  aspectHi_mean <- apply(aspectCompHi[-1,], MARGIN = 2, FUN=mean)
  aspectHi_min <- apply(aspectCompHi[-1,], MARGIN = 2, min)
  aspectHi_max <- apply(aspectCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=aspectComp18[1,], y=aspect18_mean,
       ylab=NA, cex.axis = 1.5,
       type="l", lwd=2,
       ylim=c(0.9e-2,0.05),
       log="y")
  polygon(x=c(aspectComp18[1,],rev(aspectComp18[1,])),
          y=c(aspect18_min,
              rev(aspect18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=aspectComp19[1,], y=aspect19_mean,
        lwd=2, col = "blue")
  polygon(x=c(aspectComp19[1,],rev(aspectComp19[1,])),
          y=c(aspect19_min,
              rev(aspect19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=aspectComp20[1,], y=aspect20_mean,
        lwd=2,
        col="red")
  polygon(x=c(aspectComp20[1,],rev(aspectComp20[1,])),
          y=c(aspect20_min,
              rev(aspect20_max)),
          border=NA, col=adjustcolor("red",0.2))

  abline(v=c(0,90,180,270,360), lty=2)
  
  # lines(x=aspectCompLo[1,], y=aspectLo_mean,
  #       lwd=2,
  #       col="orange")
  # polygon(x=c(aspectCompLo[1,],rev(aspectCompLo[1,])),
  #         y=c(aspectLo_min,
  #             rev(aspectLo_max)),
  #         border=NA, col=adjustcolor("orange",0.2))
  # 
  # lines(x=aspectCompHi[1,], y=aspectHi_mean,
  #       lwd=2,
  #       col="green")
  # polygon(x=c(aspectCompHi[1,],rev(aspectCompHi[1,])),
  #         y=c(aspectHi_min,
  #             rev(aspectHi_max)),
  #         border=NA, col=adjustcolor("green",0.2))
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  par(las=1, mar=c(3,5,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 2,
       xlim=range(aspectComp18[1,]),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5)
  lines(dens19, col = "blue", lwd = 2)
  lines(dens20, col="red", lwd = 2)
  
  
  #save(aspectComp18, aspectComp19, aspectComp20, aspectCompLo, aspectCompHi, file="prelimAspectComp.RData")
  
#### Soil type by year ####
  
  # Make a vector of unique soil types
  soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                          PropGap = NA,
                          nGap = NA,
                          AreaSampled = NA)
  
  # 2017 to 2018
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::mask(dchm17to18,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- raster::mask(gaps17to18,soil[soil$SOIL==soilTypes$Soil[i],])
    values_areai <- raster::getValues(areai)
    values_gapsi <- raster::getValues(gapsi)
    soilTypes$PropGap[i] <- length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$nGap[i] <- length(unique(gapsi@data@values[!is.na(gapsi@data@values)]))/length(values_areai[!is.na(values_areai)])
    soilTypes$AreaSampled[i] <-length(values_areai[!is.na(values_areai)])/(10000)
  }
  soilTypes17to18 <- soilTypes
  
  # 2018 to 2019  
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::mask(dchm18to19,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- raster::mask(gaps18to19,soil[soil$SOIL==soilTypes$Soil[i],])
    values_areai <- raster::getValues(areai)
    values_gapsi <- raster::getValues(gapsi)
    soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$nGap[i] <- length(unique(gapsi@data@values[!is.na(gapsi@data@values)]))/length(values_areai[!is.na(values_areai)])
    soilTypes$AreaSampled[i] <-length(values_areai[!is.na(values_areai)])/(10000)
  }
  soilTypes18to19 <- soilTypes
  
  # 2019 to 2020  
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::mask(dchm19to20,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- raster::mask(gaps19to20,soil[soil$SOIL==soilTypes$Soil[i],])
    values_areai <- raster::getValues(areai)
    values_gapsi <- raster::getValues(gapsi)
    soilTypes$PropGap[i] <-  length(values_gapsi[!is.na(values_gapsi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$nGap[i] <- length(unique(gapsi@data@values[!is.na(gapsi@data@values)]))/length(values_areai[!is.na(values_areai)])
    soilTypes$AreaSampled[i] <-length(values_areai[!is.na(values_areai)])/(10000)
  }
  soilTypes19to20 <- soilTypes
  # correct for longer sampling interval
  soilTypes19to20$PropGap <- soilTypes19to20$PropGap/(13/12)
  
  # Make a vector of unique soil types
  soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                          PropLo = NA,
                          PropHi = NA)
  
  # 2009 lo and hi height areas
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::mask(dem09,soil[soil$SOIL==soilTypes$Soil[i],])
    lo_i <- raster::mask(lo09,soil[soil$SOIL==soilTypes$Soil[i],])
    hi_i <- raster::mask(hi09,soil[soil$SOIL==soilTypes$Soil[i],])
    
    values_areai <- raster::getValues(areai)
    values_loi <- raster::getValues(lo_i)
    values_hii <- raster::getValues(hi_i)
    soilTypes$PropLo[i] <-  length(values_loi[!is.na(values_loi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$PropHi[i] <-  length(values_hii[!is.na(values_hii)])/length(values_areai[!is.na(values_areai)])
    
  }
  PropLoHi_All <- soilTypes
  
  # 2009 lo and hi height areas -- Old Growth
  soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                          PropLo_Old = NA,
                          PropHi_Old = NA,
                          Area_Old = NA)
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::mask(dem09,soil[soil$SOIL==soilTypes$Soil[i],])
    lo_i <- raster::mask(lo09,soil[soil$SOIL==soilTypes$Soil[i],])
    hi_i <- raster::mask(hi09,soil[soil$SOIL==soilTypes$Soil[i],])
    
    areai <- raster::mask(areai,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
    lo_i <- raster::mask(lo_i,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
    hi_i <- raster::mask(hi_i,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
    
    values_areai <- raster::getValues(areai)
    values_loi <- raster::getValues(lo_i)
    values_hii <- raster::getValues(hi_i)
    soilTypes$PropLo_Old[i] <-  length(values_loi[!is.na(values_loi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$PropHi_Old[i] <-  length(values_hii[!is.na(values_hii)])/length(values_areai[!is.na(values_areai)])
    
    soilTypes$Area_Old[i] <- length(values_areai[!is.na(values_areai)])/(10000)
  }
  PropLoHi_Old <- soilTypes
  
  # 2009 lo and hi height areas -- Secondary
  soilTypes <- data.frame(Soil = as.character(unique(soil@data$SOIL)),
                          PropLo_Sec = NA,
                          PropHi_Sec = NA,
                          Area_Sec = NA)
  for(i in 1:dim(soilTypes)[1]){
    areai <- raster::mask(dem09,soil[soil$SOIL==soilTypes$Soil[i],])
    lo_i <- raster::mask(lo09,soil[soil$SOIL==soilTypes$Soil[i],])
    hi_i <- raster::mask(hi09,soil[soil$SOIL==soilTypes$Soil[i],])
    
    areai <- raster::mask(areai,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
    lo_i <- raster::mask(lo_i,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
    hi_i <- raster::mask(hi_i,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
    
    values_areai <- raster::getValues(areai)
    values_loi <- raster::getValues(lo_i)
    values_hii <- raster::getValues(hi_i)
    soilTypes$PropLo_Sec[i] <-  length(values_loi[!is.na(values_loi)])/length(values_areai[!is.na(values_areai)])
    soilTypes$PropHi_Sec[i] <-  length(values_hii[!is.na(values_hii)])/length(values_areai[!is.na(values_areai)])
    
    soilTypes$Area_Sec[i] <- length(values_areai[!is.na(values_areai)])/(10000)
    
  }
  PropLoHi_Sec <- soilTypes
  
  
  names(soilTypes17to18) <- c("Soil","PropGap17to18","nGap17to18","Area17to18")
  names(soilTypes18to19) <- c("Soil","PropGap18to19","nGap18to19","Area18to19")
  names(soilTypes19to20) <- c("Soil","PropGap19to20","nGap19to20","Area19to20")
  
  soilTypes <- merge(soilTypes17to18,soilTypes18to19)
  soilTypes <- merge(soilTypes, soilTypes19to20)
  soilTypes <- merge(soilTypes, PropLoHi_All)
  soilTypes <- merge(soilTypes, PropLoHi_Old)
  soilTypes <- merge(soilTypes, PropLoHi_Sec)
  
  minSz <- 50
  
  soilTypes$Good18 <- F
  soilTypes[soilTypes$Area17to18>minSz,"Good18"] <- T
  soilTypes$Good19 <- F
  soilTypes[soilTypes$Area18to19>minSz,"Good19"] <- T
  soilTypes$Good20 <- F
  soilTypes[soilTypes$Area19to20>minSz,"Good20"] <- T

  
  

  par(mfrow=c(1,2),
      mar=c(4,4,1,1))
  plot(PropGap18to19~PropGap17to18,
       data=soilTypes[soilTypes$Good18==T & soilTypes$Good19==T,],
       pch=20,
       xlim=c(0.01,0.04),
       ylim=c(0.01,0.04),
       ylab="2018-2019",
       xlab="2017-2018")
  abline(a=0,b=1)
  summary(lm(PropGap18to19~PropGap17to18,
             data=soilTypes[soilTypes$Good18==T & soilTypes$Good19==T,],))

  
  plot(PropGap19to20~PropGap18to19,
       data=soilTypes[soilTypes$Good19==T & soilTypes$Good20==T,],
       pch=20,
       xlim=c(0.01,0.04),
       ylim=c(0.01,0.04),
       ylab="2019-2020",
       xlab="2018-2019")
  abline(a=0,b=1)
  summary(lm(PropGap19to20~PropGap18to19,
             data=soilTypes[soilTypes$Good19==T & soilTypes$Good20==T,]))

  
  plot(PropGap19to20~PropGap17to18, data=soilTypes[soilTypes$Good==T,],
       pch=20,
       xlim=c(0.01,0.04),
       ylim=c(0.01,0.04),
       ylab="2019-2020",
       xlab="2017-2018")
  abline(a=0,b=1)
  summary(lm(PropGap19to20~PropGap17to18,data=soilTypes[soilTypes$Good==T,]))
  summary(lm(nGap19to20~nGap17to18,data=soilTypes[soilTypes$Good==T,]))
  
  
#### Fairchild soil: distribution in gaps vs non gaps for best smoothing scale: slope ####
  endYrs <- c("2018","2019","2020")
  
  slopeComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[1]],
                           allLayer = allLayerList[[1]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Fairchild",])
  
  slopeComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[2]],
                           allLayer = allLayerList[[2]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Fairchild",])
  
  slopeComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[3]],
                           allLayer = allLayerList[[3]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Fairchild",])
  
  slopeCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryLoCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Fairchild",])
  
  slopeCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryHiCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Fairchild",])
  
  #load("prelimSlopeComp_Fairchild.RData")
  
  slope18_mean <- apply(slopeComp18[-1,], MARGIN = 2, mean)
  slope18_min <- apply(slopeComp18[-1,], MARGIN = 2, min)
  slope18_max <- apply(slopeComp18[-1,], MARGIN = 2, max)
  
  slope19_mean <- apply(slopeComp19[-1,], MARGIN = 2, mean)
  slope19_min <- apply(slopeComp19[-1,], MARGIN = 2, min)
  slope19_max <- apply(slopeComp19[-1,], MARGIN = 2, max)
  
  slope20_mean <- apply(slopeComp20[-1,], MARGIN = 2, mean)*(12/13)
  slope20_min <- apply(slopeComp20[-1,], MARGIN = 2, min)*(12/13)
  slope20_max <- apply(slopeComp20[-1,], MARGIN = 2, max)*(12/13)
  
  slopeLo_mean <- apply(slopeCompLo[-1,], MARGIN = 2, mean)
  slopeLo_min <- apply(slopeCompLo[-1,], MARGIN = 2, min)
  slopeLo_max <- apply(slopeCompLo[-1,], MARGIN = 2, max)
  
  slopeHi_mean <- apply(slopeCompHi[-1,], MARGIN = 2, mean)
  slopeHi_min <- apply(slopeCompHi[-1,], MARGIN = 2, min)
  slopeHi_max <- apply(slopeCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=slopeComp18[1,], y=slope18_mean,
       ylab=NA,cex.axis=1.5,
       type="l", lwd=2,
       ylim=c(0.7e-2,0.04),
       log="y")
  polygon(x=c(slopeComp18[1,],rev(slopeComp18[1,])),
          y=c(slope18_min,
              rev(slope18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=slopeComp19[1,], y=slope19_mean,
        lwd=2, col = "blue")
  polygon(x=c(slopeComp19[1,],rev(slopeComp19[1,])),
          y=c(slope19_min,
              rev(slope19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=slopeComp20[1,], y=slope20_mean,
        lwd=2,
        col="red")
  polygon(x=c(slopeComp20[1,],rev(slopeComp20[1,])),
          y=c(slope20_min,
              rev(slope20_max)),
          border=NA, col=adjustcolor("red",0.2))

  
  lines(x=slopeCompLo[1,], y=slopeLo_mean,
        lwd=2,
        col="orange")
  polygon(x=c(slopeCompLo[1,],rev(slopeCompLo[1,])),
          y=c(slopeLo_min,
              rev(slopeLo_max)),
          border=NA, col=adjustcolor("orange",0.2))

  lines(x=slopeCompHi[1,], y=slopeHi_mean,
        lwd=2,
        col="green")
  polygon(x=c(slopeCompHi[1,],rev(slopeCompHi[1,])),
          y=c(slopeHi_min,
              rev(slopeHi_max)),
          border=NA, col=adjustcolor("green",0.2))
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::mask(topoVals, soil[soil$SOIL=="Fairchild",])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18b <- density(topoVals[!is.na(all18)],na.rm=T)
  dens19b <- density(topoVals[!is.na(all19)],na.rm=T)
  dens20b <- density(topoVals[!is.na(all20)],na.rm=T)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 3,
       xlim=range(slopeComp18[1,]),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5,
       col = adjustcolor("black",0.5))
  lines(dens19, col = adjustcolor("black",0.5), lwd = 3)
  lines(dens20, col= adjustcolor("black",0.5), lwd = 3)
  lines(dens18b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens19b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens20b, col = adjustcolor("green",0.5), lwd = 3)
  
  
  
  
  
  #save(slopeComp18, slopeComp19, slopeComp20, slopeCompLo, slopeCompHi, file="prelimSlopeComp_Fairchild.RData")
  
#### Fairchild soil: distribution in gaps vs non gaps for best smoothing scale: curvature ####

  curvComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Fairchild",])
  
  curvComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[2]],
                          allLayer = allLayerList[[2]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Fairchild",])
  
  curvComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[3]],
                          allLayer = allLayerList[[3]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Fairchild",])
  
  curvCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryLoCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Fairchild",])
  
  curvCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryHiCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Fairchild",])
  
  #load("prelimCurvComp_Fairchild.RData")
  
  curv18_mean <- apply(curvComp18[-1,], MARGIN = 2, FUN=mean)
  curv18_min <- apply(curvComp18[-1,], MARGIN = 2, min)
  curv18_max <- apply(curvComp18[-1,], MARGIN = 2, max)
  
  curv19_mean <- apply(curvComp19[-1,], MARGIN = 2, FUN=mean)
  curv19_min <- apply(curvComp19[-1,], MARGIN = 2, min)
  curv19_max <- apply(curvComp19[-1,], MARGIN = 2, max)
  
  curv20_mean <- apply(curvComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  curv20_min <- apply(curvComp20[-1,], MARGIN = 2, min)*(12/13)
  curv20_max <- apply(curvComp20[-1,], MARGIN = 2, max)*(12/13)
  
  curvLo_mean <- apply(curvCompLo[-1,], MARGIN = 2, FUN=mean)
  curvLo_min <- apply(curvCompLo[-1,], MARGIN = 2, min)
  curvLo_max <- apply(curvCompLo[-1,], MARGIN = 2, max)
  
  curvHi_mean <- apply(curvCompHi[-1,], MARGIN = 2, FUN=mean)
  curvHi_min <- apply(curvCompHi[-1,], MARGIN = 2, min)
  curvHi_max <- apply(curvCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=curvComp18[1,], y=curv18_mean,
       ylab=NA, cex.axis = 1.5,
       type="l", lwd=2,
       ylim=c(0.5e-2,0.04),
       log="y")
  
  polygon(x=c(curvComp18[1,],rev(curvComp18[1,])),
          y=c(curv18_min,
              rev(curv18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=curvComp19[1,], y=curv19_mean,
        lwd=2, col = "blue")
  polygon(x=c(curvComp19[1,],rev(curvComp19[1,])),
          y=c(curv19_min,
              rev(curv19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=curvComp20[1,], y=curv20_mean,
        lwd=2,
        col="red")
  polygon(x=c(curvComp20[1,],rev(curvComp20[1,])),
          y=c(curv20_min,
              rev(curv20_max)),
          border=NA, col=adjustcolor("red",0.2))

  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::mask(topoVals, soil[soil$SOIL=="Fairchild",])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18b <- density(topoVals[!is.na(all18)],na.rm=T)
  dens19b <- density(topoVals[!is.na(all19)],na.rm=T)
  dens20b <- density(topoVals[!is.na(all20)],na.rm=T)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 3,
       xlim=range(curvComp18[1,]),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5,
       col = adjustcolor("black",0.5))
  lines(dens19, col = adjustcolor("black",0.5), lwd = 3)
  lines(dens20, col= adjustcolor("black",0.5), lwd = 3)
  lines(dens18b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens19b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens20b, col = adjustcolor("green",0.5), lwd = 3)
  
  
  
  #save(curvComp18, curvComp19, curvComp20, curvCompLo, curvCompHi, file="prelimCurvComp_Fairchild.RData")
  
#### Fairchild soil: distribution in gaps vs non gaps for best smoothing scale: aspect ####
  endYrs <- c("2019","2019","2020")
  
  aspectComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = gapLayerList[[1]],
                            allLayer = allLayerList[[1]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Fairchild",])
  
  aspectComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = gapLayerList[[2]],
                            allLayer = allLayerList[[2]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Fairchild",])
  
  aspectComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = gapLayerList[[3]],
                            allLayer = allLayerList[[3]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Fairchild",])
  
  aspectCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = raster::raster("binaryLoCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_39.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Fairchild",])
  
  aspectCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = raster::raster("binaryHiCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_39.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Fairchild",])
  
  #load("prelimAspectComp_Fairchild.RData")
  
  aspect18_mean <- apply(aspectComp18[-1,], MARGIN = 2, FUN=mean)
  aspect18_min <- apply(aspectComp18[-1,], MARGIN = 2, min)
  aspect18_max <- apply(aspectComp18[-1,], MARGIN = 2, max)
  
  aspect19_mean <- apply(aspectComp19[-1,], MARGIN = 2, FUN=mean)
  aspect19_min <- apply(aspectComp19[-1,], MARGIN = 2, min)
  aspect19_max <- apply(aspectComp19[-1,], MARGIN = 2, max)
  
  aspect20_mean <- apply(aspectComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  aspect20_min <- apply(aspectComp20[-1,], MARGIN = 2, min)*(12/13)
  aspect20_max <- apply(aspectComp20[-1,], MARGIN = 2, max)*(12/13)
  
  aspectLo_mean <- apply(aspectCompLo[-1,], MARGIN = 2, FUN=mean)
  aspectLo_min <- apply(aspectCompLo[-1,], MARGIN = 2, min)
  aspectLo_max <- apply(aspectCompLo[-1,], MARGIN = 2, max)
  
  aspectHi_mean <- apply(aspectCompHi[-1,], MARGIN = 2, FUN=mean)
  aspectHi_min <- apply(aspectCompHi[-1,], MARGIN = 2, min)
  aspectHi_max <- apply(aspectCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=aspectComp18[1,], y=aspect18_mean,
       ylab=NA, cex.axis = 1.5,
       type="l", lwd=2,
       ylim=c(0.5e-2,0.05),
       log="y")
  polygon(x=c(aspectComp18[1,],rev(aspectComp18[1,])),
          y=c(aspect18_min,
              rev(aspect18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=aspectComp19[1,], y=aspect19_mean,
        lwd=2, col = "blue")
  polygon(x=c(aspectComp19[1,],rev(aspectComp19[1,])),
          y=c(aspect19_min,
              rev(aspect19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=aspectComp20[1,], y=aspect20_mean,
        lwd=2,
        col="red")
  polygon(x=c(aspectComp20[1,],rev(aspectComp20[1,])),
          y=c(aspect20_min,
              rev(aspect20_max)),
          border=NA, col=adjustcolor("red",0.2))
  
  abline(v=c(0,90,180,270,360), lty=2)
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::mask(topoVals, soil[soil$SOIL=="Fairchild",])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18b <- density(topoVals[!is.na(all18)],na.rm=T)
  dens19b <- density(topoVals[!is.na(all19)],na.rm=T)
  dens20b <- density(topoVals[!is.na(all20)],na.rm=T)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 3,
       xlim=range(aspectComp18[1,]),
       ylim=range(dens20b$y),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5,
       col = adjustcolor("black",0.5))
  lines(dens19, col = adjustcolor("black",0.5), lwd = 3)
  lines(dens20, col= adjustcolor("black",0.5), lwd = 3)
  lines(dens18b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens19b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens20b, col = adjustcolor("green",0.5), lwd = 3)
  
  #save(aspectComp18, aspectComp19, aspectComp20, aspectCompLo, aspectCompHi, file="prelimAspectComp_Fairchild.RData")
  
#### Zetek soil: distribution in gaps vs non gaps for best smoothing scale: slope ####

  slopeComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[1]],
                           allLayer = allLayerList[[1]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Zetek",])
  
  slopeComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[2]],
                           allLayer = allLayerList[[2]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Zetek",])
  
  slopeComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[3]],
                           allLayer = allLayerList[[3]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Zetek",])
  
  slopeCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryLoCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Zetek",])
  
  slopeCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryHiCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = soil[soil$SOIL=="Zetek",])
  
  #load("prelimSlopeComp_Zetek.RData")
  
  slope18_mean <- apply(slopeComp18[-1,], MARGIN = 2, mean)
  slope18_min <- apply(slopeComp18[-1,], MARGIN = 2, min)
  slope18_max <- apply(slopeComp18[-1,], MARGIN = 2, max)
  
  slope19_mean <- apply(slopeComp19[-1,], MARGIN = 2, mean)
  slope19_min <- apply(slopeComp19[-1,], MARGIN = 2, min)
  slope19_max <- apply(slopeComp19[-1,], MARGIN = 2, max)
  
  slope20_mean <- apply(slopeComp20[-1,], MARGIN = 2, mean)*(12/13)
  slope20_min <- apply(slopeComp20[-1,], MARGIN = 2, min)*(12/13)
  slope20_max <- apply(slopeComp20[-1,], MARGIN = 2, max)*(12/13)
  
  slopeLo_mean <- apply(slopeCompLo[-1,], MARGIN = 2, mean)
  slopeLo_min <- apply(slopeCompLo[-1,], MARGIN = 2, min)
  slopeLo_max <- apply(slopeCompLo[-1,], MARGIN = 2, max)
  
  slopeHi_mean <- apply(slopeCompHi[-1,], MARGIN = 2, mean)
  slopeHi_min <- apply(slopeCompHi[-1,], MARGIN = 2, min)
  slopeHi_max <- apply(slopeCompHi[-1,], MARGIN = 2, max)
  
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=slopeComp18[1,], y=slope18_mean,
       ylab=NA,cex.axis=1.5,
       xlim=c(0,32),
       type="l", lwd=2,
       ylim=c(0.7e-2,0.051),
       log="y")
  polygon(x=c(slopeComp18[1,],rev(slopeComp18[1,])),
          y=c(slope18_min,
              rev(slope18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=slopeComp19[1,], y=slope19_mean,
        lwd=2, col = "blue")
  polygon(x=c(slopeComp19[1,],rev(slopeComp19[1,])),
          y=c(slope19_min,
              rev(slope19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=slopeComp20[1,], y=slope20_mean,
        lwd=2,
        col="red")
  polygon(x=c(slopeComp20[1,],rev(slopeComp20[1,])),
          y=c(slope20_min,
              rev(slope20_max)),
          border=NA, col=adjustcolor("red",0.2))
  
  
  lines(x=slopeCompLo[1,], y=slopeLo_mean,
        lwd=2,
        col="orange")
  polygon(x=c(slopeCompLo[1,],rev(slopeCompLo[1,])),
          y=c(slopeLo_min,
              rev(slopeLo_max)),
          border=NA, col=adjustcolor("orange",0.2))
  
  lines(x=slopeCompHi[1,], y=slopeHi_mean,
        lwd=2,
        col="green")
  polygon(x=c(slopeCompHi[1,],rev(slopeCompHi[1,])),
          y=c(slopeHi_min,
              rev(slopeHi_max)),
          border=NA, col=adjustcolor("green",0.2))
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::mask(topoVals, soil[soil$SOIL=="Zetek",])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18b <- density(topoVals[!is.na(all18)],na.rm=T)
  dens19b <- density(topoVals[!is.na(all19)],na.rm=T)
  dens20b <- density(topoVals[!is.na(all20)],na.rm=T)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 3,
       xlim=range(slopeComp18[1,]),
       ylim=range(dens18b$y),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5,
       col = adjustcolor("black",0.5))
  lines(dens19, col = adjustcolor("black",0.5), lwd = 3)
  lines(dens20, col= adjustcolor("black",0.5), lwd = 3)
  lines(dens18b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens19b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens20b, col = adjustcolor("green",0.5), lwd = 3)
  
  
  #save(slopeComp18, slopeComp19, slopeComp20, slopeCompLo, slopeCompHi, file="prelimSlopeComp_Zetek.RData")
  
#### Zetek soil: distribution in gaps vs non gaps for best smoothing scale: curvature ####
  
  curvComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Zetek",])
  
  curvComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[2]],
                          allLayer = allLayerList[[2]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Zetek",])
  
  curvComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[3]],
                          allLayer = allLayerList[[3]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Zetek",])
  
  curvCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryLoCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Zetek",])
  
  curvCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryHiCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = soil[soil$SOIL=="Zetek",])
  
  #load("prelimCurvComp_Zetek.RData")
  
  curv18_mean <- apply(curvComp18[-1,], MARGIN = 2, FUN=mean)
  curv18_min <- apply(curvComp18[-1,], MARGIN = 2, min)
  curv18_max <- apply(curvComp18[-1,], MARGIN = 2, max)
  
  curv19_mean <- apply(curvComp19[-1,], MARGIN = 2, FUN=mean)
  curv19_min <- apply(curvComp19[-1,], MARGIN = 2, min)
  curv19_max <- apply(curvComp19[-1,], MARGIN = 2, max)
  
  curv20_mean <- apply(curvComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  curv20_min <- apply(curvComp20[-1,], MARGIN = 2, min)*(12/13)
  curv20_max <- apply(curvComp20[-1,], MARGIN = 2, max)*(12/13)
  
  curvLo_mean <- apply(curvCompLo[-1,], MARGIN = 2, FUN=mean)
  curvLo_min <- apply(curvCompLo[-1,], MARGIN = 2, min)
  curvLo_max <- apply(curvCompLo[-1,], MARGIN = 2, max)
  
  curvHi_mean <- apply(curvCompHi[-1,], MARGIN = 2, FUN=mean)
  curvHi_min <- apply(curvCompHi[-1,], MARGIN = 2, min)
  curvHi_max <- apply(curvCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=curvComp18[1,], y=curv18_mean,
       ylab=NA, cex.axis = 1.5,
       xlim=c(-2.7,3),
       type="l", lwd=2,
       ylim=c(0.6e-2,0.05),
       log="y")
  
  polygon(x=c(curvComp18[1,],rev(curvComp18[1,])),
          y=c(curv18_min,
              rev(curv18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=curvComp19[1,], y=curv19_mean,
        lwd=2, col = "blue")
  polygon(x=c(curvComp19[1,],rev(curvComp19[1,])),
          y=c(curv19_min,
              rev(curv19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=curvComp20[1,], y=curv20_mean,
        lwd=2,
        col="red")
  polygon(x=c(curvComp20[1,],rev(curvComp20[1,])),
          y=c(curv20_min,
              rev(curv20_max)),
          border=NA, col=adjustcolor("red",0.2))
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::mask(topoVals, soil[soil$SOIL=="Zetek",])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18b <- density(topoVals[!is.na(all18)],na.rm=T)
  dens19b <- density(topoVals[!is.na(all19)],na.rm=T)
  dens20b <- density(topoVals[!is.na(all20)],na.rm=T)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 3,
       xlim=range(curvComp18[1,]),
       ylim=range(dens18b$y),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5,
       col = adjustcolor("black",0.5))
  lines(dens19, col = adjustcolor("black",0.5), lwd = 3)
  lines(dens20, col= adjustcolor("black",0.5), lwd = 3)
  lines(dens18b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens19b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens20b, col = adjustcolor("green",0.5), lwd = 3)
  
  #save(curvComp18, curvComp19, curvComp20, curvCompLo, curvCompHi, file="prelimCurvComp_Zetek.RData")

#### Zetek soil: distribution in gaps vs non gaps for best smoothing scale: aspect ####
  aspectComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = gapLayerList[[1]],
                            allLayer = allLayerList[[1]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Zetek",])
  
  aspectComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = gapLayerList[[2]],
                            allLayer = allLayerList[[2]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Zetek",])
  
  aspectComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = gapLayerList[[3]],
                            allLayer = allLayerList[[3]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Zetek",])
  
  aspectCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = raster::raster("binaryLoCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_39.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Zetek",])
  
  aspectCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif",
                            gapLayer = raster::raster("binaryHiCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_39.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = soil[soil$SOIL=="Zetek",])
  
  #load("prelimAspectComp_Zetek.RData")
  
  aspect18_mean <- apply(aspectComp18[-1,], MARGIN = 2, FUN=mean)
  aspect18_min <- apply(aspectComp18[-1,], MARGIN = 2, min)
  aspect18_max <- apply(aspectComp18[-1,], MARGIN = 2, max)
  
  aspect19_mean <- apply(aspectComp19[-1,], MARGIN = 2, FUN=mean)
  aspect19_min <- apply(aspectComp19[-1,], MARGIN = 2, min)
  aspect19_max <- apply(aspectComp19[-1,], MARGIN = 2, max)
  
  aspect20_mean <- apply(aspectComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  aspect20_min <- apply(aspectComp20[-1,], MARGIN = 2, min)*(12/13)
  aspect20_max <- apply(aspectComp20[-1,], MARGIN = 2, max)*(12/13)
  
  aspectLo_mean <- apply(aspectCompLo[-1,], MARGIN = 2, FUN=mean)
  aspectLo_min <- apply(aspectCompLo[-1,], MARGIN = 2, min)
  aspectLo_max <- apply(aspectCompLo[-1,], MARGIN = 2, max)
  
  aspectHi_mean <- apply(aspectCompHi[-1,], MARGIN = 2, FUN=mean)
  aspectHi_min <- apply(aspectCompHi[-1,], MARGIN = 2, min)
  aspectHi_max <- apply(aspectCompHi[-1,], MARGIN = 2, max)
  
  par(mar=c(3,5,1,1),oma=c(0,0,0,0),las=1)
  plot(x=aspectComp18[1,], y=aspect18_mean,
       ylab=NA, cex.axis = 1.5,
       type="l", lwd=2,
       ylim=c(0.3e-2,0.11),
       log="y")
  polygon(x=c(aspectComp18[1,],rev(aspectComp18[1,])),
          y=c(aspect18_min,
              rev(aspect18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=aspectComp19[1,], y=aspect19_mean,
        lwd=2, col = "blue")
  polygon(x=c(aspectComp19[1,],rev(aspectComp19[1,])),
          y=c(aspect19_min,
              rev(aspect19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=aspectComp20[1,], y=aspect20_mean,
        lwd=2,
        col="red")
  polygon(x=c(aspectComp20[1,],rev(aspectComp20[1,])),
          y=c(aspect20_min,
              rev(aspect20_max)),
          border=NA, col=adjustcolor("red",0.2))
  
  abline(v=c(0,90,180,270,360), lty=2)
  
  # Make topography density plot
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18 <- density(topoVals[!is.na(all18)])
  dens19 <- density(topoVals[!is.na(all19)])
  dens20 <- density(topoVals[!is.na(all20)])
  
  topoVals <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_39.tif")
  topoVals <- raster::crop(topoVals, allLayerList[[1]])
  topoVals <- raster::mask(topoVals, soil[soil$SOIL=="Zetek",])
  topoVals <- raster::values(topoVals)
  all18 <- raster::values(allLayerList[[1]])
  all19 <- raster::values(allLayerList[[2]])
  all20 <- raster::values(allLayerList[[3]])
  dens18b <- density(topoVals[!is.na(all18)],na.rm=T)
  dens19b <- density(topoVals[!is.na(all19)],na.rm=T)
  dens20b <- density(topoVals[!is.na(all20)],na.rm=T)
  
  par(las=1, mar=c(3,4,1,1),oma=c(0,0,0,0))
  plot(dens18, lwd = 3,
       xlim=range(aspectComp18[1,]),
       ylim=range(dens18b$y),
       bty="n",
       main=NA,
       xlab=NA,ylab=NA,
       cex.axis=1.5,
       col = adjustcolor("black",0.5))
  lines(dens19, col = adjustcolor("black",0.5), lwd = 3)
  lines(dens20, col= adjustcolor("black",0.5), lwd = 3)
  lines(dens18b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens19b, col = adjustcolor("green",0.5), lwd = 3)
  lines(dens20b, col = adjustcolor("green",0.5), lwd = 3)
  
  
  #save(aspectComp18, aspectComp19, aspectComp20, aspectCompLo, aspectCompHi, file="prelimAspectComp_Zetek.RData")
  
#### Young forest: distribution in gaps vs non gaps for best smoothing scale: slope ####
  
  slopeComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[1]],
                           allLayer = allLayerList[[1]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[2]],
                           allLayer = allLayerList[[2]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[3]],
                           allLayer = allLayerList[[3]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryLoCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryHiCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  #load("prelimSlopeComp_Secondary.RData")
  
  slope18_mean <- apply(slopeComp18[-1,], MARGIN = 2, mean)
  slope18_min <- apply(slopeComp18[-1,], MARGIN = 2, min)
  slope18_max <- apply(slopeComp18[-1,], MARGIN = 2, max)
  
  slope19_mean <- apply(slopeComp19[-1,], MARGIN = 2, mean)
  slope19_min <- apply(slopeComp19[-1,], MARGIN = 2, min)
  slope19_max <- apply(slopeComp19[-1,], MARGIN = 2, max)
  
  slope20_mean <- apply(slopeComp20[-1,], MARGIN = 2, mean)
  slope20_min <- apply(slopeComp20[-1,], MARGIN = 2, min)
  slope20_max <- apply(slopeComp20[-1,], MARGIN = 2, max)
  
  slopeLo_mean <- apply(slopeCompLo[-1,], MARGIN = 2, mean)
  slopeLo_min <- apply(slopeCompLo[-1,], MARGIN = 2, min)
  slopeLo_max <- apply(slopeCompLo[-1,], MARGIN = 2, max)
  
  slopeHi_mean <- apply(slopeCompHi[-1,], MARGIN = 2, mean)
  slopeHi_min <- apply(slopeCompHi[-1,], MARGIN = 2, min)
  slopeHi_max <- apply(slopeCompHi[-1,], MARGIN = 2, max)
  
  
  plot(x=slopeComp18[1,], y=slope18_mean,
       type="l", lwd=2,
       ylim=c(1e-2,0.06),
       log="y")
  polygon(x=c(slopeComp18[1,],rev(slopeComp18[1,])),
          y=c(slope18_min,
              rev(slope18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=slopeComp19[1,], y=slope19_mean,
        lwd=2, col = "blue")
  polygon(x=c(slopeComp19[1,],rev(slopeComp19[1,])),
          y=c(slope19_min,
              rev(slope19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=slopeComp20[1,], y=slope20_mean,
        lwd=2,
        col="red")
  polygon(x=c(slopeComp20[1,],rev(slopeComp20[1,])),
          y=c(slope20_min,
              rev(slope20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=gapPct17to18, col="black", lty=2, lwd=2)
  abline(h=gapPct18to19, col="blue", lty=2, lwd=2)
  abline(h=gapPct19to20, col="red", lty=2, lwd=2)
  
  lines(x=slopeCompLo[1,], y=slopeLo_mean,
        lwd=2,
        col="orange")
  polygon(x=c(slopeCompLo[1,],rev(slopeCompLo[1,])),
          y=c(slopeLo_min,
              rev(slopeLo_max)),
          border=NA, col=adjustcolor("orange",0.2))
  
  lines(x=slopeCompHi[1,], y=slopeHi_mean,
        lwd=2,
        col="green")
  polygon(x=c(slopeCompHi[1,],rev(slopeCompHi[1,])),
          y=c(slopeHi_min,
              rev(slopeHi_max)),
          border=NA, col=adjustcolor("green",0.2))
  
  #save(slopeComp18, slopeComp19, slopeComp20, slopeCompLo, slopeCompHi, file="prelimSlopeComp_Secondary.RData")
  
#### Young forest: distribution in gaps vs non gaps for best smoothing scale: curvature ####
  curvComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[2]],
                          allLayer = allLayerList[[2]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[3]],
                          allLayer = allLayerList[[3]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryLoCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryHiCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  #load("prelimCurvComp_OldGrowth.RData")
  
  curv18_mean <- apply(curvComp18[-1,], MARGIN = 2, FUN=mean)
  curv18_min <- apply(curvComp18[-1,], MARGIN = 2, min)
  curv18_max <- apply(curvComp18[-1,], MARGIN = 2, max)
  
  curv19_mean <- apply(curvComp19[-1,], MARGIN = 2, FUN=mean)
  curv19_min <- apply(curvComp19[-1,], MARGIN = 2, min)
  curv19_max <- apply(curvComp19[-1,], MARGIN = 2, max)
  
  curv20_mean <- apply(curvComp20[-1,], MARGIN = 2, FUN=mean)
  curv20_min <- apply(curvComp20[-1,], MARGIN = 2, min)
  curv20_max <- apply(curvComp20[-1,], MARGIN = 2, max)
  
  curvLo_mean <- apply(curvCompLo[-1,], MARGIN = 2, FUN=mean)
  curvLo_min <- apply(curvCompLo[-1,], MARGIN = 2, min)
  curvLo_max <- apply(curvCompLo[-1,], MARGIN = 2, max)
  
  curvHi_mean <- apply(curvCompHi[-1,], MARGIN = 2, FUN=mean)
  curvHi_min <- apply(curvCompHi[-1,], MARGIN = 2, min)
  curvHi_max <- apply(curvCompHi[-1,], MARGIN = 2, max)
  
  plot(x=curvComp18[1,], y=curv18_mean,
       type="l", lwd=2,
       ylim=c(1e-2,0.05),
       log="y")
  polygon(x=c(curvComp18[1,],rev(curvComp18[1,])),
          y=c(curv18_min,
              rev(curv18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=curvComp19[1,], y=curv19_mean,
        lwd=2, col = "blue")
  polygon(x=c(curvComp19[1,],rev(curvComp19[1,])),
          y=c(curv19_min,
              rev(curv19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=curvComp20[1,], y=curv20_mean,
        lwd=2,
        col="red")
  polygon(x=c(curvComp20[1,],rev(curvComp20[1,])),
          y=c(curv20_min,
              rev(curv20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=1,lty=2)
  
  #load("prelimCurvComp_Secondary.RData")
  
  curv18_mean <- apply(curvComp18[-1,], MARGIN = 2, FUN=mean)
  curv18_min <- apply(curvComp18[-1,], MARGIN = 2, min)
  curv18_max <- apply(curvComp18[-1,], MARGIN = 2, max)
  
  curv19_mean <- apply(curvComp19[-1,], MARGIN = 2, FUN=mean)
  curv19_min <- apply(curvComp19[-1,], MARGIN = 2, min)
  curv19_max <- apply(curvComp19[-1,], MARGIN = 2, max)
  
  curv20_mean <- apply(curvComp20[-1,], MARGIN = 2, FUN=mean)
  curv20_min <- apply(curvComp20[-1,], MARGIN = 2, min)
  curv20_max <- apply(curvComp20[-1,], MARGIN = 2, max)
  
  curvLo_mean <- apply(curvCompLo[-1,], MARGIN = 2, FUN=mean)
  curvLo_min <- apply(curvCompLo[-1,], MARGIN = 2, min)
  curvLo_max <- apply(curvCompLo[-1,], MARGIN = 2, max)
  
  curvHi_mean <- apply(curvCompHi[-1,], MARGIN = 2, FUN=mean)
  curvHi_min <- apply(curvCompHi[-1,], MARGIN = 2, min)
  curvHi_max <- apply(curvCompHi[-1,], MARGIN = 2, max)
  
  plot(x=curvComp18[1,], y=curv18_mean,
       type="l", lwd=2,
       ylim=c(0,0.1))
  polygon(x=c(curvComp18[1,],rev(curvComp18[1,])),
          y=c(curv18_min,
              rev(curv18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=curvComp19[1,], y=curv19_mean,
        lwd=2, col = "blue")
  polygon(x=c(curvComp19[1,],rev(curvComp19[1,])),
          y=c(curv19_min,
              rev(curv19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=curvComp20[1,], y=curv20_mean,
        lwd=2,
        col="red")
  polygon(x=c(curvComp20[1,],rev(curvComp20[1,])),
          y=c(curv20_min,
              rev(curv20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=1,lty=2)
  
  
  #save(curvComp18, curvComp19, curvComp20, curvCompLo, curvCompHi, file="prelimCurvComp_Secondary.RData")
  
#### Young forest: distribution in gaps vs non gaps for best smoothing scale: aspect ####
  aspectComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = gapLayerList[[1]],
                            allLayer = allLayerList[[1]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = gapLayerList[[2]],
                            allLayer = allLayerList[[2]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = gapLayerList[[3]],
                            allLayer = allLayerList[[3]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = raster::raster("binaryLoCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = raster::raster("binaryHiCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  #load("prelimAspectComp_Secondary.RData")
  
  aspect18_mean <- apply(aspectComp18[-1,], MARGIN = 2, FUN=mean)
  aspect18_min <- apply(aspectComp18[-1,], MARGIN = 2, min)
  aspect18_max <- apply(aspectComp18[-1,], MARGIN = 2, max)
  
  aspect19_mean <- apply(aspectComp19[-1,], MARGIN = 2, FUN=mean)
  aspect19_min <- apply(aspectComp19[-1,], MARGIN = 2, min)
  aspect19_max <- apply(aspectComp19[-1,], MARGIN = 2, max)
  
  aspect20_mean <- apply(aspectComp20[-1,], MARGIN = 2, FUN=mean)
  aspect20_min <- apply(aspectComp20[-1,], MARGIN = 2, min)
  aspect20_max <- apply(aspectComp20[-1,], MARGIN = 2, max)
  
  aspectLo_mean <- apply(aspectCompLo[-1,], MARGIN = 2, FUN=mean)
  aspectLo_min <- apply(aspectCompLo[-1,], MARGIN = 2, min)
  aspectLo_max <- apply(aspectCompLo[-1,], MARGIN = 2, max)
  
  aspectHi_mean <- apply(aspectCompHi[-1,], MARGIN = 2, FUN=mean)
  aspectHi_min <- apply(aspectCompHi[-1,], MARGIN = 2, min)
  aspectHi_max <- apply(aspectCompHi[-1,], MARGIN = 2, max)
  
  plot(x=aspectComp18[1,], y=aspect18_mean,
       type="l", lwd=2,
       ylim=c(0.7e-2,0.04),
       log="y")
  polygon(x=c(aspectComp18[1,],rev(aspectComp18[1,])),
          y=c(aspect18_min,
              rev(aspect18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=aspectComp19[1,], y=aspect19_mean,
        lwd=2, col = "blue")
  polygon(x=c(aspectComp19[1,],rev(aspectComp19[1,])),
          y=c(aspect19_min,
              rev(aspect19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=aspectComp20[1,], y=aspect20_mean,
        lwd=2,
        col="red")
  polygon(x=c(aspectComp20[1,],rev(aspectComp20[1,])),
          y=c(aspect20_min,
              rev(aspect20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=1,lty=2)
  abline(v=c(0,90,180,270,360), lty=2)
  
  #save(aspectComp18, aspectComp19, aspectComp20, aspectCompLo, aspectCompHi, file="prelimAspectComp_Secondary.RData")
  
#### Old growth forest: distribution in gaps vs non gaps for best smoothing scale: slope ####
  
  slopeComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[1]],
                           allLayer = allLayerList[[1]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[2]],
                           allLayer = allLayerList[[2]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = gapLayerList[[3]],
                           allLayer = allLayerList[[3]],
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryLoCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  slopeCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif",
                           gapLayer = raster::raster("binaryHiCanopy.tif"),
                           allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif"),
                           sampleDir = "D:/BCI_Spatial/BlockSamples/",
                           nBins = 50,
                           mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  #load("prelimSlopeComp_OldGrowth.RData")
  
  slope18_mean <- apply(slopeComp18[-1,], MARGIN = 2, mean)
  slope18_min <- apply(slopeComp18[-1,], MARGIN = 2, min)
  slope18_max <- apply(slopeComp18[-1,], MARGIN = 2, max)
  
  slope19_mean <- apply(slopeComp19[-1,], MARGIN = 2, mean)
  slope19_min <- apply(slopeComp19[-1,], MARGIN = 2, min)
  slope19_max <- apply(slopeComp19[-1,], MARGIN = 2, max)
  
  slope20_mean <- apply(slopeComp20[-1,], MARGIN = 2, mean)*(12/13)
  slope20_min <- apply(slopeComp20[-1,], MARGIN = 2, min)*(12/13)
  slope20_max <- apply(slopeComp20[-1,], MARGIN = 2, max)*(12/13)
  
  slopeLo_mean <- apply(slopeCompLo[-1,], MARGIN = 2, mean)
  slopeLo_min <- apply(slopeCompLo[-1,], MARGIN = 2, min)
  slopeLo_max <- apply(slopeCompLo[-1,], MARGIN = 2, max)
  
  slopeHi_mean <- apply(slopeCompHi[-1,], MARGIN = 2, mean)
  slopeHi_min <- apply(slopeCompHi[-1,], MARGIN = 2, min)
  slopeHi_max <- apply(slopeCompHi[-1,], MARGIN = 2, max)
  
  
  plot(x=slopeComp18[1,], y=slope18_mean,
       type="l", lwd=2,
       ylim=c(1e-2,0.06),
       log="y")
  polygon(x=c(slopeComp18[1,],rev(slopeComp18[1,])),
          y=c(slope18_min,
              rev(slope18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=slopeComp19[1,], y=slope19_mean,
        lwd=2, col = "blue")
  polygon(x=c(slopeComp19[1,],rev(slopeComp19[1,])),
          y=c(slope19_min,
              rev(slope19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=slopeComp20[1,], y=slope20_mean,
        lwd=2,
        col="red")
  polygon(x=c(slopeComp20[1,],rev(slopeComp20[1,])),
          y=c(slope20_min,
              rev(slope20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=gapPct17to18, col="black", lty=2, lwd=2)
  abline(h=gapPct18to19, col="blue", lty=2, lwd=2)
  abline(h=gapPct19to20, col="red", lty=2, lwd=2)
  
  lines(x=slopeCompLo[1,], y=slopeLo_mean,
        lwd=2,
        col="orange")
  polygon(x=c(slopeCompLo[1,],rev(slopeCompLo[1,])),
          y=c(slopeLo_min,
              rev(slopeLo_max)),
          border=NA, col=adjustcolor("orange",0.2))
  
  lines(x=slopeCompHi[1,], y=slopeHi_mean,
        lwd=2,
        col="green")
  polygon(x=c(slopeCompHi[1,],rev(slopeCompHi[1,])),
          y=c(slopeHi_min,
              rev(slopeHi_max)),
          border=NA, col=adjustcolor("green",0.2))
  
  #save(slopeComp18, slopeComp19, slopeComp20, slopeCompLo, slopeCompHi, file="prelimSlopeComp_OldGrowth.RData")
  
#### Old growth forest: distribution in gaps vs non gaps for best smoothing scale: curvature ####
  curvComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[1]],
                          allLayer = allLayerList[[1]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[2]],
                          allLayer = allLayerList[[2]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = gapLayerList[[3]],
                          allLayer = allLayerList[[3]],
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryLoCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  curvCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif",
                          gapLayer = raster::raster("binaryHiCanopy.tif"),
                          allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                          sampleDir = "D:/BCI_Spatial/BlockSamples/",
                          nBins = 50,
                          mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  #load("prelimCurvComp_OldGrowth.RData")
  
  curv18_mean <- apply(curvComp18[-1,], MARGIN = 2, FUN=mean)
  curv18_min <- apply(curvComp18[-1,], MARGIN = 2, min)
  curv18_max <- apply(curvComp18[-1,], MARGIN = 2, max)
  
  curv19_mean <- apply(curvComp19[-1,], MARGIN = 2, FUN=mean)
  curv19_min <- apply(curvComp19[-1,], MARGIN = 2, min)
  curv19_max <- apply(curvComp19[-1,], MARGIN = 2, max)
  
  curv20_mean <- apply(curvComp20[-1,], MARGIN = 2, FUN=mean)*(12/13)
  curv20_min <- apply(curvComp20[-1,], MARGIN = 2, min)*(12/13)
  curv20_max <- apply(curvComp20[-1,], MARGIN = 2, max)*(12/13)
  
  curvLo_mean <- apply(curvCompLo[-1,], MARGIN = 2, FUN=mean)
  curvLo_min <- apply(curvCompLo[-1,], MARGIN = 2, min)
  curvLo_max <- apply(curvCompLo[-1,], MARGIN = 2, max)
  
  curvHi_mean <- apply(curvCompHi[-1,], MARGIN = 2, FUN=mean)
  curvHi_min <- apply(curvCompHi[-1,], MARGIN = 2, min)
  curvHi_max <- apply(curvCompHi[-1,], MARGIN = 2, max)
  
  plot(x=curvComp18[1,], y=curv18_mean,
       type="l", lwd=2,
       ylim=c(1e-2,0.05),
       log="y")
  polygon(x=c(curvComp18[1,],rev(curvComp18[1,])),
          y=c(curv18_min,
              rev(curv18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=curvComp19[1,], y=curv19_mean,
        lwd=2, col = "blue")
  polygon(x=c(curvComp19[1,],rev(curvComp19[1,])),
          y=c(curv19_min,
              rev(curv19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=curvComp20[1,], y=curv20_mean,
        lwd=2,
        col="red")
  polygon(x=c(curvComp20[1,],rev(curvComp20[1,])),
          y=c(curv20_min,
              rev(curv20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=1,lty=2)
  
  
  #save(curvComp18, curvComp19, curvComp20, curvCompLo, curvCompHi, file="prelimCurvComp_OldGrowth.RData")
  
#### Old growth forest: distribution in gaps vs non gaps for best smoothing scale: aspect ####
  aspectComp18 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = gapLayerList[[1]],
                            allLayer = allLayerList[[1]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectComp19 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = gapLayerList[[2]],
                            allLayer = allLayerList[[2]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectComp20 <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = gapLayerList[[3]],
                            allLayer = allLayerList[[3]],
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectCompLo <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = raster::raster("binaryLoCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  aspectCompHi <- ratioCalc(topoLayer = "D:/BCI_Spatial/BCI_Topo/Aspect_smooth_15.tif",
                            gapLayer = raster::raster("binaryHiCanopy.tif"),
                            allLayer = raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_15.tif"),
                            sampleDir = "D:/BCI_Spatial/BlockSamples/",
                            nBins = 50,
                            mask = buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  
  #load("prelimAspectComp_OldGrowth.RData")
  
  aspect18_mean <- apply(aspectComp18[-1,], MARGIN = 2, FUN=mean)
  aspect18_min <- apply(aspectComp18[-1,], MARGIN = 2, min)
  aspect18_max <- apply(aspectComp18[-1,], MARGIN = 2, max)
  
  aspect19_mean <- apply(aspectComp19[-1,], MARGIN = 2, FUN=mean)
  aspect19_min <- apply(aspectComp19[-1,], MARGIN = 2, min)
  aspect19_max <- apply(aspectComp19[-1,], MARGIN = 2, max)
  
  aspect20_mean <- apply(aspectComp20[-1,], MARGIN = 2, FUN=mean)
  aspect20_min <- apply(aspectComp20[-1,], MARGIN = 2, min)
  aspect20_max <- apply(aspectComp20[-1,], MARGIN = 2, max)
  
  aspectLo_mean <- apply(aspectCompLo[-1,], MARGIN = 2, FUN=mean)
  aspectLo_min <- apply(aspectCompLo[-1,], MARGIN = 2, min)
  aspectLo_max <- apply(aspectCompLo[-1,], MARGIN = 2, max)
  
  aspectHi_mean <- apply(aspectCompHi[-1,], MARGIN = 2, FUN=mean)
  aspectHi_min <- apply(aspectCompHi[-1,], MARGIN = 2, min)
  aspectHi_max <- apply(aspectCompHi[-1,], MARGIN = 2, max)
  
  plot(x=aspectComp18[1,], y=aspect18_mean,
       type="l", lwd=2,
       ylim=c(1e-2,0.06),
       log="y")
  polygon(x=c(aspectComp18[1,],rev(aspectComp18[1,])),
          y=c(aspect18_min,
              rev(aspect18_max)),
          border=NA, col=adjustcolor("black",0.2))
  
  lines(x=aspectComp19[1,], y=aspect19_mean,
        lwd=2, col = "blue")
  polygon(x=c(aspectComp19[1,],rev(aspectComp19[1,])),
          y=c(aspect19_min,
              rev(aspect19_max)),
          border=NA, col=adjustcolor("blue",0.2))
  
  lines(x=aspectComp20[1,], y=aspect20_mean,
        lwd=2,
        col="red")
  polygon(x=c(aspectComp20[1,],rev(aspectComp20[1,])),
          y=c(aspect20_min,
              rev(aspect20_max)),
          border=NA, col=adjustcolor("red",0.2))
  abline(h=1,lty=2)
  abline(v=c(0,90,180,270,360), lty=2)
  
  
  
  #save(aspectComp18, aspectComp19, aspectComp20, aspectCompLo, aspectCompHi, file="prelimAspectComp_OldGrowth.RData")
  
#### Average gap formation rates in space ####
  soilTypes$avgRate <- (1/2)*(soilTypes$PropGap18to19+soilTypes$PropGap19to20)

  manual.col <- colorRampPalette(c("#b5b5b5","#57471d"))
  soilTypes$col <- manual.col(length(unique(soilTypes$avgRate)))
  
  soil$plotCol <- "pink"
  soil$avgRate <- NA
  for(i in 1:dim(soilTypes)[1]){
    if(soilTypes$Good19[i]==T & soilTypes$Good18[i]==T){
      soil[soil$SOIL==soilTypes$Soil[i],"plotCol"] <- soilTypes$col[i]
      soil[soil$SOIL==soilTypes$Soil[i],"avgRate"] <- soilTypes$avgRate[i]
      
    }
  }
  
  soilRaster <- raster::rasterize(x=soil, y=chm17, field="avgRate")
  raster::plot(soilRaster,
               col= manual.col(length(unique(soilTypes$avgRate))),
               box = F, axes=F)
  

  
  # Do soil type patterns match low and hi forest patterns?
  par(las=1,mar=c(3,5,1,1),oma=c(0,0,0,0))
  plot(PropLo~avgRate, data = soilTypes[soilTypes$Good19==T,],
       pch=20,cex=2,
       xlab=NA,ylab=NA,cex.axis=1.5)
  summary(lm(PropLo~avgRate, data = soilTypes[soilTypes$Good19==T,]))
  
  plot(PropHi~avgRate, data = soilTypes[soilTypes$Good19==T,],
       pch=20,cex=2,
       xlab=NA,ylab=NA,cex.axis=1.5)
  summary(lm(PropHi~avgRate, data = soilTypes[soilTypes$Good19==T,]))
  
  # Only old growth
  par(las=1,mar=c(3,5,1,1),oma=c(0,0,0,0))
  plot(PropLo_Old~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Old>minSz,],
       pch=20,
       xlab=NA,ylab=NA,cex.axis=1.5)
  summary(lm(PropLo_Old~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Old>minSz,],))
  
  plot(PropHi_Old~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Old>minSz,],
       pch=20, cex = 2,
       xlab=NA,ylab=NA,cex.axis=1.5)
  summary(lm(PropHi_Old~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Old>minSz,],))
  
  # Only secondary
  par(las=1,mar=c(3,5,1,1),oma=c(0,0,0,0))
  plot(PropLo_Sec~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Sec>minSz,],
       pch=20,
       xlab=NA,ylab=NA,cex.axis=1.5)
  summary(lm(PropLo_Sec~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Sec>minSz,],))
  
  plot(PropHi_Sec~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Sec>minSz,],
       pch=20,
       xlab=NA,ylab=NA,cex.axis=1.5)
  summary(lm(PropHi_Sec~avgRate, data = soilTypes[soilTypes$Good19==T & soilTypes$Area_Sec>minSz,],))
  
  
#### Correlations between average formation rates and topography ####
  
  # Slope
  
  soilTypes$avgSlope <- NA
  for(i in 1:dim(soilTypes)[1]){
    slopei <- raster::crop(slope,soil[soil$SOIL==soilTypes$Soil[i],])
    values_slopei <- raster::getValues(slopei)
    soilTypes$avgSlope[i] <- mean(values_slopei)
  }
  
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  plot(avgRate~avgSlope, data = soilTypes[soilTypes$Good==T,],
       pch=20,
       xlab="Average slope", ylab="Average gap formation rate")
  summary(lm(avgRate~avgSlope, data = soilTypes[soilTypes$Good==T,]))
  
  # Aspect
  # East-west (value of 1 is facing east, value of -1 is facing due west)
  # North-south (value of 1 is facing north, value of -1 is facing due south)
  soilTypes$avgEastWest <- NA
  soilTypes$avgNorthSouth <- NA
  
  for(i in 1:dim(soilTypes)[1]){
    aspecti <- raster::crop(aspect,soil[soil$SOIL==soilTypes$Soil[i],])
    
    valuesi <- raster::getValues(aspecti)
    valuesi <- cos(valuesi*(pi/180))
    soilTypes$avgEastWest[i] <- mean(valuesi)
    
    valuesi <- raster::getValues(aspecti)
    valuesi <- sin(valuesi*(pi/180))
    soilTypes$avgNorthSouth[i] <- mean(valuesi)
    
  }
  
  par(mfrow=c(1,2),mar=c(4,4,1,1))
  plot(avgRate~avgEastWest, data = soilTypes[soilTypes$Good==T,],
       pch=20,
       xlab="Average East-West aspect", ylab="Average gap formation rate")
  summary(lm(avgRate~avgEastWest, data = soilTypes[soilTypes$Good==T,]))
  
  plot(avgRate~avgNorthSouth, data = soilTypes[soilTypes$Good==T,],
       pch=20,
       xlab="Average North-South aspect", ylab=NA)
  summary(lm(avgRate~avgNorthSouth, data = soilTypes[soilTypes$Good==T,]))
  
  
#### Calculate slope and easting for each gap ####
  IDs18 <- raster::unique(gaps17to18)
  IDs19 <- raster::unique(gaps18to19)
  IDs20 <- raster::unique(gaps19to20)
  
  gapTopo <- data.frame(gapID = c(IDs18, IDs19, IDs20),
                        year = c(rep(2018, length(IDs18)),
                                 rep(2019, length(IDs19)),
                                 rep(2020, length(IDs20))),
                        gapArea = NA,
                        slope = NA,
                        easting = NA,
                        northing = NA)
  
  gapVals17to18 <- raster::values(gaps17to18)
  gapVals18to19 <- raster::values(gaps18to19)
  gapVals19to20 <- raster::values(gaps19to20)
  slopeVals <- raster::values(slope)
  eastVals <- raster::values (east)
  northVals <- raster::values(north)
  
  
  for(i in 1:dim(gapTopo)[1]){
    if(gapTopo$year[i]==2018){
      # calculate gap area
      gapTopo$gapArea[i] <- length(gapVals17to18[gapVals17to18==gapTopo$gapID[i]
                                                 & !is.na(gapVals17to18)])
      
      # calculate average slope in gap
      gapTopo$slope[i] <- mean(slopeVals[gapVals17to18==gapTopo$gapID[i]
                                         & !is.na(gapVals17to18)])
      
      # calculate average easting in gap
      gapTopo$easting[i] <- mean(eastVals[gapVals17to18==gapTopo$gapID[i]
                                          & !is.na(gapVals17to18)])
      
      # calcuate average northing in gap
      gapTopo$northing[i] <- mean(northVals[gapVals17to18==gapTopo$gapID[i]
                                            & !is.na(gapVals17to18)])
    }
    
    if(gapTopo$year[i]==2019){
      # calculate gap area
      gapTopo$gapArea[i] <- length(gapVals18to19[gapVals18to19==gapTopo$gapID[i]
                                                 & !is.na(gapVals18to19)])
      
      # calculate average slope in gap
      gapTopo$slope[i] <- mean(slopeVals[gapVals18to19==gapTopo$gapID[i]
                                         & !is.na(gapVals18to19)])
      
      # calculate average easting in gap
      gapTopo$easting[i] <- mean(eastVals[gapVals18to19==gapTopo$gapID[i]
                                          & !is.na(gapVals18to19)])
      
      # calcuate average northing in gap
      gapTopo$northing[i] <- mean(northVals[gapVals18to19==gapTopo$gapID[i]
                                            & !is.na(gapVals18to19)])
    }
    
    if(gapTopo$year[i]==2020){
      # calculate gap area
      gapTopo$gapArea[i] <- length(gapVals19to20[gapVals19to20==gapTopo$gapID[i]
                                                 & !is.na(gapVals19to20)])
      
      # calculate average slope in gap
      gapTopo$slope[i] <- mean(slopeVals[gapVals19to20==gapTopo$gapID[i]
                                         & !is.na(gapVals19to20)])
      
      # calculate average easting in gap
      gapTopo$easting[i] <- mean(eastVals[gapVals19to20==gapTopo$gapID[i]
                                          & !is.na(gapVals19to20)])
      
      # calcuate average northing in gap
      gapTopo$northing[i] <- mean(northVals[gapVals19to20==gapTopo$gapID[i]
                                            & !is.na(gapVals19to20)])
    }
    print(i)
    
  }
  
  write.csv(gapTopo, file="gapTopoStats.csv", row.names = F)
  
  
  
#### Set up functions for MCMC routine ####
  # Load data
  gapData <- read.csv("gapTopoStats.csv")
  
  
  # library(Rlab)
  require(VGAM)
  require(poweRlaw)
  
  # define minimum gap size
  minGap <- 10
  
  calcul_lambda <- function(theta, var)
  {  
    pred =1+exp(theta[1]+theta[2]*var)  #Univariate model cf (eq.8) of paper
    return(pred)
  }
  
  fct <- function(x)
    return(dpldis(x[1],minGap,x[2],log=T))
  
  loglikelihood_pareto= function(data, xmin =minGap ,theta,var){
    lambda_pred = calcul_lambda(theta,var)
    d_sim = cbind(data,lambda_pred)
    temp = apply(d_sim,1,fct)
    LL= as.matrix(temp)
    ll_area_tot <- colSums(LL)
    return(ll_area_tot)
  }
  
#### Do MCMC univariate: slope ####
  # Make a vector of gap areas (must convert to integer if not already)
  d <- gapData[gapData$gapArea>=minGap & !is.na(gapData$slope),"gapArea"]
  
  # Transform slope data
  p <- sqrt(gapData[gapData$gapArea>=minGap & !is.na(gapData$slope),"slope"])
  var <- c(scale(p,center= TRUE,scale= TRUE))
  
  #________________________________
  #       START with model paramter (THETA)
  # At this step all parameter must be adjusted, 
  #_________________________________ 
  theta <- c(0.5,0.4)
  nb_param <- length(theta)
  step <- c(0.05,0.05)
  count= 0
  nb_it <- 1000 #number of iterations
  nb_sample <- 2    
  m_prior <- 0.5
  sd_prior <- 500
  theta <- matrix(theta,nrow=nb_param,ncol = 2,byrow=F)
  v_e1 <- var
  save_res  = "TRUE" # TRUE or FALSE  # If we want to save results
  
  #_________________________________ 
  #      COMPUTE LILKELIHOO_PARETO for the fist elements of theta and var     
  #_________________________________ 
  llp = loglikelihood_pareto(data = d  ,xmin = minGap ,theta,var)
  llp_res=c(llp, rep(0, nb_it))
  
  #________________________________
  #      Prepare OUTPUT             
  #________________________________
  accept <- rep(0,nb_param)
  res <- array(dim=c(nb_it+1,nb_param,nb_sample))
  res[1,,] <- theta
  
  #________________________________
  #   CHAINE DE Metropolis_Hastings      
  #________________________________
  for(compt in 1:nb_it)
  {
    
    ordre <- sample(1 : nb_param)
    
    for (i in ordre)
    {
      
      theta_new <- theta #candidat
      theta_new[i,] <- rnorm(2,theta[i,],step[i])
      
      llp_new = loglikelihood_pareto(d,xmin=minGap,theta=theta_new, var = var)
      
      logr = llp_new - llp +dnorm(theta[i,],theta_new[i,],step[i],log=T) - 
        dnorm(theta_new[i,],theta[i,],step[i],log=T) +
        dnorm(theta_new[i,], m_prior, sd_prior, log=T) - 
        dnorm(theta[i,], m_prior, sd_prior, log=T)
      r <- exp(logr)
      
      nb_ok <- which(!is.na(r))
      r_ok <- na.omit(r)
      
      accepte <- nb_ok[(runif(length(r_ok),min=0,max=1) < r_ok)]
      llp_res[compt+1]=llp_new
      theta[,accepte] <- theta_new[,accepte]
      llp <- llp_new
    }
    
    res[compt+1,,] <- theta
    print(compt)
  }
  
  #________________________________
  #      Plot results 
  #________________________________
  windows(10,10)
  par(mfrow=c(2,2),mar=c(3,2,2,2))
  
  theta_est = res[,,1]
  plot(theta_est[100:1000,1],col="blue",type = "l",main ="Theta 0" )
  plot(theta_est[100:1000,2],col="magenta",type = "l",main = "Theta 1")
  theta.1 = theta_est[100:1000,1]
  plot(density(theta.1))
  theta.2 = theta_est[100:1000,2]
  plot(density(theta.2),
       main = "Slope effect")
  abline(v=0,lty=2)
  
  #_____________________________
  # Save result
  #_____________________________
  if (save_res == TRUE) {
    write.csv(res,file="slope_MCMC.csv")
  }
  
  
  
#### Do MCMC univariate: easting ####
  # Make a vector of gap areas (must convert to integer if not already)
  d <- gapData[gapData$gapArea>=minGap & !is.na(gapData$easting),"gapArea"]
  
  # Transform slope data
  p <- gapData[gapData$gapArea>=minGap & !is.na(gapData$easting),"easting"]
  var <- c(scale(p,center= TRUE,scale= TRUE))
  
  #________________________________
  #       START with model paramter (THETA)
  # At this step all parameter must be adjusted, 
  #_________________________________ 
  theta <- c(0.5,0.4)
  nb_param <- length(theta)
  step <- c(0.05,0.05)
  count= 0
  nb_it <- 1000 #number of iterations
  nb_sample <- 2    
  m_prior <- 0.5
  sd_prior <- 500
  theta <- matrix(theta,nrow=nb_param,ncol = 2,byrow=F)
  v_e1 <- var
  save_res  = "TRUE" # TRUE or FALSE  # If we want to save results
  
  #_________________________________ 
  #      COMPUTE LILKELIHOO_PARETO for the fist elements of theta and var     
  #_________________________________ 
  llp = loglikelihood_pareto(data = d  ,xmin = minGap ,theta,var)
  llp_res=c(llp, rep(0, nb_it))
  
  #________________________________
  #      Prepare OUTPUT             
  #________________________________
  accept <- rep(0,nb_param)
  res <- array(dim=c(nb_it+1,nb_param,nb_sample))
  res[1,,] <- theta
  
  #________________________________
  #   CHAINE DE Metropolis_Hastings      
  #________________________________
  for(compt in 1:nb_it)
  {
    
    ordre <- sample(1 : nb_param)
    
    for (i in ordre)
    {
      
      theta_new <- theta #candidat
      theta_new[i,] <- rnorm(2,theta[i,],step[i])
      
      llp_new = loglikelihood_pareto(d,xmin=minGap,theta=theta_new, var = var)
      
      logr = llp_new - llp +dnorm(theta[i,],theta_new[i,],step[i],log=T) - 
        dnorm(theta_new[i,],theta[i,],step[i],log=T) +
        dnorm(theta_new[i,], m_prior, sd_prior, log=T) - 
        dnorm(theta[i,], m_prior, sd_prior, log=T)
      r <- exp(logr)
      
      nb_ok <- which(!is.na(r))
      r_ok <- na.omit(r)
      
      accepte <- nb_ok[(runif(length(r_ok),min=0,max=1) < r_ok)]
      llp_res[compt+1]=llp_new
      theta[,accepte] <- theta_new[,accepte]
      llp <- llp_new
    }
    
    res[compt+1,,] <- theta
    print(compt)
  }
  
  #________________________________
  #      Plot results 
  #________________________________
  windows(10,10)
  par(mfrow=c(2,2),mar=c(3,2,2,2))
  
  theta_est = res[,,1]
  plot(theta_est[,1],col="blue",type = "l",main ="Theta 0" )
  plot(theta_est[,2],col="magenta",type = "l",main = "Theta 1")
  theta.1 = theta_est[100:1000,1]
  plot(density(theta.1))
  theta.2 = theta_est[100:1000,2]
  plot(density(theta.2),
       main = "Easting effect")
  abline(v=0,lty=2)
  
  #_____________________________
  # Save result
  #_____________________________
  if (save_res == TRUE) {
    write.csv(res,file="easting_MCMC.csv")
  }
  
  
  
  
#### Do MCMC univariate: northing ####
  # Make a vector of gap areas (must convert to integer if not already)
  d <- gapData[gapData$gapArea>=minGap & !is.na(gapData$northing),"gapArea"]
  
  # Transform slope data
  p <- gapData[gapData$gapArea>=minGap & !is.na(gapData$northing),"northing"]
  var <- c(scale(p,center= TRUE,scale= TRUE))
  
  #________________________________
  #       START with model paramter (THETA)
  # At this step all parameter must be adjusted, 
  #_________________________________ 
  theta <- c(0.5,0.4)
  nb_param <- length(theta)
  step <- c(0.05,0.05)
  count= 0
  nb_it <- 1000 #number of iterations
  nb_sample <- 2    
  m_prior <- 0.5
  sd_prior <- 500
  theta <- matrix(theta,nrow=nb_param,ncol = 2,byrow=F)
  v_e1 <- var
  save_res  = "TRUE" # TRUE or FALSE  # If we want to save results
  
  #_________________________________ 
  #      COMPUTE LILKELIHOO_PARETO for the fist elements of theta and var     
  #_________________________________ 
  llp = loglikelihood_pareto(data = d  ,xmin = minGap ,theta,var)
  llp_res=c(llp, rep(0, nb_it))
  
  #________________________________
  #      Prepare OUTPUT             
  #________________________________
  accept <- rep(0,nb_param)
  res <- array(dim=c(nb_it+1,nb_param,nb_sample))
  res[1,,] <- theta
  
  #________________________________
  #   CHAINE DE Metropolis_Hastings      
  #________________________________
  for(compt in 1:nb_it)
  {
    
    ordre <- sample(1 : nb_param)
    
    for (i in ordre)
    {
      
      theta_new <- theta #candidat
      theta_new[i,] <- rnorm(2,theta[i,],step[i])
      
      llp_new = loglikelihood_pareto(d,xmin=minGap,theta=theta_new, var = var)
      
      logr = llp_new - llp +dnorm(theta[i,],theta_new[i,],step[i],log=T) - 
        dnorm(theta_new[i,],theta[i,],step[i],log=T) +
        dnorm(theta_new[i,], m_prior, sd_prior, log=T) - 
        dnorm(theta[i,], m_prior, sd_prior, log=T)
      r <- exp(logr)
      
      nb_ok <- which(!is.na(r))
      r_ok <- na.omit(r)
      
      accepte <- nb_ok[(runif(length(r_ok),min=0,max=1) < r_ok)]
      llp_res[compt+1]=llp_new
      theta[,accepte] <- theta_new[,accepte]
      llp <- llp_new
    }
    
    res[compt+1,,] <- theta
    print(compt)
  }
  
  #________________________________
  #      Plot results 
  #________________________________
  windows(10,10)
  par(mfrow=c(2,2),mar=c(3,2,2,2))
  
  theta_est = res[,,1]
  plot(theta_est[,1],col="blue",type = "l",main ="Theta 0" )
  plot(theta_est[,2],col="magenta",type = "l",main = "Theta 1")
  theta.1 = theta_est[100:1000,1]
  plot(density(theta.1))
  theta.2 = theta_est[100:1000,2]
  plot(density(theta.2),
       main = "Northing effect")
  abline(v=0,lty=2)
  
  #_____________________________
  # Save result
  #_____________________________
  if (save_res == TRUE) {
    write.csv(res,file="northing_MCMC.csv")
  }
  
#### Forest age ####
  age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  
  # Make new age class
  age@data$AgeClass <- "Other"
  age@data$AgeClass[age@data$Mascaro_Co == "> 400"] <- "OldGrowth"
  age@data$AgeClass[age@data$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
  
  oldGrowth <- raster::aggregate(age[age@data$AgeClass=="OldGrowth",])
  secondary <- raster::aggregate(age[age@data$AgeClass=="Secondary",])
  
  # 2017-2018
  areaOld <- raster::crop(dchm17to18,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  gapsOld <- raster::crop(gaps17to18,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  values_areaOld <- raster::getValues(areaOld)
  values_gapsOld <- raster::getValues(gapsOld)
  length(values_gapsOld[!is.na(values_gapsOld)])/length(values_areaOld[!is.na(values_areaOld)])
  
  areaSec <- raster::crop(dchm17to18,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  gapsSec <- raster::crop(gaps17to18,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  values_areaSec <- raster::getValues(areaSec)
  values_gapsSec <- raster::getValues(gapsSec)
  length(values_gapsSec[!is.na(values_gapsSec)])/length(values_areaSec[!is.na(values_areaSec)])
  
  # 2018-2019
  areaOld <- raster::crop(dchm18to19,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  gapsOld <- raster::crop(gaps18to19,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  values_areaOld <- raster::getValues(areaOld)
  values_gapsOld <- raster::getValues(gapsOld)
  length(values_gapsOld[!is.na(values_gapsOld)])/length(values_areaOld[!is.na(values_areaOld)])
  
  areaSec <- raster::crop(dchm18to19,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  gapsSec <- raster::crop(gaps18to19,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  values_areaSec <- raster::getValues(areaSec)
  values_gapsSec <- raster::getValues(gapsSec)
  length(values_gapsSec[!is.na(values_gapsSec)])/length(values_areaSec[!is.na(values_areaSec)])
  
  # 2019-2020
  areaOld <- raster::crop(dchm19to20,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  gapsOld <- raster::crop(gaps19to20,buffer[buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  values_areaOld <- raster::getValues(areaOld)
  values_gapsOld <- raster::getValues(gapsOld)
  (12/13)*length(values_gapsOld[!is.na(values_gapsOld)])/length(values_areaOld[!is.na(values_areaOld)])
  
  areaSec <- raster::crop(dchm19to20,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  gapsSec <- raster::crop(gaps19to20,buffer[!buffer$Mascaro_Co=="> 400" & !is.na(buffer$Mascaro_Co),])
  values_areaSec <- raster::getValues(areaSec)
  values_gapsSec <- raster::getValues(gapsSec)
  (12/13)*length(values_gapsSec[!is.na(values_gapsSec)])/length(values_areaSec[!is.na(values_areaSec)])
  
  
#### GRAVEYARD ####    
  # Get soil type grap frequency from lidar
  soilTypes$PropGapStanding <- NA
  dsm <- raster::raster("~/Desktop/Postdoc/BCI lidar 2009/BCI_DSM_2009.tif")
  dsm <- raster::crop(dsm, raster::extent(dCHM))
  dsm_mask <- dsm; dsm_mask[is.na(dCHM_binary)] <- NA
  dem <- raster::raster("~/Desktop/Postdoc/BCI lidar 2009/DEM_clipMCH.tif")
  dem <- raster::resample(dem, dsm)
  dem <- raster::crop(dem, raster::extent(dCHM))
  chm <- dsm-dem
  chm[is.na(dCHM_binary)] <- NA
  
  for(i in 1:dim(soilTypes)[1]){
    soili <- raster::crop(chm,soil[soil$SOIL==soilTypes$Soil[i],])
    gapsi <- ForestGapR::getForestGaps(soili, threshold = 2)
    
    valuesi <- raster::getValues(soili)
    valuesi <- valuesi[!is.na(valuesi)]
    
    gapsArea<- sum(ForestGapR::GapStats(gapsi,soili)$gap_area)
    
    soilTypes$PropGapStanding[i] <- gapsArea/length(valuesi)
  }
  
  plot(PropGap~PropGapStanding, data = soilTypes, pch=20,
       xlab="Proportion standing gaps in 2009",
       ylab = "Proportion new gaps 2018-2019")
  summary(lm(PropGap~PropGapStanding, data = soilTypes))
  
  
  