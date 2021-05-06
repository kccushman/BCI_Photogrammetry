#### Read data ####
  
  # Read validation gap shapefiles from Raquel's 50 ha plot analysis
    valGaps <- rgdal::readOGR("Raquel50haData/gaps_2014_2019_all/gaps_2014_2019_all.shp")
    
    # Convert date to date class
    valGaps$dateClass <- as.Date(valGaps$date)
    
    # Merge with height drop data and only keep gaps with height drop >=5 m
    valHt <- read.csv("gaps1419_centroid_join_heightdrop.csv")
    names(valHt)[4] <- "htDrop"
    
    valGaps$htDrop <- NA
    for(i in 1:nrow(valGaps@data)){
      ht_i <- valHt[valHt$area_m2==valGaps@data$area_m2[i],"htDrop"]
      if(length(ht_i)==1){
        valGaps$htDrop[i] <- ht_i
      }
    }
    
    # Split into years that correspond to whole island flights
    date15 <- as.Date("2015-07-1")
    date18 <- as.Date("2018-06-20")

    valGaps15to18 <- valGaps[valGaps$dateClass>date15 & valGaps$dateClass<date18,]

  # Read canopy height and canopy height change rasters
    chm09 <- raster::raster("CHM_2009_QAQC.tif")
    chm15 <- raster::raster("CHM_2015_QAQC_tin.tif")
    chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
    chm20 <- raster::raster("CHM_2020_QAQC_tin.tif")
    
    dchm15to18 <- raster::raster("dCHM15to18_tin.tif")

  # Read new gap shapefiles
    gaps15to18 <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
    
  # Read plot outline  
    plotShp <- rgdal::readOGR("D:/BCI_Spatial/BCI50ha/BCI_50ha.shp")
    plotShp <- sp::spTransform(plotShp, sp::proj4string(gaps15to18))
    
#### Remove gaps that cross plot border ####
    
  # my data
      gaps15to18 <- raster::crop(gaps15to18,plotShp,
                                 snap = "in")
      
      gaps15to18$area2 <- NA
      for(i in 1:length(gaps15to18)){
        gaps15to18$area2[i] <- raster::area(gaps15to18[i,])
      }
      
     # gaps15to18 <- gaps15to18[gaps15to18$area2>=25,]
  
  # Raquel's data
      valGaps15to18 <- rgeos::gBuffer(valGaps15to18, width=0, byid=T)
      
      valGaps15to18 <- raster::crop(valGaps15to18,plotShp,
                                 snap = "in")
      
      valGaps15to18$area2 <- NA
      for(i in 1:length(valGaps15to18)){
        valGaps15to18$area2[i] <- raster::area(valGaps15to18[i,])
      }
      
     # valGaps15to18 <- valGaps15to18[valGaps15to18$area2>=25,]
      

#### Make initial plot ####
    
    # 2015 to 2018
    raster::plot(dchm15to18,
                 ext = raster::extent(valGaps))
    raster::plot(plotShp, add=T)
    raster::plot(valGaps15to18, col = NA, 
                 border = "black", lwd = 1, add=T)
    raster::plot(gaps15to18, col = NA, 
                 border = "red", lty=2, add=T)
    
#### Calculate precision and recall for various threshold values ####
    
  # Try alternate circularity metric
  gaps15to18$ratioCirc <- 4*pi*gaps15to18$area/(gaps15to18$perimeter^2)
    
    
  thresholdData <- data.frame(Ratio = seq(0.3,1.3,0.1),
                              Precision_n = NA,
                              Recall_n = NA,
                              Precision_area = NA,
                              Recall_area = NA)
    
 
       
  for(i in 1:nrow(thresholdData)){
    
    # Only keep data above threshold
    data_i <- gaps15to18[gaps15to18$ratio>=thresholdData$Ratio[i], ]
    
      # Calculate precision: how many of my gaps are also in Raquel's data?
      data_i$observed <- NA
      for(j in 1:length(data_i)){
        data_i$observed[j] <- rgeos::gIntersects(data_i[j,], valGaps15to18)
      }
      
      # Proportion of gaps recalled
  
      thresholdData[i,"Precision_n"] <- length(data_i[data_i$observed==T,])/length(data_i)
      
      # Proportion of gap area (based on T/F) recalled
      thresholdData[i,"Precision_area"]  <- sum(data_i$area[data_i$observed==T])/sum(data_i$area)
      
      
      # Calculate recall: how many of Raquel's gaps are also in my data?
        valGaps15to18$observed <- NA
        for(j in 1:length(valGaps15to18)){
          valGaps15to18$observed[j] <- rgeos::gIntersects(valGaps15to18[j,], data_i)
        }
        
        # Proportion of gaps recalled
        thresholdData[i,"Recall_n"] <- length(valGaps15to18[valGaps15to18$observed==T,])/length(valGaps15to18)
        
        # Proportion of gap area (based on T/F) recalled
        thresholdData[i,"Recall_area"] <- sum(valGaps15to18$area_m2[valGaps15to18$observed==T])/sum(valGaps15to18$area_m2)
        
      print(i)
    
  }  
  
  altThresholdData <- data.frame(Ratio = seq(0.05,0.2,0.01),
                              Precision_n = NA,
                              Recall_n = NA,
                              Precision_area = NA,
                              Recall_area = NA)
  
  
  
  for(i in 1:nrow(altThresholdData)){
    
    # Only keep data above threshold
    data_i <- gaps15to18[gaps15to18$ratioCirc>=altThresholdData$Ratio[i], ]
    
    # Calculate precision: how many of my gaps are also in Raquel's data?
    data_i$observed <- NA
    for(j in 1:length(data_i)){
      data_i$observed[j] <- rgeos::gIntersects(data_i[j,], valGaps15to18)
    }
    
    # Proportion of gaps recalled
    
    altThresholdData[i,"Precision_n"] <- length(data_i[data_i$observed==T,])/length(data_i)
    
    # Proportion of gap area (based on T/F) recalled
    altThresholdData[i,"Precision_area"]  <- sum(data_i$area[data_i$observed==T])/sum(data_i$area)
    
    
    # Calculate recall: how many of Raquel's gaps are also in my data?
    valGaps15to18$observed <- NA
    for(j in 1:length(valGaps15to18)){
      valGaps15to18$observed[j] <- rgeos::gIntersects(valGaps15to18[j,], data_i)
    }
    
    # Proportion of gaps recalled
    altThresholdData[i,"Recall_n"] <- length(valGaps15to18[valGaps15to18$observed==T,])/length(valGaps15to18)
    
    # Proportion of gap area (based on T/F) recalled
    altThresholdData[i,"Recall_area"] <- sum(valGaps15to18$area_m2[valGaps15to18$observed==T])/sum(valGaps15to18$area_m2)
    
    print(i)
    
  }  
      
      # Precision (black) how many of my gaps are in Raquel's data
      # Recall (red) how many of Raquel's gaps are in my data
      
    par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
      plot(Precision_n~Ratio, data=thresholdData,
           ylim=c(0,1),
           ylab="Precision/recall",
           xlab="Ratio value",
           type = "l")
      lines(Precision_area~Ratio, data=thresholdData,
            lty=2)
      lines(Recall_n~Ratio, data=thresholdData,
            col="red",
            lty=1)
      lines(Recall_area~Ratio, data=thresholdData,
            col="red",
            lty=2)
      
      par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
      plot(Precision_n~Ratio, data=altThresholdData ,
           ylim=c(0,1),
           ylab="Precision/recall",
           xlab="Ratio value",
           type = "l")
      lines(Precision_area~Ratio, data=altThresholdData ,
            lty=2)
      lines(Recall_n~Ratio, data=altThresholdData ,
            col="red",
            lty=1)
      lines(Recall_area~Ratio, data=altThresholdData ,
            col="red",
            lty=2)
      
    
    
    
    

      
  
  # Save shapefiles to look in ArcGIS
      rgdal::writeOGR(valGaps15to18,
                      dsn = "Monthly_50haGaps_2015-2018",
                      layer = "MonthlyGaps", 
                      driver = "ESRI Shapefile")
      
      rgdal::writeOGR(gaps15to18,
                      dsn = "Annual_50haGaps_2015-2018",
                      layer = "AnnualGaps", 
                      driver = "ESRI Shapefile")
      
#### Crop Raquel's images around my gaps that she doesn't see ####
      # Calculate precision: how many of my gaps are also in Raquel's data?
      gaps15to18$observed <- NA
      for(j in 1:length(gaps15to18)){
        gaps15to18$observed[j] <- rgeos::gIntersects(gaps15to18[j,], valGaps15to18)
      }
          
  missingGaps <- gaps15to18[gaps15to18$observed==F,]    
  missingGaps <- missingGaps[order(-missingGaps$area),]
  
  imageList <- list.files("C:/Users/cushmank/Desktop/RaquelOrthoAligned/", full.names = T, pattern = ".tif")
  
  
  # Save Raquel's images
  for(i in 1:10){
    # find extent of gap to extract from images
      extent_i <- raster::extent(missingGaps[i,])
    
    # Make a new folder to store results
      dir.create(paste0("C:/Users/cushmank/Desktop/RaquelOrthoAligned/MissingGap_",i))
      
    # Loop through all dates
      for(j in 1:length(imageList)){
        image_j <- raster::brick(imageList[j])
        date_j <- strsplit(strsplit(imageList[j], split="/")[[1]][6], split=".tif")[[1]]
        
        image_crop <- raster::crop(image_j,
                                   raster::extent(extent_i@xmin - 20,
                                                  extent_i@xmax + 20,
                                                  extent_i@ymin - 20,
                                                  extent_i@ymax + 20))
        
        raster::writeRaster(image_crop, file = paste0("C:/Users/cushmank/Desktop/RaquelOrthoAligned/MissingGap_",i,"/",date_j,"_",i,".tiff"))
      }
  }
  
  
  
  # Make a plot of my raster heights
  for(i in 1:20){
    
    # find extent of gap to extract from images
    extent_i <- raster::extent(missingGaps[i,])
    
    # Correct
    chm09i <- raster::crop(chm09,
                           raster::extent(extent_i@xmin - 10,
                                          extent_i@xmax + 10,
                                          extent_i@ymin - 10,
                                          extent_i@ymax + 10))
    chm15i <- raster::crop(chm15,
                           raster::extent(extent_i@xmin - 10,
                                          extent_i@xmax + 10,
                                          extent_i@ymin - 10,
                                          extent_i@ymax + 10))
    chm18i <- raster::crop(chm18,
                           raster::extent(extent_i@xmin - 10,
                                          extent_i@xmax + 10,
                                          extent_i@ymin - 10,
                                          extent_i@ymax + 10))
    
    # try some thresholds
    vals09i <- raster::values(chm09i)
    vals15i <- raster::values(chm15i)
    vals18i <- raster::values(chm18i)
    
    toChange <- which((vals15i-vals09i)>3 & (vals18i-vals15i)< -1)
    newVals <- vals09i[toChange]
    chm15c <- chm15i; raster::values(chm15c)[toChange] <- newVals
    
    
    colBrks <- c(seq(0,30,2),50)
    
    par(mfrow=c(2,3), mar=c(1,1,1,1), oma=c(0,0,1,6))
    raster::plot(chm09i,
                 breaks=colBrks,
                 col=rev(terrain.colors(length(colBrks))),
                 bty="n",box=F,xaxt="n",yaxt="n",legend=F,
                 main = "2009 canopy height (m)")
    raster::plot(missingGaps[i,], border="red",add=T)
    
    raster::plot(chm15i,
                 breaks=colBrks,
                 col=rev(terrain.colors(length(colBrks))),
                 bty="n",box=F,xaxt="n",yaxt="n",legend=F,
                 main = "2015 canopy height (m)")
    raster::plot(missingGaps[i,], border="red",add=T)
    
    raster::plot(chm18i,
                 breaks=colBrks,
                 col=rev(terrain.colors(length(colBrks))),
                 bty="n",box=F,xaxt="n",yaxt="n",
                 main = "2018 canopy height (m)")
    raster::plot(missingGaps[i,], border="red",add=T)
    
    raster::plot((chm15i-chm09i),
                 breaks=seq(-30,30,3),
                 col=rev(rainbow(21)),
                 bty="n",box=F,xaxt="n",yaxt="n", legend=F,
                 main = "2009-15 ht change (m)")
    raster::plot(missingGaps[i,], border="red",add=T)
    
    raster::plot((chm18i-chm15i),
                 breaks=seq(-30,30,3),
                 col=rev(rainbow(21)),
                 bty="n",box=F,xaxt="n",yaxt="n", legend=T,
                 main = "2015-18 ht change (m)")
    raster::plot(missingGaps[i,], border="red",add=T)
    
  }
  
  
  rgdal::writeOGR(missingGaps[c(3,8),],
                  dsn = "GapsToPlot",
                  layer = "GapsToPlot", 
                  driver = "ESRI Shapefile")
      

  
#### plot height value histograms for the 50ha plot ####
  vals09 <- raster::values(raster::crop(chm09,plotShp))
  dens09 <- density(vals09[!is.na(vals09)])
  
  vals15 <- raster::values(raster::crop(chm15,plotShp))
  dens15 <- density(vals15[!is.na(vals15)])
  
  vals18 <- raster::values(raster::crop(chm18,plotShp))
  dens18 <- density(vals18[!is.na(vals18)])
  
  par(mfrow=c(1,1), mar=c(4,4,1,1))
  plot(dens15,
       xlim=c(0,50),
       main = "Canopy height distribution",
       col="red",lwd=2,
       xlab="canopy height (m)")
  lines(dens09,col="black",lwd=2)
  lines(dens18,col="blue",lwd=2)
  legend(x=0,y=0.055,
         c("2009 (lidar)","2015 (photogram.)","2018 (photogram.)"),
         col=c("black","red","blue"),
         bty="n",
         lwd=2)
  
#### Sandbox: figure out good metric for correcting data ####
  chm09i <- raster::crop(chm09,
                         raster::extent(extent_i@xmin - 20,
                                        extent_i@xmax + 20,
                                        extent_i@ymin - 20,
                                        extent_i@ymax + 20))
  chm15i <- raster::crop(chm15,
                         raster::extent(extent_i@xmin - 20,
                                        extent_i@xmax + 20,
                                        extent_i@ymin - 20,
                                        extent_i@ymax + 20))
  chm18i <- raster::crop(chm18,
                         raster::extent(extent_i@xmin - 20,
                                        extent_i@xmax + 20,
                                        extent_i@ymin - 20,
                                        extent_i@ymax + 20))

  # try some thresholds
  vals09i <- raster::values(chm09i)
  vals15i <- raster::values(chm15i)
  vals18i <- raster::values(chm18i)
  
  # WORKS toChange <- which(vals09i<35 & (vals15i-vals09i)>5 & vals18i<35)
  toChange <- which((vals15i-vals09i)>3 & (vals18i-vals15i)< -1)
  newVals <- vals09i[toChange]
  chm15c <- chm15i; raster::values(chm15c)[toChange] <- newVals
                         
  # plot
  par(mfrow=c(2,3), mar=c(1,1,1,1), oma=c(0,0,1,6))
  raster::plot(raster::crop(chm09,
                            raster::extent(extent_i@xmin - 20,
                                           extent_i@xmax + 20,
                                           extent_i@ymin - 20,
                                           extent_i@ymax + 20)),
               breaks=colBrks,
               col=rev(terrain.colors(length(colBrks))),
               bty="n",box=F,xaxt="n",yaxt="n",legend=F)
  raster::plot(missingGaps[i,], border="red",add=T)
  
  raster::plot(raster::crop(chm15,
                            raster::extent(extent_i@xmin - 20,
                                           extent_i@xmax + 20,
                                           extent_i@ymin - 20,
                                           extent_i@ymax + 20)),
               breaks=colBrks,
               col=rev(terrain.colors(length(colBrks))),
               bty="n",box=F,xaxt="n",yaxt="n",legend=F)
  raster::plot(missingGaps[i,], border="red",add=T)
  
  raster::plot(raster::crop(chm18,
                            raster::extent(extent_i@xmin - 20,
                                           extent_i@xmax + 20,
                                           extent_i@ymin - 20,
                                           extent_i@ymax + 20)),
               breaks=colBrks,
               col=rev(terrain.colors(length(colBrks))),
               bty="n",box=F,xaxt="n",yaxt="n")
  raster::plot(missingGaps[i,], border="red",add=T)
  
  raster::plot(chm15c,
               breaks=colBrks,
               col=rev(terrain.colors(length(colBrks))),
               bty="n",box=F,xaxt="n",yaxt="n")
  raster::plot(missingGaps[i,], border="red",add=T)
  
  raster::plot(chm18i-chm15c,
               breaks=c(-50,-5,seq(-3,45,2)),
               col=c(viridis::viridis(30)),
               bty="n",box=F,xaxt="n",yaxt="n")
  raster::plot(missingGaps[i,], border="red",add=T)
  
  raster::plot(chm18i-chm15c,
               breaks=c(-50,-8,seq(-3,45,2)),
               col=c("red",viridis::viridis(30)),
               bty="n",box=F,xaxt="n",yaxt="n")
  raster::plot(missingGaps[i,], border="red",add=T)
  
  
#### Correct 50 ha plot and test against Raquel's data ####
  
  # make a new plot polygon with buffer
  plotBuff <- rgeos::gBuffer(plotShp, width = 20)
  
  chm09p <- raster::crop(chm09, plotBuff); chm09p <- raster::mask(chm09p, plotBuff)
  chm15p <- raster::crop(chm15, plotBuff); chm15p <- raster::mask(chm15p, plotBuff)
  chm18p <- raster::crop(chm18, plotBuff); chm18p <- raster::mask(chm18p, plotBuff)
  
  vals09p <- raster::values(chm09p)
  vals15p <- raster::values(chm15p)
  vals18p <- raster::values(chm18p)
  
  
  toChange <- which((vals15p-vals09p) > 3 & (vals18p-vals15p) < -1)
  newVals <- vals09p[toChange]
  chm15pc <- chm15p
  raster::values(chm15pc)[toChange] <- newVals
  
  dchm15to18c <- chm18p-chm15pc
  
  # raster::plot(dchm15to18c,breaks=c(-50,-5,seq(-3,45,2)),
  #              col=c("red",viridis::viridis(30)))
  # raster::plot(chm18p-chm15p,breaks=c(-50,-8,seq(-3,45,2)),
  #              col=c("red",viridis::viridis(30)))
  
  # Define gaps
  
    # Make new function
    getForestGaps <- function (chm_layer, threshold = 10, size = c(1, 10^4)) 
      {
        chm_layer[chm_layer > threshold] <- NA
        chm_layer[chm_layer <= threshold] <- 1
        gaps <- raster::clump(chm_layer, directions = 4, gap = FALSE)
        rcl <- raster::freq(gaps)
        rcl[, 2] <- rcl[, 2] * raster::res(chm_layer)[1]^2
        rcl <- cbind(rcl[, 1], rcl)
        z <- raster::reclassify(gaps, rcl = rcl, right = NA)
        z[is.na(gaps)] <- NA
        gaps[z > size[2]] <- NA
        gaps[z < size[1]] <- NA
        gaps <- raster::clump(gaps, directions = 4, gap = FALSE)
        names(gaps) <- "gaps"
        return(gaps)
      }
  
    # Do function
    plotGaps <- getForestGaps(dchm15to18c, threshold = -5, size=c(25,10^10))
    
    # Create a Spatial Polygon Data Frame object, where each polygon is a gap
    plotGaps <- ForestGapR::GapSPDF(plotGaps)
    
    # Calculate the area and perimeter from each gap object
    plotGaps@data$area <- NA
    plotGaps@data$perimeter <- NA
    for(i in 1:length(plotGaps)){
      plotGaps[plotGaps$gap_id==i,"area"] <- raster::area(plotGaps[plotGaps$gap_id==i,])
      
      perims_j <- c()
      for(j in 1:length(plotGaps[plotGaps$gap_id==i,]@polygons[[1]]@Polygons)){
        
        coordsj <- plotGaps[plotGaps$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
        
        lengths_k <- c()
        for(k in 2:dim(coordsj)[1]){
          lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
          
        }
        perims_j[j] <- sum(lengths_k)
        
        if(plotGaps[plotGaps$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
          perims_j[j] <- 0
        }
        
      }
      plotGaps[plotGaps$gap_id==i,"perimeter"] <- sum(perims_j)
      
    }
    
    # Calculate the ratio of area to perimeter
    plotGaps@data$ratio <- plotGaps@data$area/plotGaps@data$perimeter
    
    # Try alternate circularity metric
    plotGaps$ratioCirc <- 4*pi*plotGaps$area/(plotGaps$perimeter^2)
    

    # Remove gaps that cross plot border
      plotGaps <- raster::crop(plotGaps,plotShp,
                                 snap = "in")
      
      plotGaps$area2 <- NA
      for(i in 1:length(plotGaps)){
        plotGaps$area2[i] <- raster::area(plotGaps[i,])
      }
      
      plotGaps <- plotGaps[plotGaps$area2>=25,]
      names(plotGaps)[1:2] <- c("X1","X2")
    

    # Read Raquel's data that I have checked  
      valRef <- rgdal::readOGR("toCheckRaquel/toCheckRaquel.shp")
      
      valGaps15to18$vslChck <- 0
      for(i in 1:nrow(valGaps15to18@data)){
        if(!is.na(valGaps15to18@data$htDrop[i])){
          vslChck <- valRef[valRef$htDrop==valGaps15to18@data$htDrop[i],"vslChck"]
          if(length(vslChck)==1){
            valGaps15to18@data$vslChck[i] <- vslChck$vslChck
          }
        }
      }

      
    # Merge my data with visual check results
      
      valNew <- rgdal::readOGR("toCheckKC_2/toCheckKC_2.shp")
      plotGaps$vslChck <- NA
      for(i in 1:nrow(plotGaps@data)){
        vslChck <- valNew[valNew$X1==plotGaps@data$X1[i],"vslChck"]
        if(length(vslChck)==1){
          plotGaps@data$vslChck[i] <- vslChck$vslChck
        }
      }
      
    # Which of Raquel's gaps do I see?
      valGaps15to18$observed <- NA
      for(j in 1:length(valGaps15to18)){
        valGaps15to18$observed[j] <- rgeos::gIntersects(valGaps15to18[j,], plotGaps)
      }
      
      
    # Only keep Raquel's gaps with height drop for the same height and area values
    valGapsSz <- valGaps15to18[valGaps15to18$area_m2 >= 25,]
    valGapsHt <- valGaps15to18[valGaps15to18$htDrop <= -5 & !is.na(valGaps15to18$htDrop) ,]
    
    
    # Look for overlap--which of my gaps are in Raquel's data?
    plotGaps$observedAll <- NA
    for(j in 1:length(plotGaps)){
      plotGaps$observedAll[j] <- rgeos::gIntersects(plotGaps[j,], valGaps15to18)
    }
    
    # What about if you ignore Raquel's size threshold?
    plotGaps$observedArea <- NA
    for(j in 1:length(plotGaps)){
      plotGaps$observedArea[j] <- rgeos::gIntersects(plotGaps[j,], valGapsSz)
    }
    
    # What about if you ignore Raquel's ht drop threshold?
    plotGaps$observedHt <- NA
    for(j in 1:length(plotGaps)){
      plotGaps$observedHt[j] <- rgeos::gIntersects(plotGaps[j,], valGapsHt)
    }

#### Make shapefiles of gaps to check ####
  
    # Proportion of gaps recalled
    length(valGaps2[valGaps2$observed==T,])/length(valGaps2)
    
    # Proportion of gap area (based on T/F) recalled
    sum(valGaps2$area_m2[valGaps2$observed==T])/sum(valGaps2$area_m2)
    
    # Look at classification of checked gaps
    hist(valGaps2@data[valGaps2$observed==F,"vslChck"],breaks=seq(0.5,10.5,1))
    hist(toCheckKC@data[toCheckKC$observed==F,"vslChck"],
         breaks=seq(0.5,10.5,1))
    
    # what proportion of my gap area is false positive? just 3% now
    sum(toCheckKC@data[toCheckKC$observed==F & toCheckKC$vslChck %in% c(2,10),"area"])/sum(toCheckKC@data[,"area"])
    
    # what proportion of my gap area is in decaying trees? also about 3%
    sum(toCheckKC@data[toCheckKC$observed==F & toCheckKC$vslChck %in% c(7,9),"area"])/sum(toCheckKC@data[,"area"])
    
    # what proportion of my gap area is in other missed trees? 1%
    sum(toCheckKC@data[toCheckKC$observed==F & toCheckKC$vslChck %in% c(1),"area"])/sum(toCheckKC@data[,"area"])
    
    # what proportion is due to mismatch in space, size, or time? 4%
    sum(toCheckKC@data[toCheckKC$observed==F & toCheckKC$vslChck %in% c(5,6,8),"area"])/sum(toCheckKC@data[,"area"])
    

    # Missed gaps tend to be smaller and shorter than all gaps
    mean(valGaps2@data[valGaps2$observed==F & valGaps2$vslChck==1,"area_m2"])
    mean(valGaps2@data[valGaps2$observed==T,"area_m2"])
    mean(valGaps2@data[valGaps2$observed==F & valGaps2$vslChck==1,"htDrop"])
    mean(valGaps2@data[valGaps2$observed==T,"htDrop"])
    
    toCheckRaquel <- valGaps2
    
    # Add field to be filled in for type of error
    toCheckKC$visualCheck <- 0
    toCheckRaquel$visualCheck <- 0
    
    # Write shapefiles
    plotGaps@data$vslChck[is.na(plotGaps@data$vslChck)] <- 0
    rgdal::writeOGR(plotGaps,
                    dsn = "toCheckKC_2",
                    layer = "toCheckKC_2", 
                    driver = "ESRI Shapefile")

    
    rgdal::writeOGR(toCheckRaquel,
                    dsn = "toCheckRaquel",
                    layer = "toCheckRaquel", 
                    driver = "ESRI Shapefile")
    
#### Plot errors -- monthly data ####
    
    valGapsSzHt <- valGaps15to18[valGaps15to18$area2 >= 25 & valGaps15to18$htDrop <= -5 & !is.na(valGaps15to18$htDrop),]
    # make a jiggered date field to plot
    valGapsSzHt$datePlot <-  valGapsSzHt$dateClass + rnorm(mean = 0, sd = 7, n = nrow(valGapsSzHt))
    
    
    par(mar=c(4,4,1,1))
    
    # Observed with ht and area agreement
    plot(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observed==T,],
         pch=20,
         col=adjustcolor("grey",0.6),
         log="y",
         ylab="Gap area (m2)",
         xlab = "Date of gap observation (jiggered for plotting)",
         main = "Monthly gaps (> 25 m2 and with height drop > 5 m)")
    
    # BLUE: likely size discrepancy (potential overestimation in monthly data)
      points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==c(3,4),],
             pch=20,
             col=adjustcolor("lightblue",0.6))

    # RED: false positive in monthly data
      # Total false positive
      points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==2,],
             pch=20,
             col=adjustcolor("red",0.6))
      
    # PURPLE: mismatch in space
      # Spatial alignment mismatch
      points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==5,],
             pch=20,
             col=adjustcolor("purple",0.6))

    # GREEN: monthly data are correct but missed by me
      # Unclear why Raquel misses this gap
    points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==1,],
           pch=20,
           col=adjustcolor("green",0.6))

    legend(x=as.Date("2017-01-01"),
           y=650,
           bty="n",
           c("Correctly observed in annual data",
             "Likely size discrepancy (monthly overestimate?)",
             "False positive",
             "Not observed in annual data and correct",
             "Spatial misalignment with annual data"),
           col=adjustcolor(c("grey","lightblue","red","green","purple"),c(0.6)),
          pch=20,
          pt.cex=1.5,
          cex=0.8)
    
    sum(valGapsSzHt$area2[valGapsSzHt$observed==F])/sum(valGapsSzHt$area2)
    sum(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck==2])/sum(valGapsSzHt$area2)
    sum(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck %in% c(3,4)])/sum(valGapsSzHt$area2)
    sum(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck %in% c(1)])/sum(valGapsSzHt$area2)
    sum(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck %in% c(5)])/sum(valGapsSzHt$area2)
    
    length(valGapsSzHt$area2[valGapsSzHt$observed==F])/length(valGapsSzHt$area2)
    length(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck==2])/length(valGapsSzHt$area2)
    length(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck %in% c(3,4)])/length(valGapsSzHt$area2)
    length(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck %in% c(1)])/length(valGapsSzHt$area2)
    length(valGapsSzHt$area2[valGapsSzHt$observed==F  & valGapsSzHt$vslChck %in% c(5)])/length(valGapsSzHt$area2)
    
    
    par(mar=c(4,4,1,1))
    
    # Observed with ht and area agreement
    plot(htDrop~datePlot, data=valGapsSzHt[valGapsSzHt$observed==T,],
         pch=20,
         col=adjustcolor("grey",0.6),
         ylab="Gap height drop (m)",
         xlab = "Date of gap observation (jiggered for plotting)")
    
    # BLUE: likely size discrepancy (potential overestimation in monthly data)
    points(htDrop~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==c(3,4),],
           pch=20,
           col=adjustcolor("lightblue",0.6))
    
    # RED: false positive in monthly data
    # Total false positive
    points(htDrop~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==2,],
           pch=20,
           col=adjustcolor("red",0.6))
    
    # PURPLE: mismatch in space
    # Spatial alignment mismatch
    points(htDrop~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==5,],
           pch=20,
           col=adjustcolor("purple",0.6))
    
    # GREEN: monthly data are correct but missed by me
    # Unclear why Raquel misses this gap
    points(htDrop~datePlot, data=valGapsSzHt[valGapsSzHt$observed==F & valGapsSzHt$vslChck==1,],
           pch=20,
           col=adjustcolor("green",0.6))
    
    legend(x=as.Date("2017-01-01"),
           y=-23,
           bty="n",
           c("Correctly observed in annual data",
             "Likely size discrepancy (monthly overestimate?)",
             "False positive",
             "Not observed in annual data and correct",
             "Spatial misalignment with annual data"),
           col=adjustcolor(c("grey","lightblue","red","green","purple"),c(0.6)),
           pch=20,
           pt.cex=1.5,
           cex=0.8)

#### Plot errors -- island-wide  data ####
    par(mar=c(4,4,1,1))
    
    # Observed with ht and area agreement
    plot(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==T & plotGaps$observedHt==T,],
         ylim=c(0,0.7),
         pch=20,
         col=adjustcolor("grey",0.6),
         log="x",
         xlab="Gap area (m2)",
         ylab = "Gap circularity (4*pi*area/perim^2)",
         main = "Island-wide gaps")
    
    # BLUE: observed by Raquel's data when contraints are not imposed
    # Observed by Raquel's data are less than 25m2
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==F & plotGaps$observedHt==T,],
           pch=20,
           col=adjustcolor("lightblue",0.6))
    # Observed by Raquel's data have ht drop less than 5 m
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==T & plotGaps$observedHt==F & !(plotGaps$vslChck==9),],
           pch=20,
           col=adjustcolor("darkblue",0.6))
    # Observed by Raquel's data are less than 25m2 AND have ht drop less than 5m (NONE)
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==F & plotGaps$observedHt==F & !(plotGaps$vslChck==9),],
           pch=20,
           col=adjustcolor("blue",0.5))
    
    # RED: false positive in my data
    # Total false positive
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==2,],
           pch=20,
           col=adjustcolor("red",0.6))
    # Overestimation of 
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==10,],
           pch=20,
           col=adjustcolor("brown1",0.6))
    
    # PURPLE: mismatch in space/time
    # Spatial alignment mismatch
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==5,],
           pch=20,
           col=adjustcolor("purple",0.6))
    # Temporal mismatch (i.e. I see beginning of slowly dying tree and Raquel records the later fall)
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck %in% c(6,8),],
           pch=20,
           col=adjustcolor("violet",0.6))
    
    # GREEN: my data are correct
    # Unclear why Raquel misses this gap
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==1,],
           pch=20,
           col=adjustcolor("green",0.6))
    # Decaying tree: total miss
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==7,],
           pch=20,
           col=adjustcolor("limegreen",0.6))
    # Decaying tree: partial miss
    points(ratioCirc~area2, data=plotGaps@data[plotGaps$vslChck==9,],
           pch=20,
           col=adjustcolor("lightgreen",0.6))
    
    legend(x=100,
           y=0.72,
           bty="n",
           c("Correctly observed in monthly data",
             "Observed in monthly data but monthly area < 25 m2",
             "Observed in monthly data but monthly ht drop < 5 m",
             "False positive: total error",
             "False positive: size overestimation",
             "Not observed in monthly data but correct: discrete disturbance",
             "Not observed in monthly data but correct: decaying tree totally missed",
             "Not observed in monthly data but correct: decaying tree partially missed",
             "Spatial misalignment with montly data",
             "Temporal misalignment with montly data (decaying tree)"),
           col=adjustcolor(c("grey","lightblue","darkblue","red","brown1","green","limegreen","lightgreen","purple","violet"),c(0.6)),
           pch=20,
           pt.cex=1.5,
           cex=0.8)
    
    
    # Observed with ht and area agreement
    plot(ratio~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==T & plotGaps$observedHt==T,],
         pch=20,
         col=adjustcolor("grey",0.6),
         log="x",
         xlab="Gap area (m2)",
         ylab = "Old shape metric (area/perim)")
    
    # BLUE: observed by Raquel's data when contraints are not imposed
    # Observed by Raquel's data are less than 25m2
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==F & plotGaps$observedHt==T,],
           pch=20,
           col=adjustcolor("lightblue",0.6))
    # Observed by Raquel's data have ht drop less than 5 m
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==T & plotGaps$observedHt==F & !(plotGaps$vslChck==9),],
           pch=20,
           col=adjustcolor("darkblue",0.6))
    # Observed by Raquel's data are less than 25m2 AND have ht drop less than 5m (NONE)
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==T & plotGaps$observedArea==F & plotGaps$observedHt==F & !(plotGaps$vslChck==9),],
           pch=20,
           col=adjustcolor("blue",0.5))
    
    # RED: false positive in my data
    # Total false positive
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==2,],
           pch=20,
           col=adjustcolor("red",0.6))
    # Overestimation of 
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==10,],
           pch=20,
           col=adjustcolor("brown1",0.6))
    
    # PURPLE: mismatch in space/time
    # Spatial alignment mismatch
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==5,],
           pch=20,
           col=adjustcolor("purple",0.6))
    # Temporal mismatch (i.e. I see beginning of slowly dying tree and Raquel records the later fall)
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck %in% c(6,8),],
           pch=20,
           col=adjustcolor("violet",0.6))
    
    # GREEN: my data are correct
    # Unclear why Raquel misses this gap
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==1,],
           pch=20,
           col=adjustcolor("green",0.6))
    # Decaying tree: total miss
    points(ratio~area2, data=plotGaps@data[plotGaps$observedAll==F & plotGaps$vslChck==7,],
           pch=20,
           col=adjustcolor("limegreen",0.6))
    # Decaying tree: partial miss
    points(ratio~area2, data=plotGaps@data[plotGaps$vslChck==9,],
           pch=20,
           col=adjustcolor("lightgreen",0.6))
    
    legend(x=24,
           y=4,
           bty="n",
           c("Correctly observed in monthly data",
             "Observed in monthly data but monthly area < 25 m2",
             "Observed in monthly data but monthly ht drop < 5 m",
             "False positive: total error",
             "False positive: size overestimation",
             "Not observed in monthly data but correct: discrete disturbance",
             "Not observed in monthly data but correct: decaying tree totally missed",
             "Not observed in monthly data but correct: decaying tree partially missed",
             "Spatial misalignment with monthly data",
             "Temporal misalignment with monthly data (decaying tree)"),
           col=adjustcolor(c("grey","lightblue","darkblue","red","brown1","green","limegreen","lightgreen","purple","violet"),c(0.6)),
           pch=20,
           pt.cex=1.5,
           cex=0.8)
    
  # By area  
    # Total area not observed
    sum(plotGaps$area2[plotGaps$observedAll==F])/sum(plotGaps$area2)
    
    sum(plotGaps$area2[plotGaps$observedAll==T & plotGaps$observedArea==T & plotGaps$observedHt==T])/sum(plotGaps$area2)
    sum(plotGaps$area2[plotGaps$observedAll==T])/sum(plotGaps$area2)
    
    # False positive area
    sum(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(2,10)])/sum(plotGaps$area2)
    # Area missing from monthly data
    sum(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(1)])/sum(plotGaps$area2)
    sum(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(7,9)])/sum(plotGaps$area2)
    # Spatial or temporal mismatch
    sum(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(5,6,8)])/sum(plotGaps$area2)
    
  # By number of gaps  
    # Total area not observed
    length(plotGaps$area2[plotGaps$observedAll==F])/length(plotGaps$area2)
    
    length(plotGaps$area2[plotGaps$observedAll==T & plotGaps$observedArea==T & plotGaps$observedHt==T])/length(plotGaps$area2)
    length(plotGaps$area2[plotGaps$observedAll==T])/length(plotGaps$area2)
    
    # False positive area
    length(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(2,10)])/length(plotGaps$area2)
    # Area missing from monthly data
    length(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(1,7,9)])/length(plotGaps$area2)
    # Spatial or temporal mismatch
    length(plotGaps$area2[plotGaps$observedAll==F & plotGaps$vslChck %in% c(5,6,8)])/length(plotGaps$area2)
    