#### find height correction for 2015 data that best matched hi-res data ####
  
  # Read plot outline  
  plotShp <- rgdal::readOGR("Data_Ancillary/BCI50ha/BCI_50ha.shp")
  plotShp <- sp::spTransform(plotShp, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

  # Load uncorrected 2015 canopy height raster
  chm15 <- raster::raster("Data_HeightRasters/CHM_2015_QAQC_wBias.tif")
  
  # Load 2008, 2018 canopy height rasters
  chm09 <- raster::raster("Data_HeightRasters/CHM_2009_QAQC.tif")
  chm18 <- raster::raster("Data_HeightRasters/CHM_2018_QAQC.tif")
  
  # Crop canopy height rasters to 50 ha plot and pull values
  vals09 <- raster::values(raster::crop(chm09,plotShp))
  dens09 <- density(vals09[!is.na(vals09)])
  
  vals15 <- raster::values(raster::crop(chm15,plotShp))
  dens15 <- density(vals15[!is.na(vals15)])
  
  vals18 <- raster::values(raster::crop(chm18,plotShp))
  dens18 <- density(vals18[!is.na(vals18)])
  
  # Monthly data from close in time
  dsm15_50ha <- raster::raster("Data_HeightRasters/DSM_50haPlot_20150629_geo.tif")
  
    # Aggregate 50 ha plot data to get close to correct resolution taking highest value
    dsm15_50ha <- raster::aggregate(dsm15_50ha,14,fun=max)
    
    # Crop to plot
    dsm15_50ha <- raster::crop(dsm15_50ha, plotShp)
    
    # Read in BCI DEM to get height above ground
    dem <- raster::raster("Data_HeightRasters/LidarDEM_BCI.tif")
    dem <- raster::crop(dem, plotShp)
    
    # Resample so both match
    dsm15_50ha <- raster::resample(dsm15_50ha,raster::crop(chm09,plotShp))
    dem <- raster::resample(dem,raster::crop(chm09,plotShp))
    
    # Subtract to get CHM
    chm15_50ha <- dsm15_50ha-dem
    
    # Get density
    vals15_50ha <- raster::values(chm15_50ha)
    
    # Calculate offset (because hi-res data are not vertically aligned, assume that talles points are the same height)
    offset_50ha <- max(vals15_50ha,na.rm=T) -max(vals15,na.rm=T)
    vals15_50ha <- vals15_50ha - offset_50ha
    dens15_50ha <- density(vals15_50ha[!is.na(vals15_50ha)])
  
  # Do sensitivity analysis to find best parameters to use
    correctionParams <- data.frame(par1 = rep(seq(2,10,0.5),length(seq(1,10,0.5))),
                                   par2 = rep(seq(1,10,0.5),each=length(seq(2,10,0.5))),
                                   stat = NA)
    
    for(i in 1:nrow(correctionParams)){
      toChange <- which((vals15-vals09)>correctionParams$par1[i] & (vals18-vals15) < (-1*correctionParams$par2[i]))
      newVals <- vals09[toChange]
      vals15c <- vals15; vals15c[toChange] <- newVals
      
      correctionParams$stat[i] <- kruskal.test(list(vals15c,vals15_50ha))$statistic
    }
    
    #par1 <- correctionParams[which(correctionParams$stat==min(correctionParams$stat)),"par1"]
    #par2 <- correctionParams[which(correctionParams$stat==min(correctionParams$stat)),"par2"]
    

    # BEST: par1 = 5 and par2 = 1
   
    
  # PLOT RESULTS   
    par(mfrow=c(1,2),mar=c(4,5,2,4))
    
  # a. plot sensitivity analysis results  
    htThreshKruskall <- matrix(data = correctionParams$stat,
                           nrow = length(seq(2,10,0.5)),
                           ncol = length(seq(1,10,0.5)),
                           byrow = F,
                           dimnames = list(seq(2,10,0.5),seq(1,10,0.5)))
    htThreshKruskall <- htThreshKruskall-min(htThreshKruskall)
    
    plotMin <- min(htThreshKruskall)
    plotMax <- max(htThreshKruskall)
    plotBreaks <- c(0,1,2,4,8,16,24,seq(100,1200,100))
    
    library(plot.matrix)
    plot(htThreshKruskall, breaks=plotBreaks,
         col = rev(wesanderson::wes_palette("Zissou1", length(plotBreaks)-1, type = "continuous")),
         ylab = expression("Height increase '09-'15 (m)"),
         xlab = expression("Height decrease '15-'18 (m)"),
         main = expression(Delta~"rank sum"))
    mtext("a", side=3, outer=F, adj = -0.1, line=1.5, cex=1.5)
    
    # Distribution of canopy heights with correction
      toChange <- which((vals15-vals09)> par1 & (vals18-vals15)< (-1*par2))
      newVals <- vals09[toChange]
      vals15c <- vals15; vals15c[toChange] <- newVals
      dens15c <- density(vals15c[!is.na(vals15c)])
      
    plot(dens15,
         xlim=c(0,65),
         main = NA,
         col="red",lwd=2,
         xlab="Canopy height (m)")
    lines(dens15c, col="red", lwd=2, lty=2)
    lines(dens15_50ha, col="black", lwd=2, lty=1)
    # lines(dens09,col="black",lwd=2)
    # lines(dens18,col="blue",lwd=2)
    # legend(x=30,y=0.055,
    #        c("2009 (lidar)","2015 (photogram.)","2018 (photogram.)"),
    #        col=c("black","red","blue"),
    #        bty="n",
    #        lwd=2,
    #        cex=0.9)
    legend(x=30,y=0.055,
           c("Island-wide: original",
             "Island-wide: corrected",
             "50 ha plot: hi-resolution"),
           col=c("red","red","black"),
           lty=c(1,2,1),
           bty="n",
           lwd=2,
           cex=1)
    mtext("b", side=3, outer=F, adj = -0.1, line=1.5, cex=1.5)
  