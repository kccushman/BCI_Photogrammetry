#### Load data and results ####
library(INLA)
load("INLA/INLA_prelim_40m_tin.RData")
load("INLA/INLA_fullModelResult.RData")
load("INLA/INLA_fullModelResult_separate.RData")

  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)


# smoothScales <- c(1,2,3,4,6,8,12,16,24,32,48,64)
# curvScale <- 8
# slopeScale <- 24


#### Calculate predicted values and residuals ####
  
## Fixed and random effects  
  
  # Get fitted values and reorder to match original data frame
  bci.gapsAll$pred <- model_full$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  bci.gapsAll$resd <- bci.gapsAll$gapPropCens - bci.gapsAll$pred
  
## Fixed effects only  
  
    # Intercept
    bci.gapsAll$fix_int <- model_full$summary.fixed$mean[model_full$names.fixed=="(Intercept)"]
    
    # Curvature
    bci.gapsAll$fix_C <- bci.gapsAll$Sc_curvMean_8*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_curvMean_8"])
    
    
    # Slope
    bci.gapsAll$fix_S <- bci.gapsAll$Sc_slopeMean_24*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_slopeMean_24"])
    bci.gapsAll$fix_S2 <- bci.gapsAll$Sc_slopeMean_24_sq*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_slopeMean_24_sq"])
    
    # Height above drainage
    bci.gapsAll$fix_H <- bci.gapsAll$Sc_drainMean*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_drainMean"])
    bci.gapsAll$fix_H2 <- bci.gapsAll$Sc_drainMean_sq*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_drainMean_sq"])
    
    # # Soil form
    # bci.gapsAll$fix_soilForm <- 0
    # 
    # bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="MottledHeavyClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormMottledHeavyClay"]
    # bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="PaleSwellingClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormPaleSwellingClay"]
    # bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="RedLightClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormRedLightClay"]
    # 
    
    # Soil parent material
    bci.gapsAll$fix_soilParent <- 0
    
    bci.gapsAll[!is.na(bci.gapsAll$soilParent) & bci.gapsAll$soilParent=="Andesite","fix_soilParent"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilParentAndesite"]
    bci.gapsAll[!is.na(bci.gapsAll$soilParent) & bci.gapsAll$soilParent=="CaimitoMarineSedimentary","fix_soilParent"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilParentCaimitoMarineSedimentary"]
    bci.gapsAll[!is.na(bci.gapsAll$soilParent) & bci.gapsAll$soilParent=="CaimitoVolcanic","fix_soilParent"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilParentCaimitoVolcanic"]
    
    # Year 
    bci.gapsAll$fix_year <- 0
    bci.gapsAll[bci.gapsAll$Year=="2018","fix_year"] <- model_full$summary.fixed$mean[model_full$names.fixed=="Year2018"]
    
    # Age 
    bci.gapsAll$fix_age <- 0
    bci.gapsAll[bci.gapsAll$age=="Secondary" & !(is.na(bci.gapsAll$age)),"fix_age"] <- model_full$summary.fixed$mean[model_full$names.fixed=="ageSecondary"]
    
    
    # All
    bci.gapsAll$fix_sum <- bci.gapsAll$fix_int + bci.gapsAll$fix_C + bci.gapsAll$fix_S + bci.gapsAll$fix_S2 + bci.gapsAll$fix_H + bci.gapsAll$fix_H2 +bci.gapsAll$fix_soilParent + bci.gapsAll$fix_age + bci.gapsAll$fix_year
    
    # Predicted value with just fixed effects
    bci.gapsAll$fix_pred <- exp(bci.gapsAll$fix_sum)/(1 + exp(bci.gapsAll$fix_sum))
    
    bci.gapsAll$fix_resd <- bci.gapsAll$gapPropCens - bci.gapsAll$fix_pred
    
  
  # load buffer for BCI
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs") 
  
  fixRaster18 <- raster::raster(x = matrix(data = bci.gapsAll[bci.gapsAll$Year=="2018","fix_pred"],
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps18)@xmin,
                                xmx = raster::extent(bci.gaps18)@xmax,
                                ymn = raster::extent(bci.gaps18)@ymin,
                                ymx = raster::extent(bci.gaps18)@ymax)
  fixRaster18 <- raster::mask(fixRaster18, buffer)
  
  
  fixRaster20 <- raster::raster(x = matrix(data = bci.gapsAll[bci.gapsAll$Year=="2020","fix_pred"],
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps20)@xmin,
                                xmx = raster::extent(bci.gaps20)@xmax,
                                ymn = raster::extent(bci.gaps20)@ymin,
                                ymx = raster::extent(bci.gaps20)@ymax)
  fixRaster20 <- raster::mask(fixRaster20, buffer)
  
#### Look at model residuals ####
  
pdf(file = "residualsPlots.pdf", onefile=T, height=9, width=6.5)
  
# CURVATURE 
  
  par(mfrow=c(4,1), las=1, mar=c(3,3,1,1), oma=c(2,2,1,1))
  
  # Both years, fixed and random effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_8"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    plot(resd~curvMean_8, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, fixed + random effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_8 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_8 <= curvRMSE$max[i],"resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
         pch=19,
         col = "orange")
    
  # Both years, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_8"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    plot(fix_resd~curvMean_8, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, only fixed effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_8 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_8 <= curvRMSE$max[i],"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    
  # 2015-2018, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"curvMean_8"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    plot(fix_resd~curvMean_8, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2015-2018, only fixed effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_8 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_8 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    
  # 2018-2020, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"curvMean_8"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    plot(fix_resd~curvMean_8, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2018-2020, only fixed effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_8 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_8 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    
    mtext("Curvature (LaPlacian convexity)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
# SLOPE 
    
    # Both years, fixed and random effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_24"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    plot(resd~slopeMean_24, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, fixed + random effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_24 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_24 <= slopeRMSE$max[i],"resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    
    # Both years, only fixed effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_24"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    plot(fix_resd~slopeMean_24, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, only fixed effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_24 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_24 <= slopeRMSE$max[i],"fix_resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    
    # 2015-2018, only fixed effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"slopeMean_24"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    plot(fix_resd~slopeMean_24, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2015-2018, only fixed effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_24 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_24 <= slopeRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    
    # 2018-2020, only fixed effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"slopeMean_24"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    plot(fix_resd~slopeMean_24, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2018-2020, only fixed effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_24 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_24 <= slopeRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    
    mtext("Slope (degrees)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
# HEIGHT ABOVE DRAINAGE 
    
    # Both years, fixed and random effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    plot(resd~drainMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, fixed + random effects")
    abline(h=0,col="red")
    lines(x = drainSpline$x, y= drainSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                           max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i],"resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    
    # Both years, only fixed effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    plot(fix_resd~drainMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, only fixed effects")
    abline(h=0,col="red")
    lines(x = drainSpline$x, y= drainSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                           max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i],"fix_resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    
    # 2015-2018, only fixed effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    plot(fix_resd~drainMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2015-2018, only fixed effects")
    abline(h=0,col="red")
    lines(x = drainSpline$x, y= drainSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                           max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    
    # 2018-2020, only fixed effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    plot(fix_resd~drainMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2018-2020, only fixed effects")
    abline(h=0,col="red")
    lines(x = drainSpline$x, y= drainSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                           max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    
    mtext("HAND (m)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
# ASPECT
    
    
    # Both years, fixed and random effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    plot(resd~aspectMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, fixed + random effects")
    abline(h=0,col="red")
    lines(x = aspectSpline$x, y= aspectSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    
    # Both years, only fixed effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    plot(fix_resd~aspectMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, only fixed effects")
    abline(h=0,col="red")
    lines(x = aspectSpline$x, y= aspectSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"fix_resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    
    # 2015-2018, only fixed effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    plot(fix_resd~aspectMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2015-2018, only fixed effects")
    abline(h=0,col="red")
    lines(x = aspectSpline$x, y= aspectSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    
    # 2018-2020, only fixed effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    plot(fix_resd~aspectMean, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2018-2020, only fixed effects")
    abline(h=0,col="red")
    lines(x = aspectSpline$x, y= aspectSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")  
    
    mtext("Aspect (degrees)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
 
dev.off() 

#### Make plots of observed, average predicted, and SD predicted values for each year ####

  allMeans <-  model_full$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
#  allMeans[is.na(bci.gapsAll$gapPropCens)] <- NA
  
  allSds <-  model_full$summary.fitted.values$sd[order(bci.gapsAll_Order$Order)]
#  allSds[is.na(bci.gapsAll$gapPropCens)] <- NA
  
  bci.gaps18$predictedMean <- allMeans[1:(nCellX*nCellY)]
  bci.gaps20$predictedMean <- allMeans[(1 + nCellX*nCellY):(2*nCellX*nCellY)]

  bci.gaps18$predictedSD <- allSds[1:(nCellX*nCellY)]
  bci.gaps20$predictedSD <- allSds[(1 + nCellX*nCellY):(2*nCellX*nCellY)]


  # convert to a raster object
  
  # 2018
  predMeanRaster18 <- raster::raster(x = matrix(data = bci.gaps18$predictedMean,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps18)@xmin,
                                     xmx = raster::extent(bci.gaps18)@xmax,
                                     ymn = raster::extent(bci.gaps18)@ymin,
                                     ymx = raster::extent(bci.gaps18)@ymax)
  predMeanRaster18 <- raster::mask(predMeanRaster18, buffer)
  predMeanRaster18[predMeanRaster18<0] <- NA
  
  
  obsRaster18 <- raster::raster(x = matrix(data = bci.gaps18$gapPropCens,
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps18)@xmin,
                                xmx = raster::extent(bci.gaps18)@xmax,
                                ymn = raster::extent(bci.gaps18)@ymin,
                                ymx = raster::extent(bci.gaps18)@ymax)
  obsRaster18 <- raster::mask(obsRaster18, buffer)
  
  
  
  predSdRaster18 <- raster::raster(x = matrix(data = bci.gaps18$predictedSD,
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps18)@xmin,
                                   xmx = raster::extent(bci.gaps18)@xmax,
                                   ymn = raster::extent(bci.gaps18)@ymin,
                                   ymx = raster::extent(bci.gaps18)@ymax)
  predSdRaster18[is.na(predMeanRaster18)] <- NA
  
  
  # 2020
  predMeanRaster20 <- raster::raster(x = matrix(data = bci.gaps20$predictedMean,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps20)@xmin,
                                     xmx = raster::extent(bci.gaps20)@xmax,
                                     ymn = raster::extent(bci.gaps20)@ymin,
                                     ymx = raster::extent(bci.gaps20)@ymax)
  predMeanRaster20 <- raster::mask(predMeanRaster20, buffer)
  predMeanRaster20[predMeanRaster20<0] <- NA
  
  
  obsRaster20 <- raster::raster(x = matrix(data = bci.gaps20$gapPropCens,
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps20)@xmin,
                                xmx = raster::extent(bci.gaps20)@xmax,
                                ymn = raster::extent(bci.gaps20)@ymin,
                                ymx = raster::extent(bci.gaps20)@ymax)
  obsRaster20 <- raster::mask(obsRaster20, buffer)
  
  
  
  predSdRaster20 <- raster::raster(x = matrix(data = bci.gaps20$predictedSD,
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps20)@xmin,
                                   xmx = raster::extent(bci.gaps20)@xmax,
                                   ymn = raster::extent(bci.gaps20)@ymin,
                                   ymx = raster::extent(bci.gaps20)@ymax)
  predSdRaster20[is.na(predMeanRaster20)] <- NA
  
  # Plot rasters
  colBreaks <- seq(0,0.50,0.05)

  
  # 2015-2108
  jpeg(height=3000,width=2000,file = "Figure S10. Raster patterns.jpg")
  par(mfrow=c(4,3))
  raster::plot(obsRaster18,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Observed")
  raster::plot(buffer, add=T)
  
  raster::plot(predMeanRaster18,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Predicted (mean)")
  raster::plot(buffer, add=T)
  
  raster::plot(predSdRaster18,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(seq(0,0.1,0.01))),
               breaks = seq(0,0.1,0.01),
               main = "Predicted (SD)")
  raster::plot(buffer, add=T)
  mtext("2015 - 2018", outer=F, side = 3, line = -1.5)
  
  # 2018-2020
  raster::plot(obsRaster20,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Observed")
  raster::plot(buffer, add=T)
  
  raster::plot(predMeanRaster20,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Predicted (mean)")
  raster::plot(buffer, add=T)
  
  raster::plot(predSdRaster20,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(seq(0,0.1,0.01))),
               breaks = seq(0,0.1,0.01),
               main = "Predicted (SD)")
  raster::plot(buffer, add=T)
  mtext("2018 - 2020", outer=F, side = 3, line = -1.5)
  
dev.off()

#### Figure 5: Make plots of average spatial pattern across all years ####

  avgStackFix <- raster::stack(fixRaster18,fixRaster20)
  
  avgPredictedRaster <- raster::calc(avgStackFix, mean, na.rm=T)
  avgPredictedRaster[avgPredictedRaster==0] <- NA
  
  avgStackObs <- raster::stack(obsRaster18,obsRaster20)
  avgObservedRaster <- raster::calc(avgStackObs, mean, na.rm=T)
  avgObservedRaster[avgObservedRaster==0] <- NA


  # avgBreaks <- c(seq(0,0.12,0.005))
  # raster::plot(sdPredictedRaster,
  #              col=viridis::viridis(50),
  #              main = "Average spatial pattern -- SD predicted")
  
  par(mfrow=c(1,1), mar=c(2,2,1,2), oma=c(0,0,0,0))
  raster::plot(avgPredictedRaster,
               col=viridis::viridis(50),
               bty="n", box=F, xaxt="n", yaxt="n",
               main = "All years: fixed effects only")  
  raster::plot(predMeanRaster18,
               col=viridis::viridis(50),
               bty="n", box=F, xaxt="n", yaxt="n",
               main = "Predicted (fixed + random) 2015-2018")  
  raster::plot(predMeanRaster20,
               col=viridis::viridis(50),
               bty="n", box=F, xaxt="n", yaxt="n",
               main = "Predicted (fixed + random) 2018-2020")  
  
#### Look at  model R2 ####

  par(las=1, mfrow=c(1,1))
  plot(gapPropCens~pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       col = adjustcolor("black", 0.1),
       ylab = "Observed disturbance proportion",
       xlab = "Predicted disturbance proportion",
       pch=19)
  abline(a=0,b=1,col="red")
  
  # Both years
    # All observations, fixed and random effects
    summary(lm(gapPropCens~pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),]))
    # Just fixed effects
    summary(lm(gapPropCens~fix_pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),]))
  # 2015-2018
    # All observations, fixed and random effects
    summary(lm(gapPropCens~pred, data = bci.gapsAll[bci.gapsAll$Year=="2018" & !is.na(bci.gapsAll$gapPropCens),]))
    # Just fixed effects
    summary(lm(gapPropCens~fix_pred, data = bci.gapsAll[bci.gapsAll$Year=="2018" & !is.na(bci.gapsAll$gapPropCens),]))
  # 2018-2020
    # All observations, fixed and random effects
    summary(lm(gapPropCens~pred, data = bci.gapsAll[bci.gapsAll$Year=="2020" & !is.na(bci.gapsAll$gapPropCens),]))
    # Just fixed effects
    summary(lm(gapPropCens~fix_pred, data = bci.gapsAll[bci.gapsAll$Year=="2020" & !is.na(bci.gapsAll$gapPropCens),]))
  
#### Figure 7: Low canopy area and average spatial pattern ####
  
  blockData <- read.csv("bootstrapBlocks.csv")
  lo09 <- raster::raster("binaryLoCanopy.tif")

  blockData$propFix <- NA
  blockData$propObs <- NA
  blockData$propLow <- NA
  
  for(i in 1:nrow(blockData)){
    gapFix <- raster::crop(avgPredictedRaster, raster::extent(blockData$xmin[i],
                                                       blockData$xmax[i],
                                                       blockData$ymin[i],
                                                       blockData$ymax[i]))
    blockData$propFix[i] <- mean(gapFix@data@values, na.rm=T)
    
    gapObs <- raster::crop(avgObservedRaster, raster::extent(blockData$xmin[i],
                                                              blockData$xmax[i],
                                                              blockData$ymin[i],
                                                              blockData$ymax[i]))
    blockData$propObs[i] <- mean(gapObs@data@values, na.rm=T)
    
    loSub <- raster::crop(lo09, raster::extent(blockData$xmin[i],
                                                       blockData$xmax[i],
                                                       blockData$ymin[i],
                                                       blockData$ymax[i]))
    blockData$propLow[i] <- length(loSub@data@values[!is.na(loSub@data@values) & loSub@data@values==0])/length(loSub@data@values[!is.na(loSub@data@values) & loSub@data@values==1])
    print(i)
  }
  
  par(mfrow=c(1,2),mar=c(4,4,1,1))
  plot(propLow~propObs, data = blockData[blockData$area>100000,],
       log = "xy",
       xlim=c(0.01,0.055),
       ylim=c(0.02,0.25),
       pch = 19, 
       xlab = "Average observed disturbance frequency",
       ylab = "Proportion of area with canopy < 10 m")
  text("a", x = 0.01, y = 0.25)
  plot(propLow~propFix, data = blockData[blockData$area>100000,],
       pch = 19, 
       xlim=c(0.005,0.018),
       ylim=c(0.02,0.25),
       log = "xy",
       xlab = "Average predicted disturbance frequency",
       ylab = NA)
  text("b", x = 0.005, y = 0.25)
  
  summary(lm(log(propLow)~log(propObs), data = blockData[blockData$area>100000,]))
  summary(lm(log(propLow)~log(propFix), data = blockData[blockData$area>100000,]))
  
  
#### Figure 6: Fixed effects sizes ####  
  fixedResults <- model_full$summary.fixed
  
  # separated by years
  fixedResults1 <- model_full1$summary.fixed
  fixedResults2 <- model_full2$summary.fixed
  
  
  
  fixedNames <- c("Curvature (linear)", "Slope (linear)", "Slope (quadratic)",                    
                  "Height above drainage (linear)","Height above drainage (quadratic)",
                  "Soil parent: Andesite", "Soil parent: Caimito marine", "Soil parent: Caimito volcanic",
                  "Forest age: Secondary",
                  "Year: 2015-2018")

  par(mfrow=c(1,1), mar=c(4,1,1,1))
  plot(x =fixedResults[2:nrow(fixedResults),"mean"],
       y = nrow(fixedResults):2,
       xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(0,0.6),
       ylim=c(2,11),
       pch = 19, 
       cex = 1,
       xlab = "Fixed effect",
       ylab = NA, yaxt = "n")
  arrows(x0 = fixedResults[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults[2:nrow(fixedResults),"0.975quant"],
         y0 = nrow(fixedResults):2,
         y1 = nrow(fixedResults):2,
         angle = 90, code=3,
         length = 0.05)
  
  points(x = fixedResults1[2:nrow(fixedResults),"mean"],
       y = nrow(fixedResults):2 - 0.2,
       pch = 19, 
       cex = 1,
       col = "blue")
  arrows(x0 = fixedResults1[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults1[2:nrow(fixedResults),"0.975quant"],
         y0 = nrow(fixedResults):2- 0.2,
         y1 = nrow(fixedResults):2- 0.2,
         angle = 90, code=3,
         length = 0.05,
         col="blue")
  
  points(x = fixedResults2[2:nrow(fixedResults),"mean"],
         y = nrow(fixedResults):2 - 0.4,
         pch = 19, 
         cex = 1,
         col = "lightblue")
  arrows(x0 = fixedResults2[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults2[2:nrow(fixedResults),"0.975quant"],
         y0 = nrow(fixedResults):2- 0.4,
         y1 = nrow(fixedResults):2- 0.4,
         angle = 90, code=3,
         length = 0.05,
         col="lightblue")
  
  
  abline(v=0, lty=2)
  text(fixedNames,
        x = fixedResults[2:nrow(fixedResults),"0.975quant"]+0.005,
       y = nrow(fixedResults):2,
       pos = 4)
  
  
  
#### Figure 4: Raster plots of landscape predictors ####
  
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
  # Read plot outline  
  plotShp <- rgdal::readOGR("D:/BCI_Spatial/BCI50ha/BCI_50ha.shp")
  plotShp <- sp::spTransform(plotShp, sp::proj4string(buffer))
  
  curvRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_8.tif")
  slopeRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_16.tif")
  drainRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/distAboveStream_1000.tif")
  
  curvRaster <- raster::mask(curvRaster, buffer)
  slopeRaster <- raster::mask(slopeRaster, buffer)
  drainRaster <- raster::mask(drainRaster, buffer)
  

  
  par(mar=c(1,1,1,1))
  raster::plot(curvRaster, col = viridis::cividis(128),
               bty="n", box=F, yaxt="n", xaxt="n",
               main = "Curvature")
  raster::plot(plotShp,add=T, lwd=2, border="white")
  
  raster::plot(slopeRaster, col = viridis::plasma(128),
               bty="n", box=F, yaxt="n", xaxt="n",
               main = "Slope")
  raster::plot(plotShp,add=T, lwd=2, border="white")
  
  raster::plot(drainRaster, col = viridis::viridis(128),
               bty="n", box=F, yaxt="n", xaxt="n",
               main = "Height above drainage")
  raster::plot(plotShp,add=T, lwd=2, border="white")
  
  soil$formCol = NA
    soil[soil$SoilForm=="BrownFineLoam","formCol"] = wesanderson::wes_palette("Rushmore1",5)[1]
    soil[soil$SoilForm=="PaleSwellingClay","formCol"] = wesanderson::wes_palette("Rushmore1",5)[3]
    soil[soil$SoilForm=="MottledHeavyClay","formCol"] = wesanderson::wes_palette("Rushmore1",5)[4]
    soil[soil$SoilForm=="RedLightClay","formCol"] = wesanderson::wes_palette("Rushmore1",5)[5]
    
  soil$parentCol = NA
    soil[soil$SoilParent=="CaimitoVolcanic","parentCol"] = wesanderson::wes_palette("Chevalier1",4)[1]
    soil[soil$SoilParent=="Andesite","parentCol"] = wesanderson::wes_palette("Chevalier1",4)[2]
    soil[soil$SoilParent=="CaimitoMarineSedimentary","parentCol"] = wesanderson::wes_palette("Chevalier1",4)[3]
    soil[soil$SoilParent=="Bohio","parentCol"] = wesanderson::wes_palette("Chevalier1",4)[4]
    
    soil <- sp::spTransform(soil, sp::proj4string(buffer))
    
  raster::plot(soil,
               main = "Soil form",
               col = soil$formCol,
               border="NA")
  raster::plot(plotShp,add=T, lwd=2, border="white")
  
  raster::plot(soil,
               main = "Soil parent material",
               col = soil$parentCol,
               border="NA")
  raster::plot(plotShp,add=T, lwd=2, border="white")
  
  age$ageCol = NA
    age[age$AgeClass=="OldGrowth","ageCol"] <- wesanderson::wes_palette("Moonrise2",4)[1]
    age[age$AgeClass=="Secondary","ageCol"] <- wesanderson::wes_palette("Moonrise2",4)[2]
    age[age$AgeClass=="Other","ageCol"] <- wesanderson::wes_palette("Moonrise2",4)[3]
  
  age <- raster::crop(age,buffer)
  raster::plot(age ,
               main = "Forest age",
               col = age$ageCol,
               border="NA")
  raster::plot(plotShp,add=T, lwd=2, border="white")
  
  
  
  
#### Fig??. Proportion of short canopy in each cell ####
  load("INLA/INLA_prelim_40m_tin.RData")
  
  
  chm15 <- raster::raster("CHM_2015_QAQC_tin.tif")
  chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
  chm20 <- raster::raster("CHM_2020_QAQC_tin.tif")
  
  
  propShort <- function(gapPoly, gapLayer, baseLayer, cellSz, nX, nY){
    
    shortVals <- data.frame(shortProp = rep(NA, length(gapPoly)))
    
    baseCrop <- raster::crop(baseLayer, raster::extent(gapPoly))
    baseCropNew <- raster::values(baseCrop)
    baseCropNew[!is.na(baseCropNew)] <- 1
    raster::values(baseCrop) <- baseCropNew
    
    shortCrop <- raster::crop(gapLayer, raster::extent(gapPoly))
    shortCropNew <- raster::values(shortCrop)
    shortCropNew[!is.na(shortCropNew) & (shortCropNew>=10 | shortCropNew<5)] <- 0
    shortCropNew[!is.na(shortCropNew) & (shortCropNew<10 & shortCropNew>=5)] <- 1

    raster::values(shortCrop) <- shortCropNew
    
    resampleBase <- raster::aggregate(baseCrop, cellSz, fun = sum)
    resampleShort <- raster::aggregate(shortCrop, cellSz, fun = sum)
    
    valuesBase <- raster::values(resampleBase)
    valuesShort <- raster::values(resampleShort)
    valuesShort[is.na(valuesShort)] <- 0
    
    # reorder to proper order for INLA object
    shortVals$shortProp <- valuesShort[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]/valuesBase[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]
    
    return(shortVals)
  }
  
  
  
  # 2015-2018
  shortArea <- propShort(gapPoly = bci.gaps18,
                     gapLayer = chm15,
                     baseLayer = d15to18tall,
                     cellSz = cellSize,
                     nX = nCellX,
                     nY = nCellY)
  bci.gaps18$shortProp <- shortArea$shortProp
  
  plot(gapProp~shortProp, data = bci.gaps18[!is.na(bci.gaps18$areaObs) & bci.gaps18$areaObs>(0.5*cellSize*cellSize/10000),],
       pch=19,
       col=adjustcolor("black",alpha.f = 0.1))
  summary(lm(gapProp~shortProp, data = bci.gaps18[!is.na(bci.gaps18$areaObs) & bci.gaps18$areaObs>(0.5*cellSize*cellSize/10000),]))
  
  
  # 2018-2020
  shortArea <- propShort(gapPoly = bci.gaps20,
                         gapLayer = chm18,
                         baseLayer = d18to20tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
  bci.gaps20$shortProp <- shortArea$shortProp
  
  plot(gapProp~shortProp, data = bci.gaps20[!is.na(bci.gaps20$areaObs) & bci.gaps20$areaObs>(0.5*cellSize*cellSize/10000),],
       pch=19,
       col=adjustcolor("black",alpha.f = 0.1))
  summary(lm(gapProp~shortProp, data = bci.gaps20[!is.na(bci.gaps20$areaObs) & bci.gaps20$areaObs>(0.5*cellSize*cellSize/10000),]))
  
  
  # how much was lo canopy in lidar data?
  dsm09 <- raster::raster("DSM_2009.tif")
  dem09 <- raster::raster("D:/BCI_Spatial/BCI_Topo/LidarDEM_BCI.tif")
  dem09 <- raster::crop(dem09,dsm09)
  raster::crs(dem09) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
  raster::crs(dsm09) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
  dem09 <- raster::resample(dem09,dsm09)
  
  chm09 <- dsm09-dem09
  chm09 <- raster::mask(chm09, buffer)
  chm09 <- raster::mask(chm09, ageUse)  
  
  chm09vals <- raster::values(chm09)
  length(chm09vals[chm09vals<15 & !is.na(chm09vals)])/ length(chm09vals[!is.na(chm09vals)])
  
  chm15vals <- raster::values(chm15)
  length(chm15vals[chm15vals<15 & !is.na(chm15vals)])/ length(chm15vals[!is.na(chm15vals)])
  length(chm15vals[is.na(chm15vals)])/length(chm15vals)
  
  
  chm18vals <- raster::values(chm18)
  length(chm18vals[chm18vals<15 & !is.na(chm18vals)])/ length(chm18vals[!is.na(chm18vals)])
  length(chm18vals[is.na(chm18vals)])/length(chm18vals)
  
  chm20vals <- raster::values(chm20)
  length(chm20vals[chm20vals<15 & !is.na(chm20vals)])/ length(chm20vals[!is.na(chm20vals)])
  length(chm20vals[is.na(chm20vals)])/length(chm20vals)
  
#### Fig S12. Fixed effects ####

  slopeAll <- bci.gapsAll$fix_S + bci.gapsAll$fix_S2
  htAll <- bci.gapsAll$fix_H + bci.gapsAll$fix_H2
  curvAll <- bci.gapsAll$fix_C
  
  par(mfrow=c(2,3), mar=c(0,4,0,1), oma=c(5,1,3,1))
  
  hist(bci.gapsAll$curvMean_8[!is.na(bci.gapsAll$gapPropCens)],
       xlim=range(bci.gapsAll$curvMean_8[!is.na(bci.gapsAll$gapPropCens)]),
       border="white",col="black",
       ylim=c(0,6000),
       main=NA, xaxt="n")
  
  hist(bci.gapsAll$slopeMean_24[!is.na(bci.gapsAll$gapPropCens)],
       xlim=range(bci.gapsAll$slopeMean_24[!is.na(bci.gapsAll$gapPropCens)]),
       border="white",col="black",
       ylim=c(0,6000),
       main=NA, xaxt="n", yaxt="n", ylab=NA)

  hist(bci.gapsAll$drainMean[!is.na(bci.gapsAll$gapPropCens)],
       xlim=range(bci.gapsAll$drainMean[!is.na(bci.gapsAll$gapPropCens)]),
       border="white",col="black",
       ylim=c(0,6000),
       main=NA, xaxt="n", yaxt="n", ylab=NA)
  
  plot(y = curvAll[order(bci.gapsAll$curvMean_8)],
       x = bci.gapsAll$curvMean_8[order(bci.gapsAll$curvMean_8)],
       xlim=range(bci.gapsAll$curvMean_8[!is.na(bci.gapsAll$gapPropCens)]),
       type= "l", lwd=2,
       ylim=c(-0.4,0.3), 
       ylab = "Fixed effect",
       col="black",
       cex=1.5)
  text("a", x = -2.7, y = 0.3)
  mtext("Curvature (LaPlacian convexity)",side=1,outer=F, line=3, cex = 0.8)
  
  
  plot(y = slopeAll[order(bci.gapsAll$slopeMean_24)],
       x = bci.gapsAll$slopeMean_24[order(bci.gapsAll$slopeMean_24)],       
       xlim=range(bci.gapsAll$slopeMean_24[!is.na(bci.gapsAll$gapPropCens)]),
       type= "l", lwd=2,
       ylim=c(-0.4,0.3),
       ylab = NA,
       col="black",
       cex=1.5)  
  text("b", x = 2, y = 0.3)
  mtext("Slope (degree)",side=1,outer=F, line=3, cex = 0.8)
  
  plot(y = htAll[order(bci.gapsAll$drainMean)],
       x = bci.gapsAll$drainMean[order(bci.gapsAll$drainMean)],
       xlim=range(bci.gapsAll$drainMean[!is.na(bci.gapsAll$gapPropCens)]),
       type= "l", lwd=2,
       pch=20, ylim=c(-0.4,0.3),
       col="black",
       ylab = NA,
       cex=1.5)
  text("c", x = 0, y = 0.3)
  mtext("Height above drainage (m)",side=1,outer=F, line=3, cex = 0.8)
  
  

  
#### MEETING PREP: look at distribution of random vs. fixed effects ####  
  # summed fixed effects
    hist(bci.gapsAll$fix_sum[!is.na(bci.gapsAll$gapPropCens)])
  # random effects
    hist(model_full$summary.random[[1]]$mean[!is.na(bci.gapsAll_Order$gapPropCens)],
         breaks=seq(-3,4.6,0.2),
         col="black",border="white",
         xlab="Random effects",main=NA)
  
  # total predicted values
    hist(bci.gapsAll$pred[!is.na(bci.gapsAll$gapPropCens)])
  # predicted values with just fixed effects  
    hist(bci.gapsAll$fix_pred[!is.na(bci.gapsAll$gapPropCens)])

#### MEETING PREP: aggregate and look at R2 ####          
    
  #Create rasters where all non-NA values are 1
    predMeanRaster18n <- predMeanRaster18
    predMeanRaster18n[!is.na(raster::values(predMeanRaster18n))] <- 1
    fixRaster18n <- fixRaster18
    fixRaster18n[!is.na(raster::values(fixRaster18n))] <- 1
    predMeanRaster20n <- predMeanRaster20
    predMeanRaster20n[!is.na(raster::values(predMeanRaster20n))] <- 1
    fixRaster20n <- fixRaster20
    fixRaster20n[!is.na(raster::values(fixRaster20n))] <- 1
    
  # Create data frame to store results  
    agResults <- data.frame(agBy = 1:20,
                            fixR2_18 = NA,
                            allR2_18 = NA,
                            fixR2_20 = NA,
                            allR2_20 = NA)
    
    for(i in 1:nrow(agResults)){
      
      agPred18_all <- raster::aggregate(predMeanRaster18, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred18_all_n <- raster::aggregate(predMeanRaster18n, fact = agResults$agBy[i], fun = sum, na.rm=T)
      agPred18_fix <- raster::aggregate(fixRaster18, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred18_fix_n <- raster::aggregate(fixRaster18n, fact = agResults$agBy[i], fun = sum, na.rm=T)
      agObs18 <- raster::aggregate(obsRaster18, fact = agResults$agBy[i], fun = mean, na.rm=T)
      
      minObs <- 0.9*agResults$agBy[i]^2
      useObs <- raster::values(agPred18_all_n)>minObs
      
      agResults$fixR2_18[i] <- summary(lm(raster::values(agPred18_fix)[useObs]~raster::values(agObs18)[useObs]))$r.squared
      agResults$allR2_18[i] <- summary(lm(raster::values(agPred18_all)[useObs]~raster::values(agObs18)[useObs]))$r.squared
      
      agPred20_all <- raster::aggregate(predMeanRaster20, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred20_all_n <- raster::aggregate(predMeanRaster20n, fact = agResults$agBy[i], fun = sum, na.rm=T)
      agPred20_fix <- raster::aggregate(fixRaster20, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred20_fix_n <- raster::aggregate(fixRaster20n, fact = agResults$agBy[i], fun = sum, na.rm=T)
      agObs20 <- raster::aggregate(obsRaster20, fact = agResults$agBy[i], fun = mean, na.rm=T)
      
      minObs <- 0.5*agResults$agBy[i]^2
      useObs <- raster::values(agPred20_all_n)>minObs
      
      agResults$fixR2_20[i] <- summary(lm(raster::values(agPred20_fix)[useObs]~raster::values(agObs20)[useObs]))$r.squared
      agResults$allR2_20[i] <- summary(lm(raster::values(agPred20_all)[useObs]~raster::values(agObs20)[useObs]))$r.squared
      
    }
    
    # Plot results
    
    par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,1,1,1), las=1)
    
    xVals <- (agResults$agBy*40)^2/10000
    
    plot(x = xVals,
         y = agResults$fixR2_18,
         ylim=range(agResults[,c("fixR2_18","fixR2_20")]),
         xlab = NA,
         ylab = "Observed variation explained (R^2)",
         main = "Fixed effects only",
         type = "l",
         col = "blue",
         lwd=2)
    lines(x = xVals,
         y = agResults$fixR2_20,
         col="lightblue",
         lwd=2)
    
    plot(x = xVals,
         y = agResults$allR2_18,
         ylim=range(agResults[,c("allR2_18","allR2_20")]),
         xlab = NA,
         ylab = NA,
         main = "Fixed + random effects",
         type = "l",
         col = "blue",
         lwd=2)
    lines(x = xVals,
          y = agResults$allR2_20,
          col="lightblue",
          lwd=2)
    
    mtext("Spatial resolution (ha)", side=1, outer=T, line=-1)
