#### Load data and results ####
library(INLA)
load("INLA/INLA_prelim_40m_tin.RData")
load("INLA/INLA_fullModelResult.RData")
load("INLA/INLA_fullModelResult_separate.RData")
load("INLA/INLA_fullModelResult_noLargeGaps.RData")
load("INLA/INLA_fullModelResult_initialHt.RData")

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

  # Colors used in plots
  col18 <- "blue"
  col20 <- "#d95f02"
  colBro <- wesanderson::wes_palette("Rushmore1",5)[1]
  colMot <- wesanderson::wes_palette("Rushmore1",5)[4]
  colPal <- wesanderson::wes_palette("Rushmore1",5)[3]
  colRed <- wesanderson::wes_palette("Rushmore1",5)[5]
  colAnd <- wesanderson::wes_palette("Chevalier1",4)[2]
  colBoh <- wesanderson::wes_palette("Chevalier1",4)[4]
  colMar <- wesanderson::wes_palette("Chevalier1",4)[3]
  colVol <- wesanderson::wes_palette("Chevalier1",4)[1]
  colOld <- wesanderson::wes_palette("Moonrise2",4)[1]
  colSec <- wesanderson::wes_palette("Moonrise2",4)[2]
  
# smoothScales <- c(1,2,3,4,6,8,12,16,24,32,48,64)
# curvScale <- 2
# slopeScale <- 16


#### Calculate predicted values and residuals ####
  
## Fixed and random effects  
  
  # Get fitted values and reorder to match original data frame
  bci.gapsAll$pred <- model_full$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  bci.gapsAll$resd <- bci.gapsAll$gapPropCens - bci.gapsAll$pred
  
## Fixed effects only  
  
    # Intercept
    bci.gapsAll$fix_int <- model_full$summary.fixed$mean[model_full$names.fixed=="(Intercept)"]
    
    # Curvature
    bci.gapsAll$fix_C <- bci.gapsAll$Sc_curvMean_2*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_curvMean_2"])
    
    
    # Slope
    bci.gapsAll$fix_S <- bci.gapsAll$Sc_slopeMean_16*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_slopeMean_16"])
    bci.gapsAll$fix_S2 <- bci.gapsAll$Sc_slopeMean_16_Sq*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_slopeMean_16_Sq"])
    
    # Height above drainage
    bci.gapsAll$fix_H <- bci.gapsAll$Sc_drainMean*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_drainMean"])
    bci.gapsAll$fix_H2 <- bci.gapsAll$Sc_drainMean_Sq*(model_full$summary.fixed$mean[model_full$names.fixed=="Sc_drainMean_Sq"])
    
    # Soil form
    bci.gapsAll$fix_soilForm <- 0

    bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="MottledHeavyClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormMottledHeavyClay"]
    bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="PaleSwellingClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormPaleSwellingClay"]
    bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="RedLightClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormRedLightClay"]

    
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
    bci.gapsAll$fix_sum <- bci.gapsAll$fix_int + bci.gapsAll$fix_C + bci.gapsAll$fix_S + bci.gapsAll$fix_S2 + bci.gapsAll$fix_H + bci.gapsAll$fix_H2 + bci.gapsAll$fix_soilParent + bci.gapsAll$fix_soilForm + bci.gapsAll$fix_age + bci.gapsAll$fix_year
    
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
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    plot(resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
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
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i],"resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
         pch=19,
         col = "orange")
    
  # Both years, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    plot(fix_resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
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
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i],"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    
  # 2015-2018, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    plot(fix_resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
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
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    
  # 2018-2020, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    plot(fix_resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
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
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i]
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
  
  # # Plot rasters
  # colBreaks <- seq(0,0.50,0.05)
  # 
  # 
#   # 2015-2108
#   jpeg(height=3000,width=2000,file = "Figure S10. Raster patterns.jpg")
#   par(mfrow=c(4,3))
#   raster::plot(obsRaster18,
#                bty = "n", box = F, xaxt="n", yaxt="n",
#                col = viridis::viridis(length(colBreaks)),
#                breaks = colBreaks,
#                main = "Observed")
#   raster::plot(buffer, add=T)
#   
#   raster::plot(predMeanRaster18,
#                bty = "n", box = F, xaxt="n", yaxt="n",
#                col = viridis::viridis(length(colBreaks)),
#                breaks = colBreaks,
#                main = "Predicted (mean)")
#   raster::plot(buffer, add=T)
#   
#   raster::plot(predSdRaster18,
#                bty = "n", box = F, xaxt="n", yaxt="n",
#                col = viridis::viridis(length(seq(0,0.1,0.01))),
#                breaks = seq(0,0.1,0.01),
#                main = "Predicted (SD)")
#   raster::plot(buffer, add=T)
#   mtext("2015 - 2018", outer=F, side = 3, line = -1.5)
#   
#   # 2018-2020
#   raster::plot(obsRaster20,
#                bty = "n", box = F, xaxt="n", yaxt="n",
#                col = viridis::viridis(length(colBreaks)),
#                breaks = colBreaks,
#                main = "Observed")
#   raster::plot(buffer, add=T)
#   
#   raster::plot(predMeanRaster20,
#                bty = "n", box = F, xaxt="n", yaxt="n",
#                col = viridis::viridis(length(colBreaks)),
#                breaks = colBreaks,
#                main = "Predicted (mean)")
#   raster::plot(buffer, add=T)
#   
#   raster::plot(predSdRaster20,
#                bty = "n", box = F, xaxt="n", yaxt="n",
#                col = viridis::viridis(length(seq(0,0.1,0.01))),
#                breaks = seq(0,0.1,0.01),
#                main = "Predicted (SD)")
#   raster::plot(buffer, add=T)
#   mtext("2018 - 2020", outer=F, side = 3, line = -1.5)
#   
# dev.off()

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
               col=viridis::cividis(50),
               bty="n", box=F, xaxt="n", yaxt="n")  
  raster::plot(predMeanRaster18,
               col=viridis::viridis(50),
               breaks = seq(range(raster::values(predMeanRaster20),na.rm=T)[1],
                            range(raster::values(predMeanRaster20),na.rm=T)[2],
                            length.out = 50),
               bty="n", box=F, xaxt="n", yaxt="n",
               main = "Predicted (fixed + random) 2015-2018")  
  raster::plot(predMeanRaster20,
               col=viridis::viridis(50),
               breaks = seq(range(raster::values(predMeanRaster20),na.rm=T)[1],
                            range(raster::values(predMeanRaster20),na.rm=T)[2],
                            length.out = 50),
               bty="n", box=F, xaxt="n", yaxt="n",
               main = "Predicted (fixed + random) 2018-2020")  
  
#### Plot S#: Predicted vs observed values at 40 m resolution ####
summary(lm(gapPropCens~pred, data=bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),]))
  
  smoothAll <- loess.smooth(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"pred"],
                            y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],
                            degree=1, span = 1.2)
  smooth18 <- loess.smooth(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2018","pred"],
                            y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2018","gapPropCens"],
                            degree=1, span = 1.2)
  smooth20 <- loess.smooth(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2020","pred"],
                            y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2020","gapPropCens"],
                            degree=1, span = 1.2)
  
  
  par(las=1, mfrow=c(1,3), mar=c(3,3,1,0),oma=c(1,3,1,1))
  plot(gapPropCens~pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       col = adjustcolor("black", 0.08),
       main = "All years",
       ylab = NA,
       xlab = NA,
       cex = 0.5,
       pch=19)
  lines(smoothAll, col="red",lwd=2, lty=2)
  abline(a=0,b=1,col="red")
  mtext("Observed disturbance proportion", side=2,outer=T, las=0)
  mtext("Predicted disturbance proportion (fixed + random effects)", side=1,outer=T)

  
  plot(gapPropCens~pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year=="2018",],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       col = adjustcolor(col18, 0.08),
       main = "2015-2018",
       ylab = NA,
       xlab = NA,
       yaxt="n",
       cex = 0.5,
       pch=19)
  lines(smooth18, col="red", lwd=2, lty=2)
  abline(a=0,b=1,col="red")
  
  plot(gapPropCens~pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year=="2020",],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"],na.rm=T),
       col = adjustcolor(col20, 0.08),
       main = "2015-2018",
       ylab = NA,
       xlab = NA,
       yaxt="n",
       cex = 0.5,
       pch=19)
  lines(smooth20, col="red", lwd=2, lty=2)
  abline(a=0,b=1,col="red")
  
  
  smoothAllb <- loess.smooth(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_pred"],
                            y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"],
                            degree=1, span = 1.2)
  smooth18b <- loess.smooth(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2018","fix_pred"],
                           y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2018","resd"],
                           degree=1, span = 1.2)
  smooth20b <- loess.smooth(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2020","fix_pred"],
                           y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens)& bci.gapsAll$Year=="2020","resd"],
                           degree=1, span = 1.2)
  
  par(las=1, mfrow=c(1,3), mar=c(3,3,1,0),oma=c(1,3,1,1))
  plot(resd~fix_pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_pred"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"],na.rm=T),
       col = adjustcolor("black", 0.08),
       #main = "All years",
       ylab = NA,
       xlab = NA,
       cex = 0.5,
       pch=19)
  lines(smoothAllb, col="red",lwd=2, lty=1)
  mtext("Predicted disturbance distribubtion (fixed effects only)", side=1,outer=T)
  mtext("Residual", side=2,outer=T, las=0)

  
  plot(resd~fix_pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year=="2018",],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_pred"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"],na.rm=T),
       col = adjustcolor(col18, 0.08),
       #main = "2015-2018",
       ylab = NA,
       xlab = NA,
       yaxt="n",
       cex = 0.5,
       pch=19)
  lines(smooth18b, col="red", lwd=2, lty=1)

  plot(resd~fix_pred, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year=="2020",],
       xlim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_pred"],na.rm=T),
       ylim=range(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"],na.rm=T),
       col = adjustcolor(col20, 0.08),
       #main = "2015-2018",
       ylab = NA,
       xlab = NA,
       yaxt="n",
       cex = 0.5,
       pch=19)
  lines(smooth20b, col="red", lwd=2, lty=1)

#### Figure 4c: Fixed effects sizes ####  
  fixedResults <- model_full$summary.fixed
  
  # separated by years
  fixedResults1 <- model_full1$summary.fixed
  fixedResults2 <- model_full2$summary.fixed
  
  # main model and second interval without large gaps
  fixedResults_alt <- model_full_alt$summary.fixed
  fixedResults2_alt <- model_full2_alt$summary.fixed
  
  
  # VERSION 1: main figure with only full model
  
  fixedNames <- c("Curvature (linear)", "Slope (linear)", "Slope (quadratic)",                    
                  "Height above drainage (linear)","Height above drainage (quadratic)",
                  "Soil parent: Andesite", "Soil parent: Caimito marine", "Soil parent: Caimito volcanic",
                  "Soil form: Heavy clay", "Soil form: Swelling clay", "Soil form: Light clay",
                  "Forest age: Secondary",
                  "Year: 2015-2018")
  
  par(mfrow=c(1,1), mar=c(4,1,1,1))
  plot(x = fixedResults[2:nrow(fixedResults),"mean"],
       y = nrow(fixedResults):2,
       xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(0,0.6),
       ylim=c(2,nrow(fixedResults)),
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
  
  abline(v=0, lty=2)
  text(fixedNames,
       x = fixedResults[2:nrow(fixedResults),"0.975quant"]+0.005,
       y = nrow(fixedResults):2,
       pos = 4)
  
  # VERSION 2: separate results for each year

  
  fixedNames <- c("Curvature (linear)", "Slope (linear)", "Slope (quadratic)",                    
                  "Height above drainage (linear)","Height above drainage (quadratic)",
                  "Soil parent: Andesite", "Soil parent: Caimito marine", "Soil parent: Caimito volcanic",
                  "Soil form: Heavy clay", "Soil form: Swelling clay", "Soil form: Light clay",
                  "Forest age: Secondary",
                  "Year: 2015-2018")

  par(mfrow=c(1,1), mar=c(4,1,1,1))
  plot(x = fixedResults[2:nrow(fixedResults),"mean"],
       y = 1.2*nrow(fixedResults):2,
       xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(-0.4,0.6),
       ylim=c(2,1.2*nrow(fixedResults)),
       xlab=NA,
       pch = 19, 
       cex = 1,
       ylab = NA, yaxt = "n")
  arrows(x0 = fixedResults[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2,
         y1 = 1.2*nrow(fixedResults):2,
         angle = 90, code=3,
         length = 0.05)
  
  points(x = fixedResults1[2:nrow(fixedResults),"mean"],
       y = 1.2*nrow(fixedResults):2 - 0.2,
       pch = 19, 
       cex = 1,
       col = col18)
  arrows(x0 = fixedResults1[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults1[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2- 0.2,
         y1 = 1.2*nrow(fixedResults):2- 0.2,
         angle = 90, code=3,
         length = 0.05,
         col = col18)
  
  points(x = fixedResults2[2:nrow(fixedResults),"mean"],
         y = 1.2*nrow(fixedResults):2 - 0.4,
         pch = 19, 
         cex = 1,
         col = col20)
  arrows(x0 = fixedResults2[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults2[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2- 0.4,
         y1 = 1.2*nrow(fixedResults):2- 0.4,
         angle = 90, code=3,
         length = 0.05,
         col = col20)
  
  points(x = fixedResults2_alt[2:nrow(fixedResults),"mean"],
         y = 1.2*nrow(fixedResults):2 - 0.6,
         pch = 1, 
         cex = 1,
         col = col20)
  arrows(x0 = fixedResults2_alt[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults2_alt[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2- 0.6,
         y1 = 1.2*nrow(fixedResults):2- 0.6,
         angle = 90, code=3,
         length = 0.05,
         lty=2,
         col = col20)
  
  
  abline(v=0, lty=2)
  text(fixedNames,
      x = fixedResults[2:nrow(fixedResults),"0.975quant"]+0.03,
       y = 1.2*nrow(fixedResults):2,
       pos = 4)
  
  legend(x = -1.3,
         y = 6,
         c("Full model",
           "2015-2018 only",
           "2018-2020 only",
           "2018-2020 only (no blowdown)"),
         bty="n",
         col = c("black", col18, col20, col20),
         pch=c(19,19,19,1))
  text("c", x = -1.25, y = 1.2*nrow(fixedResults))
  mtext("Fixed effect", side=1, outer=F, line=2)
  
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
  
#### Fig S12. Fixed topographic effects ####

  slopeAll <- bci.gapsAll$fix_S + bci.gapsAll$fix_S2
  htAll <- bci.gapsAll$fix_H + bci.gapsAll$fix_H2
  curvAll <- bci.gapsAll$fix_C
  
  # calculate mean abs contribution of slope, HAND, curvature to linear predictor
  mean(abs(curvAll[!is.na(bci.gapsAll$gapPropCens)]))
  mean(abs(slopeAll[!is.na(bci.gapsAll$gapPropCens)]))
  mean(abs(htAll[!is.na(bci.gapsAll$gapPropCens)]))
  
  # Make a simplified version for plotting
  
    # Find the relationship between each metric and its scaled (and squared) value
      slopeSq <- bci.gapsAll$slopeMean_16^2
      slopeScaleLinear <- lm(Sc_slopeMean_16~slopeMean_16, data = bci.gapsAll)
      slopeScaleQuad <- lm(bci.gapsAll$Sc_slopeMean_16_Sq~slopeSq)
      
      curveScaleLinear <- lm(Sc_curvMean_2~curvMean_2, data = bci.gapsAll)
      
      drainSq <- bci.gapsAll$drainMean^2
      drainScaleLinear <- lm(Sc_drainMean~drainMean, data = bci.gapsAll)
      drainScaleQuad <- lm(bci.gapsAll$Sc_drainMean_Sq~drainSq)
    
    # Find the range of curvature, slope, and HAND; choose new values for plotting
      slopeRange <- seq(range(bci.gapsAll$slopeMean_16)[1],
                        range(bci.gapsAll$slopeMean_16)[2],
                        length.out=100)
      curvRange <- seq(range(bci.gapsAll$curvMean_2)[1],
                       range(bci.gapsAll$curvMean_2)[2],
                       length.out=100)
      drainRange <- seq(range(bci.gapsAll$drainMean)[1],
                        range(bci.gapsAll$drainMean)[2],
                        length.out=100)
    
    # Generate new y values for ranges from linear models
      curvVals <- model_full$summary.fixed[2,"mean"]*predict(curveScaleLinear, data.frame(curvMean_2 = curvRange))
      
      slopeVals <- (model_full$summary.fixed[3,"mean"]*predict(slopeScaleLinear, data.frame(slopeMean_16 = slopeRange)) 
                    + model_full$summary.fixed[4,"mean"]*predict(slopeScaleQuad, data.frame(slopeSq = slopeRange^2)))
      
      drainVals <- (model_full$summary.fixed[5,"mean"]*predict(drainScaleLinear, data.frame(drainMean = drainRange)) 
                    + model_full$summary.fixed[6,"mean"]*predict(drainScaleQuad, data.frame(drainSq = drainRange^2)))
        
  
  par(mfrow=c(2,3), mar=c(0,4,0,1), oma=c(5,1,3,1), las=1)
  
  hist(bci.gapsAll$curvMean_2[!is.na(bci.gapsAll$gapPropCens)],
       xlim=range(bci.gapsAll$curvMean_2[!is.na(bci.gapsAll$gapPropCens)]),
       border="white",col="black",
       ylim=c(0,6000),
       main=NA, xaxt="n")
  
  hist(bci.gapsAll$slopeMean_16[!is.na(bci.gapsAll$gapPropCens)],
       xlim=range(bci.gapsAll$slopeMean_16[!is.na(bci.gapsAll$gapPropCens)]),
       border="white",col="black",
       ylim=c(0,6000),
       main=NA, xaxt="n", yaxt="n", ylab=NA)

  hist(bci.gapsAll$drainMean[!is.na(bci.gapsAll$gapPropCens)],
       xlim=range(bci.gapsAll$drainMean[!is.na(bci.gapsAll$gapPropCens)]),
       border="white",col="black",
       ylim=c(0,6000),
       main=NA, xaxt="n", yaxt="n", ylab=NA)
  
  plot(y = curvVals,
       x = curvRange, 
       xlim=range(bci.gapsAll$curvMean_2[!is.na(bci.gapsAll$gapPropCens)]),
       type= "l", lwd=2,
       ylim=c(-0.4,0.35), 
       ylab = "Fixed effect",
       col="black",
       cex=1.5)
  text("a", x = -3.9, y = 0.35)
  mtext("Curvature (LaPlacian convexity)",side=1,outer=F, line=3, cex = 0.8)
  
  
  plot(y = slopeVals,
       x = slopeRange,       
       xlim=range(bci.gapsAll$slopeMean_16[!is.na(bci.gapsAll$gapPropCens)]),
       type= "l", lwd=2,
       ylim=c(-0.4,0.35),
       ylab = NA,
       col="black",
       cex=1.5)  
  text("b", x = 2, y = 0.35)
  mtext("Slope (degree)",side=1,outer=F, line=3, cex = 0.8)
  
  plot(y = drainVals,
       x = drainRange,
       xlim=range(bci.gapsAll$drainMean[!is.na(bci.gapsAll$gapPropCens)]),
       type= "l", lwd=2,
       pch=20, ylim=c(-0.4,0.35),
       col="black",
       ylab = NA,
       cex=1.5)
  text("c", x = 0, y = 0.35)
  mtext("Height above drainage (m)",side=1,outer=F, line=3, cex = 0.8)
  
  
  
#### Figure 7: aggregate and look at R2 ####          
    
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
                            fixR_18 = NA,
                            fixR_18lo = NA,
                            fixR_18hi = NA,
                            allR_18 = NA,
                            allR_18lo = NA,
                            allR_18hi = NA,
                            fixR_20 = NA,
                            fixR_20lo = NA,
                            fixR_20hi = NA,
                            allR_20 = NA,
                            allR_20lo = NA,
                            allR_20hi = NA)
    
    for(i in 1:nrow(agResults)){
      
      # Get a matrix of observed cells
      minObs <- 0.75*agResults$agBy[i]^2
      
      agPred18_all <- raster::aggregate(predMeanRaster18, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred18_fix <- raster::aggregate(fixRaster18, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agObs18 <- raster::aggregate(obsRaster18, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred18_all_n <- raster::aggregate(predMeanRaster18n, fact = agResults$agBy[i], fun = sum, na.rm=T)
      
      # Get matrix of values for each aggregated raster
      agPred18_all_mat <- raster::values(agPred18_all, format = "matrix")
      agPred18_fix_mat <- raster::values(agPred18_fix, format = "matrix")
      agObs18_mat <- raster::values(agObs18, format = "matrix")
      nMat <- raster::values(agPred18_all_n, format = "matrix")
      
        # matrix of cells that need to be combined
        combineMat <- which(nMat<minObs, arr.ind = T)
        
        while(nrow(combineMat)>0){
        
          # get surrounding cells for each matrix
            rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
            rows <- rows[which(rows>0 & rows <= nrow(nMat))]
            cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
            cols <- cols[which(cols>0 & cols <= ncol(nMat))]
            
            agPred18_all_mat_j <- agPred18_all_mat[rows,cols]
            agPred18_fix_mat_j <- agPred18_fix_mat[rows,cols]
            agObs18_mat_j <- agObs18_mat[rows,cols]
            nMat_j <- nMat[rows,cols]
              # set cell of interest to 0
              nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- NaN
          
          # Only proceed if there are non-NA neighboring cells
          if(length(c(nMat_j)[!is.na(c(nMat_j))])>0){  
              
            # find neighboring cell with most nearby observations
              new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
              
            
            # replace value for new cell with a weighted mean of all other values
              agPred18_all_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agPred18_all_mat_j[new_j[1],new_j[2]],agPred18_all_mat[combineMat[1,1],combineMat[1,2]]),
                                                                     w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
              agPred18_fix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agPred18_fix_mat_j[new_j[1],new_j[2]],agPred18_fix_mat[combineMat[1,1],combineMat[1,2]]),
                                                                     w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
              agObs18_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs18_mat_j[new_j[1],new_j[2]],agObs18_mat[combineMat[1,1],combineMat[1,2]]),
                                                                     w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
              nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
              
            # replace neighborhood in original matrices
              agPred18_all_mat[rows,cols] <- agPred18_all_mat_j
              agPred18_fix_mat[rows,cols] <- agPred18_fix_mat_j
              agObs18_mat[rows,cols] <- agObs18_mat_j
              nMat[rows,cols] <- nMat_j
          }
              
          # Replace original observation with NaN
            agPred18_all_mat[combineMat[1,1],combineMat[1,2]] <- NaN
            agPred18_fix_mat[combineMat[1,1],combineMat[1,2]] <- NaN
            agObs18_mat[combineMat[1,1],combineMat[1,2]] <- NaN
            nMat[combineMat[1,1],combineMat[1,2]] <- NaN
            
          # Redefine combineMat
            combineMat <- which(nMat<minObs, arr.ind = T)
        }
        
      agResults$fixR_18[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$estimate
      agResults$fixR_18lo[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[1]
      agResults$fixR_18hi[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[2]
      agResults$allR_18[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_all_mat))$estimate
      agResults$allR_18lo[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_all_mat))$conf.int[1]
      agResults$allR_18hi[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_all_mat))$conf.int[2]
      
      
      agPred20_all <- raster::aggregate(predMeanRaster20, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred20_fix <- raster::aggregate(fixRaster20, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agObs20 <- raster::aggregate(obsRaster20, fact = agResults$agBy[i], fun = mean, na.rm=T)
      agPred20_all_n <- raster::aggregate(predMeanRaster20n, fact = agResults$agBy[i], fun = sum, na.rm=T)
      
      # Get matrix of values for each aggregated raster
      agPred20_all_mat <- raster::values(agPred20_all, format = "matrix")
      agPred20_fix_mat <- raster::values(agPred20_fix, format = "matrix")
      agObs20_mat <- raster::values(agObs20, format = "matrix")
      nMat <- raster::values(agPred20_all_n, format = "matrix")
      
      # matrix of cells that need to be combined
      combineMat <- which(nMat<minObs, arr.ind = T)
      
      while(nrow(combineMat)>0){
        
        # get surrounding cells for each matrix
        rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
        rows <- rows[which(rows>0 & rows <= nrow(nMat))]
        cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
        cols <- cols[which(cols>0 & cols <= ncol(nMat))]
        
        agPred20_all_mat_j <- agPred20_all_mat[rows,cols]
        agPred20_fix_mat_j <- agPred20_fix_mat[rows,cols]
        agObs20_mat_j <- agObs20_mat[rows,cols]
        nMat_j <- nMat[rows,cols]
        # set cell of interest to 0
        nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- NaN
        
        # Only proceed if there are non-NA neighboring cells
        if(length(c(nMat_j)[!is.na(c(nMat_j))])>0){
          
          # find neighboring cell with most nearby observations
          new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
          
          
          # replace value for new cell with a weighted mean of all other values
          agPred20_all_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agPred20_all_mat_j[new_j[1],new_j[2]],agPred20_all_mat[combineMat[1,1],combineMat[1,2]]),
                                                                 w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agPred20_fix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agPred20_fix_mat_j[new_j[1],new_j[2]],agPred20_fix_mat[combineMat[1,1],combineMat[1,2]]),
                                                                 w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agObs20_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs20_mat_j[new_j[1],new_j[2]],agObs20_mat[combineMat[1,1],combineMat[1,2]]),
                                                            w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
          
          # replace neighborhood in original matrices
          agPred20_all_mat[rows,cols] <- agPred20_all_mat_j
          agPred20_fix_mat[rows,cols] <- agPred20_fix_mat_j
          agObs20_mat[rows,cols] <- agObs20_mat_j
          nMat[rows,cols] <- nMat_j
        }
        
        # Replace original observation with NaN
        agPred20_all_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agPred20_fix_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agObs20_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        nMat[combineMat[1,1],combineMat[1,2]] <- NaN
      
        # Redefine combineMat
        combineMat <- which(nMat<minObs, arr.ind = T)
      }
      
      agResults$fixR_20[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$estimate
      agResults$fixR_20lo[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[1]
      agResults$fixR_20hi[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[2]
      agResults$allR_20[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_all_mat))$estimate
      agResults$allR_20lo[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_all_mat))$conf.int[1]
      agResults$allR_20hi[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_all_mat))$conf.int[2]
      
    }
    
    # Plot results
    
    par(mfrow=c(1,2), mar=c(4,1,1,0), oma=c(1,4,1,1), las=1)
    
    xVals <- (agResults$agBy*40)^2/10000
    
    plot(x = xVals,
         y = agResults$fixR_18,
         log="x",
         ylim=c(0,1.02),
         xlab = NA,
         ylab = NA,
         type = "l",
         col = col18,
         lwd=2)
    text("a. Fixed effects only", x = 0.15, y = 1.02, adj=0)
    lines(x = xVals,
         y = agResults$fixR_20,
         col=col20,
         lwd=2)
    polygon(x = c(xVals,rev(xVals),xVals[1]),
            y = c(agResults$fixR_18lo,rev(agResults$fixR_18hi),agResults$fixR_18lo[1]),
            col=adjustcolor(col18,0.15),
            border=NA)
    polygon(x = c(xVals,rev(xVals),xVals[1]),
            y = c(agResults$fixR_20lo,rev(agResults$fixR_20hi),agResults$fixR_20lo[1]),
            col=adjustcolor(col20,0.15),
            border=NA)
    legend(c("2015-2018",
             "2018-2020"),
           x=0.15,y=0.95,
           col=c(col18,col20),
           bty="n",
           lwd=2)
    
    plot(x = xVals,
         y = agResults$allR_18,
         ylim=c(0,1.02),
         log="x",
         yaxt="n",
         xlab = NA,
         ylab = NA,
         type = "l",
         col = col18,
         lwd=2)
    text("b. Fixed + random effects", x = 0.15, y = 1.02, adj=0)
    lines(x = xVals,
          y = agResults$allR_20,
          col=col20,
          lwd=2)
    polygon(x = c(xVals,rev(xVals),xVals[1]),
            y = c(agResults$allR_18lo,rev(agResults$allR_18hi),agResults$allR_18lo[1]),
            col=adjustcolor(col18,0.15),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(agResults$allR_20lo,rev(agResults$allR_20hi)),
            col=adjustcolor(col20,0.15),
            border=NA)

    
    mtext("Spatial grain (ha)", side=1, outer=T, line=-2)
    mtext("Pearson correlation (r)", side=2, outer=T, line=1.5, las=0)

#### Define forest age and soil type polygons ####
    # Forest age polygon
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    age$AgeClass <- "Other"
    age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
    age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    ageUse <- age[!(age$AgeClass=="Other"),]
    
    # Soil type polygon  
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(age))
    
    # Define parent material and soil form from soil class
    soil$SoilParent <- NA
    soil[soil$SOIL=="AVA", c("SoilParent")] <- c("Andesite")
    soil[soil$SOIL=="Barbour", c("SoilParent")] <- c("CaimitoVolcanic")
    soil[soil$SOIL=="Fairchild",c("SoilParent")] <- c("Bohio")
    soil[soil$SOIL=="Gross",c("SoilParent")] <- c("Bohio")
    soil[soil$SOIL=="Harvard",c("SoilParent")] <- c("CaimitoVolcanic")
    soil[soil$SOIL=="Hood",c("SoilParent")] <- c("CaimitoVolcanic")
    soil[soil$SOIL=="Lake",c("SoilParent")] <- c("Andesite")
    soil[soil$SOIL=="Lutz",c("SoilParent")] <- c("CaimitoMarineSedimentary")
    soil[soil$SOIL=="Marron",c("SoilParent")] <- c("Andesite")
    soil[soil$SOIL=="Poacher",c("SoilParent")] <- c("CaimitoMarineSedimentary")
    soil[soil$SOIL=="Standley",c("SoilParent")] <- c("Bohio")
    soil[soil$SOIL=="Wetmore",c("SoilParent")] <- c("CaimitoMarineSedimentary")
    soil[soil$SOIL=="Zetek",c("SoilParent")] <- c("CaimitoMarineSedimentary")
    
    soil$SoilForm <- NA
    soil[soil$SOIL=="AVA", c("SoilForm")] <- c("RedLightClay")
    soil[soil$SOIL=="Barbour", c("SoilForm")] <- c("PaleSwellingClay")
    soil[soil$SOIL=="Fairchild",c("SoilForm")] <- c("RedLightClay")
    soil[soil$SOIL=="Gross",c("SoilForm")] <- c("PaleSwellingClay")
    soil[soil$SOIL=="Harvard",c("SoilForm")] <- c("RedLightClay")
    soil[soil$SOIL=="Hood",c("SoilForm")] <- c("BrownFineLoam")
    soil[soil$SOIL=="Lake",c("SoilForm")] <- c("PaleSwellingClay")
    soil[soil$SOIL=="Lutz",c("SoilForm")] <- c("MottledHeavyClay")
    soil[soil$SOIL=="Marron",c("SoilForm")] <- c("BrownFineLoam")
    soil[soil$SOIL=="Poacher",c("SoilForm")] <- c("RedLightClay")
    soil[soil$SOIL=="Standley",c("SoilForm")] <- c("BrownFineLoam")
    soil[soil$SOIL=="Wetmore",c("SoilForm")] <- c("BrownFineLoam")
    soil[soil$SOIL=="Zetek",c("SoilForm")] <- c("PaleSwellingClay")
    
#### Figure 6ab: Comparison of 2009 lidar data and average spatial pattern ####
    
# Raster of low canopy area in 2009
  lo09 <- raster::raster("binaryLoCanopy.tif") # pixels with value 1 are => 10 m height, 0 are < 10 m
# Raster of canopy height in 2009
  chm09 <- raster::raster("CHM_2009_QAQC.tif") 
    
      
# Make raster of overall sampling effort
  avgRasterN <- 0.5*raster::calc(raster::stack(predMeanRaster18n,predMeanRaster20n),sum, na.rm=T)
    
# Make raster of low canopy proportion using the same sampling as the INLA analysis
  loCrop <- raster::crop(lo09, raster::extent(bci.gaps18))
  rasterAll09 <- loCrop
  rasterAll09[!is.na(raster::values(rasterAll09))] <- 1
  resampleAll09 <- raster::aggregate(rasterAll09, cellSize, fun= sum)
  resampleHi09 <- raster::aggregate(loCrop, cellSize, fun= sum)
  rasterLo09 <- (resampleAll09-resampleHi09)/resampleAll09

# Make raster of mean canopy height using the same sampling as the INLA analysis
  chmCrop <- raster::crop(chm09, raster::extent(bci.gaps18))
  rasterchm09 <- raster::aggregate(chmCrop, cellSize, fun= mean, na.rm=T)
  
# Create data frame to store results  
  loResults <- data.frame(agBy = 1:20,
                          fixRlo = NA,
                          fixRlo_lo = NA,
                          fixRlo_hi = NA,
                          obsRlo = NA,
                          obsRlo_lo = NA,
                          obsRlo_hi = NA,
                          fixRchm = NA,
                          fixRchm_lo = NA,
                          fixRchm_hi = NA,
                          obsRchm = NA,
                          obsRchm_lo = NA,
                          obsRchm_hi = NA,
                          N = NA)
  
  for(i in 1:nrow(loResults)){
    
    # Get a matrix of observed cells
    minObs <- 0.75*loResults$agBy[i]^2
    
    # Aggregate rasters
      # Aggregate low canopy area raster
      agLo09 <- raster::aggregate(rasterLo09, fact = loResults$agBy[i], fun = mean, na.rm=T)
      # Aggregate mean canopy height raster
      agChm09 <- raster::aggregate(rasterchm09, fact = loResults$agBy[i], fun = mean, na.rm=T)
      # Aggregate average fixed effects raster
      agFix <- raster::aggregate(avgPredictedRaster, fact = loResults$agBy[i], fun = mean, na.rm=T)
      # Aggregate average observed disturbance raster
      agObs <- raster::aggregate(avgObservedRaster, fact = loResults$agBy[i], fun = mean, na.rm=T)
      # Aggregate sampling effort
      agN <- raster::aggregate(avgRasterN, fact = loResults$agBy[i], fun = sum, na.rm=T)
    
    # Get matrix of values for each aggregated raster
      agLo09_mat <- raster::values(agLo09, format = "matrix")
      agChm09_mat <- raster::values(agChm09, format = "matrix")
      agFix_mat <- raster::values(agFix, format = "matrix")
      agObs_mat <- raster::values(agObs, format = "matrix")
      nMat <- raster::values(agN, format = "matrix")
    
    # matrix of cells that need to be combined
    combineMat <- which(nMat<minObs & nMat>0, arr.ind = T)
    
    while(nrow(combineMat)>0){
      
      # get surrounding cells for each matrix
      rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
      rows <- rows[which(rows>0 & rows <= nrow(nMat))]
      cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
      cols <- cols[which(cols>0 & cols <= ncol(nMat))]
      
      agLo09_mat_j <- agLo09_mat[rows,cols]
      agChm09_mat_j <- agChm09_mat[rows,cols]
      agFix_mat_j <- agFix_mat[rows,cols]
      agObs_mat_j <- agObs_mat[rows,cols]
      nMat_j <- nMat[rows,cols]
      
      # set cell of interest to 0
      nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- 0
      
      # Only proceed if there are non-NA neighboring cells
      if(length(c(nMat_j)[c(nMat_j)>0])>0){  
        
        # find neighboring cell with most nearby observations
        new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
        
        # replace value for new cell with a weighted mean of all other values
        agLo09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agLo09_mat_j[new_j[1],new_j[2]],agLo09_mat[combineMat[1,1],combineMat[1,2]]),
                                                               w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agChm09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agChm09_mat_j[new_j[1],new_j[2]],agChm09_mat[combineMat[1,1],combineMat[1,2]]),
                                                               w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agFix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agFix_mat_j[new_j[1],new_j[2]],agFix_mat[combineMat[1,1],combineMat[1,2]]),
                                                             w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agObs_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs_mat_j[new_j[1],new_j[2]],agObs_mat[combineMat[1,1],combineMat[1,2]]),
                                                             w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))

        nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
        
        # replace neighborhood in original matrices
        agLo09_mat[rows,cols] <- agLo09_mat_j
        agChm09_mat[rows,cols] <- agChm09_mat_j
        agFix_mat[rows,cols] <- agFix_mat_j
        agObs_mat[rows,cols] <- agObs_mat_j
        nMat[rows,cols] <- nMat_j
      }
      
      # Replace original observation with NaN
      agLo09_mat[combineMat[1,1],combineMat[1,2]] <- NaN
      agChm09_mat[combineMat[1,1],combineMat[1,2]] <- NaN
      agFix_mat[combineMat[1,1],combineMat[1,2]] <- NaN
      agObs_mat[combineMat[1,1],combineMat[1,2]] <- NaN
      nMat[combineMat[1,1],combineMat[1,2]] <- 0
      
      # Redefine combineMat
      combineMat <- which(nMat<minObs & nMat >0, arr.ind = T)
    }
    
    valsKeep <- which(c(nMat)>minObs)
    loResults$fixRlo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$estimate
    loResults$fixRlo_lo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[1]
    loResults$fixRlo_hi[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[2]
    loResults$obsRlo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$estimate
    loResults$obsRlo_lo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[1]
    loResults$obsRlo_hi[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[2]
    
    loResults$fixRchm[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$estimate
    loResults$fixRchm_lo[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat[valsKeep]))$conf.int[1]
    loResults$fixRchm_hi[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[2]
    loResults$obsRchm[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$estimate
    loResults$obsRchm_lo[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[1]
    loResults$obsRchm_hi[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[2]
    
  }
  
  
  # Plot results
  
  # Significance with aggregation scale    
    par(mfrow=c(1,2), mar=c(4,1,1,1), oma=c(2,4,1,1), las=1)
    
    xVals <- ((loResults$agBy*40)^2)/10000
    
    plot(x = xVals,
         y = loResults$fixRlo,
         log = "x",
         ylim=c(-1,1.02),
         xlab = NA,
         ylab = NA,
         type = "l",
         col = "blue",
         lwd=2)
    par(las=0)
    mtext(expression("Pearson correlation (r)"),
          side=2, outer=T, line=2)
    par(las=1)
    text("a. Proportion of low canopy area",
         x = 0.15,
         y = 1.02, adj=0)
    abline(h=0,lty=2)
    lines(x = xVals,
          y = loResults$obsRlo,
          col="lightblue",
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$fixRlo_lo,rev(loResults$fixRlo_hi)),
            col=adjustcolor("blue",0.15),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$obsRlo_lo,rev(loResults$obsRlo_hi)),
            col=adjustcolor("lightblue",0.4),
            border=NA)
    legend(x = 0.15,
           y = 0,
           c("Predicted frequency (fixed effects)",
             "Average observed frequency"),
           col=c("blue","lightblue"),
           lty=c(1,1),
           lwd=2,
           bty="n")
    
    plot(x = xVals,
         y = loResults$fixRchm,
         log = "x",
         yaxt="n",
         ylim=c(-1,1.02),
         xlab = NA,
         ylab = NA,
         type = "l",
         col = "blue",
         lwd=2)
    par(las=0)
    abline(h=0,lty=2)

    par(las=1)
    lines(x = xVals,
          y = loResults$obsRchm,
          col="lightblue",
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$fixRchm_lo,rev(loResults$fixRchm_hi)),
            col=adjustcolor("blue",0.15),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$obsRchm_lo,rev(loResults$obsRchm_hi)),
            col=adjustcolor("lightblue",0.4),
            border=NA)
    
  

    text("b. Mean canopy height",
         x = 0.15,
         y = 1.02, adj=0)

    
    mtext("Spatial resolution (ha)", side=1, outer=T, line=-1)
    
    
    
#### Figures S23: Plot correlations with best scales ####
    
    # Get forest age and soil type values for each pixel
    
    # Scale 1
    agScale_1 <- 5

    # Aggregate low canopy area raster
    agLo09_1 <- raster::aggregate(rasterLo09, fact = agScale_1, fun = mean, na.rm=T)
    # Aggregate mean canopy height raster
    agChm09_1 <- raster::aggregate(rasterchm09, fact = agScale_1, fun = mean, na.rm=T)
    # Aggregate average fixed effects raster
    agFix_1 <- raster::aggregate(avgPredictedRaster, fact = agScale_1, fun = mean, na.rm=T)
    # Aggregate average observed disturance raster
    agObs_1 <- raster::aggregate(avgObservedRaster, fact = agScale_1, fun = mean, na.rm=T)
    # Aggregate sampling effort
    agN_1 <- raster::aggregate(avgRasterN, fact = agScale_1, fun = sum, na.rm=T)
    # Find observations with 90% of values present
    minObs_1 <- 0.75*loResults$agBy[agScale_1]^2
    useObs_1 <- raster::values(agN_1)>minObs_1
    
    # Make forest age raster
    agAge_1 <- agLo09_1
    agAge_1[ageUse[ageUse$AgeClass=="OldGrowth",]] <- 1
    agAge_1[ageUse[ageUse$AgeClass=="Secondary",]] <- 2
    
    # Make soil parent material raster
    agParent_1 <- agLo09_1
    agParent_1[soil[soil$SoilParent=="CaimitoVolcanic",]] <- 1
    agParent_1[soil[soil$SoilParent=="Andesite",]] <- 2
    agParent_1[soil[soil$SoilParent=="Bohio",]] <- 3
    agParent_1[soil[soil$SoilParent=="CaimitoMarineSedimentary",]] <- 4
    
    # Make soil form raster
    agForm_1 <- agLo09_1
    agForm_1[soil[soil$SoilForm=="RedLightClay",]] <- 1
    agForm_1[soil[soil$SoilForm=="PaleSwellingClay",]] <- 2
    agForm_1[soil[soil$SoilForm=="BrownFineLoam",]] <- 3
    agForm_1[soil[soil$SoilForm=="MottledHeavyClay",]] <- 4
    
    # Combine all into data frame
    bestAgResults_1 <- data.frame(agLo09 = raster::values(agLo09_1),
                                agChm09 = raster::values(agChm09_1),
                                agFix = raster::values(agFix_1),
                                agObs = raster::values(agObs_1),
                                agAge = raster::values(agAge_1),
                                agParent = raster::values(agParent_1),
                                agForm = raster::values(agForm_1),
                                agN = raster::values(agN_1))
    # Only keep rows with enough observations
    bestAgResults_1 <- bestAgResults_1[useObs_1,]
    # Rename age, parent material, soil form codes
    bestAgResults_1[bestAgResults_1$agAge==1,"agAge"] <- "OldGrowth"
    bestAgResults_1[bestAgResults_1$agAge==2,"agAge"] <- "Secondary"
    bestAgResults_1[bestAgResults_1$agParent==1,"agParent"] <- "CaimitoVolcanic"
    bestAgResults_1[bestAgResults_1$agParent==2,"agParent"] <- "Andesite"
    bestAgResults_1[bestAgResults_1$agParent==3,"agParent"] <- "Bohio"
    bestAgResults_1[bestAgResults_1$agParent==4,"agParent"] <- "CaimitoMarineSedimentary"
    bestAgResults_1[bestAgResults_1$agForm==1,"agForm"] <- "RedLightClay"
    bestAgResults_1[bestAgResults_1$agForm==2,"agForm"] <- "PaleSwellingClay"
    bestAgResults_1[bestAgResults_1$agForm==3,"agForm"] <- "BrownFineLoam"
    bestAgResults_1[bestAgResults_1$agForm==4,"agForm"] <- "MottledHeavyClay"
    
    # Scale 2
    agScale_2 <- 8
    
    # Aggregate low canopy area raster
    agLo09_2 <- raster::aggregate(rasterLo09, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate mean canopy height raster
    agChm09_2 <- raster::aggregate(rasterchm09, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate average fixed effects raster
    agFix_2 <- raster::aggregate(avgPredictedRaster, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate average observed disturance raster
    agObs_2 <- raster::aggregate(avgObservedRaster, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate sampling effort
    agN_2 <- raster::aggregate(avgRasterN, fact = agScale_2, fun = sum, na.rm=T)
    # Find observations with 90% of values present
    minObs_2 <- 0.75*loResults$agBy[agScale_2]^2
    useObs_2 <- raster::values(agN_2)>minObs_2
    
    # Make forest age raster
    agAge_2 <- agLo09_2
    agAge_2[ageUse[ageUse$AgeClass=="OldGrowth",]] <- 1
    agAge_2[ageUse[ageUse$AgeClass=="Secondary",]] <- 2
    
    # Make soil parent material raster
    agParent_2 <- agLo09_2
    agParent_2[soil[soil$SoilParent=="CaimitoVolcanic",]] <- 1
    agParent_2[soil[soil$SoilParent=="Andesite",]] <- 2
    agParent_2[soil[soil$SoilParent=="Bohio",]] <- 3
    agParent_2[soil[soil$SoilParent=="CaimitoMarineSedimentary",]] <- 4
    
    # Make soil form raster
    agForm_2 <- agLo09_2
    agForm_2[soil[soil$SoilForm=="RedLightClay",]] <- 1
    agForm_2[soil[soil$SoilForm=="PaleSwellingClay",]] <- 2
    agForm_2[soil[soil$SoilForm=="BrownFineLoam",]] <- 3
    agForm_2[soil[soil$SoilForm=="MottledHeavyClay",]] <- 4
    
    # Combine all into data frame
    bestAgResults_2 <- data.frame(agLo09 = raster::values(agLo09_2),
                                  agChm09 = raster::values(agChm09_2),
                                  agFix = raster::values(agFix_2),
                                  agObs = raster::values(agObs_2),
                                  agAge = raster::values(agAge_2),
                                  agParent = raster::values(agParent_2),
                                  agForm = raster::values(agForm_2),
                                  agN = raster::values(agN_2))
    # Only keep rows with enough observations
    bestAgResults_2 <- bestAgResults_2[useObs_2,]
    # Rename age, parent material, soil form codes
    bestAgResults_2[bestAgResults_2$agAge==1,"agAge"] <- "OldGrowth"
    bestAgResults_2[bestAgResults_2$agAge==2,"agAge"] <- "Secondary"
    bestAgResults_2[bestAgResults_2$agParent==1,"agParent"] <- "CaimitoVolcanic"
    bestAgResults_2[bestAgResults_2$agParent==2,"agParent"] <- "Andesite"
    bestAgResults_2[bestAgResults_2$agParent==3,"agParent"] <- "Bohio"
    bestAgResults_2[bestAgResults_2$agParent==4,"agParent"] <- "CaimitoMarineSedimentary"
    bestAgResults_2[bestAgResults_2$agForm==1,"agForm"] <- "RedLightClay"
    bestAgResults_2[bestAgResults_2$agForm==2,"agForm"] <- "PaleSwellingClay"
    bestAgResults_2[bestAgResults_2$agForm==3,"agForm"] <- "BrownFineLoam"
    bestAgResults_2[bestAgResults_2$agForm==4,"agForm"] <- "MottledHeavyClay"
    
    # Scale 3
    agScale_3 <- 18
    
    # Aggregate low canopy area raster
    agLo09_3 <- raster::aggregate(rasterLo09, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate mean canopy height raster
    agChm09_3 <- raster::aggregate(rasterchm09, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate average fixed effects raster
    agFix_3 <- raster::aggregate(avgPredictedRaster, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate average observed disturance raster
    agObs_3 <- raster::aggregate(avgObservedRaster, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate sampling effort
    agN_3 <- raster::aggregate(avgRasterN, fact = agScale_3, fun = sum, na.rm=T)
    # Find observations with 90% of values present
    minObs_3 <- 0.75*loResults$agBy[agScale_3]^2
    useObs_3 <- raster::values(agN_3)>minObs_3
    
    # Make forest age raster
    agAge_3 <- agLo09_3
    agAge_3[ageUse[ageUse$AgeClass=="OldGrowth",]] <- 1
    agAge_3[ageUse[ageUse$AgeClass=="Secondary",]] <- 2
    
    # Make soil parent material raster
    agParent_3 <- agLo09_3
    agParent_3[soil[soil$SoilParent=="CaimitoVolcanic",]] <- 1
    agParent_3[soil[soil$SoilParent=="Andesite",]] <- 2
    agParent_3[soil[soil$SoilParent=="Bohio",]] <- 3
    agParent_3[soil[soil$SoilParent=="CaimitoMarineSedimentary",]] <- 4
    
    # Make soil form raster
    agForm_3 <- agLo09_3
    agForm_3[soil[soil$SoilForm=="RedLightClay",]] <- 1
    agForm_3[soil[soil$SoilForm=="PaleSwellingClay",]] <- 2
    agForm_3[soil[soil$SoilForm=="BrownFineLoam",]] <- 3
    agForm_3[soil[soil$SoilForm=="MottledHeavyClay",]] <- 4
    
    # Combine all into data frame
    bestAgResults_3 <- data.frame(agLo09 = raster::values(agLo09_3),
                                  agChm09 = raster::values(agChm09_3),
                                  agFix = raster::values(agFix_3),
                                  agObs = raster::values(agObs_3),
                                  agAge = raster::values(agAge_3),
                                  agParent = raster::values(agParent_3),
                                  agForm = raster::values(agForm_3),
                                  agN = raster::values(agN_3))
    # Only keep rows with enough observations
    bestAgResults_3 <- bestAgResults_3[useObs_3,]
    # Rename age, parent material, soil form codes
    bestAgResults_3[bestAgResults_3$agAge==1,"agAge"] <- "OldGrowth"
    bestAgResults_3[bestAgResults_3$agAge==2,"agAge"] <- "Secondary"
    bestAgResults_3[bestAgResults_3$agParent==1,"agParent"] <- "CaimitoVolcanic"
    bestAgResults_3[bestAgResults_3$agParent==2,"agParent"] <- "Andesite"
    bestAgResults_3[bestAgResults_3$agParent==3,"agParent"] <- "Bohio"
    bestAgResults_3[bestAgResults_3$agParent==4,"agParent"] <- "CaimitoMarineSedimentary"
    bestAgResults_3[bestAgResults_3$agForm==1,"agForm"] <- "RedLightClay"
    bestAgResults_3[bestAgResults_3$agForm==2,"agForm"] <- "PaleSwellingClay"
    bestAgResults_3[bestAgResults_3$agForm==3,"agForm"] <- "BrownFineLoam"
    bestAgResults_3[bestAgResults_3$agForm==4,"agForm"] <- "MottledHeavyClay"
    
    
  # Plot 1: Low canopy with predicted values from fixed effects
    PtCex_1 <- 0.8
    PtCex_2 <- 1
    PtCex_3 <- 2
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "4.00 ha resolution",
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFix")]),
           y=max(bestAgResults_1[,c("agLo09")]),
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
  
    plot(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agFix")]),
         ylim=range(bestAgResults_2[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "10.24 ha resolution",
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "51.84 ha resolution",
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # By soil parent material
    plot(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFix")]),
           y=max(bestAgResults_1[,c("agLo09")]),
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agFix")]),
         ylim=range(bestAgResults_2[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(agLo09~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFix")]),
           y=max(bestAgResults_1[,c("agLo09")]),
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agFix")]),
         ylim=range(bestAgResults_2[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(agLo09~agFix, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(agLo09~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Predicted disturbance frequency (fixed effects only)",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Low canopy area in 2009", side=2, outer=T, line=1.5)
    par(las=1)
    
    # Plot 2: Low canopy vs observed values
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "2.56 ha resolution",
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agObs")]),
           y=max(bestAgResults_1[,c("agLo09")]),
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agObs")]),
         ylim=range(bestAgResults_2[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "10.24 ha resolution",
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "64.00 ha resolution",
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # By soil parent material
    plot(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agObs")]),
           y=max(bestAgResults_1[,c("agLo09")]),
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agObs")]),
         ylim=range(bestAgResults_2[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(agLo09~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agObs")]),
           y=max(bestAgResults_1[,c("agLo09")]),
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agObs")]),
         ylim=range(bestAgResults_2[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(agLo09~agObs, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agLo09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(agLo09~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Average observed disturbance frequency",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Low canopy area in 2009", side=2, outer=T, line=1.5)
    par(las=1)
    
    # Plot 3: Mean canopy height with predicted values from fixed effects
    PtCex_1 <- 0.8
    PtCex_2 <- 1
    PtCex_3 <- 2
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "2.56 ha resolution",
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFix")]),
           y=max(bestAgResults_1[,c("agChm09")]),
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agFix")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "10.24 ha resolution",
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "64.00 ha resolution",
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # By soil parent material
    plot(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFix")]),
           y=max(bestAgResults_1[,c("agChm09")]),
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agFix")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFix")]),
           y=max(bestAgResults_1[,c("agChm09")]),
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agFix")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(agChm09~agFix, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Predicted disturbance frequency (fixed effects only)",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Mean canopy height in 2009", side=2, outer=T, line=1.5)
    par(las=1)
    
    # Plot 3: Mean canopy height vs observed values
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "2.56 ha resolution",
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agObs")]),
           y=max(bestAgResults_1[,c("agChm09")]),
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agObs")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "10.24 ha resolution",
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = "64.00 ha resolution",
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # By soil parent material
    plot(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agObs")]),
           y=max(bestAgResults_1[,c("agChm09")]),
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agObs")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agObs")]),
           y=max(bestAgResults_1[,c("agChm09")]),
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agObs")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(agChm09~agObs, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Average observed disturbance frequency",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Mean canopy height in 2009", side=2, outer=T, line=1.5)
    par(las=1)   
    
#### Figure 6: relationship between disturbance and standing canopy height at 2 spatial scales ####
    
    bestAgResults_1$agLo09Pct <- 100*bestAgResults_1$agLo09
    bestAgResults_2$agLo09Pct <- 100*bestAgResults_2$agLo09
    bestAgResults_3$agLo09Pct <- 100*bestAgResults_3$agLo09
    
    ## Be sure to run section above first!
    axisCex <- 1.2
    
    par(mfrow=c(4,2), mar=c(2,3,1,1), oma=c(4,4,1,1))
    
    plot(agLo09Pct~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agLo09Pct")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("a", x=0.0075, y = 25, cex = axisCex)
    mtext("Low canopy area (%)", side=2, outer=F, line=4, las=0, at = -0.1)
    
    points(agLo09Pct~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    
    
    plot(agLo09Pct~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agLo09Pct")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         yaxt="n",
         pch=19,
         cex=PtCex_1,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("b", x=0.0145, y = 25, cex = axisCex)
    
    points(agLo09Pct~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.06,
           y=28,
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agLo09Pct~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agLo09Pct")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("c", x=0.0082, y = 9.5, cex = axisCex)
    points(agLo09Pct~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    plot(agLo09Pct~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agLo09Pct")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         yaxt="n",
         pch=19,
         cex=PtCex_3,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("d", x=0.0145, y = 9.5, cex = axisCex)
    points(agLo09Pct~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    plot(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agFix")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("e", x=0.0075, y = 30, cex = axisCex)
    mtext("Mean canopy height (m)", side=2, outer=F, line=4, las=0, at=10)
    points(agChm09~agFix, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)

    
    plot(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agObs")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         yaxt="n",
         pch=19,
         cex=PtCex_1,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("f", x=0.01, y = 30, cex = axisCex)
    
    points(agChm09~agObs, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    
    plot(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agFix")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("g", x=0.0082, y = 26.5, cex = axisCex)
    points(agChm09~agFix, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    mtext("Predicted disturbance rate",
          side=1, outer=F, line=2.5)  
    
    plot(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agObs")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         yaxt="n",
         pch=19,
         cex=PtCex_3,
         cex.axis = axisCex,
         xlab = NA, ylab=NA)
    text("h", x=0.0145, y = 26.5, cex = axisCex)
    
    points(agChm09~agObs, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Observed disturbance rate",
          side=1, outer=F, line=2.5)  
    

#### Figure 6cd: Comparison of 2009 lidar data and average spatial pattern ####
    
    # Make alternate versions of all raster layers--only old growth forests
    rasterLo09b <- raster::mask(rasterLo09, ageUse[ageUse$AgeClass=="OldGrowth",])
    rasterchm09b <- raster::mask(rasterchm09, ageUse[ageUse$AgeClass=="OldGrowth",])
    avgPredictedRasterb <- raster::mask(avgPredictedRaster, ageUse[ageUse$AgeClass=="OldGrowth",])
    avgObservedRasterb <- raster::mask(avgObservedRaster, ageUse[ageUse$AgeClass=="OldGrowth",])
    avgRasterNb <- raster::mask(avgRasterN, ageUse[ageUse$AgeClass=="OldGrowth",])
    
    
    # Create data frame to store results  
    loResults_b <- data.frame(agBy = 1:20,
                            fixRlo = NA,
                            fixRlo_lo = NA,
                            fixRlo_hi = NA,
                            obsRlo = NA,
                            obsRlo_lo = NA,
                            obsRlo_hi = NA,
                            fixRchm = NA,
                            fixRchm_lo = NA,
                            fixRchm_hi = NA,
                            obsRchm = NA,
                            obsRchm_lo = NA,
                            obsRchm_hi = NA,
                            N = NA)
    
    for(i in 1:nrow(loResults_b)){
      
      # Get a matrix of observed cells
      minObs <- 0.75*loResults_b$agBy[i]^2
      
      # Aggregate rasters
      # Aggregate low canopy area raster
      agLo09 <- raster::aggregate(rasterLo09b, fact = loResults_b$agBy[i], fun = mean, na.rm=T)
      # Aggregate mean canopy height raster
      agChm09 <- raster::aggregate(rasterchm09b, fact = loResults_b$agBy[i], fun = mean, na.rm=T)
      # Aggregate average fixed effects raster
      agFix <- raster::aggregate(avgPredictedRasterb, fact = loResults_b$agBy[i], fun = mean, na.rm=T)
      # Aggregate average observed disturbance raster
      agObs <- raster::aggregate(avgObservedRasterb, fact = loResults_b$agBy[i], fun = mean, na.rm=T)
      # Aggregate sampling effort
      agN <- raster::aggregate(avgRasterNb, fact = loResults_b$agBy[i], fun = sum, na.rm=T)
      
      # Get matrix of values for each aggregated raster
      agLo09_mat <- raster::values(agLo09, format = "matrix")
      agChm09_mat <- raster::values(agChm09, format = "matrix")
      agFix_mat <- raster::values(agFix, format = "matrix")
      agObs_mat <- raster::values(agObs, format = "matrix")
      nMat <- raster::values(agN, format = "matrix")
      
      # matrix of cells that need to be combined
      combineMat <- which(nMat<minObs & nMat>0, arr.ind = T)
      
      while(nrow(combineMat)>0){
        
        # get surrounding cells for each matrix
        rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
        rows <- rows[which(rows>0 & rows <= nrow(nMat))]
        cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
        cols <- cols[which(cols>0 & cols <= ncol(nMat))]
        
        agLo09_mat_j <- agLo09_mat[rows,cols]
        agChm09_mat_j <- agChm09_mat[rows,cols]
        agFix_mat_j <- agFix_mat[rows,cols]
        agObs_mat_j <- agObs_mat[rows,cols]
        nMat_j <- nMat[rows,cols]
        
        # set cell of interest to 0
        nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- 0
        
        # Only proceed if there are non-NA neighboring cells
        if(length(c(nMat_j)[c(nMat_j)>0])>0){  
          
          # find neighboring cell with most nearby observations
          new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
          
          # replace value for new cell with a weighted mean of all other values
          agLo09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agLo09_mat_j[new_j[1],new_j[2]],agLo09_mat[combineMat[1,1],combineMat[1,2]]),
                                                           w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agChm09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agChm09_mat_j[new_j[1],new_j[2]],agChm09_mat[combineMat[1,1],combineMat[1,2]]),
                                                            w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agFix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agFix_mat_j[new_j[1],new_j[2]],agFix_mat[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agObs_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs_mat_j[new_j[1],new_j[2]],agObs_mat[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          
          nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
          
          # replace neighborhood in original matrices
          agLo09_mat[rows,cols] <- agLo09_mat_j
          agChm09_mat[rows,cols] <- agChm09_mat_j
          agFix_mat[rows,cols] <- agFix_mat_j
          agObs_mat[rows,cols] <- agObs_mat_j
          nMat[rows,cols] <- nMat_j
        }
        
        # Replace original observation with NaN
        agLo09_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agChm09_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agFix_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agObs_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        nMat[combineMat[1,1],combineMat[1,2]] <- 0
        
        # Redefine combineMat
        combineMat <- which(nMat<minObs & nMat >0, arr.ind = T)
      }
      
      valsKeep <- which(c(nMat)>minObs)
      loResults_b$fixRlo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$estimate
      loResults_b$fixRlo_lo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[1]
      loResults_b$fixRlo_hi[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[2]
      loResults_b$obsRlo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$estimate
      loResults_b$obsRlo_lo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[1]
      loResults_b$obsRlo_hi[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[2]
      
      loResults_b$fixRchm[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$estimate
      loResults_b$fixRchm_lo[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat[valsKeep]))$conf.int[1]
      loResults_b$fixRchm_hi[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[2]
      loResults_b$obsRchm[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$estimate
      loResults_b$obsRchm_lo[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[1]
      loResults_b$obsRchm_hi[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[2]
      
    }
    
    # Make alternate versions of all raster layers--only secondary forests
    rasterLo09c <- raster::mask(rasterLo09, ageUse[ageUse$AgeClass=="Secondary",])
    rasterchm09c <- raster::mask(rasterchm09, ageUse[ageUse$AgeClass=="Secondary",])
    avgPredictedRasterc <- raster::mask(avgPredictedRaster, ageUse[ageUse$AgeClass=="Secondary",])
    avgObservedRasterc <- raster::mask(avgObservedRaster, ageUse[ageUse$AgeClass=="Secondary",])
    avgRasterNc <- raster::mask(avgRasterN, ageUse[ageUse$AgeClass=="Secondary",])
    
    
    # Create data frame to store results  
    loResults_c <- data.frame(agBy = 1:20,
                              fixRlo = NA,
                              fixRlo_lo = NA,
                              fixRlo_hi = NA,
                              obsRlo = NA,
                              obsRlo_lo = NA,
                              obsRlo_hi = NA,
                              fixRchm = NA,
                              fixRchm_lo = NA,
                              fixRchm_hi = NA,
                              obsRchm = NA,
                              obsRchm_lo = NA,
                              obsRchm_hi = NA,
                              N = NA)
    
    for(i in 1:nrow(loResults_c)){
      
      # Get a matrix of observed cells
      minObs <- 0.75*loResults_c$agBy[i]^2
      
      # Aggregate rasters
      # Aggregate low canopy area raster
      agLo09 <- raster::aggregate(rasterLo09c, fact = loResults_c$agBy[i], fun = mean, na.rm=T)
      # Aggregate mean canopy height raster
      agChm09 <- raster::aggregate(rasterchm09c, fact = loResults_c$agBy[i], fun = mean, na.rm=T)
      # Aggregate average fixed effects raster
      agFix <- raster::aggregate(avgPredictedRasterc, fact = loResults_c$agBy[i], fun = mean, na.rm=T)
      # Aggregate average observed disturbance raster
      agObs <- raster::aggregate(avgObservedRasterc, fact = loResults_c$agBy[i], fun = mean, na.rm=T)
      # Aggregate sampling effort
      agN <- raster::aggregate(avgRasterNc, fact = loResults_c$agBy[i], fun = sum, na.rm=T)
      
      # Get matrix of values for each aggregated raster
      agLo09_mat <- raster::values(agLo09, format = "matrix")
      agChm09_mat <- raster::values(agChm09, format = "matrix")
      agFix_mat <- raster::values(agFix, format = "matrix")
      agObs_mat <- raster::values(agObs, format = "matrix")
      nMat <- raster::values(agN, format = "matrix")
      
      # matrix of cells that need to be combined
      combineMat <- which(nMat<minObs & nMat>0, arr.ind = T)
      
      while(nrow(combineMat)>0){
        
        # get surrounding cells for each matrix
        rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
        rows <- rows[which(rows>0 & rows <= nrow(nMat))]
        cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
        cols <- cols[which(cols>0 & cols <= ncol(nMat))]
        
        agLo09_mat_j <- agLo09_mat[rows,cols]
        agChm09_mat_j <- agChm09_mat[rows,cols]
        agFix_mat_j <- agFix_mat[rows,cols]
        agObs_mat_j <- agObs_mat[rows,cols]
        nMat_j <- nMat[rows,cols]
        
        # set cell of interest to 0
        nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- 0
        
        # Only proceed if there are non-NA neighboring cells
        if(length(c(nMat_j)[c(nMat_j)>0])>0){  
          
          # find neighboring cell with most nearby observations
          new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
          
          # replace value for new cell with a weighted mean of all other values
          agLo09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agLo09_mat_j[new_j[1],new_j[2]],agLo09_mat[combineMat[1,1],combineMat[1,2]]),
                                                           w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agChm09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agChm09_mat_j[new_j[1],new_j[2]],agChm09_mat[combineMat[1,1],combineMat[1,2]]),
                                                            w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agFix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agFix_mat_j[new_j[1],new_j[2]],agFix_mat[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          agObs_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs_mat_j[new_j[1],new_j[2]],agObs_mat[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
          
          nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
          
          # replace neighborhood in original matrices
          agLo09_mat[rows,cols] <- agLo09_mat_j
          agChm09_mat[rows,cols] <- agChm09_mat_j
          agFix_mat[rows,cols] <- agFix_mat_j
          agObs_mat[rows,cols] <- agObs_mat_j
          nMat[rows,cols] <- nMat_j
        }
        
        # Replace original observation with NaN
        agLo09_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agChm09_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agFix_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        agObs_mat[combineMat[1,1],combineMat[1,2]] <- NaN
        nMat[combineMat[1,1],combineMat[1,2]] <- 0
        
        # Redefine combineMat
        combineMat <- which(nMat<minObs & nMat >0, arr.ind = T)
      }
      
      valsKeep <- which(c(nMat)>minObs)
      loResults_c$fixRlo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$estimate
      loResults_c$fixRlo_lo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[1]
      loResults_c$fixRlo_hi[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[2]
      loResults_c$obsRlo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$estimate
      loResults_c$obsRlo_lo[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[1]
      loResults_c$obsRlo_hi[i] <- cor.test(x = c(agLo09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[2]
      
      loResults_c$fixRchm[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$estimate
      loResults_c$fixRchm_lo[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat[valsKeep]))$conf.int[1]
      loResults_c$fixRchm_hi[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agFix_mat)[valsKeep])$conf.int[2]
      loResults_c$obsRchm[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$estimate
      loResults_c$obsRchm_lo[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[1]
      loResults_c$obsRchm_hi[i] <- cor.test(x = c(agChm09_mat)[valsKeep], y = c(agObs_mat)[valsKeep])$conf.int[2]
    }
    
    
#### Figure 6: plot all ####
    
    # Plot results
    
    # Significance with aggregation scale    
    axisCex <- 1.2
    par(mfrow=c(3,2), mar=c(1,1,0,0), oma=c(5,4,3,1), las=1)
    
    xVals <- ((loResults$agBy*40)^2)/10000
    
    plot(x = xVals,
         y = loResults$fixRlo,
         ylim=c(-1,1.02),
         xlab = NA,
         ylab = NA,
         log = "x",
         type = "l",
         col = "black",
         cex.axis = axisCex,
         xaxt="n",
         lwd=2)
    par(las=0)
    text("All forest", x=0.15, y=-0.8, cex=1.5, adj=0)
    mtext("Proportion of low canopy area", side=3, outer=F)
    mtext(expression("Pearson correlation (r)"),
          side=2, outer=T, line=2)
    par(las=1)
    text("a",
         x = 0.15,
         y = 1, adj=0, cex=1.5)
    abline(h=0,lty=2)
    lines(x = xVals,
          y = loResults$obsRlo,
          col="grey",
          lty=2,
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$fixRlo_lo,rev(loResults$fixRlo_hi)),
            col=adjustcolor("black",0.25),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$obsRlo_lo,rev(loResults$obsRlo_hi)),
            col=adjustcolor("grey",0.35),
            border=NA)
    legend(x = 0.15,
           y = 0,
           c("Predicted frequency (fixed effects)",
             "Average observed frequency"),
           col=c("black","grey"),
           cex = 1.2,
           lty=c(1,2),
           lwd=2,
           bty="n")
    
    plot(x = xVals,
         y = loResults$fixRchm,
         log = "x",
         yaxt="n",
         ylim=c(-1,1.02),
         xlab = NA,
         ylab = NA,
         type = "l",
         cex.axis = axisCex,
         col = "black",
         xaxt="n",
         lwd=2)
    abline(h=0,lty=2)
    mtext("Mean canopy height", side=3, outer=F)
    
    lines(x = xVals,
          y = loResults$obsRchm,
          col="grey",
          lty=2,
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$fixRchm_lo,rev(loResults$fixRchm_hi)),
            col=adjustcolor("black",0.25),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults$obsRchm_lo,rev(loResults$obsRchm_hi)),
            col=adjustcolor("grey",0.35),
            border=NA)
    text("b",
         x = 0.15,
         y = 1, adj=0, cex=1.5)
    
    
    plot(x = xVals,
         y = loResults_b$fixRlo,
         cex.axis = axisCex,
         ylim=c(-1,1.02),
         log = "x",
         xlab = NA,
         ylab = NA,
         type = "l",
         xaxt="n",
         lty=1,
         col = colOld,
         lwd=2)
    par(las=0)
    text("Old growth forest", x=0.15, y=-0.8, cex=1.5, adj=0, col=colOld)
    par(las=1)
    text("c",
         x = 0.15,
         y = 1, adj=0, cex=1.5)
    abline(h=0,lty=2)
    lines(x = xVals,
          y = loResults_b$obsRlo,
          lty=2,
          col="grey",
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_b$fixRlo_lo,rev(loResults_b$fixRlo_hi)),
            col=adjustcolor(colOld,0.35),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_b$obsRlo_lo,rev(loResults_b$obsRlo_hi)),
            col=adjustcolor("grey",0.35),
            border=NA)

    
    plot(x = xVals,
         y = loResults_b$fixRchm,
         log = "x",
         xaxt="n",
         yaxt="n",
         ylim=c(-1,1.02),
         lty=1,
         xlab = NA,
         ylab = NA,
         type = "l",
         col = colOld,
         lwd=2)
    par(las=0)
    abline(h=0,lty=2)
    
    par(las=1)
    lines(x = xVals,
          y = loResults_b$obsRchm,
          col="grey",
          lty=2,
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_b$fixRchm_lo,rev(loResults_b$fixRchm_hi)),
            col=adjustcolor(colOld,0.35),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_b$obsRchm_lo,rev(loResults_b$obsRchm_hi)),
            col=adjustcolor("grey",0.35),
            border=NA)
    
    text("d",
         x = 0.15,
         y = 1, adj=0, cex=1.5)
    
    
    plot(x = xVals,
         y = loResults_c$fixRlo,
         cex.axis = axisCex,
         log = "x",
         ylim=c(-1,1.02),
         xlab = NA,
         ylab = NA,
         type = "l",
         lty=1,
         col = colSec,
         lwd=2)
    par(las=0)
    text("Secondary forest", x=0.15, y=-0.8, cex=1.5, adj=0, col=colSec)
    
    par(las=1)
    text("e",
         x = 0.15,
         y = 1, adj=0, cex=1.5)
    abline(h=0,lty=2)
    lines(x = xVals,
          y = loResults_c$obsRlo,
          lty=2,
          col="grey",
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_c$fixRlo_lo,rev(loResults_c$fixRlo_hi)),
            col=adjustcolor(colSec, 0.35),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_c$obsRlo_lo,rev(loResults_c$obsRlo_hi)),
            col=adjustcolor("grey", 0.35),
            border=NA)
    
    
    plot(x = xVals,
         y = loResults_c$fixRchm,
         cex.axis = axisCex,
         log = "x",
         yaxt="n",
         ylim=c(-1,1.02),
         lty=1,
         xlab = NA,
         ylab = NA,
         type = "l",
         col = colSec,
         lwd=2)
    par(las=0)
    abline(h=0,lty=2)
    
    par(las=1)
    lines(x = xVals,
          y = loResults_c$obsRchm,
          col="grey",
          lty=2,
          lwd=2)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_c$fixRchm_lo,rev(loResults_c$fixRchm_hi)),
            col=adjustcolor(colSec, 0.35),
            border=NA)
    polygon(x = c(xVals,rev(xVals)),
            y = c(loResults_c$obsRchm_lo,rev(loResults_c$obsRchm_hi)),
            col=adjustcolor("grey", 0.35),
            border=NA)
    
    text("f",
         x = 0.15,
         y = 1, adj=0, cex=1.5)
    
    mtext("Spatial resolution (ha)", side=1, outer=T, line=2)
    
    
    
    
#### RESIVION: Fixed effects from model with initial height ####
    
    # Original model
    fixedResults <- model_full$summary.fixed
    
    # Initial height models
    fixedResults_Ht <- model_full_ht$summary.fixed
    fixedResults_logHt <- model_full_htlog$summary.fixed
    
    fixedNames <- c("Curvature (linear)", "Slope (linear)", "Slope (quadratic)",                    
                    "Height above drainage (linear)","Height above drainage (quadratic)",
                    "Soil parent: Andesite", "Soil parent: Caimito marine", "Soil parent: Caimito volcanic",
                    "Soil form: Heavy clay", "Soil form: Swelling clay", "Soil form: Light clay",
                    "Forest age: Secondary",
                    "Year: 2015-2018", "Initial canopy height")
    
    par(mfrow=c(1,1), mar=c(4,1,1,1))
    plot(x = fixedResults[2:nrow(fixedResults),"mean"],
         y = 1.2*nrow(fixedResults_Ht):3,
         xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(-0.4,0.6),
         ylim=c(2,1.2*nrow(fixedResults_Ht)),
         pch = 19, 
         cex = 1,
         xlab = "Fixed effect",
         ylab = NA, yaxt = "n")
    arrows(x0 = fixedResults[2:nrow(fixedResults),"0.025quant"],
           x1 = fixedResults[2:nrow(fixedResults),"0.975quant"],
           y0 = 1.2*nrow(fixedResults_Ht):3,
           y1 = 1.2*nrow(fixedResults_Ht):3,
           angle = 90, code=3,
           length = 0.05)
    
    points(x = fixedResults_Ht[2:nrow(fixedResults_Ht),"mean"],
           y = 1.2*nrow(fixedResults_Ht):2 - 0.2,
           pch = 19, 
           cex = 1,
           col = "grey")
    arrows(x0 = fixedResults_Ht[2:nrow(fixedResults_Ht),"0.025quant"],
           x1 = fixedResults_Ht[2:nrow(fixedResults_Ht),"0.975quant"],
           y0 = 1.2*nrow(fixedResults_Ht):2- 0.2,
           y1 = 1.2*nrow(fixedResults_Ht):2- 0.2,
           angle = 90, code=3,
           length = 0.05,
           col = "grey")
    
    points(x = fixedResults_logHt[2:nrow(fixedResults_logHt),"mean"],
           y = 1.2*nrow(fixedResults_logHt):2 - 0.4,
           pch = 19, 
           cex = 1,
           col = "red")
    arrows(x0 = fixedResults_logHt[2:nrow(fixedResults_logHt),"0.025quant"],
           x1 = fixedResults_logHt[2:nrow(fixedResults_logHt),"0.975quant"],
           y0 = 1.2*nrow(fixedResults_logHt):2- 0.4,
           y1 = 1.2*nrow(fixedResults_logHt):2- 0.4,
           angle = 90, code=3,
           length = 0.05,
           col = "red")
    
    
    abline(v=0, lty=2)
    text(fixedNames,
         x = fixedResults_Ht[2:nrow(fixedResults_Ht),"0.975quant"]+0.05,
         y = 1.2*nrow(fixedResults_Ht):2,
         pos = 4)
    
    legend(x = -1.2,
           y = 1.2*nrow(fixedResults_Ht),
           c("Original model",
             "Canopy height",
             "Log canopy height"),
           bty="n",
           col = c("black", "grey", "red"),
           pch=c(19,19,19,1))
    
#### REVISION: Calculate predicted values and residuals from model with initial height ####
    
    ## Fixed and random effects  
    
    # Get fitted values and reorder to match original data frame
    bci.gapsAll$pred_Ht <- model_full_htlog$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
    bci.gapsAll$resd_Ht <- bci.gapsAll$gapPropCens - bci.gapsAll$predHt
    
    ## Fixed effects only  
    
    # Intercept
    bci.gapsAll$fix_int_Ht <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="(Intercept)"]
    
    # Curvature
    bci.gapsAll$fix_C_Ht <- bci.gapsAll$Sc_curvMean_2*(model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Sc_curvMean_2"])
    
    # Slope
    bci.gapsAll$fix_S_Ht <- bci.gapsAll$Sc_slopeMean_16*(model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Sc_slopeMean_16"])
    bci.gapsAll$fix_S2_Ht <- bci.gapsAll$Sc_slopeMean_16_Sq*(model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Sc_slopeMean_16_Sq"])
    
    # Height above drainage
    bci.gapsAll$fix_H_Ht <- bci.gapsAll$Sc_drainMean*(model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Sc_drainMean"])
    bci.gapsAll$fix_H2_Ht <- bci.gapsAll$Sc_drainMean_Sq*(model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Sc_drainMean_Sq"])
    
    # Soil form
    bci.gapsAll$fix_soilForm_Ht <- 0
    
    bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="MottledHeavyClay","fix_soilForm_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="soilFormMottledHeavyClay"]
    bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="PaleSwellingClay","fix_soilForm_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="soilFormPaleSwellingClay"]
    bci.gapsAll[!is.na(bci.gapsAll$soilForm) & bci.gapsAll$soilForm=="RedLightClay","fix_soilForm_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="soilFormRedLightClay"]
    
    
    # Soil parent material
    bci.gapsAll$fix_soilParent_Ht <- 0
    
    bci.gapsAll[!is.na(bci.gapsAll$soilParent) & bci.gapsAll$soilParent=="Andesite","fix_soilParent_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="soilParentAndesite"]
    bci.gapsAll[!is.na(bci.gapsAll$soilParent) & bci.gapsAll$soilParent=="CaimitoMarineSedimentary","fix_soilParent_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="soilParentCaimitoMarineSedimentary"]
    bci.gapsAll[!is.na(bci.gapsAll$soilParent) & bci.gapsAll$soilParent=="CaimitoVolcanic","fix_soilParent_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="soilParentCaimitoVolcanic"]
    
    # Year 
    bci.gapsAll$fix_year_Ht <- 0
    bci.gapsAll[bci.gapsAll$Year=="2018","fix_year_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Year2018"]
    
    # Age 
    bci.gapsAll$fix_age_Ht <- 0
    bci.gapsAll[bci.gapsAll$age=="Secondary" & !(is.na(bci.gapsAll$age)),"fix_age_Ht"] <- model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="ageSecondary"]
    
    # Initial canopy height
    bci.gapsAll$fix_Ht <- bci.gapsAll$Sc_initialHtlog*(model_full_htlog$summary.fixed$mean[model_full_htlog$names.fixed=="Sc_initialHtlog"])
    
    # All
    bci.gapsAll$fix_sum_Ht <- bci.gapsAll$fix_int_Ht + bci.gapsAll$fix_C_Ht + bci.gapsAll$fix_S_Ht + bci.gapsAll$fix_S2_Ht + bci.gapsAll$fix_H_Ht + bci.gapsAll$fix_H2_Ht + bci.gapsAll$fix_soilParent_Ht + bci.gapsAll$fix_soilForm_Ht + bci.gapsAll$fix_age_Ht + bci.gapsAll$fix_year_Ht + bci.gapsAll$fix_Ht
    
    # Predicted value with just fixed effects
    bci.gapsAll$fix_pred_Ht <- exp(bci.gapsAll$fix_sum_Ht)/(1 + exp(bci.gapsAll$fix_sum_Ht))
    
    bci.gapsAll$fix_resd_Ht <- bci.gapsAll$gapPropCens - bci.gapsAll$fix_pred_Ht
    
    
   
#### REVISION: Look at model residuals from model with initial height ####
    
    pdf(file = "residualsPlots_logHeight.pdf", onefile=T, height=9, width=6.5)
    
    # CURVATURE 
    
    par(mfrow=c(4,1), las=1, mar=c(3,3,1,1), oma=c(2,2,1,1))
    
    # Both years, fixed and random effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    curvSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd_Ht"], df=30)
    
    plot(resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, fixed + random effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = curvSpline_Ht$x, y= curvSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i],"resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i],"resd_Ht"]
      curvRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
   
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=curvRMSE,
           pch=1,
           col = "purple")
    
    # Both years, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    curvSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd_Ht"], df=30)
    
    plot(fix_resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, only fixed effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = curvSpline_Ht$x, y= curvSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i],"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i],"fix_resd_Ht"]
      curvRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=curvRMSE,
           pch=1,
           col = "purple")
    
    # 2015-2018, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    curvSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd_Ht"], df=30)
    
    plot(fix_resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2015-2018, only fixed effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = curvSpline_Ht$x, y= curvSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd_Ht"]
      curvRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=curvRMSE,
           pch=1,
           col = "purple")
    
    # 2018-2020, only fixed effects
    curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    curvSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"curvMean_2"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd_Ht"], df=30)
    
    plot(fix_resd~curvMean_2, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2018-2020, only fixed effects")
    abline(h=0,col="red")
    lines(x = curvSpline$x, y= curvSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = curvSpline_Ht$x, y= curvSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                           max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
    for(i in 1:nrow(curvRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$curvMean_2 > curvRMSE$min[i]
                           & bci.gapsAll$curvMean_2 <= curvRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd_Ht"]
      curvRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=curvRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=curvRMSE,
           pch=1,
           col = "purple")
    
    mtext("Curvature (LaPlacian convexity)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
    # SLOPE 
    
    # Both years, fixed and random effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    slopeSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd_Ht"], df=30)
    
    plot(resd~slopeMean_16, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, fixed + random effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = slopeSpline_Ht$x, y= slopeSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i],"resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i],"resd_Ht"]
      slopeRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=slopeRMSE,
           pch=1,
           col = "purple")
    
    # Both years, only fixed effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    slopeSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd_Ht"], df=30)
    
    plot(fix_resd~slopeMean_16, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "All years, only fixed effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = slopeSpline_Ht$x, y= slopeSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i],"fix_resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i],"fix_resd_Ht"]
      slopeRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=slopeRMSE,
           pch=1,
           col = "purple")
    
    # 2015-2018, only fixed effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    slopeSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd_Ht"], df=30)
    
    plot(fix_resd~slopeMean_16, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2015-2018, only fixed effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = slopeSpline_Ht$x, y= slopeSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd_Ht"]
      slopeRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=slopeRMSE,
           pch=1,
           col = "purple")
    
    # 2018-2020, only fixed effects
    slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    slopeSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"slopeMean_16"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd_Ht"], df=30)
    
    plot(fix_resd~slopeMean_16, data = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,],
         ylim=c(-0.02,0.05),
         pch=20,
         col = adjustcolor("black", 0.05),
         cex = 0.6,
         ylab = NA,
         xlab = NA,
         main = "2018-2020, only fixed effects")
    abline(h=0,col="red")
    lines(x = slopeSpline$x, y= slopeSpline$y, col=adjustcolor("orange",0.7),lty=1,lwd=2)
    lines(x = slopeSpline_Ht$x, y= slopeSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                           max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA,
                           RMSE_Ht = NA)
    slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
    for(i in 1:nrow(slopeRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$slopeMean_16 > slopeRMSE$min[i]
                           & bci.gapsAll$slopeMean_16 <= slopeRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd_Ht"]
      slopeRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=slopeRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=slopeRMSE,
           pch=1,
           col = "purple")
    
    mtext("Slope (degrees)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
    # HEIGHT ABOVE DRAINAGE 
    
    # Both years, fixed and random effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    drainSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd_Ht"], df=30)
    
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
    lines(x = drainSpline_Ht$x, y= drainSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                            max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i],"resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i],"resd_Ht"]
      drainRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=drainRMSE,
           pch=1,
           col = "purple")
    
    # Both years, only fixed effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    drainSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd_Ht"], df=30)
    
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
    lines(x = drainSpline_Ht$x, y= drainSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                            max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i],"fix_resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i],"fix_resd_Ht"]
      drainRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=drainRMSE,
           pch=1,
           col = "purple")
    
    # 2015-2018, only fixed effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    drainSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd_Ht"], df=30)
    
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
    lines(x = drainSpline_Ht$x, y= drainSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                            max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd_Ht"]
      drainRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=drainRMSE,
           pch=1,
           col = "purple")
    
    # 2018-2020, only fixed effects
    drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    drainSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"drainMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd_Ht"], df=30)
    
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
    lines(x = drainSpline_Ht$x, y= drainSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                            max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
    for(i in 1:nrow(drainRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$drainMean > drainRMSE$min[i]
                           & bci.gapsAll$drainMean <= drainRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd_Ht"]
      drainRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=drainRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=drainRMSE,
           pch=1,
           col = "purple")
    
    mtext("HAND (m)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
    # ASPECT
    
    # Both years, fixed and random effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd"], df=30)
    aspectSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"resd_Ht"], df=30)
    
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
    lines(x = aspectSpline_Ht$x, y= aspectSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"resd_Ht"]
      aspectRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=aspectRMSE,
           pch=1,
           col = "purple")
    
    # Both years, only fixed effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd"], df=30)
    aspectSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"fix_resd_Ht"], df=30)
    
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
    lines(x = aspectSpline_Ht$x, y= aspectSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"fix_resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"fix_resd_Ht"]
      aspectRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=aspectRMSE,
           pch=1,
           col = "purple")
    
    # 2015-2018, only fixed effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd"], df=30)
    aspectSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2018,"fix_resd_Ht"], df=30)
    
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
    lines(x = aspectSpline_Ht$x, y= aspectSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i]
                           & bci.gapsAll$Year==2018,"fix_resd_Ht"]
      aspectRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=aspectRMSE,
           pch=1,
           col = "purple")
    
    # 2018-2020, only fixed effects
    aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd"], df=30)
    aspectSpline_Ht <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"aspectMean"], y = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) & bci.gapsAll$Year==2020,"fix_resd_Ht"], df=30)
    
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
    lines(x = aspectSpline_Ht$x, y= aspectSpline_Ht$y, col=adjustcolor("purple",0.7),lty=3,lwd=2)
    
    aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                            max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                            mid = NA,
                            RMSE = NA,
                            RMSE_Ht = NA)
    aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
    for(i in 1:nrow(aspectRMSE)){
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd"]
      aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
      
      res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                           & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                           & bci.gapsAll$aspectMean <= aspectRMSE$max[i]
                           & bci.gapsAll$Year==2020,"fix_resd_Ht"]
      aspectRMSE$RMSE_Ht[i] <- sqrt(mean(res_i^2))
    }
    
    points(RMSE~mid, data=aspectRMSE,
           pch=19,
           col = "orange")
    points(RMSE_Ht~mid, data=aspectRMSE,
           pch=1,
           col = "purple")
    
    mtext("Aspect (degrees)", side=1, outer=T)
    par(las=0)
    mtext("Residual value", side=2, outer=T)
    par(las=1)
    
    dev.off() 
