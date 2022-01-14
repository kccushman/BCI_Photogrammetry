#### Load data and results ####
library(INLA)
load("Data_INLA/INLA_40m.RData")
load("Data_INLA/INLA_fullModelResult.RData")
load("Data_INLA/INLA_fullModelResult_separate.RData")
load("Data_INLA/INLA_fullModelResult_noLargeGaps.RData")
load("Data_INLA/INLA_fullModelResult_initialHt.RData")

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
  
  # load buffer for BCI
  buffer <- rgdal::readOGR("Data_Ancillary/BCI_Outline_Minus25/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
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
    
    # Predicted value with just fixed effects (adjusting for asymmetry of logit function)
    
      # Vector of random effects, in right order
      all_random <- model_full$summary.random[[1]]$mean[order(bci.gapsAll_Order$Order)]
      # Keep only non-NA values
      all_random[is.na(bci.gapsAll$age)] <- NA
      all_random <- all_random[!is.na(all_random)]
      # define inverse logit function
      invLogit <- function(x){
        return(exp(x)/(1+exp(x)))
      }
      
      bci.gapsAll$fix_pred <- NA
      for(i in 1:nrow(bci.gapsAll)){
        if(!is.na(bci.gapsAll$fix_sum[i])){
          bci.gapsAll$fix_pred[i] <- mean(invLogit(bci.gapsAll$fix_sum[i]+all_random))
        }
      }
  
    # bci.gapsAll$fix_pred <- exp(bci.gapsAll$fix_sum)/(1 + exp(bci.gapsAll$fix_sum))
    
    bci.gapsAll$fix_resd <- bci.gapsAll$gapPropCens - bci.gapsAll$fix_pred
  
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
  
#### Make plots of observed, average predicted, and SD predicted values ####

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

  ## Rasters of average spatial pattern across all years

    avgStackFix <- raster::stack(fixRaster18,fixRaster20)
    avgPredictedRaster <- raster::calc(avgStackFix, mean, na.rm=T)
    raster::writeRaster(avgPredictedRaster, file = "Data_INLA/avgPredictedFix.tif")
  
    avgStackObs <- raster::stack(obsRaster18,obsRaster20)
    avgObservedRaster <- raster::calc(avgStackObs, mean, na.rm=T)
    avgObservedRaster[avgObservedRaster==0] <- NA




#### Figure 4ab: aggregate and look at correlation with observed rates ####          
    
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
    
    # NOTE: for first iteration at original spatial grain
    # this will produce a warning saying "nothing to aggregate"
    
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
  cexLab <- 0.6
  cexAx <- 0.5
    
  pdf("Figure4ab.pdf", width=4.33, height=2)
  
    par(mfrow=c(1,2), mar=c(1,1,0,0), oma=c(2,2,0.5,1), las=1)
    
    xVals <- (agResults$agBy*40)^2/10000
    
    plot(x = xVals,
         y = agResults$fixR_18,
         log="x",
         ylim=c(0,1.02),
         xlab = NA,
         ylab = NA,
         type = "l",
         col = col18,
         lwd=2,
         cex.axis = cexAx)
    text("a. Fixed effects only", x = 0.15, y = 1.02, adj=0, cex=cexAx)
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
           lwd=2,
           cex=cexLab)
    
    plot(x = xVals,
         y = agResults$allR_18,
         ylim=c(0,1.02),
         log="x",
         yaxt="n",
         xlab = NA,
         ylab = NA,
         type = "l",
         col = col18,
         lwd=2,
         cex.axis = cexAx)
    text("b. Fixed + random effects", x = 0.15, y = 1.02, adj=0, cex=cexAx)
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
    
    
    mtext("Spatial grain (ha)", side=1, outer=T, line=0.5, cex=cexLab)
    mtext(bquote("Pearson correlation ("~italic(rho)~")"), side=2, outer=T, line=0.5, las=0, cex=cexLab)
    
  dev.off()  
    
#### Figure 4c: Fixed effects sizes ####  
  fixedResults <- model_full$summary.fixed
  
  # separated by years
  fixedResults1 <- model_full1$summary.fixed
  fixedResults2 <- model_full2$summary.fixed
  
  # main model and second interval without large gaps
  fixedResults_alt <- model_full_alt$summary.fixed
  fixedResults2_alt <- model_full2_alt$summary.fixed
  
  
  fixedNames <- c("Curvature (linear)", "Slope (linear)", "Slope (quadratic)",                    
                  "Height above drainage (linear)","Height above drainage (quadratic)",
                  "Soil parent: Andesite", "Soil parent: Caimito marine", "Soil parent: Caimito volcanic",
                  "Soil form: Mottled heavy clay", "Soil form: Pale swelling clay", "Soil form: Red light clay",
                  "Forest age: Secondary",
                  "Year: 2015-2018")
# Plot
  cexPt <- 0.4
  lengthArr <- 0.01
  
pdf("Figure4c.pdf", width=4.33, height=3)
  
  par(mfrow=c(1,1), mar=c(1,1,0,0), oma=c(2,2,0.5,1))
  
  plot(x = fixedResults[2:nrow(fixedResults),"mean"],
       y = 1.2*nrow(fixedResults):2,
       xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(-0.4,0.6),
       ylim=c(2,1.2*nrow(fixedResults)),
       xlab=NA,
       pch = 19, 
       cex = cexPt,
       cex.axis = cexAx,
       col = adjustcolor("black",0.99),
       ylab = NA, yaxt = "n")
  arrows(x0 = fixedResults[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2,
         y1 = 1.2*nrow(fixedResults):2,
         angle = 90, code=3,
         length = lengthArr)
  
  points(x = fixedResults1[2:nrow(fixedResults),"mean"],
       y = 1.2*nrow(fixedResults):2 - 0.2,
       pch = 19, 
       cex = cexPt,
       col = adjustcolor(col18,0.99))
  arrows(x0 = fixedResults1[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults1[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2- 0.2,
         y1 = 1.2*nrow(fixedResults):2- 0.2,
         angle = 90, code=3,
         length = lengthArr,
         col = col18)
  
  points(x = fixedResults2[2:nrow(fixedResults),"mean"],
         y = 1.2*nrow(fixedResults):2 - 0.4,
         pch = 19, 
         cex = cexPt,
         col = adjustcolor(col20,0.99))
  arrows(x0 = fixedResults2[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults2[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2- 0.4,
         y1 = 1.2*nrow(fixedResults):2- 0.4,
         angle = 90, code=3,
         length = lengthArr,
         col = col20)
  
  points(x = fixedResults2_alt[2:nrow(fixedResults),"mean"],
         y = 1.2*nrow(fixedResults):2 - 0.6,
         pch = 1, 
         cex = cexPt,
         col = adjustcolor(col20,0.99))
  arrows(x0 = fixedResults2_alt[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults2_alt[2:nrow(fixedResults),"0.975quant"],
         y0 = 1.2*nrow(fixedResults):2- 0.6,
         y1 = 1.2*nrow(fixedResults):2- 0.6,
         angle = 90, code=3,
         length = lengthArr,
         lty=2,
         col = col20)
  
  
  abline(v=0, lty=2)
  text(fixedNames,
      x = fixedResults[2:nrow(fixedResults),"0.975quant"]+0.03,
       y = 1.2*nrow(fixedResults):2,
       pos = 4, cex = cexAx)
  
  legend(x = -1.3,
         y = 5.5,
         c("Full model",
           "2015-2018 only",
           "2018-2020 only",
           "2018-2020 only (no blowdown)"),
         bty="n",
         col = adjustcolor(c("black", col18, col20, col20),0.99),
         pch=c(19,19,19,1),
         cex=cexLab)
  text("c", x = -1.25, y = 1.2*nrow(fixedResults), cex=cexAx)
  mtext("Fixed effect", side=1, outer=F, line=2, cex=cexLab)
dev.off()  

#### Define forest age and soil type polygons ####
    # Forest age polygon
    age <- rgdal::readOGR("Data_Ancillary/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    age$AgeClass <- "Other"
    age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
    age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    ageUse <- age[!(age$AgeClass=="Other"),]
    
    # Soil type polygon  
    soil <- rgdal::readOGR("Data_Ancillary/BCI_Soils/BCI_Soils.shp")
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
    
#### Figure 6: Comparison of 2009 lidar data and average spatial pattern (log-log scales) ####
    
    #Create rasters where all non-NA values are 1
    predMeanRaster18n <- predMeanRaster18
    predMeanRaster18n[!is.na(raster::values(predMeanRaster18n))] <- 1
    fixRaster18n <- fixRaster18
    fixRaster18n[!is.na(raster::values(fixRaster18n))] <- 1
    predMeanRaster20n <- predMeanRaster20
    predMeanRaster20n[!is.na(raster::values(predMeanRaster20n))] <- 1
    fixRaster20n <- fixRaster20
    fixRaster20n[!is.na(raster::values(fixRaster20n))] <- 1
    
    # Raster of low canopy area in 2009
    lo09 <- raster::raster("Data_HeightRasters/loCanopy2009.tif") # pixels with value 1 are => 10 m height, 0 are < 10 m
    # Raster of canopy height in 2009
    chm09 <- raster::raster("Data_HeightRasters/CHM_2009_QAQC.tif") 
    
    
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
      
      # Replace 0 values for log transformation
      agLo09_mat[agLo09_mat==0 & !is.na(agLo09_mat)] <- 1/3200
      
      loResults$fixRlo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$estimate
      loResults$fixRlo_lo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults$fixRlo_hi[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults$obsRlo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$estimate
      loResults$obsRlo_lo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults$obsRlo_hi[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[2]
      
      loResults$fixRchm[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$estimate
      loResults$fixRchm_lo[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults$fixRchm_hi[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults$obsRchm[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$estimate
      loResults$obsRchm_lo[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults$obsRchm_hi[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[2]
      
    }
    
    
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
      
      # Replace 0 values for log transformation
      agLo09_mat[agLo09_mat==0 & !is.na(agLo09_mat)] <- 1/3200
      
      loResults_b$fixRlo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$estimate
      loResults_b$fixRlo_lo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults_b$fixRlo_hi[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_b$obsRlo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$estimate
      loResults_b$obsRlo_lo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_b$obsRlo_hi[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[2]
      
      loResults_b$fixRchm[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$estimate
      loResults_b$fixRchm_lo[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat[valsKeep])))$conf.int[1]
      loResults_b$fixRchm_hi[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_b$obsRchm[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$estimate
      loResults_b$obsRchm_lo[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_b$obsRchm_hi[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[2]
      
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
      
      # Replace 0 values for log transformation
      agLo09_mat[agLo09_mat==0 & !is.na(agLo09_mat)] <- 1/3200
      
      valsKeep <- which(c(nMat)>minObs)
      loResults_c$fixRlo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$estimate
      loResults_c$fixRlo_lo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults_c$fixRlo_hi[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_c$obsRlo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$estimate
      loResults_c$obsRlo_lo[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_c$obsRlo_hi[i] <- cor.test(x = log(c(agLo09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[2]
      
      loResults_c$fixRchm[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$estimate
      loResults_c$fixRchm_lo[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat[valsKeep])))$conf.int[1]
      loResults_c$fixRchm_hi[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_c$obsRchm[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$estimate
      loResults_c$obsRchm_lo[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_c$obsRchm_hi[i] <- cor.test(x = log(c(agChm09_mat)[valsKeep]), y = log(c(agObs_mat)[valsKeep]))$conf.int[2]
    }
    
# Plot results
    
  # Significance with aggregation scale    
  axisCex <- 1.2
  labCex <- 0.8
  
  
  pdf("Figure_6.pdf", width=4.33, height = 5)
  
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
    text("All forest", x=0.15, y=-0.8, cex=axisCex, adj=0)
    mtext("Proportion of low canopy area", side=3, outer=F, cex = labCex)
    mtext(bquote("Pearson correlation ("~italic(rho)~")"), side=2, outer=T, line=2, cex = labCex)
    par(las=1)
    text("a",
         x = 0.15,
         y = 0.95, adj=0, cex=axisCex)
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
           c("Predicted rate (fixed effects)",
             "Observed rate"),
           col=c("black","grey"),
           cex = 0.9,
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
    mtext("Mean canopy height", side=3, outer=F, cex = labCex)
    
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
         y = 0.95, adj=0, cex=axisCex)
    
    
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
         y = 0.95, adj=0, cex=axisCex)
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
         y = 0.95, adj=0, cex=axisCex)
    
    
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
         y = 0.95, adj=0, cex=axisCex)
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
         y = 0.95, adj=0, cex=axisCex)
    
    mtext("Spatial resolution (ha)", side=1, outer=T, line=2, cex = labCex)
  dev.off()  
    
    
#### Figure S14: Look at alternate models soil, forest age, and topography separately ####
  
  # SOIL ONLY
  
  #1. Load these alternate models    
  load("Data_INLA/INLA_ModelResult_soilsOnly.RData")
  
  #2. First, calculate alternate predicted values with and without random effects
  
  # Duplicate 'bci.gapsAll' objects for alternate model
  bci.gapsSoil <- bci.gapsAll
  
  # Get fitted values and reorder to match original data frame
  bci.gapsSoil$pred <- model_soilsOnly$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  
  ## Fixed effects only
  
  # Intercept
  bci.gapsSoil$fix_int <- model_soilsOnly$summary.fixed$mean[model_soilsOnly$names.fixed=="(Intercept)"]
  
  # Soil form
  bci.gapsSoil$fix_soilForm <- 0
  
  bci.gapsSoil[!is.na(bci.gapsSoil$soilForm) & bci.gapsSoil$soilForm=="MottledHeavyClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormMottledHeavyClay"]
  bci.gapsSoil[!is.na(bci.gapsSoil$soilForm) & bci.gapsSoil$soilForm=="PaleSwellingClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormPaleSwellingClay"]
  bci.gapsSoil[!is.na(bci.gapsSoil$soilForm) & bci.gapsSoil$soilForm=="RedLightClay","fix_soilForm"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilFormRedLightClay"]
  
  
  # Soil parent material
  bci.gapsSoil$fix_soilParent <- 0
  
  bci.gapsSoil[!is.na(bci.gapsSoil$soilParent) & bci.gapsSoil$soilParent=="Andesite","fix_soilParent"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilParentAndesite"]
  bci.gapsSoil[!is.na(bci.gapsSoil$soilParent) & bci.gapsSoil$soilParent=="CaimitoMarineSedimentary","fix_soilParent"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilParentCaimitoMarineSedimentary"]
  bci.gapsSoil[!is.na(bci.gapsSoil$soilParent) & bci.gapsSoil$soilParent=="CaimitoVolcanic","fix_soilParent"] <- model_full$summary.fixed$mean[model_full$names.fixed=="soilParentCaimitoVolcanic"]
  
  # Year 
  bci.gapsSoil$fix_year <- 0
  bci.gapsSoil[bci.gapsSoil$Year=="2018","fix_year"] <- model_soilsOnly$summary.fixed$mean[model_soilsOnly$names.fixed=="Year2018"]
  
  # All
  bci.gapsSoil$fix_sum <- bci.gapsSoil$fix_int + bci.gapsSoil$fix_soilForm + bci.gapsSoil$fix_soilParent + bci.gapsSoil$fix_year
  
  # Predicted value with just fixed effects
  bci.gapsSoil$fix_pred <- exp(bci.gapsSoil$fix_sum)/(1 + exp(bci.gapsSoil$fix_sum))
  
  # Make predicted raster (fixed and random effects)
  predMeanRaster18Soil <- raster::raster(x = matrix(data = bci.gapsSoil[bci.gapsSoil$Year=="2018","pred"],
                                                    nrow = nCellY,
                                                    ncol = nCellX,
                                                    byrow = F),
                                         xmn = raster::extent(bci.gaps18)@xmin,
                                         xmx = raster::extent(bci.gaps18)@xmax,
                                         ymn = raster::extent(bci.gaps18)@ymin,
                                         ymx = raster::extent(bci.gaps18)@ymax)
  predMeanRaster18Soil <- raster::mask(predMeanRaster18Soil, buffer)
  
  
  predMeanRaster20Soil <- raster::raster(x = matrix(data = bci.gapsSoil[bci.gapsSoil$Year=="2020","pred"],
                                                    nrow = nCellY,
                                                    ncol = nCellX,
                                                    byrow = F),
                                         xmn = raster::extent(bci.gaps20)@xmin,
                                         xmx = raster::extent(bci.gaps20)@xmax,
                                         ymn = raster::extent(bci.gaps20)@ymin,
                                         ymx = raster::extent(bci.gaps20)@ymax)
  predMeanRaster20Soil <- raster::mask(predMeanRaster20Soil, buffer)
  
  # Make predicted raster (fixed effects only)
  fixRaster18Soil <- raster::raster(x = matrix(data = bci.gapsSoil[bci.gapsSoil$Year=="2018","fix_pred"],
                                               nrow = nCellY,
                                               ncol = nCellX,
                                               byrow = F),
                                    xmn = raster::extent(bci.gaps18)@xmin,
                                    xmx = raster::extent(bci.gaps18)@xmax,
                                    ymn = raster::extent(bci.gaps18)@ymin,
                                    ymx = raster::extent(bci.gaps18)@ymax)
  fixRaster18Soil <- raster::mask(fixRaster18Soil, buffer)
  
  
  fixRaster20Soil <- raster::raster(x = matrix(data = bci.gapsSoil[bci.gapsSoil$Year=="2020","fix_pred"],
                                               nrow = nCellY,
                                               ncol = nCellX,
                                               byrow = F),
                                    xmn = raster::extent(bci.gaps20)@xmin,
                                    xmx = raster::extent(bci.gaps20)@xmax,
                                    ymn = raster::extent(bci.gaps20)@ymin,
                                    ymx = raster::extent(bci.gaps20)@ymax)
  fixRaster20Soil <- raster::mask(fixRaster20Soil, buffer)
  
  
  #3. Aggregate and look at correlations
  
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
  agResultsSoil <- data.frame(agBy = 1:20,
                              fixR_18 = NA,
                              fixR_18lo = NA,
                              fixR_18hi = NA,
                              fixR_20 = NA,
                              fixR_20lo = NA,
                              fixR_20hi = NA)
  
  for(i in 1:nrow(agResultsSoil)){
    
    # Get a matrix of observed cells
    minObs <- 0.75*agResultsSoil$agBy[i]^2
    
    agPred18_all <- raster::aggregate(predMeanRaster18Soil, fact = agResultsSoil$agBy[i], fun = mean, na.rm=T)
    agPred18_fix <- raster::aggregate(fixRaster18Soil, fact = agResultsSoil$agBy[i], fun = mean, na.rm=T)
    agObs18 <- raster::aggregate(obsRaster18, fact = agResultsSoil$agBy[i], fun = mean, na.rm=T)
    agPred18_all_n <- raster::aggregate(predMeanRaster18n, fact = agResultsSoil$agBy[i], fun = sum, na.rm=T)
    
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
    
    agResultsSoil$fixR_18[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$estimate
    agResultsSoil$fixR_18lo[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[1]
    agResultsSoil$fixR_18hi[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[2]
    
    agPred20_all <- raster::aggregate(predMeanRaster20Soil, fact = agResultsSoil$agBy[i], fun = mean, na.rm=T)
    agPred20_fix <- raster::aggregate(fixRaster20Soil, fact = agResultsSoil$agBy[i], fun = mean, na.rm=T)
    agObs20 <- raster::aggregate(obsRaster20, fact = agResultsSoil$agBy[i], fun = mean, na.rm=T)
    agPred20_all_n <- raster::aggregate(predMeanRaster20n, fact = agResultsSoil$agBy[i], fun = sum, na.rm=T)
    
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
    
    agResultsSoil$fixR_20[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$estimate
    agResultsSoil$fixR_20lo[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[1]
    agResultsSoil$fixR_20hi[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[2]
    
  }
  
  # TOPOGRAPHY ONLY
  
  #1. Load these alternate models    
  load("Data_INLA/INLA_ModelResult_topoOnly.RData")
  
  #2. First, calculate alternate predicted values with and without random effects
  
  # Duplicate 'bci.gapsAll' objects for alternate model
  bci.gapsTopo <- bci.gapsAll
  
  # Get fitted values and reorder to match original data frame
  bci.gapsTopo$pred <- model_topoOnly$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  
  ## Fixed effects only
  
  # Intercept
  bci.gapsTopo$fix_int <- model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="(Intercept)"]
  
  # Curvature
  bci.gapsTopo$fix_C <- bci.gapsTopo$Sc_curvMean_2*(model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="Sc_curvMean_2"])
  
  # Slope
  bci.gapsTopo$fix_S <- bci.gapsTopo$Sc_slopeMean_16*(model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="Sc_slopeMean_16"])
  bci.gapsTopo$fix_S2 <- bci.gapsTopo$Sc_slopeMean_16_Sq*(model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="Sc_slopeMean_16_Sq"])
  
  # Height above drainage
  bci.gapsTopo$fix_H <- bci.gapsTopo$Sc_drainMean*(model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="Sc_drainMean"])
  bci.gapsTopo$fix_H2 <- bci.gapsTopo$Sc_drainMean_Sq*(model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="Sc_drainMean_Sq"])
  
  # Year 
  bci.gapsTopo$fix_year <- 0
  bci.gapsTopo[bci.gapsTopo$Year=="2018","fix_year"] <- model_topoOnly$summary.fixed$mean[model_topoOnly$names.fixed=="Year2018"]
  
  # All
  bci.gapsTopo$fix_sum <- bci.gapsTopo$fix_int + bci.gapsTopo$fix_C + bci.gapsTopo$fix_S + bci.gapsTopo$fix_S2 + bci.gapsTopo$fix_H + bci.gapsTopo$fix_H2 + bci.gapsTopo$fix_year
  
  # Predicted value with just fixed effects
  # Vector of random effects, in right order
  all_random <- model_topoOnly$summary.random[[1]]$mean[order(bci.gapsAll_Order$Order)]
  # Keep only non-NA values
  all_random[is.na(bci.gapsAll$age)] <- NA
  all_random <- all_random[!is.na(all_random)]
  
  bci.gapsTopo$fix_pred <- NA
  for(i in 1:nrow(bci.gapsAll)){
    if(!is.na(bci.gapsTopo$fix_sum[i])){
      bci.gapsTopo$fix_pred[i] <- mean(exp(bci.gapsTopo$fix_sum[i] + all_random)/(1 + exp(bci.gapsTopo$fix_sum[i] + all_random)))
    }
  }
  
  # Make predicted raster (fixed and random effects)
  predMeanRaster18Topo <- raster::raster(x = matrix(data = bci.gapsTopo[bci.gapsTopo$Year=="2018","pred"],
                                                    nrow = nCellY,
                                                    ncol = nCellX,
                                                    byrow = F),
                                         xmn = raster::extent(bci.gaps18)@xmin,
                                         xmx = raster::extent(bci.gaps18)@xmax,
                                         ymn = raster::extent(bci.gaps18)@ymin,
                                         ymx = raster::extent(bci.gaps18)@ymax)
  predMeanRaster18Topo <- raster::mask(predMeanRaster18Topo, buffer)
  
  
  predMeanRaster20Topo <- raster::raster(x = matrix(data = bci.gapsTopo[bci.gapsTopo$Year=="2020","pred"],
                                                    nrow = nCellY,
                                                    ncol = nCellX,
                                                    byrow = F),
                                         xmn = raster::extent(bci.gaps20)@xmin,
                                         xmx = raster::extent(bci.gaps20)@xmax,
                                         ymn = raster::extent(bci.gaps20)@ymin,
                                         ymx = raster::extent(bci.gaps20)@ymax)
  predMeanRaster20Topo <- raster::mask(predMeanRaster20Topo, buffer)
  
  # Make predicted raster (fixed effects only)
  fixRaster18Topo <- raster::raster(x = matrix(data = bci.gapsTopo[bci.gapsTopo$Year=="2018","fix_pred"],
                                               nrow = nCellY,
                                               ncol = nCellX,
                                               byrow = F),
                                    xmn = raster::extent(bci.gaps18)@xmin,
                                    xmx = raster::extent(bci.gaps18)@xmax,
                                    ymn = raster::extent(bci.gaps18)@ymin,
                                    ymx = raster::extent(bci.gaps18)@ymax)
  fixRaster18Topo <- raster::mask(fixRaster18Topo, buffer)
  
  
  fixRaster20Topo <- raster::raster(x = matrix(data = bci.gapsTopo[bci.gapsTopo$Year=="2020","fix_pred"],
                                               nrow = nCellY,
                                               ncol = nCellX,
                                               byrow = F),
                                    xmn = raster::extent(bci.gaps20)@xmin,
                                    xmx = raster::extent(bci.gaps20)@xmax,
                                    ymn = raster::extent(bci.gaps20)@ymin,
                                    ymx = raster::extent(bci.gaps20)@ymax)
  fixRaster20Topo <- raster::mask(fixRaster20Topo, buffer)
  
  
  #3. Aggregate and look at correlations
  
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
  agResultsTopo <- data.frame(agBy = 1:20,
                              fixR_18 = NA,
                              fixR_18lo = NA,
                              fixR_18hi = NA,
                              fixR_20 = NA,
                              fixR_20lo = NA,
                              fixR_20hi = NA)
  
  for(i in 1:nrow(agResultsTopo)){
    
    # Get a matrix of observed cells
    minObs <- 0.75*agResultsTopo$agBy[i]^2
    
    agPred18_all <- raster::aggregate(predMeanRaster18Topo, fact = agResultsTopo$agBy[i], fun = mean, na.rm=T)
    agPred18_fix <- raster::aggregate(fixRaster18Topo, fact = agResultsTopo$agBy[i], fun = mean, na.rm=T)
    agObs18 <- raster::aggregate(obsRaster18, fact = agResultsTopo$agBy[i], fun = mean, na.rm=T)
    agPred18_all_n <- raster::aggregate(predMeanRaster18n, fact = agResultsTopo$agBy[i], fun = sum, na.rm=T)
    
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
    
    agResultsTopo$fixR_18[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$estimate
    agResultsTopo$fixR_18lo[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[1]
    agResultsTopo$fixR_18hi[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[2]
    
    agPred20_all <- raster::aggregate(predMeanRaster20Topo, fact = agResultsTopo$agBy[i], fun = mean, na.rm=T)
    agPred20_fix <- raster::aggregate(fixRaster20Topo, fact = agResultsTopo$agBy[i], fun = mean, na.rm=T)
    agObs20 <- raster::aggregate(obsRaster20, fact = agResultsTopo$agBy[i], fun = mean, na.rm=T)
    agPred20_all_n <- raster::aggregate(predMeanRaster20n, fact = agResultsTopo$agBy[i], fun = sum, na.rm=T)
    
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
    
    agResultsTopo$fixR_20[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$estimate
    agResultsTopo$fixR_20lo[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[1]
    agResultsTopo$fixR_20hi[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[2]
  }        
  
  # AGE ONLY
  
  #1. Load these alternate models    
  load("Data_INLA/INLA_ModelResult_ageOnly.RData")
  
  #2. First, calculate alternate predicted values with and without random effects
  
  # Duplicate 'bci.gapsAll' objects for alternate model
  bci.gapsAge <- bci.gapsAll
  
  # Get fitted values and reorder to match original data frame
  bci.gapsAge$pred <- model_ageOnly$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  
  ## Fixed effects only
  
  # Intercept
  bci.gapsAge$fix_int <- model_ageOnly$summary.fixed$mean[model_ageOnly$names.fixed=="(Intercept)"]
  
  
  # Year 
  bci.gapsAge$fix_year <- 0
  bci.gapsAge[bci.gapsAge$Year=="2018","fix_year"] <- model_ageOnly$summary.fixed$mean[model_ageOnly$names.fixed=="Year2018"]
  
  # Age 
  bci.gapsAge$fix_age <- 0
  bci.gapsAge[bci.gapsAge$age=="Secondary" & !(is.na(bci.gapsAge$age)),"fix_age"] <- model_ageOnly$summary.fixed$mean[model_ageOnly$names.fixed=="ageSecondary"]
  
  # All
  bci.gapsAge$fix_sum <- bci.gapsAge$fix_int +bci.gapsAge$fix_age + bci.gapsAge$fix_year
  
  # Predicted value with just fixed effects
  # Vector of random effects, in right order
  all_random <- model_ageOnly$summary.random[[1]]$mean[order(bci.gapsAll_Order$Order)]
  # Keep only non-NA values
  all_random[is.na(bci.gapsAll$age)] <- NA
  all_random <- all_random[!is.na(all_random)]
  
  bci.gapsAge$fix_pred <- NA
  for(i in 1:nrow(bci.gapsAll)){
    if(!is.na(bci.gapsAge$fix_sum[i])){
      bci.gapsAge$fix_pred[i] <- mean(exp(bci.gapsAge$fix_sum[i] + all_random)/(1 + exp(bci.gapsAge$fix_sum[i] + all_random)))
    }
  }
  
  # Make predicted raster (fixed and random effects)
  predMeanRaster18Age <- raster::raster(x = matrix(data = bci.gapsAge[bci.gapsAge$Year=="2018","pred"],
                                                   nrow = nCellY,
                                                   ncol = nCellX,
                                                   byrow = F),
                                        xmn = raster::extent(bci.gaps18)@xmin,
                                        xmx = raster::extent(bci.gaps18)@xmax,
                                        ymn = raster::extent(bci.gaps18)@ymin,
                                        ymx = raster::extent(bci.gaps18)@ymax)
  predMeanRaster18Age <- raster::mask(predMeanRaster18Age, buffer)
  
  
  predMeanRaster20Age <- raster::raster(x = matrix(data = bci.gapsAge[bci.gapsAge$Year=="2020","pred"],
                                                   nrow = nCellY,
                                                   ncol = nCellX,
                                                   byrow = F),
                                        xmn = raster::extent(bci.gaps20)@xmin,
                                        xmx = raster::extent(bci.gaps20)@xmax,
                                        ymn = raster::extent(bci.gaps20)@ymin,
                                        ymx = raster::extent(bci.gaps20)@ymax)
  predMeanRaster20Age <- raster::mask(predMeanRaster20Age, buffer)
  
  # Make predicted raster (fixed effects only)
  fixRaster18Age <- raster::raster(x = matrix(data = bci.gapsAge[bci.gapsAge$Year=="2018","fix_pred"],
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps18)@xmin,
                                   xmx = raster::extent(bci.gaps18)@xmax,
                                   ymn = raster::extent(bci.gaps18)@ymin,
                                   ymx = raster::extent(bci.gaps18)@ymax)
  fixRaster18Age <- raster::mask(fixRaster18Age, buffer)
  
  
  fixRaster20Age <- raster::raster(x = matrix(data = bci.gapsAge[bci.gapsAge$Year=="2020","fix_pred"],
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps20)@xmin,
                                   xmx = raster::extent(bci.gaps20)@xmax,
                                   ymn = raster::extent(bci.gaps20)@ymin,
                                   ymx = raster::extent(bci.gaps20)@ymax)
  fixRaster20Age <- raster::mask(fixRaster20Age, buffer)
  
  
  #3. Aggregate and look at correlations
  
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
  agResultsAge <- data.frame(agBy = 1:20,
                             fixR_18 = NA,
                             fixR_18lo = NA,
                             fixR_18hi = NA,
                             fixR_20 = NA,
                             fixR_20lo = NA,
                             fixR_20hi = NA)
  
  for(i in 1:nrow(agResultsAge)){
    
    # Get a matrix of observed cells
    minObs <- 0.75*agResultsAge$agBy[i]^2
    
    agPred18_all <- raster::aggregate(predMeanRaster18Age, fact = agResultsAge$agBy[i], fun = mean, na.rm=T)
    agPred18_fix <- raster::aggregate(fixRaster18Age, fact = agResultsAge$agBy[i], fun = mean, na.rm=T)
    agObs18 <- raster::aggregate(obsRaster18, fact = agResultsAge$agBy[i], fun = mean, na.rm=T)
    agPred18_all_n <- raster::aggregate(predMeanRaster18n, fact = agResultsAge$agBy[i], fun = sum, na.rm=T)
    
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
    
    agResultsAge$fixR_18[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$estimate
    agResultsAge$fixR_18lo[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[1]
    agResultsAge$fixR_18hi[i] <- cor.test(x = c(agObs18_mat), y = c(agPred18_fix_mat))$conf.int[2]
    
    agPred20_all <- raster::aggregate(predMeanRaster20Age, fact = agResultsAge$agBy[i], fun = mean, na.rm=T)
    agPred20_fix <- raster::aggregate(fixRaster20Age, fact = agResultsAge$agBy[i], fun = mean, na.rm=T)
    agObs20 <- raster::aggregate(obsRaster20, fact = agResultsAge$agBy[i], fun = mean, na.rm=T)
    agPred20_all_n <- raster::aggregate(predMeanRaster20n, fact = agResultsAge$agBy[i], fun = sum, na.rm=T)
    
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
    
    agResultsAge$fixR_20[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$estimate
    agResultsAge$fixR_20lo[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[1]
    agResultsAge$fixR_20hi[i] <- cor.test(x = c(agObs20_mat), y = c(agPred20_fix_mat))$conf.int[2]
  }        
  
  ## Plot results compared to original model
  
  # Take average of both time intervals
  agResults$fixR_mean <- 0.5*(agResults$fixR_18 + agResults$fixR_20)
  agResults$fixR_lo <- 0.5*(agResults$fixR_18lo + agResults$fixR_18lo)
  agResults$fixR_hi <- 0.5*(agResults$fixR_18hi + agResults$fixR_20hi)
  agResultsSoil$fixR_mean <- 0.5*(agResultsSoil$fixR_18 + agResultsSoil$fixR_20)
  agResultsSoil$fixR_lo <- 0.5*(agResultsSoil$fixR_18lo + agResultsSoil$fixR_18lo)
  agResultsSoil$fixR_hi <- 0.5*(agResultsSoil$fixR_18hi + agResultsSoil$fixR_20hi)
  agResultsTopo$fixR_mean <- 0.5*(agResultsTopo$fixR_18 + agResultsTopo$fixR_20)
  agResultsTopo$fixR_lo <- 0.5*(agResultsTopo$fixR_18lo + agResultsTopo$fixR_18lo)
  agResultsTopo$fixR_hi <- 0.5*(agResultsTopo$fixR_18hi + agResultsTopo$fixR_20hi)
  agResultsAge$fixR_mean <- 0.5*(agResultsAge$fixR_18 + agResultsAge$fixR_20)
  agResultsAge$fixR_lo <- 0.5*(agResultsAge$fixR_18lo + agResultsAge$fixR_18lo)
  agResultsAge$fixR_hi <- 0.5*(agResultsAge$fixR_18hi + agResultsAge$fixR_20hi)
  
  
  par(mfrow=c(1,1), mar=c(4,1,1,0), oma=c(1,4,1,1), las=1)
  
  xVals <- (agResults$agBy*40)^2/10000
  
  plot(x = xVals,
       y = agResults$fixR_mean,
       log="x",
       ylim=c(0,1.02),
       xlab = NA,
       ylab = NA,
       type = "l",
       col = "black",
       lwd=2)
  lines(x = xVals,
        y = agResultsSoil$fixR_mean,
        col="brown",
        lwd=2, lty=3)
  lines(x = xVals,
        y = agResultsTopo$fixR_mean,
        col="darkgreen",
        lwd=2, lty=3)
  lines(x = xVals,
        y = agResultsAge$fixR_mean,
        col="grey",
        lwd=2, lty=3)
  
  legend(c("Original model",
           "Soil only",
           "Forest age only",
           "Topography only"),
         x=0.15,y=1,
         lty=c(1,3,3,3),
         col=c("black","brown","grey","darkgreen"),
         bty="n",
         lwd=2)
  abline(h=0,lty=2)
  
  mtext("Spatial grain (ha)", side=1, outer=T, line=-2)
  mtext(bquote("Pearson correlation ("~italic(rho)~")"), side=2, outer=T, line=1.5, las=0)
  
  
  mean(agResults$fixR_mean^2-agResultsSoil$fixR_mean^2)
  mean(agResults$fixR_mean^2-agResultsAge$fixR_mean^2)
  mean(agResults$fixR_mean^2-agResultsTopo$fixR_mean^2)
  
#### Figures S15-18: Plot correlations at 3 different spatial grains ####
    
    # Get forest age and soil type values for each pixel
    
    # Scale 1
    agScale_1 <- 3
    
    # Get a matrix of observed cells
    minObs <- 0.75*agScale_1^2
    
    # Aggregate rasters
      # Aggregate low canopy area raster
      agLo09 <- raster::aggregate(rasterLo09, fact = agScale_1, fun = mean, na.rm=T)
      # Aggregate mean canopy height raster
      agChm09 <- raster::aggregate(rasterchm09, fact = agScale_1, fun = mean, na.rm=T)
      # Aggregate average fixed effects raster
      agFix <- raster::aggregate(avgPredictedRaster, fact = agScale_1, fun = mean, na.rm=T)
      # Aggregate average observed disturbance raster
      agObs <- raster::aggregate(avgObservedRaster, fact = agScale_1, fun = mean, na.rm=T)
      # Aggregate sampling effort
      agN <- raster::aggregate(avgRasterN, fact = agScale_1, fun = sum, na.rm=T)
      
      # Get matrix of values for each aggregated raster
      agLo09_mat_1 <- raster::values(agLo09, format = "matrix")
      agChm09_mat_1 <- raster::values(agChm09, format = "matrix")
      agFix_mat_1 <- raster::values(agFix, format = "matrix")
      agObs_mat_1 <- raster::values(agObs, format = "matrix")
      nMat <- raster::values(agN, format = "matrix")
      
    # matrix of cells that need to be combined
    combineMat <- which(nMat<minObs & nMat>0, arr.ind = T)
    
    while(nrow(combineMat)>0){
      
      # get surrounding cells for each matrix
      rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
      rows <- rows[which(rows>0 & rows <= nrow(nMat))]
      cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
      cols <- cols[which(cols>0 & cols <= ncol(nMat))]
      
      agLo09_mat_j <- agLo09_mat_1[rows,cols]
      agChm09_mat_j <- agChm09_mat_1[rows,cols]
      agFix_mat_j <- agFix_mat_1[rows,cols]
      agObs_mat_j <- agObs_mat_1[rows,cols]
      nMat_j <- nMat[rows,cols]
      
      # set cell of interest to 0
      nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- 0
      
      # Only proceed if there are non-NA neighboring cells
      if(length(c(nMat_j)[c(nMat_j)>0])>0){  
        
        # find neighboring cell with most nearby observations
        new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
        
        # replace value for new cell with a weighted mean of all other values
        agLo09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agLo09_mat_j[new_j[1],new_j[2]],agLo09_mat_1[combineMat[1,1],combineMat[1,2]]),
                                                         w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agChm09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agChm09_mat_j[new_j[1],new_j[2]],agChm09_mat_1[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agFix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agFix_mat_j[new_j[1],new_j[2]],agFix_mat_1[combineMat[1,1],combineMat[1,2]]),
                                                        w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agObs_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs_mat_j[new_j[1],new_j[2]],agObs_mat_1[combineMat[1,1],combineMat[1,2]]),
                                                        w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        
        nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
        
        # replace neighborhood in original matrices
        agLo09_mat_1[rows,cols] <- agLo09_mat_j
        agChm09_mat_1[rows,cols] <- agChm09_mat_j
        agFix_mat_1[rows,cols] <- agFix_mat_j
        agObs_mat_1[rows,cols] <- agObs_mat_j
        nMat[rows,cols] <- nMat_j
      }
      
      # Replace original observation with NaN
      agLo09_mat_1[combineMat[1,1],combineMat[1,2]] <- NaN
      agChm09_mat_1[combineMat[1,1],combineMat[1,2]] <- NaN
      agFix_mat_1[combineMat[1,1],combineMat[1,2]] <- NaN
      agObs_mat_1[combineMat[1,1],combineMat[1,2]] <- NaN
      nMat[combineMat[1,1],combineMat[1,2]] <- 0
      
      # Redefine combineMat
      combineMat <- which(nMat<minObs & nMat >0, arr.ind = T)
    }
    
    valsRemove <- which(c(nMat)<minObs)
    agLo09_mat_1[valsRemove] <- NA
    agChm09_mat_1[valsRemove] <- NA
    agFix_mat_1[valsRemove] <- NA
    agObs_mat_1[valsRemove] <- NA
    
    agLo09_1 <- raster::raster(agLo09_mat_1, template=agLo09)
    agChm09_1 <- raster::raster(agChm09_mat_1, template=agChm09)
    agFix_1 <- raster::raster(agFix_mat_1, template=agFix)
    agObs_1 <- raster::raster(agObs_mat_1, template=agObs)
    
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
                                  agForm = raster::values(agForm_1))
    # Rename age, parent material, soil form codes
    bestAgResults_1[bestAgResults_1$agAge==1 & !is.na(bestAgResults_1$agAge),"agAge"] <- "OldGrowth"
    bestAgResults_1[bestAgResults_1$agAge==2 & !is.na(bestAgResults_1$agAge),"agAge"] <- "Secondary"
    bestAgResults_1[bestAgResults_1$agParent==1 & !is.na(bestAgResults_1$agParent),"agParent"] <- "CaimitoVolcanic"
    bestAgResults_1[bestAgResults_1$agParent==2 & !is.na(bestAgResults_1$agParent),"agParent"] <- "Andesite"
    bestAgResults_1[bestAgResults_1$agParent==3 & !is.na(bestAgResults_1$agParent),"agParent"] <- "Bohio"
    bestAgResults_1[bestAgResults_1$agParent==4 & !is.na(bestAgResults_1$agParent),"agParent"] <- "CaimitoMarineSedimentary"
    bestAgResults_1[bestAgResults_1$agForm==1 & !is.na(bestAgResults_1$agForm),"agForm"] <- "RedLightClay"
    bestAgResults_1[bestAgResults_1$agForm==2 & !is.na(bestAgResults_1$agForm),"agForm"] <- "PaleSwellingClay"
    bestAgResults_1[bestAgResults_1$agForm==3 & !is.na(bestAgResults_1$agForm),"agForm"] <- "BrownFineLoam"
    bestAgResults_1[bestAgResults_1$agForm==4 & !is.na(bestAgResults_1$agForm),"agForm"] <- "MottledHeavyClay"
    bestAgResults_1 <- bestAgResults_1[!is.na(bestAgResults_1$agObs),]
      
    # Scale 2
    agScale_2 <- 8
    
    # Get a matrix of observed cells
    minObs <- 0.75*agScale_2^2
    
    # Aggregate rasters
    # Aggregate low canopy area raster
    agLo09 <- raster::aggregate(rasterLo09, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate mean canopy height raster
    agChm09 <- raster::aggregate(rasterchm09, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate average fixed effects raster
    agFix <- raster::aggregate(avgPredictedRaster, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate average observed disturbance raster
    agObs <- raster::aggregate(avgObservedRaster, fact = agScale_2, fun = mean, na.rm=T)
    # Aggregate sampling effort
    agN <- raster::aggregate(avgRasterN, fact = agScale_2, fun = sum, na.rm=T)
    
    # Get matrix of values for each aggregated raster
    agLo09_mat_2 <- raster::values(agLo09, format = "matrix")
    agChm09_mat_2 <- raster::values(agChm09, format = "matrix")
    agFix_mat_2 <- raster::values(agFix, format = "matrix")
    agObs_mat_2 <- raster::values(agObs, format = "matrix")
    nMat <- raster::values(agN, format = "matrix")
    
    # matrix of cells that need to be combined
    combineMat <- which(nMat<minObs & nMat>0, arr.ind = T)
    
    while(nrow(combineMat)>0){
      
      # get surrounding cells for each matrix
      rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
      rows <- rows[which(rows>0 & rows <= nrow(nMat))]
      cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
      cols <- cols[which(cols>0 & cols <= ncol(nMat))]
      
      agLo09_mat_j <- agLo09_mat_2[rows,cols]
      agChm09_mat_j <- agChm09_mat_2[rows,cols]
      agFix_mat_j <- agFix_mat_2[rows,cols]
      agObs_mat_j <- agObs_mat_2[rows,cols]
      nMat_j <- nMat[rows,cols]
      
      # set cell of interest to 0
      nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- 0
      
      # Only proceed if there are non-NA neighboring cells
      if(length(c(nMat_j)[c(nMat_j)>0])>0){  
        
        # find neighboring cell with most nearby observations
        new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
        
        # replace value for new cell with a weighted mean of all other values
        agLo09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agLo09_mat_j[new_j[1],new_j[2]],agLo09_mat_2[combineMat[1,1],combineMat[1,2]]),
                                                         w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agChm09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agChm09_mat_j[new_j[1],new_j[2]],agChm09_mat_2[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agFix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agFix_mat_j[new_j[1],new_j[2]],agFix_mat_2[combineMat[1,1],combineMat[1,2]]),
                                                        w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agObs_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs_mat_j[new_j[1],new_j[2]],agObs_mat_2[combineMat[1,1],combineMat[1,2]]),
                                                        w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        
        nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
        
        # replace neighborhood in original matrices
        agLo09_mat_2[rows,cols] <- agLo09_mat_j
        agChm09_mat_2[rows,cols] <- agChm09_mat_j
        agFix_mat_2[rows,cols] <- agFix_mat_j
        agObs_mat_2[rows,cols] <- agObs_mat_j
        nMat[rows,cols] <- nMat_j
      }
      
      # Replace original observation with NaN
      agLo09_mat_2[combineMat[1,1],combineMat[1,2]] <- NaN
      agChm09_mat_2[combineMat[1,1],combineMat[1,2]] <- NaN
      agFix_mat_2[combineMat[1,1],combineMat[1,2]] <- NaN
      agObs_mat_2[combineMat[1,1],combineMat[1,2]] <- NaN
      nMat[combineMat[1,1],combineMat[1,2]] <- 0
      
      # Redefine combineMat
      combineMat <- which(nMat<minObs & nMat >0, arr.ind = T)
    }
    
    valsRemove <- which(c(nMat)<minObs)
    agLo09_mat_2[valsRemove] <- NA
    agChm09_mat_2[valsRemove] <- NA
    agFix_mat_2[valsRemove] <- NA
    agObs_mat_2[valsRemove] <- NA
    
    agLo09_2 <- raster::raster(agLo09_mat_2, template=agLo09)
    agChm09_2 <- raster::raster(agChm09_mat_2, template=agChm09)
    agFix_2 <- raster::raster(agFix_mat_2, template=agFix)
    agObs_2 <- raster::raster(agObs_mat_2, template=agObs)
    
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
                                  agForm = raster::values(agForm_2))
    # Rename age, parent material, soil form codes
    bestAgResults_2[bestAgResults_2$agAge==1 & !is.na(bestAgResults_2$agAge),"agAge"] <- "OldGrowth"
    bestAgResults_2[bestAgResults_2$agAge==2 & !is.na(bestAgResults_2$agAge),"agAge"] <- "Secondary"
    bestAgResults_2[bestAgResults_2$agParent==1 & !is.na(bestAgResults_2$agParent),"agParent"] <- "CaimitoVolcanic"
    bestAgResults_2[bestAgResults_2$agParent==2 & !is.na(bestAgResults_2$agParent),"agParent"] <- "Andesite"
    bestAgResults_2[bestAgResults_2$agParent==3 & !is.na(bestAgResults_2$agParent),"agParent"] <- "Bohio"
    bestAgResults_2[bestAgResults_2$agParent==4 & !is.na(bestAgResults_2$agParent),"agParent"] <- "CaimitoMarineSedimentary"
    bestAgResults_2[bestAgResults_2$agForm==1 & !is.na(bestAgResults_2$agForm),"agForm"] <- "RedLightClay"
    bestAgResults_2[bestAgResults_2$agForm==2 & !is.na(bestAgResults_2$agForm),"agForm"] <- "PaleSwellingClay"
    bestAgResults_2[bestAgResults_2$agForm==3 & !is.na(bestAgResults_2$agForm),"agForm"] <- "BrownFineLoam"
    bestAgResults_2[bestAgResults_2$agForm==4 & !is.na(bestAgResults_2$agForm),"agForm"] <- "MottledHeavyClay"
    bestAgResults_2 <- bestAgResults_2[!is.na(bestAgResults_2$agObs),]
    
    
    # Scale 3
    agScale_3 <- 18
    
    # Get a matrix of observed cells
    minObs <- 0.75*agScale_3^2
    
    # Aggregate rasters
    # Aggregate low canopy area raster
    agLo09 <- raster::aggregate(rasterLo09, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate mean canopy height raster
    agChm09 <- raster::aggregate(rasterchm09, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate average fixed effects raster
    agFix <- raster::aggregate(avgPredictedRaster, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate average observed disturbance raster
    agObs <- raster::aggregate(avgObservedRaster, fact = agScale_3, fun = mean, na.rm=T)
    # Aggregate sampling effort
    agN <- raster::aggregate(avgRasterN, fact = agScale_3, fun = sum, na.rm=T)
    
    # Get matrix of values for each aggregated raster
    agLo09_mat_3 <- raster::values(agLo09, format = "matrix")
    agChm09_mat_3 <- raster::values(agChm09, format = "matrix")
    agFix_mat_3 <- raster::values(agFix, format = "matrix")
    agObs_mat_3 <- raster::values(agObs, format = "matrix")
    nMat <- raster::values(agN, format = "matrix")
    
    # matrix of cells that need to be combined
    combineMat <- which(nMat<minObs & nMat>0, arr.ind = T)
    
    while(nrow(combineMat)>0){
      
      # get surrounding cells for each matrix
      rows <- (combineMat[1,1]-1):(combineMat[1,1]+1)
      rows <- rows[which(rows>0 & rows <= nrow(nMat))]
      cols <- (combineMat[1,2]-1):(combineMat[1,2]+1)
      cols <- cols[which(cols>0 & cols <= ncol(nMat))]
      
      agLo09_mat_j <- agLo09_mat_3[rows,cols]
      agChm09_mat_j <- agChm09_mat_3[rows,cols]
      agFix_mat_j <- agFix_mat_3[rows,cols]
      agObs_mat_j <- agObs_mat_3[rows,cols]
      nMat_j <- nMat[rows,cols]
      
      # set cell of interest to 0
      nMat_j[which(rows==combineMat[1,1]),which(cols==combineMat[1,2])] <- 0
      
      # Only proceed if there are non-NA neighboring cells
      if(length(c(nMat_j)[c(nMat_j)>0])>0){  
        
        # find neighboring cell with most nearby observations
        new_j <- which(nMat_j==max(nMat_j,na.rm=T), arr.ind = T)[1,]
        
        # replace value for new cell with a weighted mean of all other values
        agLo09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agLo09_mat_j[new_j[1],new_j[2]],agLo09_mat_3[combineMat[1,1],combineMat[1,2]]),
                                                         w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agChm09_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agChm09_mat_j[new_j[1],new_j[2]],agChm09_mat_3[combineMat[1,1],combineMat[1,2]]),
                                                          w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agFix_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agFix_mat_j[new_j[1],new_j[2]],agFix_mat_3[combineMat[1,1],combineMat[1,2]]),
                                                        w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        agObs_mat_j[new_j[1],new_j[2]] <- weighted.mean(x = c(agObs_mat_j[new_j[1],new_j[2]],agObs_mat_3[combineMat[1,1],combineMat[1,2]]),
                                                        w = c(nMat_j[new_j[1],new_j[2]],nMat[combineMat[1,1],combineMat[1,2]]))
        
        nMat_j[new_j[1],new_j[2]] <- nMat_j[new_j[1],new_j[2]] + nMat[combineMat[1,1],combineMat[1,2]]
        
        # replace neighborhood in original matrices
        agLo09_mat_3[rows,cols] <- agLo09_mat_j
        agChm09_mat_3[rows,cols] <- agChm09_mat_j
        agFix_mat_3[rows,cols] <- agFix_mat_j
        agObs_mat_3[rows,cols] <- agObs_mat_j
        nMat[rows,cols] <- nMat_j
      }
      
      # Replace original observation with NaN
      agLo09_mat_3[combineMat[1,1],combineMat[1,2]] <- NaN
      agChm09_mat_3[combineMat[1,1],combineMat[1,2]] <- NaN
      agFix_mat_3[combineMat[1,1],combineMat[1,2]] <- NaN
      agObs_mat_3[combineMat[1,1],combineMat[1,2]] <- NaN
      nMat[combineMat[1,1],combineMat[1,2]] <- 0
      
      # Redefine combineMat
      combineMat <- which(nMat<minObs & nMat >0, arr.ind = T)
    }
    
    valsRemove <- which(c(nMat)<minObs)
    agLo09_mat_3[valsRemove] <- NA
    agChm09_mat_3[valsRemove] <- NA
    agFix_mat_3[valsRemove] <- NA
    agObs_mat_3[valsRemove] <- NA
    
    agLo09_3 <- raster::raster(agLo09_mat_3, template=agLo09)
    agChm09_3 <- raster::raster(agChm09_mat_3, template=agChm09)
    agFix_3 <- raster::raster(agFix_mat_3, template=agFix)
    agObs_3 <- raster::raster(agObs_mat_3, template=agObs)
    
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
                                  agForm = raster::values(agForm_3))
    # Rename age, parent material, soil form codes
    bestAgResults_3[bestAgResults_3$agAge==1 & !is.na(bestAgResults_3$agAge),"agAge"] <- "OldGrowth"
    bestAgResults_3[bestAgResults_3$agAge==2 & !is.na(bestAgResults_3$agAge),"agAge"] <- "Secondary"
    bestAgResults_3[bestAgResults_3$agParent==1 & !is.na(bestAgResults_3$agParent),"agParent"] <- "CaimitoVolcanic"
    bestAgResults_3[bestAgResults_3$agParent==2 & !is.na(bestAgResults_3$agParent),"agParent"] <- "Andesite"
    bestAgResults_3[bestAgResults_3$agParent==3 & !is.na(bestAgResults_3$agParent),"agParent"] <- "Bohio"
    bestAgResults_3[bestAgResults_3$agParent==4 & !is.na(bestAgResults_3$agParent),"agParent"] <- "CaimitoMarineSedimentary"
    bestAgResults_3[bestAgResults_3$agForm==1 & !is.na(bestAgResults_3$agForm),"agForm"] <- "RedLightClay"
    bestAgResults_3[bestAgResults_3$agForm==2 & !is.na(bestAgResults_3$agForm),"agForm"] <- "PaleSwellingClay"
    bestAgResults_3[bestAgResults_3$agForm==3 & !is.na(bestAgResults_3$agForm),"agForm"] <- "BrownFineLoam"
    bestAgResults_3[bestAgResults_3$agForm==4 & !is.na(bestAgResults_3$agForm),"agForm"] <- "MottledHeavyClay"
    bestAgResults_3 <- bestAgResults_3[!is.na(bestAgResults_3$agObs),]
    
    
    # Plot 1: Low canopy with predicted values from fixed effects
    # Multiply values by 100 to get rates in %
    bestAgResults_1$agFixPct <- 100*bestAgResults_1$agFix
    bestAgResults_1$agObsPct <- 100*bestAgResults_1$agObs
    bestAgResults_2$agFixPct <- 100*bestAgResults_2$agFix
    bestAgResults_2$agObsPct <- 100*bestAgResults_2$agObs
    bestAgResults_3$agFixPct <- 100*bestAgResults_3$agFix
    bestAgResults_3$agObsPct <- 100*bestAgResults_3$agObs
    
    
    PtCex_1 <- 0.8
    PtCex_2 <- 1
    PtCex_3 <- 2
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agFixPct")],na.rm=T),
         ylim=c(100*min(bestAgResults_1[,c("agLo09")]),250*max(bestAgResults_1[,c("agLo09")])),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_1],"ha resolution"),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFixPct")],na.rm=T),
           y=250*max(bestAgResults_1[,c("agLo09")],na.rm=T),
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_1$agFixPct),
                                      y = log(bestAgResults_1$agLo09))$estimate,3)),
         x = 1.4, y = 0.5)
    
    plot(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agFixPct")],na.rm=T),
         ylim=100*range(bestAgResults_2[,c("agLo09")],na.rm=T),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_2],"ha resolution"),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_2$agFixPct),
                                      y = log(bestAgResults_2$agLo09))$estimate,3)),
         x = 1.4, y = 2.2)
    
    plot(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agFixPct")],na.rm=T),
         ylim=100*range(bestAgResults_3[,c("agLo09")],na.rm=T),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_3],"ha resolution"),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_3$agFixPct),
                                      y = log(bestAgResults_3$agLo09))$estimate,3)),
         x = 1.4, y = 4.5)
    
    # By soil parent material
    plot(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agFixPct")],na.rm=T),
         ylim=c(100*min(bestAgResults_1[,c("agLo09")]),250*max(bestAgResults_1[,c("agLo09")])),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFixPct")]),
           y=250*max(bestAgResults_1[,c("agLo09")]),
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agFixPct")],na.rm=T),
         ylim=100*range(bestAgResults_2[,c("agLo09")],na.rm=T),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agFixPct")],na.rm=T),
         ylim=100*range(bestAgResults_3[,c("agLo09")],na.rm=T),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agFixPct")],na.rm=T),
         ylim=c(100*min(bestAgResults_1[,c("agLo09")]),250*max(bestAgResults_1[,c("agLo09")])),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(100*agLo09~agFixPct, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=min(bestAgResults_1[,c("agFixPct")]),
           y=250*max(bestAgResults_1[,c("agLo09")]),
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agFixPct")],na.rm=T),
         ylim=100*range(bestAgResults_2[,c("agLo09")],na.rm=T),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(100*agLo09~agFixPct, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agFixPct")],na.rm=T),
         ylim=100*range(bestAgResults_3[,c("agLo09")],na.rm=T),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(100*agLo09~agFixPct, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Predicted disturbance rate (fixed effects only, % yr-1)",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Low canopy area in 2009 (%)", side=2, outer=T, line=1.5)
    par(las=1)
    
    # Plot 2: Low canopy vs observed values
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agObsPct")]),
         ylim=c(100*min(bestAgResults_1[,c("agLo09")]),250*max(bestAgResults_1[,c("agLo09")])),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_1],"ha resolution"),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.3,
           y=250*max(bestAgResults_1[,c("agLo09")]),
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_1$agObsPct),
                                      y = log(bestAgResults_1$agLo09))$estimate,3)),
         x = 10, y = 0.5)
    
    plot(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agObsPct")]),
         ylim=100*range(bestAgResults_2[,c("agLo09")]),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_2],"ha resolution"),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_2$agObsPct),
                                      y = log(bestAgResults_2$agLo09))$estimate,3)),
         x = 4, y = 2.2)
    
    plot(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agObsPct")]),
         ylim=100*range(bestAgResults_3[,c("agLo09")]),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_3],"ha resolution"),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_3$agObsPct),
                                      y = log(bestAgResults_3$agLo09))$estimate,3)),
         x = 3, y = 4.5)
    
    # By soil parent material
    plot(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agObsPct")]),
         ylim=c(100*min(bestAgResults_1[,c("agLo09")]),250*max(bestAgResults_1[,c("agLo09")])),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.3,
           y=250*max(bestAgResults_1[,c("agLo09")]),
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agObsPct")]),
         ylim=100*range(bestAgResults_2[,c("agLo09")]),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agObsPct")]),
         ylim=100*range(bestAgResults_3[,c("agLo09")]),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agObsPct")]),
         ylim=c(100*min(bestAgResults_1[,c("agLo09")]),250*max(bestAgResults_1[,c("agLo09")])),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(100*agLo09~agObsPct, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.3,
           y=250*max(bestAgResults_1[,c("agLo09")]),
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agObsPct")]),
         ylim=100*range(bestAgResults_2[,c("agLo09")]),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(100*agLo09~agObsPct, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agObsPct")]),
         ylim=100*range(bestAgResults_3[,c("agLo09")]),
         log=c("xy"),
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(100*agLo09~agObsPct, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Observed disturbance rate (% yr-1)",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Low canopy area in 2009 (%)", side=2, outer=T, line=1.5)
    par(las=1)
    
    # Plot 3: Mean canopy height with predicted values from fixed effects
    PtCex_1 <- 0.8
    PtCex_2 <- 1
    PtCex_3 <- 2
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agFixPct")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_1],"ha resolution"),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.8,
           y=18,
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_1$agFixPct),
                                      y = log(bestAgResults_1$agChm09))$estimate,3)),
         x = 1.5, y = 13)
    
    plot(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agFixPct")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_2],"ha resolution"),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_2$agFixPct),
                                      y = log(bestAgResults_2$agChm09))$estimate,3)),
         x = 1, y = 18)
    
    plot(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agFixPct")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_3],"ha resolution"),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_3$agFixPct),
                                      y = log(bestAgResults_3$agChm09))$estimate,3)),
         x = 1, y = 20)
    
    # By soil parent material
    plot(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agFixPct")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.8,
           y=18,
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agFixPct")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agFixPct")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)

    
    
    # Plot by forest age
    plot(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agFixPct")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(agChm09~agFixPct, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.8,
           y=18,
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agFixPct")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(agChm09~agFixPct, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agFixPct")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(agChm09~agFixPct, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)

    mtext("Predicted disturbance rate (fixed effects only, % yr-1)",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Mean canopy height in 2009 (m)", side=2, outer=T, line=1.5)
    par(las=1)
    
    # Plot 4: Mean canopy height vs observed values
    par(mfrow=c(3,3), mar=c(2,3,1,0), oma=c(4,4,1,1))
    
    # By soil form
    plot(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="RedLightClay",],
         xlim=range(bestAgResults_1[,c("agObsPct")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_1],"ha resolution"),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.3,
           y=18,
           c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
           col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_1$agObsPct),
                                      y = log(bestAgResults_1$agChm09))$estimate,3)),
         x = 10, y = 13)
    
    plot(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="RedLightClay",],
         xlim=range(bestAgResults_2[,c("agObsPct")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_2],"ha resolution"),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_2)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_2$agObsPct),
                                      y = log(bestAgResults_2$agChm09))$estimate,3)),
         x = 4.5, y = 18)
    
    plot(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="RedLightClay",],
         xlim=range(bestAgResults_3[,c("agObsPct")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],0.6),
         main = paste(xVals[agScale_3],"ha resolution"),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,
         xaxt="n")
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="BrownFineLoam",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="PaleSwellingClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agForm=="MottledHeavyClay",],
           col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],0.6),
           pch=19,
           cex=PtCex_3)
    text(paste0("r = ",round(cor.test(x = log(bestAgResults_3$agObsPct),
                                      y = log(bestAgResults_3$agChm09))$estimate,3)),
         x = 3, y = 20)
    
    # By soil parent material
    plot(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_1[,c("agObsPct")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.3,
           y=18,
           c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
           col=adjustcolor(c(colAnd,colBoh,colMar,colVol),1),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_2[,c("agObsPct")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoVolcanic",],
         xlim=range(bestAgResults_3[,c("agObsPct")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA,xaxt="n")
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="CaimitoMarineSedimentary",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="Andesite",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agParent=="Bohio",],
           col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],0.6),
           pch=19,
           cex=PtCex_3)
    
    # Plot by forest age
    plot(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agAge=="OldGrowth",],
         xlim=range(bestAgResults_1[,c("agObsPct")]),
         ylim=range(bestAgResults_1[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_1,
         xlab = NA, ylab=NA)
    points(agChm09~agObsPct, data = bestAgResults_1[bestAgResults_1$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_1)
    legend(x=0.3,
           y=18,
           c("Old growth","Secondary"),
           col=adjustcolor(c(colOld,colSec),0.6),
           pch=19, pt.cex=1.5,
           bty="n")
    
    plot(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agAge=="OldGrowth",],
         xlim=range(bestAgResults_2[,c("agObsPct")]),
         ylim=range(bestAgResults_2[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_2,
         xlab = NA, ylab=NA)
    points(agChm09~agObsPct, data = bestAgResults_2[bestAgResults_2$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_2)
    
    plot(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agAge=="OldGrowth",],
         xlim=range(bestAgResults_3[,c("agObsPct")]),
         ylim=range(bestAgResults_3[,c("agChm09")]),
         log="xy",
         col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],0.6),
         pch=19,
         cex=PtCex_3,
         xlab = NA, ylab=NA)
    points(agChm09~agObsPct, data = bestAgResults_3[bestAgResults_3$agAge=="Secondary",],
           col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],0.6),
           pch=19,
           cex=PtCex_3)
    
    mtext("Observed disturbance rate (% yr-1)",
          side=1, outer=T, line=1)
    par(las=0)
    mtext("Mean canopy height in 2009 (m)", side=2, outer=T, line=1.5)
    par(las=1)   
    
#### Figure S25: Comparison of 2009 lidar data and average spatial pattern (untransformed) ####
    
    #Create rasters where all non-NA values are 1
    predMeanRaster18n <- predMeanRaster18
    predMeanRaster18n[!is.na(raster::values(predMeanRaster18n))] <- 1
    fixRaster18n <- fixRaster18
    fixRaster18n[!is.na(raster::values(fixRaster18n))] <- 1
    predMeanRaster20n <- predMeanRaster20
    predMeanRaster20n[!is.na(raster::values(predMeanRaster20n))] <- 1
    fixRaster20n <- fixRaster20
    fixRaster20n[!is.na(raster::values(fixRaster20n))] <- 1
    
    # Raster of low canopy area in 2009
    lo09 <- raster::raster("Data_HeightRasters/loCanopy2009.tif") # pixels with value 1 are => 10 m height, 0 are < 10 m
    # Raster of canopy height in 2009
    chm09 <- raster::raster("Data_HeightRasters/CHM_2009_QAQC.tif") 
    
    
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
      
      # Replace 0 values for log transformation
      #agLo09_mat[agLo09_mat==0 & !is.na(agLo09_mat)] <- 1/3200
      
      loResults$fixRlo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$estimate
      loResults$fixRlo_lo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults$fixRlo_hi[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults$obsRlo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$estimate
      loResults$obsRlo_lo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults$obsRlo_hi[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[2]
      
      loResults$fixRchm[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$estimate
      loResults$fixRchm_lo[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults$fixRchm_hi[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults$obsRchm[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$estimate
      loResults$obsRchm_lo[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults$obsRchm_hi[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[2]
      
    }
    
    
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
      
      # Replace 0 values for log transformation
      #agLo09_mat[agLo09_mat==0 & !is.na(agLo09_mat)] <- 1/3200
      
      loResults_b$fixRlo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$estimate
      loResults_b$fixRlo_lo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults_b$fixRlo_hi[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_b$obsRlo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$estimate
      loResults_b$obsRlo_lo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_b$obsRlo_hi[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[2]
      
      loResults_b$fixRchm[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$estimate
      loResults_b$fixRchm_lo[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat[valsKeep])))$conf.int[1]
      loResults_b$fixRchm_hi[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_b$obsRchm[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$estimate
      loResults_b$obsRchm_lo[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_b$obsRchm_hi[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[2]
      
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
      
      # Replace 0 values for log transformation
      #agLo09_mat[agLo09_mat==0 & !is.na(agLo09_mat)] <- 1/3200
      
      valsKeep <- which(c(nMat)>minObs)
      loResults_c$fixRlo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$estimate
      loResults_c$fixRlo_lo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[1]
      loResults_c$fixRlo_hi[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_c$obsRlo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$estimate
      loResults_c$obsRlo_lo[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_c$obsRlo_hi[i] <- cor.test(x = (c(agLo09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[2]
      
      loResults_c$fixRchm[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$estimate
      loResults_c$fixRchm_lo[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat[valsKeep])))$conf.int[1]
      loResults_c$fixRchm_hi[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agFix_mat)[valsKeep]))$conf.int[2]
      loResults_c$obsRchm[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$estimate
      loResults_c$obsRchm_lo[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[1]
      loResults_c$obsRchm_hi[i] <- cor.test(x = (c(agChm09_mat)[valsKeep]), y = (c(agObs_mat)[valsKeep]))$conf.int[2]
    }
    
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
           c("Predicted disturbance rate (fixed effects)",
             "Observed disturbance rate"),
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
    
#### Figure S26: Alternate model with initial canopy height as predictor #### 
    
    fixedResults <- model_full$summary.fixed
    
    # separated by years
    fixedResultsHt <- model_full_ht$summary.fixed
    
    fixedNames <- c("Curvature (linear)", "Slope (linear)", "Slope (quadratic)",                    
                    "Height above drainage (linear)","Height above drainage (quadratic)",
                    "Soil parent: Andesite", "Soil parent: Caimito marine", "Soil parent: Caimito volcanic",
                    "Soil form: Mottled heavy clay", "Soil form: Pale swelling clay", "Soil form: Red light clay",
                    "Forest age: Secondary",
                    "Year: 2015-2018",
                    "Initial Canopy Height")
    
    par(mfrow=c(1,1), mar=c(4,1,1,1), oma=c(0,0,0,0))
    
    plot(x = fixedResults[2:nrow(fixedResults),"mean"],
         y = 1.2*nrow(fixedResultsHt):3,
         xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(-0.4,0.6),
         ylim=c(2,1.2*nrow(fixedResultsHt)),
         xlab=NA,
         pch = 19, 
         cex = 1,
         ylab = NA, yaxt = "n")
    arrows(x0 = fixedResults[2:nrow(fixedResults),"0.025quant"],
           x1 = fixedResults[2:nrow(fixedResults),"0.975quant"],
           y0 = 1.2*nrow(fixedResultsHt):3,
           y1 = 1.2*nrow(fixedResultsHt):3,
           angle = 90, code=3,
           length = 0.05)
    
    points(x = fixedResultsHt[2:nrow(fixedResultsHt),"mean"],
           y = 1.2*nrow(fixedResultsHt):2 - 0.3,
           pch = 19, 
           cex = 1,
           col = "grey")
    arrows(x0 = fixedResultsHt[2:nrow(fixedResultsHt),"0.025quant"],
           x1 = fixedResultsHt[2:nrow(fixedResultsHt),"0.975quant"],
           y0 = 1.2*nrow(fixedResultsHt):2 - 0.3,
           y1 = 1.2*nrow(fixedResultsHt):2 - 0.3,
           angle = 90, code=3,
           length = 0.05,
           col = "grey")
    
    abline(v=0, lty=2)
    text(fixedNames,
         x = c(fixedResults[2:nrow(fixedResults),"0.975quant"]+0.03,fixedResultsHt[nrow(fixedResultsHt),"0.975quant"]),
         y = 1.2*nrow(fixedResultsHt):2,
         pos = 4)
    
    legend(x = -1.3,
           y = 18,
           c("Original full model",
             "Model with canopy height"),
           bty="n",
           col = c("black", "darkgrey"),
           pch=c(19,19))
    mtext("Fixed effect", side=1, outer=F, line=2)
    