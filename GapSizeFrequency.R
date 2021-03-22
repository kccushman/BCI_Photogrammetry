# Do preliminary gap analyses

#### READ DATA ####
    
# Canopy height change rasters
  d15to18 <- raster::raster("dCHM15to18_tin.tif")
  d18to19 <- raster::raster("dCHM18to19_tin.tif")
  d19to20 <- raster::raster("dCHM19to20_tin.tif")
  d18to20 <- raster::raster("dCHM18to20_tin.tif")
  
# Canopy height change rasters where only possible gap values (> 5 m initially) are included  
  d15to18tall <- raster::raster("dCHM15to18tall_tin.tif")
  d18to19tall <- raster::raster("dCHM18to19tall_tin.tif")
  d19to20tall <- raster::raster("dCHM19to20tall_tin.tif")
  d18to20tall <- raster::raster("dCHM18to20tall_tin.tif")
  
# Gap rasters
  gaps15to18 <- raster::raster("newGaps15to18_tin.tif")
  gaps18to19 <- raster::raster("newGaps18to19_tin.tif")
  gaps19to20 <- raster::raster("newGaps19to20_tin.tif")
  gaps18to20 <- raster::raster("newGaps18to20_tin.tif")
  
# Gap polygons
  gaps15to18sp <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
  gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile_tin/gaps18to19sp.shp")
  gaps19to20sp <- rgdal::readOGR("gaps19to20_shapefile_tin/gaps19to20sp.shp")
  gaps18to20sp <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")
  
# Block values for bootstrapping
  blockData <- read.csv("bootstrapBlocks.csv")

#### AREA SAMPLED EACH YEAR ####
  
# All area
  allVals15to18 <- raster::values(d15to18)
  areaSampled15to18 <- length(allVals15to18[!is.na(allVals15to18)])/10000
  
  allVals18to19 <- raster::values(d18to19)
  areaSampled18to19 <- length(allVals18to19[!is.na(allVals18to19)])/10000
  
  allVals19to20 <- raster::values(d19to20)
  areaSampled19to20 <- length(allVals19to20[!is.na(allVals19to20)])/10000  
  
  allVals18to20 <- raster::values(d18to20)
  areaSampled18to20 <- length(allVals18to20[!is.na(allVals18to20)])/10000  
  
# Only area greater than 5 m height initially
  allVals15to18tall <- raster::values(d15to18tall)
  areaSampled15to18tall <- length(allVals15to18tall[!is.na(allVals15to18tall)])/10000
  
  allVals18to19tall <- raster::values(d18to19tall)
  areaSampled18to19tall <- length(allVals18to19tall[!is.na(allVals18to19tall)])/10000
  
  allVals19to20tall <- raster::values(d19to20tall)
  areaSampled19to20tall <- length(allVals19to20tall[!is.na(allVals19to20tall)])/10000
  
  allVals18to20tall <- raster::values(d18to20tall)
  areaSampled18to20tall <- length(allVals18to20tall[!is.na(allVals18to20tall)])/10000
    
# Total number of new canopy disturbance events
  length(c(gaps15to18sp[gaps15to18sp$use==T,]$area,
           gaps18to19sp[gaps18to19sp$use==T,]$area,
           gaps19to20sp[gaps19to20sp$use==T,]$area))
  
  length(c(gaps15to18sp[gaps15to18sp$use==T,]$area,
           gaps18to20sp[gaps18to20sp$use==T,]$area))
  
#### SUMMARY GAP STATS PER YEAR ####
  
  # Define a function to estimate the gap % and #/ha using bootstrapping, 
  # where different blocks are randomly omitted
  
  
  getGapSummary <- function(gapLayer,allLayer,bootBlocks, nBoot){
    
    # first, count the number of all pixels (with data) and gap pixels for each block
      bootBlocks$nGap <- NA
      bootBlocks$nAll <- NA
      bootBlocks$nGapID <- NA
      
      bootGapIDs <- list()
      for(i in 1:dim(bootBlocks)[1]){
        # create extent object from block limits
          blockExt <- raster::extent(as.numeric(bootBlocks[i,c("xmin","xmax","ymin","ymax")]))
          
        # crop gap and all pixel rasters to block extent
          blockGap <- raster::crop(gapLayer, blockExt)
          blockAll <- raster::crop(allLayer, blockExt)
          
        # count and save non-NA pixels in block
          valGap <- raster::values(blockGap)
          valAll <- raster::values(blockAll)
          bootBlocks$nGap[i] <- length(valGap[!is.na(valGap)])
          bootBlocks$nAll[i] <- length(valAll[!is.na(valAll)])
          
        # save unique gap IDs in a list
          bootBlocks$nGapID[i] <- length(unique(valGap[!is.na(valGap)]))
      }
      
    # Resample blocks with replacement to estimate the total 
      results <- data.frame(ID = 1:nBoot,
                            gapsPerHa = rep(NA,nBoot),
                            percentGap = rep(NA,nBoot))
      
      set.seed(1)
      for(i in 1:nBoot){
        bootSample <- sample(1:dim(bootBlocks)[1], dim(bootBlocks)[1], replace = T)
        
        # Find percent of total area in gaps
        nGapPixels <- sum(bootBlocks[bootSample,"nGap"])
        nAllPixels <- sum(bootBlocks[bootSample,"nAll"])
        results$percentGap[i] <- 100*nGapPixels/nAllPixels
        
        # Find number of gaps per ha
        nHa <- nAllPixels/10000
        nGaps <- sum(bootBlocks[bootSample,"nGapID"])
        results$gapsPerHa[i] <- nGaps/nHa
      }
    
      return(results)
  }
  
  gapSummary15to18 <- getGapSummary(gapLayer = gaps15to18,
                                    allLayer = d15to18tall,
                                    bootBlocks = blockData,
                                    nBoot = 1000)
  
  gapSummary18to19 <- getGapSummary(gapLayer = gaps18to19,
                                    allLayer = d18to19tall,
                                    bootBlocks = blockData,
                                    nBoot = 1000)
  
  gapSummary19to20 <- getGapSummary(gapLayer = gaps19to20,
                                    allLayer = d19to20tall,
                                    bootBlocks = blockData,
                                    nBoot = 1000)
  
  gapSummary18to20 <- getGapSummary(gapLayer = gaps18to20,
                                    allLayer = d18to20tall,
                                    bootBlocks = blockData,
                                    nBoot = 1000)
  
  round(mean(gapSummary15to18$gapsPerHa)/3,2)
  round(mean(gapSummary18to19$gapsPerHa),2)
  round(mean(gapSummary19to20$gapsPerHa)/(13/12),2)
  round(mean(gapSummary18to20$gapsPerHa)/(25/12),2)
  
  round(quantile(gapSummary15to18$gapsPerHa,probs = c(0.025,0.975))/3,2)
  round(quantile(gapSummary18to19$gapsPerHa,probs = c(0.025,0.975)),2)
  round(quantile(gapSummary19to20$gapsPerHa,probs = c(0.025,0.975))/(13/12),2)
  round(quantile(gapSummary18to20$gapsPerHa,probs = c(0.025,0.975))/(25/12),2)
  
  round(mean(gapSummary15to18$percentGap)/3,2)
  round(mean(gapSummary18to19$percentGap),2)
  round(mean(gapSummary19to20$percentGap)/(13/12),2)
  round(mean(gapSummary18to20$percentGap)/(25/12),2)
  
  round(quantile(gapSummary15to18$percentGap,probs = c(0.025,0.975))/3,2)
  round(quantile(gapSummary18to19$percentGap,probs = c(0.025,0.975)),2)
  round(quantile(gapSummary19to20$percentGap,probs = c(0.025,0.975))/(13/12),2)
  round(quantile(gapSummary18to20$percentGap,probs = c(0.025,0.975))/(25/12),2)
  
#### BOOTSTRAPPED SIZE FREQUENCY DISTRIBUTIONS ####
  source("makesizedistforRaquel.r")
  source("fitsizedistforRaquel.R")
  source("sizedistpart3forRaquel.R")
  
  # find max gap size
  mxSz <- quantile(c(gaps15to17sp[gaps15to17sp$use==T,]$area,
              gaps17to18sp[gaps17to18sp$use==T,]$area,
              gaps18to19sp[gaps18to19sp$use==T,]$area,
              gaps19to20sp[gaps19to20sp$use==T,]$area),0.999)
  
  
  allData15to17 <- data.frame(dbh = gaps15to17sp[gaps15to17sp$use==T,]$area,
                              quadnum = gaps15to17sp[gaps15to17sp$use==T,]$block)
  
  i <- 1
  datadir <- "SizeFreqResults/"
  site <- "BCI"
  census <- "15to17"
  szFreq15to17 <- doallbootstrapdbhfits(alldata = allData15to17,
                                        nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                        alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                        fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                        # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                        filestem="gaps15to17", # this is the beginning of the names of the output files
                                        ddiv=seq(9.5,(mxSz + 0.5),by=1))
  
  allData17to18 <- data.frame(dbh = gaps17to18sp[gaps17to18sp$use==T,]$area,
                              quadnum = gaps17to18sp[gaps17to18sp$use==T,]$block)
  datadir <- "SizeFreqResults/"
  site <- "BCI"
  census <- "17to18"
  szFreq17to18 <- doallbootstrapdbhfits(alldata = allData17to18,
                                        nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                        alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                        fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                        # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                        filestem="gaps17to18", # this is the beginning of the names of the output files
                                        ddiv=seq(9.5,(mxSz + 0.5),by=1))
  
  
  allData18to19 <- data.frame(dbh = gaps18to19sp[gaps18to19sp$use==T,]$area,
                              quadnum = gaps18to19sp[gaps18to19sp$use==T,]$block)
  datadir <- "SizeFreqResults/"
  site <- "BCI"
  census <- "18to19"
  szFreq18to19 <- doallbootstrapdbhfits(alldata = allData18to19,
                                        nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                        alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                        fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                        # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                        filestem="gaps18to19", # this is the beginning of the names of the output files
                                        ddiv=seq(9.5,(mxSz + 0.5),by=1))
  
  
  allData19to20 <- data.frame(dbh = gaps19to20sp[gaps19to20sp$use==T,]$area,
                              quadnum = gaps19to20sp[gaps19to20sp$use==T,]$block)
  i <- 1
  datadir <- "SizeFreqResults/"
  site <- "BCI"
  census <- "19to20"
  szFreq19to20 <- doallbootstrapdbhfits(alldata = allData19to20,
                                        nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                        alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                        fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                        # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                        filestem="gaps19to20", # this is the beginning of the names of the output files
                                        ddiv=seq(9.5,(mxSz + 0.5),by=1))
  
#### READ AND PLOT BOOTSTRAPPED RESULTS ####
  
  # Read distribution results
  szFreq15to17 <- read.table("SizeFreqResults/gaps15to17sizedistbsfit.txt", header = T)
  szFreq17to18 <- read.table("SizeFreqResults/gaps17to18sizedistbsfit.txt", header = T)
  szFreq18to19 <- read.table("SizeFreqResults/gaps18to19sizedistbsfit.txt", header = T)
  szFreq19to20 <- read.table("SizeFreqResults/gaps19to20sizedistbsfit.txt", header = T)
  
  
  # What distribution has the highest adjusted R2?
  
  # 2015 to 2017
  round(szFreq15to17$weibcatadjr2,2)
  round(szFreq15to17$powcatadjr2,2)
  round(szFreq15to17$expcatadjr2,2)
  
  # 2017 to 2018
  round(szFreq17to18$weibcatadjr2,2)
  round(szFreq17to18$powcatadjr2,2)
  round(szFreq17to18$expcatadjr2,2)
  
  # 2018 to 2019
  round(szFreq18to19$weibcatadjr2,2)
  round(szFreq18to19$powcatadjr2,2)
  round(szFreq18to19$expcatadjr2,2)
  
  # 2019 to 2020
  round(szFreq19to20$weibcatadjr2,2)
  round(szFreq19to20$powcatadjr2,2)
  round(szFreq19to20$expcatadjr2,2)
  
  #Weibull has highest R2
  round(szFreq15to17[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],2)
  round(szFreq17to18[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],2)
  round(szFreq18to19[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],2)
  round(szFreq19to20[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],2)
  
  round(szFreq15to17[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],1)
  round(szFreq17to18[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],1)
  round(szFreq18to19[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],1)
  round(szFreq19to20[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],1)
  
  
  # Make plot of data and best-fit lines
  
  # find max gap size
  mxSz <- max(c(gaps15to18sp[gaps15to18sp$use==T,]$area,
                gaps18to19sp[gaps18to19sp$use==T,]$area,
                gaps19to20sp[gaps19to20sp$use==T,]$area,
                gaps18to20sp[gaps18to20sp$use==T,]$area))
  
  # Divide intervals for plotting
  gapAreaRange <- c(25,mxSz)
  logRange <- log(gapAreaRange)
  brksRange_log <- seq(logRange[1],logRange[2],length.out = 18)
  brksRange <- floor(exp(brksRange_log))
  
  brksMins <- brksRange[1:length(brksRange)-1]
  brksMaxs <- brksRange[2:length(brksRange)]
  brksMids <- (brksMins + brksMaxs)/2
  
  # Make vectors of gaps of each area range, normalized by the number of "bins" combined
  makeGapVectors <- function(gapAreaVector, minarea, maxarea){
    nArea <- rep(NA,length(minarea))
    for(i in 1:length(minarea)){
      nGaps <- length(gapAreaVector[gapAreaVector>=minarea[i] & gapAreaVector<maxarea[i]])
      szRange <- maxarea[i]-minarea[i]
      nArea[i] <-nGaps/szRange
    }
    return(nArea)
  }
  
  gapSizes15to18 <- makeGapVectors(gapAreaVector = gaps15to18sp[gaps15to18sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/areaSampled15to18
  
  gapSizes18to19 <- makeGapVectors(gapAreaVector = gaps18to19sp[gaps18to19sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/areaSampled18to19
  
  gapSizes19to20 <- makeGapVectors(gapAreaVector = gaps19to20sp[gaps19to20sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/areaSampled19to20
  
  gapSizes18to20 <- makeGapVectors(gapAreaVector = gaps18to20sp[gaps18to20sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/areaSampled18to20
  
  
  # Get correct fitted y-values for each year, adjusted for truncated distributions
  xVals <- seq(gapAreaRange[1],gapAreaRange[2],length.out = 1000)
  
  getWeibEsts <- function(szFreqResults, gapSp, xVals){
    scaleWeib <- 1/(pweibull(mxSz,
                               shape= szFreqResults$weibpar1est,
                               scale = szFreqResults$weibpar2est) -
                        pweibull(10,
                                 shape= szFreqResults$weibpar1est,
                                 scale = szFreqResults$weibpar2est))
    yValsWeib <- scaleWeib*dweibull(xVals,
                                    shape= szFreqResults$weibpar1est,
                                    scale = szFreqResults$weibpar2est)*nrow(gapSp[gapSp$use==T,])
    return(yValsWeib)
  }
  
  yValsWeib17 <- getWeibEsts(szFreqResults = szFreq15to17, gapSp = gaps15to17sp, xVals)/2/areaSampled15to17
  yValsWeib18 <- getWeibEsts(szFreqResults = szFreq17to18, gapSp = gaps17to18sp, xVals)/areaSampled17to18
  yValsWeib19 <- getWeibEsts(szFreqResults = szFreq18to19, gapSp = gaps18to19sp, xVals)/areaSampled18to19
  yValsWeib20 <- getWeibEsts(szFreqResults = szFreq19to20, gapSp = gaps19to20sp, xVals)/(13/12)/areaSampled19to20
  
  # Plot 1 : three intervals
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizes15to18,gapSizes18to19,gapSizes19to20,gapSizes18to20)) + c(1e-7,0)
  
  #Define colors
  col18 <- "black"
  col19 <- "#1b9e77"
  col20 <- "#d95f02"

  par(las = 1)
  plot(x = brksMids, y = gapSizes15to18/3,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(col18,0.6),
       log = logOption,
       pch=19,
       cex.axis = 0.8,
       ylab = expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
       xlab = expression("Disturbance area (m"^"2"~")"))
  
  points(x = brksMids, y = gapSizes18to19,
         col = adjustcolor(col19,0.6), pch=19)
  points(x = brksMids, y = gapSizes19to20/(23/12),
         col = adjustcolor(col20,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("'15-'18","'18-'19","'19-'20"),
         col=adjustcolor(c(col18,col19,col20),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeib18, col = adjustcolor(col18,0.6), lwd=2)
  lines(x = xVals, y = yValsWeib19, col = adjustcolor(col19,0.6), lwd=2)
  lines(x = xVals, y = yValsWeib20, col = adjustcolor(col20,0.6), lwd=2)
  
  # Plot 2 : two intervals
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizes15to18,gapSizes18to19,gapSizes19to20,gapSizes18to20)) + c(1e-7,0)
  
  #Define colors
  col18 <- "black"
  col19 <- "#1b9e77"
  col20 <- "#d95f02"
  
  par(las = 1)
  plot(x = brksMids, y = gapSizes15to18/3,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(col18,0.6),
       log = logOption,
       pch=19,
       cex.axis = 0.8,
       ylab = expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
       xlab = expression("Disturbance area (m"^"2"~")"))
  
  points(x = brksMids, y = gapSizes18to20/(25/12),
         col = adjustcolor(col20,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("'15-'18","'18-'20"),
         col=adjustcolor(c(col18,col20),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeib18, col = adjustcolor(col18,0.6), lwd=2)
  lines(x = xVals, y = yValsWeib20, col = adjustcolor(col20,0.6), lwd=2)
  

  scaleExp20 <- 1/(pexp(mxSz,rate = szFreq19to20$exppar1est) -
                     pexp(10, rate = szFreq19to20$exppar1est))
  yValsExp20 <- scaleExp20*dexp(xVals, rate = szFreq19to20$exppar1est)*nrow(gaps19to20sp[gaps19to20sp$use==T,])
  
  yValsPow17 <- xVals^(-szFreq15to17$powpar1est)*nrow(gaps15to17sp[gaps15to17sp$use==T,])/2
  
  

  
  