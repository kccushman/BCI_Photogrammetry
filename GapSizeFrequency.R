# Do preliminary gap analyses

#### READ DATA ####
    
# Canopy height change rasters
  d15to17 <- raster::raster("dCHM15to17.tif")
  d17to18 <- raster::raster("dCHM17to18.tif")
  d18to19 <- raster::raster("dCHM18to19.tif")
  d19to20 <- raster::raster("dCHM19to20.tif")
  
# Canopy height change rastes where only possible gap values (> 5 m initially) are included  
  d15to17tall <- raster::raster("dCHM15to17tall.tif")
  d17to18tall <- raster::raster("dCHM17to18tall.tif")
  d18to19tall <- raster::raster("dCHM18to19tall.tif")
  d19to20tall <- raster::raster("dCHM19to20tall.tif")
  
# Gap rasters
  gaps15to17 <- raster::raster("newGaps15to17.tif")
  gaps17to18 <- raster::raster("newGaps17to18.tif")
  gaps18to19 <- raster::raster("newGaps18to19.tif")
  gaps19to20 <- raster::raster("newGaps19to20.tif")
  
# Gap polygons
  gaps15to17sp <- rgdal::readOGR("gaps15to17_shapefile/gaps15to17sp.shp")
  gaps17to18sp <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")
  gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")
  gaps19to20sp <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")
  
# Block values for bootstrapping
  blockData <- read.csv("bootstrapBlocks.csv")

#### AREA SAMPLED EACH YEAR ####
  
# All area
  allVals15to17 <- raster::values(d15to17)
  areaSampled15to17 <- length(allVals15to17[!is.na(allVals15to17)])/10000
  
  allVals17to18 <- raster::values(d17to18)
  areaSampled17to18 <- length(allVals17to18[!is.na(allVals17to18)])/10000
  
  allVals18to19 <- raster::values(d18to19)
  areaSampled18to19 <- length(allVals18to19[!is.na(allVals18to19)])/10000
  
  allVals19to20 <- raster::values(d19to20)
  areaSampled19to20 <- length(allVals19to20[!is.na(allVals19to20)])/10000  
  
# Only area greater than 5 m height initially
  allVals15to17tall <- raster::values(d15to17tall)
  areaSampled15to17tall <- length(allVals15to17tall[!is.na(allVals15to17tall)])/10000
  
  allVals17to18tall <- raster::values(d17to18tall)
  areaSampled17to18tall <- length(allVals17to18tall[!is.na(allVals17to18tall)])/10000
  
  allVals18to19tall <- raster::values(d18to19tall)
  areaSampled18to19tall <- length(allVals18to19tall[!is.na(allVals18to19tall)])/10000
  
  allVals19to20tall <- raster::values(d19to20tall)
  areaSampled19to20tall <- length(allVals19to20tall[!is.na(allVals19to20tall)])/10000
    
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
  
  gapSummary15to17 <- getGapSummary(gapLayer = gaps15to17,
                                    allLayer = d15to17tall,
                                    bootBlocks = blockData,
                                    nBoot = 1000)
  
  gapSummary17to18 <- getGapSummary(gapLayer = gaps17to18,
                                    allLayer = d17to18tall,
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
  
  round(mean(gapSummary15to17$gapsPerHa),2)/2
  round(mean(gapSummary17to18$gapsPerHa),2)
  round(mean(gapSummary18to19$gapsPerHa),2)
  round(mean(gapSummary19to20$gapsPerHa),2)
  round(quantile(gapSummary15to17$gapsPerHa,probs = c(0.025,0.975)),2)/2
  round(quantile(gapSummary17to18$gapsPerHa,probs = c(0.025,0.975)),2)
  round(quantile(gapSummary18to19$gapsPerHa,probs = c(0.025,0.975)),2)
  round(quantile(gapSummary19to20$gapsPerHa,probs = c(0.025,0.975)),2)
  
  round(mean(gapSummary15to17$percentGap),2)/2
  round(mean(gapSummary17to18$percentGap),2)
  round(mean(gapSummary18to19$percentGap),2)
  round(mean(gapSummary19to20$percentGap),2)
  round(quantile(gapSummary15to17$percentGap,probs = c(0.025,0.975)),2)/2
  round(quantile(gapSummary17to18$percentGap,probs = c(0.025,0.975)),2)
  round(quantile(gapSummary18to19$percentGap,probs = c(0.025,0.975)),2)
  round(quantile(gapSummary19to20$percentGap,probs = c(0.025,0.975)),2)
  
#### BOOTSTRAPPED SIZE FREQUENCY DISTRIBUTIONS ####
  source("makesizedistforRaquel.r")
  source("fitsizedistforRaquel.R")
  source("sizedistpart3forRaquel.R")
  
  # find max gap size
  mxSz <- max(c(gaps15to17sp[gaps15to17sp$use==T,]$area,
              gaps17to18sp[gaps17to18sp$use==T,]$area,
              gaps18to19sp[gaps18to19sp$use==T,]$area,
              gaps19to20sp[gaps19to20sp$use==T,]$area))
  
  
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
  
  i <- 1
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
  
  i <- 1
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
  szFreq17to18 <- read.table("SizeFreqResults/gaps17to18sizedistbsfit.txt", header = T)
  szFreq18to19 <- read.table("SizeFreqResults/gaps18to19sizedistbsfit.txt", header = T)
  szFreq19to20 <- read.table("SizeFreqResults/gaps19to20sizedistbsfit.txt", header = T)
  
  # Divide intervals for plotting
  gapAreaRange <- c(10,1000)
  logRange <- log(gapAreaRange)
  brksRange_log <- seq(logRange[1],logRange[2],length.out = 20)
  brksRange <- floor(exp(brksRange_log))
  
  brksMins <- brksRange[1:length(brksRange)-1]
  brksMaxs <- brksRange[2:length(brksRange)]
  brksMids <- (brksMins + brksMaxs)/2
  
  # Make vectors of gaps of each area range, normalized by the number of "bins" combined
  makeGapVectors <- function(gapAreaVector, minarea, maxarea){
    nArea <- rep(NA,length(minarea))
    for(i in 1:length(minarea)){
      nGaps <- length(gapAreaVector[gapAreaVector>minarea[i] & gapAreaVector<maxarea[i]])
      szRange <- maxarea[i]-minarea[i]
      nArea[i] <-nGaps/szRange
    }
    return(nArea)
  }
  
  gapSizes17to18 <- makeGapVectors(gapAreaVector = gaps17to18sp[gaps17to18sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)
  
  gapSizes18to19 <- makeGapVectors(gapAreaVector = gaps18to19sp[gaps18to19sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)
  
  gapSizes19to20 <- makeGapVectors(gapAreaVector = gaps19to20sp[gaps19to20sp$use==T,]$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)
  
  logOption <- "xy"
  #logOption <- ""
  
  
# All years together  
  plot(x = brksMids, y = gapSizes17to18,
       xlim=gapAreaRange,
       ylim=range(c(gapSizes17to18,gapSizes18to19,gapSizes19to20)),
       col = adjustcolor("black",0.5),
       log = logOption,
       pch=19,
       ylab = "n Gaps (#/m2)",
       xlab = "Gap area (m2)")
  
  points(x = brksMids, y = gapSizes18to19,
         col = adjustcolor("blue",0.5), pch=19)
  points(x = brksMids, y = gapSizes19to20,
         col = adjustcolor("red",0.5), pch=19)
  
# Look at fits of different distributions
  
  # 2017 to 2018
  szFreq17to18$weibloglike
  szFreq17to18$powloglike
  szFreq17to18$exploglike
  
  plot(x = brksMids, y = gapSizes17to18,
       xlim=gapAreaRange,
       ylim=range(c(gapSizes17to18,gapSizes18to19,gapSizes19to20)),
       col = adjustcolor("black",0.5),
       log = logOption,
       pch=19,
       ylab = "n Gaps (#/m2)",
       xlab = "Gap area (m2)")
  
  xVals <- seq(gapAreaRange[1],gapAreaRange[2],length.out = 1000)
  
  # use pweibull to correct for truncation mutliply by 1/(pwebmax - pweibwin)
  
  yValsWeib <- dweibull(xVals,
                        shape= szFreq17to18$weibpar1est,
                        scale = szFreq17to18$weibpar2est)*sum(gapSizes17to18)
  lines(x = xVals, y = yValsWeib, lwd=2)
  
  yValsExp <- dexp(xVals, rate = szFreq17to18$exppar1est)*sum(gapSizes17to18)
  lines(x = xVals, y = yValsExp, lwd=2, lty=2)
  
  yValsPow <- xVals^(1-szFreq17to18$powpar1est)*sum(gapSizes17to18)
  lines(x = xVals, y = yValsPow, lwd=2, lty=3)
  
  