# Do preliminary gap analyses

#### READ DATA ####
    
# Canopy height change rasters
  d15to18 <- raster::raster("dCHM15to18_tin.tif")
  d18to20 <- raster::raster("dCHM18to20_tin.tif")
  
# Canopy height change rasters where only possible gap values (> 5 m initially) are included  
  d15to18tall <- raster::raster("dCHM15to18tall_tin.tif")
  d18to20tall <- raster::raster("dCHM18to20tall_tin.tif")
  
# Gap rasters
  gaps15to18 <- raster::raster("newGaps15to18_tin.tif")
  gaps18to20 <- raster::raster("newGaps18to20_tin.tif")
  
# Gap polygons
  gaps15to18sp <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
  gaps18to20sp <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")
  
# Block values for bootstrapping
  blockData <- read.csv("bootstrapBlocks.csv")

#### AREA SAMPLED EACH YEAR ####
  
# All area
  allVals15to18 <- raster::values(d15to18)
  areaSampled15to18 <- length(allVals15to18[!is.na(allVals15to18)])/10000
  
  allVals18to20 <- raster::values(d18to20)
  areaSampled18to20 <- length(allVals18to20[!is.na(allVals18to20)])/10000  
  
# Only area greater than 5 m height initially
  allVals15to18tall <- raster::values(d15to18tall)
  areaSampled15to18tall <- length(allVals15to18tall[!is.na(allVals15to18tall)])/10000
  
  allVals18to20tall <- raster::values(d18to20tall)
  areaSampled18to20tall <- length(allVals18to20tall[!is.na(allVals18to20tall)])/10000
    
# Total number of new canopy disturbance events
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
  
  gapSummary18to20 <- getGapSummary(gapLayer = gaps18to20,
                                    allLayer = d18to20tall,
                                    bootBlocks = blockData,
                                    nBoot = 1000)
  
  nYr15to18 <- as.numeric(as.Date("2018-06-07") - as.Date("2015-06-26"))/365
  nYr18to20 <- as.numeric(as.Date("2020-07-31") - as.Date("2018-06-07"))/365
  
  # Mean frequency of gaps per interval--percent of area
  round(mean(gapSummary15to18$percentGap)/nYr15to18,2)
  round(mean(gapSummary18to20$percentGap)/nYr18to20,2)
  
  round(quantile(gapSummary15to18$percentGap,probs = c(0.025,0.975))/nYr15to18,2)
  round(quantile(gapSummary18to20$percentGap,probs = c(0.025,0.975))/nYr18to20,2)
  
    # % higher in first interval
    (2.57-2.04)/2.04*100
  
  # Mean frequency of gaps per interval--number of gaps
  round(mean(gapSummary15to18$gapsPerHa)/nYr15to18,2)
  round(mean(gapSummary18to20$gapsPerHa)/nYr18to20,2)
  
  round(quantile(gapSummary15to18$gapsPerHa,probs = c(0.025,0.975))/nYr15to18,2)
  round(quantile(gapSummary18to20$gapsPerHa,probs = c(0.025,0.975))/nYr18to20,2)
  
    # % higher in first interval
    (1.97-1.84)/1.84*100

  
#### BOOTSTRAPPED SIZE FREQUENCY DISTRIBUTIONS ####
  source("makesizedistforRaquel.r")
  source("fitsizedistforRaquel.R")
  source("sizedistpart3forRaquel.R")
  
  # find max gap size
  mxSz <- quantile(c(gaps15to18sp$area,
              gaps18to20sp$area),1)
  
  
  allData15to18 <- data.frame(dbh = gaps15to18sp$area,
                              quadnum = gaps15to18sp$block)
  
  i <- 1
  datadir <- "SizeFreqResults/"
  site <- "BCI"
  census <- "15to18"
  szFreq15to18 <- doallbootstrapdbhfits(alldata = allData15to18,
                                        nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                        alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                        fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                        # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                        filestem="gaps15to18", # this is the beginning of the names of the output files
                                        ddiv=seq(24.5,(mxSz + 0.5),by=1))
  
 
  allData18to20 <- data.frame(dbh = gaps18to20sp$area,
                              quadnum = gaps18to20sp$block)
  i <- 1
  datadir <- "SizeFreqResults/"
  site <- "BCI"
  census <- "18to20"
  szFreq18to20 <- doallbootstrapdbhfits(alldata = allData18to20,
                                        nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                        alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                        fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                        # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                        filestem="gaps18to20", # this is the beginning of the names of the output files
                                        ddiv=seq(24.5,(mxSz + 0.5),by=1))
  
#### READ AND PLOT BOOTSTRAPPED RESULTS ####
  
  # Read distribution results
  szFreq15to18 <- read.table("SizeFreqResults/gaps15to18sizedistbsfit.txt", header = T)
  szFreq18to20 <- read.table("SizeFreqResults/gaps18to20sizedistbsfit.txt", header = T)
  
  
  # What distribution has the lowest -log likelihood?
  
    # 2015 to 2017
    round(szFreq15to18$weibloglike,0)
    round(szFreq15to18$powloglike,0)
    round(szFreq15to18$exploglike,0)
    
    # 2018 to 2020
    round(szFreq18to20$weibloglike,0)
    round(szFreq18to20$powloglike,0)
    round(szFreq18to20$exploglike,0)
  
  # Get r-squared values for each
    
    # 2015 to 2017
    round(szFreq15to18$weibcatadjr2,3)
    round(szFreq15to18$powcatadjr2,3)
    round(szFreq15to18$expcatadjr2,3)
    
    # 2018 to 2020
    round(szFreq18to20$weibcatadjr2,3)
    round(szFreq18to20$powcatadjr2,3)
    round(szFreq18to20$expcatadjr2,3)
  
  # Weibull has the lowest -log likelihood in both intervals
  round(szFreq15to18[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],4)
  round(szFreq18to20[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],4)
  
  round(szFreq15to18[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],2)
  round(szFreq18to20[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],2)
  

  # Make plot of data and best-fit lines
  
  # find max gap size
  mxSz <- max(c(gaps15to18sp$area,
                gaps18to20sp$area))
  
  # Divide intervals for plotting
  gapAreaRange <- c(25,mxSz)
  logRange <- log(gapAreaRange)
  brksRange_log <- c(seq(logRange[1],7,length.out = 15),logRange[2])
  brksRange <- ceiling(exp(brksRange_log))
  
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
  
  gapSizes15to18 <- makeGapVectors(gapAreaVector = gaps15to18sp$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/areaSampled15to18tall
  
  gapSizes18to20 <- makeGapVectors(gapAreaVector = gaps18to20sp$area,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/areaSampled18to20tall
  
  
  # Get correct fitted y-values for each year, adjusted for truncated distributions
  xVals <- seq(gapAreaRange[1],gapAreaRange[2],length.out = 1000)
  
  getWeibEsts <- function(szFreqResults, gapSp, xVals){
    
    scaleWeib <- 1/(pweibull(mxSz,
                               shape= szFreqResults$weibpar1est,
                               scale = szFreqResults$weibpar2est) -
                        pweibull(25,
                                 shape= szFreqResults$weibpar1est,
                                 scale = szFreqResults$weibpar2est))
    yValsWeib <- scaleWeib*dweibull(xVals,
                                    shape= szFreqResults$weibpar1est,
                                    scale = szFreqResults$weibpar2est)*nrow(gapSp)
    return(yValsWeib)
  }
  
  getExpEsts <- function(szFreqResults, gapSp, xVals){
    
    
    scaleExp <- 1/(pexp(mxSz,rate = szFreqResults$exppar1est) -
                     pexp(25, rate = szFreqResults$exppar1est))
    
    yValsExp <- scaleExp*dexp(xVals, rate = szFreqResults$exppar1est, log=F)*nrow(gapSp)
    
    return(yValsExp)
  }
  
  yValsWeib18 <- getWeibEsts(szFreqResults = szFreq15to18, gapSp = gaps15to18sp, xVals)/nYr15to18/areaSampled15to18tall
  yValsWeib20 <- getWeibEsts(szFreqResults = szFreq18to20, gapSp = gaps18to20sp, xVals)/nYr18to20/areaSampled18to20tall
  
  yValsExp18 <- getExpEsts(szFreqResults = szFreq15to18, gapSp = gaps15to18sp, xVals)/nYr15to18/areaSampled15to18tall
  yValsExp20 <- getExpEsts(szFreqResults = szFreq18to20, gapSp = gaps18to20sp, xVals)/nYr18to20/areaSampled18to20tall
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizes15to18/nYr15to18,gapSizes18to20/nYr18to20)) + c(1e-7,0)
  
  #Define colors
  col18 <- "black"
  col20 <- "#d95f02"
  
  par(las = 1, mar=c(4,5,2,1))
  plot(x = brksMids, y = gapSizes15to18/nYr15to18,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(col18,0.6),
       log = logOption,
       pch=19,
       cex.axis = 0.8,
       ylab = expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
       xlab = expression("Disturbance area (m"^"2"~")"))
  
  points(x = brksMids, y = gapSizes18to20/nYr18to20,
         col = adjustcolor(col20,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("2015 - 2018","2018 - 2020"),
         col=adjustcolor(c(col18,col20),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeib18, col = adjustcolor(col18,0.6), lwd=2)
  lines(x = xVals, y = yValsWeib20, col = adjustcolor(col20,0.6), lwd=2)
  


  

  

  
  