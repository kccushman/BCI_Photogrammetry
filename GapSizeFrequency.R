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
  
# Only area greater than 10 m height initially
  allVals15to18tall <- raster::values(d15to18tall)
  areaSampled15to18tall <- length(allVals15to18tall[!is.na(allVals15to18tall)])/10000
  
  allVals18to20tall <- raster::values(d18to20tall)
  areaSampled18to20tall <- length(allVals18to20tall[!is.na(allVals18to20tall)])/10000
    
# Total number of new canopy disturbance events
  length(c(gaps15to18sp$area,
           gaps18to20sp$area))
  
# Area sampled in old growth and secondary forest [> 10 m height]
  # Forest age polygon
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    age$AgeClass <- "Other"
    age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
    age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    ageUse <- age[!(age$AgeClass=="Other"),]
  
  areaOld18 <- raster::mask(d15to18tall, ageUse[ageUse$AgeClass=="OldGrowth",])
  valsOld18 <- raster::values(areaOld18)
  areaSampledOld18 <- length(valsOld18[!is.na(valsOld18)])/10000
  
  areaSec18 <- raster::mask(d15to18tall, ageUse[ageUse$AgeClass=="Secondary",])
  valsSec18 <- raster::values(areaSec18)
  areaSampledSec18 <- length(valsSec18[!is.na(valsSec18)])/10000
  
  areaOld20 <- raster::mask(d18to20tall, ageUse[ageUse$AgeClass=="OldGrowth",])
  valsOld20 <- raster::values(areaOld20)
  areaSampledOld20 <- length(valsOld20[!is.na(valsOld20)])/10000
  
  areaSec20 <- raster::mask(d18to20tall, ageUse[ageUse$AgeClass=="Secondary",])
  valsSec20 <- raster::values(areaSec20)
  areaSampledSec20 <- length(valsSec20[!is.na(valsSec20)])/10000
  
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
          
          bootBlocks$nGap[i] <- length(valGap[!is.na(valGap) & !is.na(valAll)]) # only count gap pixels > 10 m height initially
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
  
    # % higher in second interval
    (2.05-1.79)/1.79*100
  
  # Mean frequency of gaps per interval--number of gaps
  round(mean(gapSummary15to18$gapsPerHa)/nYr15to18,2)
  round(mean(gapSummary18to20$gapsPerHa)/nYr18to20,2)
  
  round(quantile(gapSummary15to18$gapsPerHa,probs = c(0.025,0.975))/nYr15to18,2)
  round(quantile(gapSummary18to20$gapsPerHa,probs = c(0.025,0.975))/nYr18to20,2)
  
    # % higher in first interval
    (1.92-1.64)/1.64*100

  
#### BOOTSTRAPPED SIZE FREQUENCY DISTRIBUTIONS ####
  source("makesizedistforRaquel.r")
  source("fitsizedistforRaquel.R")
  source("sizedistpart3forRaquel.R")
  
  # find max gap size
  mxSz <- quantile(c(gaps15to18sp$area,
              gaps18to20sp$area),1)
  
  # find which gap sizes are NA
  gapPres <- data.frame(minSz = seq(24.5,(mxSz - 0.5),by=1),
                        maxSz = seq(25.5,(mxSz + 0.5),by=1),
                        n18 = NA,
                        n20 = NA)
  
  for(i in 1:nrow(gapPres)){
    gapPres$n18[i] <- length(gaps15to18sp$area[gaps15to18sp$area>gapPres$minSz[i] & gaps15to18sp$area<gapPres$maxSz[i]])
    gapPres$n20[i] <- length(gaps18to20sp$area[gaps18to20sp$area>gapPres$minSz[i] & gaps18to20sp$area<gapPres$maxSz[i]])
  }
  szUse <- unique(c(gapPres$minSz[(gapPres$n18>0 | gapPres$n20>0)],gapPres$maxSz[(gapPres$n18>0 | gapPres$n20>0)]))
  szUse <- szUse[order(szUse)]
  
  allData15to18 <- data.frame(dbh = gaps15to18sp$area,
                              quadnum = gaps15to18sp$block)
  
  allData18to20 <- data.frame(dbh = gaps18to20sp$area,
                              quadnum = gaps18to20sp$block)
  
#### Bootstrapped size frequency by year ####
  
    i <- 1
    datadir <- "SizeFreqResults/"
    site <- "BCI"
    census <- "15to18"
    szFreq15to18 <- doallbootstrapdbhfits(alldata = allData15to18,
                                          nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                          alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                          fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                          # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                          filestem="gaps15to18_new", # this is the beginning of the names of the output files
                                          # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                          ddiv=szUse)
    
   

    i <- 1
    datadir <- "SizeFreqResults/"
    site <- "BCI"
    census <- "18to20"
    szFreq18to20 <- doallbootstrapdbhfits(alldata = allData18to20,
                                          nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                          alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                          fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                          # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                          filestem="gaps18to20_new", # this is the beginning of the names of the output files
                                          # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                          ddiv=szUse)
    
    
#### Bootstrapped size frequency by forest age ####
    
    # Forest age polygon
      age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
      age$AgeClass <- "Other"
      age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
      age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
      ageUse <- age[!(age$AgeClass=="Other"),]
      
    # Find forest age for each gap
      gaps15to18sp$Age <- sp::over(gaps15to18sp, ageUse[,c("AgeClass")])
      gaps18to20sp$Age <- sp::over(gaps18to20sp, ageUse[,c("AgeClass")])
      
    # Run code for each age
      
      oldGrowthGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$Age=="OldGrowth"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$Age=="OldGrowth"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$Age=="OldGrowth"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$Age=="OldGrowth"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "OldGrowth"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = oldGrowthGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="gapsOldGrowth", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)
      
      
      secondaryGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$Age=="Secondary"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$Age=="Secondary"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$Age=="Secondary"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$Age=="Secondary"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "Secondary"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = secondaryGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="gapsSecondary", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)
      
#### Bootstrapped size frequency by soil parent material ####
      
    # Soil type polygon  
      soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
      soil <- sp::spTransform(soil,raster::crs(d15to18tall))
      
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
      
      # Find forest age for each gap
      gaps15to18sp$SoilParent <- sp::over(gaps15to18sp, soil[,c("SoilParent")])
      gaps18to20sp$SoilParent <- sp::over(gaps18to20sp, soil[,c("SoilParent")])
      
      # Run code for each age
      
      andesiteGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="Andesite"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="Andesite"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="Andesite"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="Andesite"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "Andesite"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = andesiteGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="gapsAndesite", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)      
      
      bohioGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="Bohio"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="Bohio"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="Bohio"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="Bohio"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "Bohio"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = bohioGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="gapsBohio", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)
      
      caimitoVolcanicGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="CaimitoVolcanic"],
                                    quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="CaimitoVolcanic"]),
                         data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="CaimitoVolcanic"],
                                    quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="CaimitoVolcanic"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "CaimitoVolcanic"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = caimitoVolcanicGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="gapsCaimitoVolcanic", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)
      
      caimitoMarineGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="CaimitoMarineSedimentary"],
                                              quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="CaimitoMarineSedimentary"]),
                                   data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="CaimitoMarineSedimentary"],
                                              quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="CaimitoMarineSedimentary"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "CaimitoMarine"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = caimitoMarineGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="gapsCaimitoMarine", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)      
#### Plot bootstrapped results by year ####
  
  # Read distribution results
  szFreq15to18 <- read.table("SizeFreqResults/gaps15to18sizedistbsfit.txt", header = T)
  szFreq18to20 <- read.table("SizeFreqResults/gaps18to20sizedistbsfit.txt", header = T)

  # What distribution has the highest likelihood?
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
    # 2015 to 2017
    round(-szFreq15to18$weibloglike-max(c(-szFreq15to18$weibloglike,-szFreq15to18$powloglike,-szFreq15to18$exploglike)),0) 
    round(-szFreq15to18$powloglike-max(c(-szFreq15to18$weibloglike,-szFreq15to18$powloglike,-szFreq15to18$exploglike)),0) 
    round(-szFreq15to18$exploglike-max(c(-szFreq15to18$weibloglike,-szFreq15to18$powloglike,-szFreq15to18$exploglike)),0) 
    
    # 2018 to 2020
    round(-szFreq18to20$weibloglike-max(c(-szFreq18to20$weibloglike,-szFreq18to20$powloglike,-szFreq18to20$exploglike)),0) 
    round(-szFreq18to20$powloglike-max(c(-szFreq18to20$weibloglike,-szFreq18to20$powloglike,-szFreq18to20$exploglike)),0) 
    round(-szFreq18to20$exploglike-max(c(-szFreq18to20$weibloglike,-szFreq18to20$powloglike,-szFreq18to20$exploglike)),0) 
  
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
  round(szFreq15to18[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreq18to20[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  
  round(szFreq15to18[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],2)
  round(szFreq18to20[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],2)
  

  # Make plot of data and best-fit lines
  
  # find max gap size
  mxSz <- max(c(gaps15to18sp$area,
                gaps18to20sp$area))
  
  # Divide intervals for plotting
  gapAreaRange <- c(25,mxSz)
  logRange <- log(gapAreaRange)
  logMid <- sum(logRange)/1.8
  brksRange_log <- c(seq(logRange[1],logMid,length.out = 20),seq(logMid,logRange[2],length.out = 6)[-1])
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
  
  # getExpEsts <- function(szFreqResults, gapSp, xVals){
  #   
  #   
  #   scaleExp <- 1/(pexp(mxSz,rate = szFreqResults$exppar1est) -
  #                    pexp(25, rate = szFreqResults$exppar1est))
  #   
  #   yValsExp <- scaleExp*dexp(xVals, rate = szFreqResults$exppar1est, log=F)*nrow(gapSp)
  #   
  #   return(yValsExp)
  # }
  
  yValsWeib18 <- getWeibEsts(szFreqResults = szFreq15to18, gapSp = gaps15to18sp, xVals)/nYr15to18/areaSampled15to18tall
  yValsWeib20 <- getWeibEsts(szFreqResults = szFreq18to20, gapSp = gaps18to20sp, xVals)/nYr18to20/areaSampled18to20tall
  
  # yValsExp18 <- getExpEsts(szFreqResults = szFreq15to18, gapSp = gaps15to18sp, xVals)/nYr15to18/areaSampled15to18tall
  # yValsExp20 <- getExpEsts(szFreqResults = szFreq18to20, gapSp = gaps18to20sp, xVals)/nYr18to20/areaSampled18to20tall
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizes15to18/nYr15to18,gapSizes18to20/nYr18to20)) + c(1e-7,0.03)
  
  #Define colors
  col18 <- "blue"
  col20 <- "#d95f02"
  
  par(las = 1, mar=c(3,5,1,3), oma=c(1,1,1,0), mfrow=c(1,2))
  plot(x = brksMids, y = gapSizes15to18/nYr15to18,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(col18,0.6),
       log = logOption,
       pch=19,
       cex.axis = 0.8,
       ylab = expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
       xlab = NA)
  mtext(expression("Disturbance area (m"^"2"~")"), side=1,outer = T)
  text("a", x = range(brksMids)[1], y = yRange[2])
  
  points(x = brksMids, y = gapSizes18to20/nYr18to20,
         col = adjustcolor(col20,0.6), pch=19)
  
  legend(x=25,y=0.00001,
         c("2015 - 2018","2018 - 2020"),
         col=adjustcolor(c(col18,col20),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeib18, col = adjustcolor(col18,0.6), lwd=2)
  lines(x = xVals, y = yValsWeib20, col = adjustcolor(col20,0.6), lwd=2)
  
  plot(x = gaps15to18sp$area,
       y = -gaps15to18sp$htDrop,
       log="xy",
       xlim=range(c(gaps15to18sp$area,gaps18to20sp$area)),
       ylim=-range(c(gaps15to18sp$htDrop,gaps18to20sp$htDrop)),
       ylab= "Average height decrease (m)",
       col = adjustcolor(col18,0.2),
       cex = 0.3,
       pch=19)
  text("b",
       x = min(c(gaps15to18sp$area,gaps18to20sp$area)),
       y = -max(c(gaps15to18sp$htDrop,gaps18to20sp$htDrop)))
  points(x = gaps18to20sp$area,
       y = -gaps18to20sp$htDrop,
       col = adjustcolor(col20,0.2),
       cex = 0.3,
       pch=19)

#### Plot bootstrapped results by forest age ####
  
  # Read distribution results
  szFreqOld <- read.table("SizeFreqResults/gapsOldGrowthsizedistbsfit.txt", header = T)
  szFreqSec <- read.table("SizeFreqResults/gapsSecondarysizedistbsfit.txt", header = T)
  
  
  # What distribution has the highest likelihood?
    # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
    
    # Old growth
    round(-szFreqOld$weibloglike-max(c(-szFreqOld$weibloglike,-szFreqOld$powloglike,-szFreqOld$exploglike)),0) 
    round(-szFreqOld$powloglike-max(c(-szFreqOld$weibloglike,-szFreqOld$powloglike,-szFreqOld$exploglike)),0) 
    round(-szFreqOld$exploglike-max(c(-szFreqOld$weibloglike,-szFreqOld$powloglike,-szFreqOld$exploglike)),0) 
    
    # Secondary
    round(-szFreqSec$weibloglike-max(c(-szFreqSec$weibloglike,-szFreqSec$powloglike,-szFreqSec$exploglike)),0) 
    round(-szFreqSec$powloglike-max(c(-szFreqSec$weibloglike,-szFreqSec$powloglike,-szFreqSec$exploglike)),0) 
    round(-szFreqSec$exploglike-max(c(-szFreqSec$weibloglike,-szFreqSec$powloglike,-szFreqSec$exploglike)),0) 
  
  # Get r-squared values for each
  
    # Old growth
    round(szFreqOld$weibcatadjr2,3)
    round(szFreqOld$powcatadjr2,3)
    round(szFreqOld$expcatadjr2,3)
    
    # Secondary
    round(szFreqSec$weibcatadjr2,3)
    round(szFreqSec$powcatadjr2,3)
    round(szFreqSec$expcatadjr2,3)
  
  # Weibull has the lowest -log likelihood in both intervals
  round(szFreqOld[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqSec[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  
  round(szFreqOld[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],2)
  round(szFreqSec[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],2)
  
  
  # Make plot of data and best-fit lines
  
  # find max gap size
  mxSz <- max(c(gaps15to18sp$area,
                gaps18to20sp$area))
  
  # Divide intervals for plotting
  gapAreaRange <- c(25,mxSz)
  logRange <- log(gapAreaRange)
  logMid <- sum(logRange)/1.8
  brksRange_log <- c(seq(logRange[1],logMid,length.out = 20),seq(logMid,logRange[2],length.out = 6)[-1])
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
  
  gapSizesOld <- makeGapVectors(gapAreaVector = oldGrowthGaps$dbh,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/(areaSampledOld18*nYr15to18 + areaSampledOld20*nYr18to20)
  
  gapSizesSec <- makeGapVectors(gapAreaVector = secondaryGaps$dbh,
                                   minarea = brksMins,
                                   maxarea = brksMaxs)/(areaSampledSec18*nYr15to18 + areaSampledSec20*nYr18to20)
  
  
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
  
  # getExpEsts <- function(szFreqResults, gapSp, xVals){
  #   
  #   
  #   scaleExp <- 1/(pexp(mxSz,rate = szFreqResults$exppar1est) -
  #                    pexp(25, rate = szFreqResults$exppar1est))
  #   
  #   yValsExp <- scaleExp*dexp(xVals, rate = szFreqResults$exppar1est, log=F)*nrow(gapSp)
  #   
  #   return(yValsExp)
  # }
  
  yValsWeibOld <- getWeibEsts(szFreqResults = szFreqOld, gapSp = oldGrowthGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledOld18+areaSampledOld20)
  yValsWeibSec <- getWeibEsts(szFreqResults = szFreqSec, gapSp = secondaryGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledSec18+areaSampledSec20)
  
  # yValsExp18 <- getExpEsts(szFreqResults = szFreqOld, gapSp = gapsOldsp, xVals)/nYrOld/areaSampledOldtall
  # yValsExp20 <- getExpEsts(szFreqResults = szFreqSec, gapSp = gapsSecsp, xVals)/nYrSec/areaSampledSectall
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizesOld,gapSizesSec)) + c(1e-7,0)
  
  #Define colors
  colOld <- "darkgreen"
  colSec <- "lightgreen"
  
  par(las = 1, mar=c(4,5,2,1))
  plot(x = brksMids, y = gapSizesOld,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(colOld,0.6),
       log = logOption,
       pch=19,
       cex.axis = 0.8,
       ylab = expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
       xlab = expression("Disturbance area (m"^"2"~")"))
  
  points(x = brksMids, y = gapSizesSec,
         col = adjustcolor(colSec,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("Old growth","Secondary"),
         col=adjustcolor(c(colOld,colSec),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeibOld, col = adjustcolor(colOld,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibSec, col = adjustcolor(colSec,0.6), lwd=2)
  
  
  
  
  
  
  
  
  