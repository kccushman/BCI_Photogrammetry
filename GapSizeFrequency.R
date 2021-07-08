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

# Area sampled by soil parent material [> 10 m height]
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
  
  areaAnd18 <- raster::mask(d15to18tall, soil[soil$SoilParent=="Andesite",])
  valsAnd18 <- raster::values(areaAnd18)
  areaSampledAnd18 <- length(valsAnd18[!is.na(valsAnd18)])/10000
  
  areaBoh18 <- raster::mask(d15to18tall, soil[soil$SoilParent=="Bohio",])
  valsBoh18 <- raster::values(areaBoh18)
  areaSampledBoh18 <- length(valsBoh18[!is.na(valsBoh18)])/10000
  
  areaMar18 <- raster::mask(d15to18tall, soil[soil$SoilParent=="CaimitoMarineSedimentary",])
  valsMar18 <- raster::values(areaMar18)
  areaSampledMar18 <- length(valsMar18[!is.na(valsMar18)])/10000
  
  areaVol18 <- raster::mask(d15to18tall, soil[soil$SoilParent=="CaimitoVolcanic",])
  valsVol18 <- raster::values(areaVol18)
  areaSampledVol18 <- length(valsVol18[!is.na(valsVol18)])/10000
  
  areaAnd20 <- raster::mask(d18to20tall, soil[soil$SoilParent=="Andesite",])
  valsAnd20 <- raster::values(areaAnd20)
  areaSampledAnd20 <- length(valsAnd20[!is.na(valsAnd20)])/10000
  
  areaBoh20 <- raster::mask(d18to20tall, soil[soil$SoilParent=="Bohio",])
  valsBoh20 <- raster::values(areaBoh20)
  areaSampledBoh20 <- length(valsBoh20[!is.na(valsBoh20)])/10000
  
  areaMar20 <- raster::mask(d18to20tall, soil[soil$SoilParent=="CaimitoMarineSedimentary",])
  valsMar20 <- raster::values(areaMar20)
  areaSampledMar20 <- length(valsMar20[!is.na(valsMar20)])/10000
  
  areaVol20 <- raster::mask(d18to20tall, soil[soil$SoilParent=="CaimitoVolcanic",])
  valsVol20 <- raster::values(areaVol20)
  areaSampledVol20 <- length(valsVol20[!is.na(valsVol20)])/10000

# Area sampled by soil form
  areaPal18 <- raster::mask(d15to18tall, soil[soil$SoilForm=="PaleSwellingClay",])
  valsPal18 <- raster::values(areaPal18)
  areaSampledPal18 <- length(valsPal18[!is.na(valsPal18)])/10000
  
  areaRed18 <- raster::mask(d15to18tall, soil[soil$SoilForm=="RedLightClay",])
  valsRed18 <- raster::values(areaRed18)
  areaSampledRed18 <- length(valsRed18[!is.na(valsRed18)])/10000
  
  areaMot18 <- raster::mask(d15to18tall, soil[soil$SoilForm=="MottledHeavyClay",])
  valsMot18 <- raster::values(areaMot18)
  areaSampledMot18 <- length(valsMot18[!is.na(valsMot18)])/10000
  
  areaBro18 <- raster::mask(d15to18tall, soil[soil$SoilForm=="BrownFineLoam",])
  valsBro18 <- raster::values(areaBro18)
  areaSampledBro18 <- length(valsBro18[!is.na(valsBro18)])/10000
  
  
  areaPal20 <- raster::mask(d18to20tall, soil[soil$SoilForm=="PaleSwellingClay",])
  valsPal20 <- raster::values(areaPal20)
  areaSampledPal20 <- length(valsPal20[!is.na(valsPal20)])/10000
  
  areaRed20 <- raster::mask(d18to20tall, soil[soil$SoilForm=="RedLightClay",])
  valsRed20 <- raster::values(areaRed20)
  areaSampledRed20 <- length(valsRed20[!is.na(valsRed20)])/10000
  
  areaMot20 <- raster::mask(d18to20tall, soil[soil$SoilForm=="MottledHeavyClay",])
  valsMot20 <- raster::values(areaMot20)
  areaSampledMot20 <- length(valsMot20[!is.na(valsMot20)])/10000
  
  areaBro20 <- raster::mask(d18to20tall, soil[soil$SoilForm=="BrownFineLoam",])
  valsBro20 <- raster::values(areaBro20)
  areaSampledBro20 <- length(valsBro20[!is.na(valsBro20)])/10000
  
  
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
      
    # Find forest age for each gap
      gaps15to18sp$Age <- sp::over(gaps15to18sp, ageUse[,c("AgeClass")])
      gaps18to20sp$Age <- sp::over(gaps18to20sp, ageUse[,c("AgeClass")])
      
    # Run code for each age
      
      oldGrowthGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$Age=="OldGrowth"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$Age=="OldGrowth"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$Age=="OldGrowth"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$Age=="OldGrowth"]))
      
      secondaryGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$Age=="Secondary"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$Age=="Secondary"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$Age=="Secondary"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$Age=="Secondary"]))
      
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
      
      # Find soil parent material for each gap
      gaps15to18sp$SoilParent <- sp::over(gaps15to18sp, soil[,c("SoilParent")])
      gaps18to20sp$SoilParent <- sp::over(gaps18to20sp, soil[,c("SoilParent")])
      
      # Run code for each soil parent material
      
      andesiteGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="Andesite"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="Andesite"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="Andesite"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="Andesite"]))
     
      bohioGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="Bohio"],
                                    quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="Bohio"]),
                         data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="Bohio"],
                                    quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="Bohio"]))
      
      caimitoVolcanicGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="CaimitoVolcanic"],
                                              quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="CaimitoVolcanic"]),
                                   data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="CaimitoVolcanic"],
                                              quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="CaimitoVolcanic"]))
      
      caimitoMarineGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilParent=="CaimitoMarineSedimentary"],
                                            quadnum = gaps15to18sp$block[gaps15to18sp$SoilParent=="CaimitoMarineSedimentary"]),
                                 data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilParent=="CaimitoMarineSedimentary"],
                                            quadnum = gaps18to20sp$block[gaps18to20sp$SoilParent=="CaimitoMarineSedimentary"]))
      
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
#### Bootstrapped size frequency by soil form ####
      
      # Find soil form for each gap
      gaps15to18sp$SoilForm <- sp::over(gaps15to18sp, soil[,c("SoilForm")])
      gaps18to20sp$SoilForm <- sp::over(gaps18to20sp, soil[,c("SoilForm")])
      
      # Run code for each soil form
      
      lightClayGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilForm=="RedLightClay"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$SoilForm=="RedLightClay"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilForm=="RedLightClay"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$SoilForm=="RedLightClay"]))
      
      swellingClayGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilForm=="PaleSwellingClay"],
                                           quadnum = gaps15to18sp$block[gaps15to18sp$SoilForm=="PaleSwellingClay"]),
                                data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilForm=="PaleSwellingClay"],
                                           quadnum = gaps18to20sp$block[gaps18to20sp$SoilForm=="PaleSwellingClay"]))
      
      heavyClayGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilForm=="MottledHeavyClay"],
                                        quadnum = gaps15to18sp$block[gaps15to18sp$SoilForm=="MottledHeavyClay"]),
                             data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilForm=="MottledHeavyClay"],
                                        quadnum = gaps18to20sp$block[gaps18to20sp$SoilForm=="MottledHeavyClay"]))
      
      fineLoamGaps <- rbind(data.frame(dbh = gaps15to18sp$area[gaps15to18sp$SoilForm=="BrownFineLoam"],
                                       quadnum = gaps15to18sp$block[gaps15to18sp$SoilForm=="BrownFineLoam"]),
                            data.frame(dbh = gaps18to20sp$area[gaps18to20sp$SoilForm=="BrownFineLoam"],
                                       quadnum = gaps18to20sp$block[gaps18to20sp$SoilForm=="BrownFineLoam"]))
      i <- 1
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "RedLightClay"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = lightClayGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="lightClayGaps", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)      
      

      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "PaleSwellingClay"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = swellingClayGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="swellingClayGaps", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)
    
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "MottledHeavyClay"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = heavyClayGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="heavyClayGaps", # this is the beginning of the names of the output files
                                            # ddiv=seq(24.5,(mxSz + 0.5),by=1) #old format
                                            ddiv=szUse)
      
      datadir <- "SizeFreqResults/"
      site <- "BCI"
      census <- "BrownFineLoam"
      szFreq15to18 <- doallbootstrapdbhfits(alldata = fineLoamGaps,
                                            nreps=1000, # number of bootstrap replicates - should be 1000 or so
                                            alpha=c(0.05,0.01), # p-values for which to report results; confidence intervals given are for 1-alpha %
                                            fitfcn=c("weib","pow","exp"), # the types of function to be fitted to the size distribution:
                                            # "weib"= 2-parameter Weibull, "pow" = power function (1 parameter), "exp"=negative exponential (1 parameter)
                                            filestem="fineLoamGaps", # this is the beginning of the names of the output files
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
  
  
  smooth18 <- loess.smooth(x = gaps15to18sp$area,
                           y = -gaps15to18sp$htDrop,
                           degree=1, span = 1.2)
  smooth20 <- loess.smooth(x = gaps18to20sp$area,
                           y = -gaps18to20sp$htDrop,
                           degree=1, span = 1.2)
  
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
  lines(smooth18, col=col18,lwd=2)
  lines(smooth20, col=col20,lwd=2)
  
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
  colOld <- wesanderson::wes_palette("Moonrise2",4)[1]
  colSec <- wesanderson::wes_palette("Moonrise2",4)[2]
  
  par(las = 1, mar=c(3,4,2,1),oma=c(2,2,1,1), mfrow=c(1,3))
  plot(x = brksMids, y = gapSizesOld,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(colOld,0.6),
       log = logOption,
       pch=19,
       cex.axis = 1,
       ylab = NA,
       xlab = NA,
       main = "Forest age")
  mtext(expression("Disturbance area (m"^"2"~")"),
        side=1, outer=T)
  par(las=0)
  mtext(expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
        side=2, outer=T)
  par(las=1)
  points(x = brksMids, y = gapSizesSec,
         col = adjustcolor(colSec,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("Old growth","Secondary"),
         col=adjustcolor(c(colOld,colSec),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeibOld, col = adjustcolor(colOld,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibSec, col = adjustcolor(colSec,0.6), lwd=2)

  
#### Plot bootstrapped results by soil parent material ####
  
  # Read distribution results
  szFreqAnd <- read.table("SizeFreqResults/gapsAndesitesizedistbsfit.txt", header = T)
  szFreqBoh <- read.table("SizeFreqResults/gapsBohiosizedistbsfit.txt", header = T)
  szFreqMar <- read.table("SizeFreqResults/gapsCaimitoMarinesizedistbsfit.txt", header = T)
  szFreqVol <- read.table("SizeFreqResults/gapsCaimitoVolcanicsizedistbsfit.txt", header = T)
  
  
  
  # What distribution has the highest likelihood?
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
  # Andesite
  round(-szFreqAnd$weibloglike-max(c(-szFreqAnd$weibloglike,-szFreqAnd$powloglike,-szFreqAnd$exploglike)),0) 
  round(-szFreqAnd$powloglike-max(c(-szFreqAnd$weibloglike,-szFreqAnd$powloglike,-szFreqAnd$exploglike)),0) 
  round(-szFreqAnd$exploglike-max(c(-szFreqAnd$weibloglike,-szFreqAnd$powloglike,-szFreqAnd$exploglike)),0) 
  
  # Bohio
  round(-szFreqBoh$weibloglike-max(c(-szFreqBoh$weibloglike,-szFreqBoh$powloglike,-szFreqBoh$exploglike)),0) 
  round(-szFreqBoh$powloglike-max(c(-szFreqBoh$weibloglike,-szFreqBoh$powloglike,-szFreqBoh$exploglike)),0) 
  round(-szFreqBoh$exploglike-max(c(-szFreqBoh$weibloglike,-szFreqBoh$powloglike,-szFreqBoh$exploglike)),0) 
  
  # Caimito marine
  round(-szFreqMar$weibloglike-max(c(-szFreqMar$weibloglike,-szFreqMar$powloglike,-szFreqMar$exploglike)),0) 
  round(-szFreqMar$powloglike-max(c(-szFreqMar$weibloglike,-szFreqMar$powloglike,-szFreqMar$exploglike)),0) 
  round(-szFreqMar$exploglike-max(c(-szFreqMar$weibloglike,-szFreqMar$powloglike,-szFreqMar$exploglike)),0) 
  
  # Caimito volcanic
  round(-szFreqVol$weibloglike-max(c(-szFreqVol$weibloglike,-szFreqVol$powloglike,-szFreqVol$exploglike)),0) 
  round(-szFreqVol$powloglike-max(c(-szFreqVol$weibloglike,-szFreqVol$powloglike,-szFreqVol$exploglike)),0) 
  round(-szFreqVol$exploglike-max(c(-szFreqVol$weibloglike,-szFreqVol$powloglike,-szFreqVol$exploglike)),0) 

  # Weibull has the lowest -log likelihood in all categories
  round(szFreqAnd[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqBoh[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqMar[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqVol[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  
  round(szFreqAnd[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  round(szFreqBoh[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  round(szFreqMar[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  round(szFreqVol[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  
  
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
  
  gapSizesAnd <- makeGapVectors(gapAreaVector = andesiteGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledAnd18*nYr15to18 + areaSampledAnd20*nYr18to20)
  gapSizesBoh <- makeGapVectors(gapAreaVector = bohioGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledBoh18*nYr15to18 + areaSampledBoh20*nYr18to20)
  gapSizesMar <- makeGapVectors(gapAreaVector = caimitoMarineGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledMar18*nYr15to18 + areaSampledMar20*nYr18to20)
  gapSizesVol <- makeGapVectors(gapAreaVector = caimitoVolcanicGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledVol18*nYr15to18 + areaSampledVol20*nYr18to20)
  
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
  
  
  yValsWeibAnd <- getWeibEsts(szFreqResults = szFreqAnd, gapSp = andesiteGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledAnd18+areaSampledAnd20)
  yValsWeibBoh <- getWeibEsts(szFreqResults = szFreqBoh, gapSp = bohioGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledBoh18+areaSampledBoh20)
  yValsWeibMar <- getWeibEsts(szFreqResults = szFreqMar, gapSp = caimitoMarineGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledMar18+areaSampledMar20)
  yValsWeibVol <- getWeibEsts(szFreqResults = szFreqVol, gapSp = caimitoVolcanicGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledVol18+areaSampledVol20)
  
  # yValsExp18 <- getExpEsts(szFreqResults = szFreqOld, gapSp = gapsOldsp, xVals)/nYrOld/areaSampledOldtall
  # yValsExp20 <- getExpEsts(szFreqResults = szFreqSec, gapSp = gapsSecsp, xVals)/nYrSec/areaSampledSectall
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizesOld,gapSizesSec)) + c(1e-7,0)
  
  #Define colors
  colAnd <- wesanderson::wes_palette("Chevalier1",4)[2]
  colBoh <- wesanderson::wes_palette("Chevalier1",4)[4]
  colMar <- wesanderson::wes_palette("Chevalier1",4)[3]
  colVol <- wesanderson::wes_palette("Chevalier1",4)[1]
  
  plot(x = brksMids, y = gapSizesAnd,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(colAnd,0.6),
       log = logOption,
       pch=19,
       cex.axis = 1,
       ylab = NA,
       xlab = NA,
       main = "Soil parent material")
  points(x = brksMids, y = gapSizesBoh,
         col = adjustcolor(colBoh,0.6), pch=19)
  points(x = brksMids, y = gapSizesMar,
         col = adjustcolor(colMar,0.6), pch=19)
  points(x = brksMids, y = gapSizesVol,
         col = adjustcolor(colVol,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
         col=adjustcolor(c(colAnd,colBoh,colMar,colVol),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeibAnd, col = adjustcolor(colAnd,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibBoh, col = adjustcolor(colBoh,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibMar, col = adjustcolor(colMar,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibVol, col = adjustcolor(colVol,0.6), lwd=2)
  
  
#### Plot bootstrapped results by soil form ####
  
  # Read distribution results
  szFreqBro <- read.table("SizeFreqResults/fineLoamGapssizedistbsfit.txt", header = T)
  szFreqMot <- read.table("SizeFreqResults/heavyClayGapssizedistbsfit.txt", header = T)
  szFreqPal <- read.table("SizeFreqResults/swellingClayGapssizedistbsfit.txt", header = T)
  szFreqRed <- read.table("SizeFreqResults/lightClayGapssizedistbsfit.txt", header = T)
  
  
  # What distribution has the highest likelihood?
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
  # Brown fine loam
  round(-szFreqBro$weibloglike-max(c(-szFreqBro$weibloglike,-szFreqBro$powloglike,-szFreqBro$exploglike)),0) 
  round(-szFreqBro$powloglike-max(c(-szFreqBro$weibloglike,-szFreqBro$powloglike,-szFreqBro$exploglike)),0) 
  round(-szFreqBro$exploglike-max(c(-szFreqBro$weibloglike,-szFreqBro$powloglike,-szFreqBro$exploglike)),0) 
  
  # Mottled heavy clay
  round(-szFreqMot$weibloglike-max(c(-szFreqMot$weibloglike,-szFreqMot$powloglike,-szFreqMot$exploglike)),0) 
  round(-szFreqMot$powloglike-max(c(-szFreqMot$weibloglike,-szFreqMot$powloglike,-szFreqMot$exploglike)),0) 
  round(-szFreqMot$exploglike-max(c(-szFreqMot$weibloglike,-szFreqMot$powloglike,-szFreqMot$exploglike)),0) 
  
  # Pale swelling clay
  round(-szFreqPal$weibloglike-max(c(-szFreqPal$weibloglike,-szFreqPal$powloglike,-szFreqPal$exploglike)),0) 
  round(-szFreqPal$powloglike-max(c(-szFreqPal$weibloglike,-szFreqPal$powloglike,-szFreqPal$exploglike)),0) 
  round(-szFreqPal$exploglike-max(c(-szFreqPal$weibloglike,-szFreqPal$powloglike,-szFreqPal$exploglike)),0) 
  
  # Red light clay
  round(-szFreqRed$weibloglike-max(c(-szFreqRed$weibloglike,-szFreqRed$powloglike,-szFreqRed$exploglike)),0) 
  round(-szFreqRed$powloglike-max(c(-szFreqRed$weibloglike,-szFreqRed$powloglike,-szFreqRed$exploglike)),0) 
  round(-szFreqRed$exploglike-max(c(-szFreqRed$weibloglike,-szFreqRed$powloglike,-szFreqRed$exploglike)),0) 
  
  # Weibull has the lowest -log likelihood in all categories
  round(szFreqBro[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqMot[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqPal[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  round(szFreqRed[,c("weibpar1est","weibpar1lo1","weibpar1hi1")],3)
  
  round(szFreqBro[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  round(szFreqMot[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  round(szFreqPal[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  round(szFreqRed[,c("weibpar2est","weibpar2lo1","weibpar2hi1")],3)
  
  
  # Make plot of best-fit lines
  
  gapSizesBro <- makeGapVectors(gapAreaVector = fineLoamGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledBro18*nYr15to18 + areaSampledBro20*nYr18to20)
  gapSizesMot <- makeGapVectors(gapAreaVector = heavyClayGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledMot18*nYr15to18 + areaSampledMot20*nYr18to20)
  gapSizesPal <- makeGapVectors(gapAreaVector = swellingClayGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledPal18*nYr15to18 + areaSampledPal20*nYr18to20)
  gapSizesRed <- makeGapVectors(gapAreaVector = lightClayGaps$dbh,
                                minarea = brksMins,
                                maxarea = brksMaxs)/(areaSampledRed18*nYr15to18 + areaSampledRed20*nYr18to20)
  
  # Get correct fitted y-values for each year, adjusted for truncated distributions
  yValsWeibBro <- getWeibEsts(szFreqResults = szFreqBro, gapSp = fineLoamGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledBro18+areaSampledBro20)
  yValsWeibMot <- getWeibEsts(szFreqResults = szFreqMot, gapSp = heavyClayGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledMot18+areaSampledMot20)
  yValsWeibPal <- getWeibEsts(szFreqResults = szFreqPal, gapSp = swellingClayGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledPal18+areaSampledPal20)
  yValsWeibRed <- getWeibEsts(szFreqResults = szFreqRed, gapSp = lightClayGaps, xVals)/(0.5*(nYr15to18+nYr18to20))/(areaSampledRed18+areaSampledRed20)
  
  # yValsExp18 <- getExpEsts(szFreqResults = szFreqOld, gapSp = gapsOldsp, xVals)/nYrOld/areaSampledOldtall
  # yValsExp20 <- getExpEsts(szFreqResults = szFreqSec, gapSp = gapsSecsp, xVals)/nYrSec/areaSampledSectall
  
  logOption <- "xy"
  
  yRange <- range(c(gapSizesOld,gapSizesSec)) + c(1e-7,0)
  
  #Define colors
  colBro <- wesanderson::wes_palette("Rushmore1",5)[1]
  colMot <- wesanderson::wes_palette("Rushmore1",5)[4]
  colPal <- wesanderson::wes_palette("Rushmore1",5)[3]
  colRed <- wesanderson::wes_palette("Rushmore1",5)[5]
  
  plot(x = brksMids, y = gapSizesBro,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(colBro,0.6),
       log = logOption,
       pch=19,
       cex.axis = 1,
       ylab = NA,
       xlab = NA,
       main = "Soil form")
  points(x = brksMids, y = gapSizesMot,
         col = adjustcolor(colMot,0.6), pch=19)
  points(x = brksMids, y = gapSizesPal,
         col = adjustcolor(colPal,0.6), pch=19)
  points(x = brksMids, y = gapSizesRed,
         col = adjustcolor(colRed,0.6), pch=19)
  
  legend(x=30,y=0.0001,
         c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
         col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
         pch=19, cex=1,
         bty="n")
  
  lines(x = xVals, y = yValsWeibBro, col = adjustcolor(colBro,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibMot, col = adjustcolor(colMot,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibPal, col = adjustcolor(colPal,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibRed, col = adjustcolor(colRed,0.6), lwd=2)
  
  