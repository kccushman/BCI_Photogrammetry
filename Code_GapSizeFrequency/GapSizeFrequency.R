#### Packages ####

# Install the following packages, if needed

# install.packages("raster") # version 3.3-13 used
# install.packages("rgdal") # version 1.5-16 used
# install.packages("sp") # version 1.4-2 used
# install.packages("viridis") # version 0.5.1 used
# install.packages("wesanderson") # version 0.3.6 used

#### READ DATA ####
    
# Canopy height change rasters
  d15to18 <- raster::raster("Data_HeightRasters/dCHM15to18.tif")
  d18to20 <- raster::raster("Data_HeightRasters/dCHM18to20.tif")
  
# Canopy height change rasters where only likely gap values (> 10 m initially) are included  
  d15to18tall <- raster::raster("Data_HeightRasters/dCHM15to18tall.tif")
  d18to20tall <- raster::raster("Data_HeightRasters/dCHM18to20tall.tif")
  
# Gap rasters
  gaps15to18 <- raster::raster("Data_HeightRasters/newGaps15to18.tif")
  gaps18to20 <- raster::raster("Data_HeightRasters/newGaps18to20.tif")
  
# Gap polygons
  gaps15to18sp <- rgdal::readOGR("Data_GapShapefiles/gaps15to18_shapefile/gaps15to18sp.shp")
  gaps18to20sp <- rgdal::readOGR("Data_GapShapefiles/gaps18to20_shapefile/gaps18to20sp.shp")
  
# Block values for bootstrapping
  blockData <- read.csv("Data_Ancillary/bootstrapBlocks.csv")
  
# Forest age polygon
  age <- rgdal::readOGR("Data_Ancillary/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  age$AgeClass <- "Other"
  age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
  age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
  ageUse <- age[!(age$AgeClass=="Other"),]
  
# Soil type polygon  
  soil <- rgdal::readOGR("Data_Ancillary/BCI_Soils/BCI_Soils.shp")
  soil <- sp::spTransform(soil,raster::crs(d15to18tall))
  
# Define soil parent material and soil form from soil class (based on Baillie soil survey)
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
  
#### CALCULATE AREA SAMPLED EACH YEAR ####
  
# All area (in ha)
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
  
# Total area of new canopy disturbances (in ha)
  sum(c(gaps15to18sp$area,
           gaps18to20sp$area))/10000
  
# Mean and median disturbance sizes
  mean(c(gaps15to18sp$area,
        gaps18to20sp$area))
  median(c(gaps15to18sp$area,
         gaps18to20sp$area))
  
  
# Area sampled in old growth and secondary forest [> 10 m height]

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
  
  # Calculate the precise numnber of years sampled in each interval
    nYr15to18 <- as.numeric(as.Date("2018-06-07") - as.Date("2015-06-26"))/365
    nYr18to20 <- as.numeric(as.Date("2020-07-31") - as.Date("2018-06-07"))/365
  
  # Mean frequency of gaps per interval--percent of area
    round(mean(gapSummary15to18$percentGap)/nYr15to18,2)
    round(mean(gapSummary18to20$percentGap)/nYr18to20,2)
  
    round(quantile(gapSummary15to18$percentGap,probs = c(0.025,0.975))/nYr15to18,2)
    round(quantile(gapSummary18to20$percentGap,probs = c(0.025,0.975))/nYr18to20,2)

  
  # Mean frequency of gaps per interval--number of gaps
  round(mean(gapSummary15to18$gapsPerHa)/nYr15to18,2)
  round(mean(gapSummary18to20$gapsPerHa)/nYr18to20,2)
  
  round(quantile(gapSummary15to18$gapsPerHa,probs = c(0.025,0.975))/nYr15to18,2)
  round(quantile(gapSummary18to20$gapsPerHa,probs = c(0.025,0.975))/nYr18to20,2)

  
#### BOOTSTRAPPED SIZE FREQUENCY DISTRIBUTIONS: SET UP ####
  
  # Read in code to fit gap size frequency distribution (from Helene)
  # NOTE: this code was originally used to fit tree size frequency distributions
      # so "dbh" is used for gap area and "quadnum" is used for bootstrapping block
  # NOTE: the files below were written by H. Muller-Landau (mullerh@si.edu)
  source("Code_GapSizeFrequency/FitSizeDist_1.R")
  source("Code_GapSizeFrequency/FitSizeDist_2.R")
  source("Code_GapSizeFrequency/FitSizeDist_3.R")
  
  # find max gap size
  mxSz <- max(c(gaps15to18sp$area, gaps18to20sp$area))
  
  # find which gap sizes are NA (omit in code to speed it up)
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
  
#### BOOTSTRAPPED MEAN AND MEDIAN GAP SIZE ####
  
  # Create a function to calculate mean and median gap size, with CIs from bootstrapping
  gapSizeSimple <- function(alldata, nreps=1000, probLevels = c(0.025,0.975)){
    
    # set seed for repeatability
    set.seed(1)
    
    # create data frame for results
    longResults <- data.frame(rep = 1:nreps,
                              nGap = NA,
                              areaGap = NA,
                              meanGap = NA,
                              medGap = NA)
    
    
    Quads <- unique(alldata$quadnum)
    nQuads <- length(Quads)
    for(i in 1:nreps){
      
      # sample from bootstrapping units with replacement
      iSample <- sample(Quads,nQuads,replace=T)
      
      # make a vector of all gaps to summarize
      iGaps <- alldata[alldata$quadnum==iSample[1],"dbh"]
      for(j in 2:nQuads){
        iGaps <- c(iGaps,alldata[alldata$quadnum==iSample[j],"dbh"])
      }
      
      # summarize gaps
      longResults[i,"nGap"] <- length(iGaps)
      longResults[i,"areaGap"] <- sum(iGaps)
      longResults[i,"meanGap"] <- mean(iGaps)
      longResults[i,"medGap"] <- median(iGaps)
    }
    
    # summary bootstraps
    shortResults <- data.frame(metric = c("lo","mean","hi"),
                              nGap = NA,
                              areaGap = NA,
                              meanGap = NA,
                              medGap = NA)
    
    shortResults[shortResults$metric %in% c("lo","hi") ,c("nGap")] <- quantile(longResults$nGap,probLevels)
    shortResults[shortResults$metric %in% c("lo","hi") ,c("areaGap")] <- quantile(longResults$areaGap,probLevels)
    shortResults[shortResults$metric %in% c("lo","hi") ,c("meanGap")] <- quantile(longResults$meanGap,probLevels)
    shortResults[shortResults$metric %in% c("lo","hi") ,c("medGap")] <- quantile(longResults$medGap,probLevels)
    
    shortResults[shortResults$metric %in% c("mean") ,c("nGap")] <- mean(longResults$nGap)
    shortResults[shortResults$metric %in% c("mean") ,c("areaGap")] <- mean(longResults$areaGap)
    shortResults[shortResults$metric %in% c("mean") ,c("meanGap")] <- mean(longResults$meanGap)
    shortResults[shortResults$metric %in% c("mean") ,c("medGap")] <- mean(longResults$medGap)
    
    # structure outputs
    allResults <- list(shortResults, longResults, nQuads)
  }
  
  
  
#### Bootstrapped size frequency by year ####
  
  # simple
    sz15to18 <- gapSizeSimple(alldata = allData15to18, nreps=1000, probLevels = c(0.025,0.975))
    sz18to20 <- gapSizeSimple(alldata = allData18to20, nreps=1000, probLevels = c(0.025,0.975))
    
  # size frequency
  
    i <- 1
    datadir <- "Results_GapSizeFrequency/"
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
    datadir <- "Results_GapSizeFrequency/"
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
      
      # simple
        szOld <- gapSizeSimple(alldata = oldGrowthGaps, nreps=1000, probLevels = c(0.025,0.975))
        szSec <- gapSizeSimple(alldata = secondaryGaps, nreps=1000, probLevels = c(0.025,0.975))
      
      # size frequency
      i <- 1
      datadir <- "Results_GapSizeFrequency/"
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
      datadir <- "Results_GapSizeFrequency/"
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
      # simple
        szAnd <- gapSizeSimple(alldata = andesiteGaps, nreps=1000, probLevels = c(0.025,0.975))
        szBoh <- gapSizeSimple(alldata = bohioGaps, nreps=1000, probLevels = c(0.025,0.975))
        szVol <- gapSizeSimple(alldata = caimitoVolcanicGaps, nreps=1000, probLevels = c(0.025,0.975))
        szMar <- gapSizeSimple(alldata = caimitoMarineGaps, nreps=1000, probLevels = c(0.025,0.975))
      
      # size frequency
       i <- 1
      datadir <- "Results_GapSizeFrequency/"
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
      datadir <- "Results_GapSizeFrequency/"
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
      datadir <- "Results_GapSizeFrequency/"
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
      datadir <- "Results_GapSizeFrequency/"
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
      
    # simple
      szRed <- gapSizeSimple(alldata = lightClayGaps, nreps=1000, probLevels = c(0.025,0.975))
      szPal <- gapSizeSimple(alldata = swellingClayGaps, nreps=1000, probLevels = c(0.025,0.975))
      szMot <- gapSizeSimple(alldata = heavyClayGaps, nreps=1000, probLevels = c(0.025,0.975))
      szBro <- gapSizeSimple(alldata = fineLoamGaps, nreps=1000, probLevels = c(0.025,0.975))
      
    # size frequency
      i <- 1
      datadir <- "Results_GapSizeFrequency/"
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
      

      datadir <- "Results_GapSizeFrequency/"
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
    
      datadir <- "Results_GapSizeFrequency/"
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
      
      datadir <- "Results_GapSizeFrequency/"
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
      
#### Figure S9: Plot bootstrapped results by year ####
      
  col18 <- "blue"
  col20 <- "#d95f02"    
  
  # Read distribution results
  szFreq15to18 <- read.table("Results_GapSizeFrequency/gaps15to18sizedistbsfit.txt", header = T)
  szFreq18to20 <- read.table("Results_GapSizeFrequency/gaps18to20sizedistbsfit.txt", header = T)

  # Calculate AIC values for each model
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
  szFreq15to18$weibAIC <-  -2*(-1*szFreq15to18$weibloglike) + 2*2
  szFreq15to18$powAIC <-  -2*(-1*szFreq15to18$powloglike) + 2*1
  szFreq15to18$expAIC <-  -2*(-1*szFreq15to18$exploglike) + 2*1
  
  szFreq18to20$weibAIC <-  -2*(-1*szFreq18to20$weibloglike) + 2*2
  szFreq18to20$powAIC <-  -2*(-1*szFreq18to20$powloglike) + 2*1
  szFreq18to20$expAIC <-  -2*(-1*szFreq18to20$exploglike) + 2*1
  
  # Find the delta AIC
  
    # 2015 to 2017
    round(szFreq15to18$weibAIC-min(c(szFreq15to18$weibAIC,szFreq15to18$powAIC,szFreq15to18$expAIC)),0) 
    round(szFreq15to18$powAIC-min(c(szFreq15to18$weibAIC,szFreq15to18$powAIC,szFreq15to18$expAIC)),0) 
    round(szFreq15to18$expAIC-min(c(szFreq15to18$weibAIC,szFreq15to18$powAIC,szFreq15to18$expAIC)),0) 
    
    # 2018 to 2020
    round(szFreq18to20$weibAIC-min(c(szFreq18to20$weibAIC,szFreq18to20$powAIC,szFreq18to20$expAIC)),0) 
    round(szFreq18to20$powAIC-min(c(szFreq18to20$weibAIC,szFreq18to20$powAIC,szFreq18to20$expAIC)),0) 
    round(szFreq18to20$expAIC-min(c(szFreq18to20$weibAIC,szFreq18to20$powAIC,szFreq18to20$expAIC)),0) 
  
  # Get r-squared values for each
    
    # 2015 to 2017
    round(szFreq15to18$weibcatadjr2,3)
    round(szFreq15to18$powcatadjr2,3)
    round(szFreq15to18$expcatadjr2,3)
    
    # 2018 to 2020
    round(szFreq18to20$weibcatadjr2,3)
    round(szFreq18to20$powcatadjr2,3)
    round(szFreq18to20$expcatadjr2,3)
  
  # Weibull has the lowest -AIC in both intervals
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
       ylab = expression("Disturbance rate (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
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
  
  
  smooth18 <- ksmooth(x = gaps15to18sp$area,
                      y = -gaps15to18sp$htDrop,
                      bandwidth = 8,
                      x.points = quantile(gaps15to18sp$area, probs = seq(0,1,0.025)))
  smooth20 <- ksmooth(x = gaps18to20sp$area,
                      y = -gaps18to20sp$htDrop,
                      bandwidth = 8,
                      x.points = quantile(gaps18to20sp$area, probs = seq(0,1,0.025)))
  
  plot(x = gaps15to18sp$area,
       y = -gaps15to18sp$htDrop,
       log="xy",
       xlim=range(c(gaps15to18sp$area,gaps18to20sp$area)),
       ylim=-range(c(gaps15to18sp$htDrop,gaps18to20sp$htDrop)),
       ylab= "Average height decrease (m)",
       col = adjustcolor(col18,0.05),
       cex = 0.3,
       pch=19)
  text("b",
       x = min(c(gaps15to18sp$area,gaps18to20sp$area)),
       y = -max(c(gaps15to18sp$htDrop,gaps18to20sp$htDrop)))
  points(x = gaps18to20sp$area,
       y = -gaps18to20sp$htDrop,
       col = adjustcolor(col20,0.05),
       cex = 0.3,
       pch=19)
  lines(smooth18, col=col18,lwd=3)
  lines(smooth20, col=col20,lwd=3)
  
#### Figure S10a: Plot bootstrapped results by forest age ####
  
  # Read distribution results
  szFreqOld <- read.table("Results_GapSizeFrequency/gapsOldGrowthsizedistbsfit.txt", header = T)
  szFreqSec <- read.table("Results_GapSizeFrequency/gapsSecondarysizedistbsfit.txt", header = T)
  
  # Calculate AIC values for each model
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
  szFreqOld$weibAIC <-  -2*(-1*szFreqOld$weibloglike) + 2*2
  szFreqOld$powAIC <-  -2*(-1*szFreqOld$powloglike) + 2*1
  szFreqOld$expAIC <-  -2*(-1*szFreqOld$exploglike) + 2*1
  
  szFreqSec$weibAIC <-  -2*(-1*szFreqSec$weibloglike) + 2*2
  szFreqSec$powAIC <-  -2*(-1*szFreqSec$powloglike) + 2*1
  szFreqSec$expAIC <-  -2*(-1*szFreqSec$exploglike) + 2*1
  
  # Find the delta AIC
  
  round(szFreqOld$weibAIC-min(c(szFreqOld$weibAIC,szFreqOld$powAIC,szFreqOld$expAIC)),0) 
  round(szFreqOld$powAIC-min(c(szFreqOld$weibAIC,szFreqOld$powAIC,szFreqOld$expAIC)),0) 
  round(szFreqOld$expAIC-min(c(szFreqOld$weibAIC,szFreqOld$powAIC,szFreqOld$expAIC)),0) 
  
  round(szFreqSec$weibAIC-min(c(szFreqSec$weibAIC,szFreqSec$powAIC,szFreqSec$expAIC)),0) 
  round(szFreqSec$powAIC-min(c(szFreqSec$weibAIC,szFreqSec$powAIC,szFreqSec$expAIC)),0) 
  round(szFreqSec$expAIC-min(c(szFreqSec$weibAIC,szFreqSec$powAIC,szFreqSec$expAIC)),0) 
  
  # Weibull has the lowest AIC in both intervals
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
  
  par(las = 1, mar=c(3,2,2,1),oma=c(2,4,1,1), mfrow=c(1,3))
  plot(x = brksMids, y = gapSizesOld,
       xlim=range(brksMids),
       ylim=yRange,
       col = adjustcolor(colOld,0.6),
       log = logOption,
       pch=19,
       cex.axis = 1.2,
       ylab = NA,
       xlab = NA,
       main = "Forest age")
  mtext(expression("Disturbance area (m"^"2"~")"),
        side=1, outer=T)
  par(las=0)
  mtext(expression("Disturbance frequency (events m"^"-2"~"yr"^"-1"~"ha"^"-1"~")"),
        side=2, outer=T, line=1.5)
  par(las=1)
  points(x = brksMids, y = gapSizesSec,
         col = adjustcolor(colSec,0.6), pch=19)
  
  legend(x=20,y=0.000005,
         c("Old growth","Secondary"),
         col=adjustcolor(c(colOld,colSec),0.6),
         pch=19, cex=1.2,
         bty="n")
  
  lines(x = xVals, y = yValsWeibOld, col = adjustcolor(colOld,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibSec, col = adjustcolor(colSec,0.6), lwd=2)

  
#### Figure S10b: Plot bootstrapped results by soil parent material ####
  
  # Read distribution results
  szFreqAnd <- read.table("Results_GapSizeFrequency/gapsAndesitesizedistbsfit.txt", header = T)
  szFreqBoh <- read.table("Results_GapSizeFrequency/gapsBohiosizedistbsfit.txt", header = T)
  szFreqMar <- read.table("Results_GapSizeFrequency/gapsCaimitoMarinesizedistbsfit.txt", header = T)
  szFreqVol <- read.table("Results_GapSizeFrequency/gapsCaimitoVolcanicsizedistbsfit.txt", header = T)
  
  # Calculate AIC values for each model
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
  szFreqAnd$weibAIC <-  -2*(-1*szFreqAnd$weibloglike) + 2*2
  szFreqAnd$powAIC <-  -2*(-1*szFreqAnd$powloglike) + 2*1
  szFreqAnd$expAIC <-  -2*(-1*szFreqAnd$exploglike) + 2*1
  
  szFreqBoh$weibAIC <-  -2*(-1*szFreqBoh$weibloglike) + 2*2
  szFreqBoh$powAIC <-  -2*(-1*szFreqBoh$powloglike) + 2*1
  szFreqBoh$expAIC <-  -2*(-1*szFreqBoh$exploglike) + 2*1
  
  szFreqMar$weibAIC <-  -2*(-1*szFreqMar$weibloglike) + 2*2
  szFreqMar$powAIC <-  -2*(-1*szFreqMar$powloglike) + 2*1
  szFreqMar$expAIC <-  -2*(-1*szFreqMar$exploglike) + 2*1
  
  szFreqVol$weibAIC <-  -2*(-1*szFreqVol$weibloglike) + 2*2
  szFreqVol$powAIC <-  -2*(-1*szFreqVol$powloglike) + 2*1
  szFreqVol$expAIC <-  -2*(-1*szFreqVol$exploglike) + 2*1
  
  
  # Find the delta AIC
  
  round(szFreqAnd$weibAIC-min(c(szFreqAnd$weibAIC,szFreqAnd$powAIC,szFreqAnd$expAIC)),0) 
  round(szFreqAnd$powAIC-min(c(szFreqAnd$weibAIC,szFreqAnd$powAIC,szFreqAnd$expAIC)),0) 
  round(szFreqAnd$expAIC-min(c(szFreqAnd$weibAIC,szFreqAnd$powAIC,szFreqAnd$expAIC)),0) 
  
  round(szFreqBoh$weibAIC-min(c(szFreqBoh$weibAIC,szFreqBoh$powAIC,szFreqBoh$expAIC)),0) 
  round(szFreqBoh$powAIC-min(c(szFreqBoh$weibAIC,szFreqBoh$powAIC,szFreqBoh$expAIC)),0) 
  round(szFreqBoh$expAIC-min(c(szFreqBoh$weibAIC,szFreqBoh$powAIC,szFreqBoh$expAIC)),0) 
  
  round(szFreqMar$weibAIC-min(c(szFreqMar$weibAIC,szFreqMar$powAIC,szFreqMar$expAIC)),0) 
  round(szFreqMar$powAIC-min(c(szFreqMar$weibAIC,szFreqMar$powAIC,szFreqMar$expAIC)),0) 
  round(szFreqMar$expAIC-min(c(szFreqMar$weibAIC,szFreqMar$powAIC,szFreqMar$expAIC)),0) 
  
  round(szFreqVol$weibAIC-min(c(szFreqVol$weibAIC,szFreqVol$powAIC,szFreqVol$expAIC)),0) 
  round(szFreqVol$powAIC-min(c(szFreqVol$weibAIC,szFreqVol$powAIC,szFreqVol$expAIC)),0) 
  round(szFreqVol$expAIC-min(c(szFreqVol$weibAIC,szFreqVol$powAIC,szFreqVol$expAIC)),0) 
  
  # Weibull has the lowest AIC in all categories
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
       yaxt="n",
       col = adjustcolor(colAnd,0.6),
       log = logOption,
       pch=19,
       cex.axis = 1.2,
       ylab = NA,
       xlab = NA,
       main = "Soil parent material")
  points(x = brksMids, y = gapSizesBoh,
         col = adjustcolor(colBoh,0.6), pch=19)
  points(x = brksMids, y = gapSizesMar,
         col = adjustcolor(colMar,0.6), pch=19)
  points(x = brksMids, y = gapSizesVol,
         col = adjustcolor(colVol,0.6), pch=19)
  
  legend(x=20,y=0.000005,
         c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
         col=adjustcolor(c(colAnd,colBoh,colMar,colVol),0.6),
         pch=19, cex=1.2,
         bty="n")
  
  lines(x = xVals, y = yValsWeibAnd, col = adjustcolor(colAnd,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibBoh, col = adjustcolor(colBoh,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibMar, col = adjustcolor(colMar,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibVol, col = adjustcolor(colVol,0.6), lwd=2)
  
  
#### Figure S10c: Plot bootstrapped results by soil form ####
  
  # Read distribution results
  szFreqBro <- read.table("Results_GapSizeFrequency/fineLoamGapssizedistbsfit.txt", header = T)
  szFreqMot <- read.table("Results_GapSizeFrequency/heavyClayGapssizedistbsfit.txt", header = T)
  szFreqPal <- read.table("Results_GapSizeFrequency/swellingClayGapssizedistbsfit.txt", header = T)
  szFreqRed <- read.table("Results_GapSizeFrequency/lightClayGapssizedistbsfit.txt", header = T)
  
  
  # Calculate AIC values for each model
  # NOTE: the log likelihood returned by the size frequency code is the negative log likelihood, so multiply by -1 again
  
  szFreqBro$weibAIC <-  -2*(-1*szFreqBro$weibloglike) + 2*2
  szFreqBro$powAIC <-  -2*(-1*szFreqBro$powloglike) + 2*1
  szFreqBro$expAIC <-  -2*(-1*szFreqBro$exploglike) + 2*1
  
  szFreqMot$weibAIC <-  -2*(-1*szFreqMot$weibloglike) + 2*2
  szFreqMot$powAIC <-  -2*(-1*szFreqMot$powloglike) + 2*1
  szFreqMot$expAIC <-  -2*(-1*szFreqMot$exploglike) + 2*1
  
  szFreqPal$weibAIC <-  -2*(-1*szFreqPal$weibloglike) + 2*2
  szFreqPal$powAIC <-  -2*(-1*szFreqPal$powloglike) + 2*1
  szFreqPal$expAIC <-  -2*(-1*szFreqPal$exploglike) + 2*1
  
  szFreqRed$weibAIC <-  -2*(-1*szFreqRed$weibloglike) + 2*2
  szFreqRed$powAIC <-  -2*(-1*szFreqRed$powloglike) + 2*1
  szFreqRed$expAIC <-  -2*(-1*szFreqRed$exploglike) + 2*1
  
  
  # Find the delta AIC
  
  round(szFreqBro$weibAIC-min(c(szFreqBro$weibAIC,szFreqBro$powAIC,szFreqBro$expAIC)),0) 
  round(szFreqBro$powAIC-min(c(szFreqBro$weibAIC,szFreqBro$powAIC,szFreqBro$expAIC)),0) 
  round(szFreqBro$expAIC-min(c(szFreqBro$weibAIC,szFreqBro$powAIC,szFreqBro$expAIC)),0) 
  
  round(szFreqMot$weibAIC-min(c(szFreqMot$weibAIC,szFreqMot$powAIC,szFreqMot$expAIC)),0) 
  round(szFreqMot$powAIC-min(c(szFreqMot$weibAIC,szFreqMot$powAIC,szFreqMot$expAIC)),0) 
  round(szFreqMot$expAIC-min(c(szFreqMot$weibAIC,szFreqMot$powAIC,szFreqMot$expAIC)),0) 
  
  round(szFreqPal$weibAIC-min(c(szFreqPal$weibAIC,szFreqPal$powAIC,szFreqPal$expAIC)),0) 
  round(szFreqPal$powAIC-min(c(szFreqPal$weibAIC,szFreqPal$powAIC,szFreqPal$expAIC)),0) 
  round(szFreqPal$expAIC-min(c(szFreqPal$weibAIC,szFreqPal$powAIC,szFreqPal$expAIC)),0) 
  
  round(szFreqRed$weibAIC-min(c(szFreqRed$weibAIC,szFreqRed$powAIC,szFreqRed$expAIC)),0) 
  round(szFreqRed$powAIC-min(c(szFreqRed$weibAIC,szFreqRed$powAIC,szFreqRed$expAIC)),0) 
  round(szFreqRed$expAIC-min(c(szFreqRed$weibAIC,szFreqRed$powAIC,szFreqRed$expAIC)),0) 
  
  
  # Weibull has the lowest AIC in all categories
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
       yaxt="n",
       pch=19,
       cex.axis = 1.2,
       ylab = NA,
       xlab = NA,
       main = "Soil form")
  points(x = brksMids, y = gapSizesMot,
         col = adjustcolor(colMot,0.6), pch=19)
  points(x = brksMids, y = gapSizesPal,
         col = adjustcolor(colPal,0.6), pch=19)
  points(x = brksMids, y = gapSizesRed,
         col = adjustcolor(colRed,0.6), pch=19)
  
  legend(x=20,y=0.000005,
         c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
         col=adjustcolor(c(colBro,colMot,colPal,colRed),0.6),
         pch=19, cex=1.2,
         bty="n")
  
  lines(x = xVals, y = yValsWeibBro, col = adjustcolor(colBro,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibMot, col = adjustcolor(colMot,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibPal, col = adjustcolor(colPal,0.6), lwd=2)
  lines(x = xVals, y = yValsWeibRed, col = adjustcolor(colRed,0.6), lwd=2)
  
  


  
#### Figure S11: Make new stacked barplots ####  

  gapBins <- c(25,50,100,200,400,800,1600,3200,6400,19200)
  
  
  # Make data frame
  gapProp <- data.frame(minSz = gapBins[1:(length(gapBins)-1)],
                        maxSz = gapBins[2:(length(gapBins))],
                        yr15to18 = NA,
                        yr18to20 = NA,
                        oldGrowth = NA,
                        secondary = NA,
                        parAnd = NA,
                        parBoh = NA,
                        parMar = NA,
                        parVol = NA,
                        formBro = NA,
                        formMot = NA,
                        formPal = NA,
                        formRed = NA)
  
  for(i in 1:(length(gapBins)-1)){
    gapProp$yr15to18[i] <- sum(allData15to18[allData15to18$dbh>=gapProp$minSz[i] & allData15to18$dbh<gapProp$maxSz[i],"dbh"])/sum(allData15to18$dbh)
    gapProp$yr18to20[i] <- sum(allData18to20[allData18to20$dbh>=gapProp$minSz[i] & allData18to20$dbh<gapProp$maxSz[i],"dbh"])/sum(allData18to20$dbh)
    # gapProp$yr15to18[i] <- sum(allData15to18[allData15to18$dbh>=gapProp$minSz[i] & allData15to18$dbh<gapProp$maxSz[i],"dbh"])/(areaSampled15to18tall)/10000
    # gapProp$yr18to20[i] <- sum(allData18to20[allData18to20$dbh>=gapProp$minSz[i] & allData18to20$dbh<gapProp$maxSz[i],"dbh"])/(areaSampled18to20tall)/10000
    
    
    gapProp$oldGrowth[i] <- sum(oldGrowthGaps[oldGrowthGaps$dbh>=gapProp$minSz[i] & oldGrowthGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(oldGrowthGaps$dbh)
    gapProp$secondary[i] <- sum(secondaryGaps[secondaryGaps$dbh>=gapProp$minSz[i] & secondaryGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(secondaryGaps$dbh)
    # gapProp$oldGrowth[i] <- sum(oldGrowthGaps[oldGrowthGaps$dbh>=gapProp$minSz[i] & oldGrowthGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledOld18+areaSampledOld20)/10000
    # gapProp$secondary[i] <- sum(secondaryGaps[secondaryGaps$dbh>=gapProp$minSz[i] & secondaryGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledSec18+areaSampledSec20)/10000
    
    
    gapProp$parAnd[i] <- sum(andesiteGaps[andesiteGaps$dbh>=gapProp$minSz[i] & andesiteGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(andesiteGaps$dbh)
    gapProp$parBoh[i] <- sum(bohioGaps[bohioGaps$dbh>=gapProp$minSz[i] & bohioGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(bohioGaps$dbh)
    gapProp$parMar[i] <- sum(caimitoMarineGaps[caimitoMarineGaps$dbh>=gapProp$minSz[i] & caimitoMarineGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(caimitoMarineGaps$dbh)
    gapProp$parVol[i] <- sum(caimitoVolcanicGaps[caimitoVolcanicGaps$dbh>=gapProp$minSz[i] & caimitoVolcanicGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(caimitoVolcanicGaps$dbh)
    # gapProp$parAnd[i] <- sum(andesiteGaps[andesiteGaps$dbh>=gapProp$minSz[i] & andesiteGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledAnd18+areaSampledAnd20)/10000
    # gapProp$parBoh[i] <- sum(bohioGaps[bohioGaps$dbh>=gapProp$minSz[i] & bohioGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledBoh18+areaSampledBoh20)/10000
    # gapProp$parMar[i] <- sum(caimitoMarineGaps[caimitoMarineGaps$dbh>=gapProp$minSz[i] & caimitoMarineGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledMar18+areaSampledMar20)/10000
    # gapProp$parVol[i] <- sum(caimitoVolcanicGaps[caimitoVolcanicGaps$dbh>=gapProp$minSz[i] & caimitoVolcanicGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledVol18+areaSampledVol20)/10000
    
    
    gapProp$formBro[i] <- sum(fineLoamGaps[fineLoamGaps$dbh>=gapProp$minSz[i] & fineLoamGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(fineLoamGaps$dbh)
    gapProp$formMot[i] <- sum(heavyClayGaps[heavyClayGaps$dbh>=gapProp$minSz[i] & heavyClayGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(heavyClayGaps$dbh)
    gapProp$formPal[i] <- sum(swellingClayGaps[swellingClayGaps$dbh>=gapProp$minSz[i] & swellingClayGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(swellingClayGaps$dbh)
    gapProp$formRed[i] <- sum(lightClayGaps[lightClayGaps$dbh>=gapProp$minSz[i] & lightClayGaps$dbh<gapProp$maxSz[i],"dbh"])/sum(lightClayGaps$dbh)
    # gapProp$formBro[i] <- sum(fineLoamGaps[fineLoamGaps$dbh>=gapProp$minSz[i] & fineLoamGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledBro18+areaSampledBro20)/10000
    # gapProp$formMot[i] <- sum(heavyClayGaps[heavyClayGaps$dbh>=gapProp$minSz[i] & heavyClayGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledMot18+areaSampledMot20)/10000
    # gapProp$formPal[i] <- sum(swellingClayGaps[swellingClayGaps$dbh>=gapProp$minSz[i] & swellingClayGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledPal18+areaSampledPal20)/10000
    # gapProp$formRed[i] <- sum(lightClayGaps[lightClayGaps$dbh>=gapProp$minSz[i] & lightClayGaps$dbh<gapProp$maxSz[i],"dbh"])/(areaSampledRed18+areaSampledRed20)/10000
  }
  
  # Make data frame
  gapPropN <- data.frame(minSz = gapBins[1:(length(gapBins)-1)],
                        maxSz = gapBins[2:(length(gapBins))],
                        yr15to18 = NA,
                        yr18to20 = NA,
                        oldGrowth = NA,
                        secondary = NA,
                        parAnd = NA,
                        parBoh = NA,
                        parMar = NA,
                        parVol = NA,
                        formBro = NA,
                        formMot = NA,
                        formPal = NA,
                        formRed = NA)
  
  for(i in 1:(length(gapBins)-1)){
    gapPropN$yr15to18[i] <- length(allData15to18[allData15to18$dbh>=gapPropN$minSz[i] & allData15to18$dbh<gapPropN$maxSz[i],"dbh"])/length(allData15to18$dbh)
    gapPropN$yr18to20[i] <- length(allData18to20[allData18to20$dbh>=gapPropN$minSz[i] & allData18to20$dbh<gapPropN$maxSz[i],"dbh"])/length(allData18to20$dbh)
    # gapPropN$yr15to18[i] <- length(allData15to18[allData15to18$dbh>=gapPropN$minSz[i] & allData15to18$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampled15to18tall)/10000
    # gapPropN$yr18to20[i] <- length(allData18to20[allData18to20$dbh>=gapPropN$minSz[i] & allData18to20$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampled18to20tall)/10000
    
    
    gapPropN$oldGrowth[i] <- length(oldGrowthGaps[oldGrowthGaps$dbh>=gapPropN$minSz[i] & oldGrowthGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(oldGrowthGaps$dbh)
    gapPropN$secondary[i] <- length(secondaryGaps[secondaryGaps$dbh>=gapPropN$minSz[i] & secondaryGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(secondaryGaps$dbh)
    # gapPropN$oldGrowth[i] <- length(oldGrowthGaps[oldGrowthGaps$dbh>=gapPropN$minSz[i] & oldGrowthGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledOld18+areaSampledOld20)/10000
    # gapPropN$secondary[i] <- length(secondaryGaps[secondaryGaps$dbh>=gapPropN$minSz[i] & secondaryGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledSec18+areaSampledSec20)/10000
    
    
    gapPropN$parAnd[i] <- length(andesiteGaps[andesiteGaps$dbh>=gapPropN$minSz[i] & andesiteGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(andesiteGaps$dbh)
    gapPropN$parBoh[i] <- length(bohioGaps[bohioGaps$dbh>=gapPropN$minSz[i] & bohioGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(bohioGaps$dbh)
    gapPropN$parMar[i] <- length(caimitoMarineGaps[caimitoMarineGaps$dbh>=gapPropN$minSz[i] & caimitoMarineGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(caimitoMarineGaps$dbh)
    gapPropN$parVol[i] <- length(caimitoVolcanicGaps[caimitoVolcanicGaps$dbh>=gapPropN$minSz[i] & caimitoVolcanicGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(caimitoVolcanicGaps$dbh)
    # gapPropN$parAnd[i] <- length(andesiteGaps[andesiteGaps$dbh>=gapPropN$minSz[i] & andesiteGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledAnd18+areaSampledAnd20)/10000
    # gapPropN$parBoh[i] <- length(bohioGaps[bohioGaps$dbh>=gapPropN$minSz[i] & bohioGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledBoh18+areaSampledBoh20)/10000
    # gapPropN$parMar[i] <- length(caimitoMarineGaps[caimitoMarineGaps$dbh>=gapPropN$minSz[i] & caimitoMarineGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledMar18+areaSampledMar20)/10000
    # gapPropN$parVol[i] <- length(caimitoVolcanicGaps[caimitoVolcanicGaps$dbh>=gapPropN$minSz[i] & caimitoVolcanicGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledVol18+areaSampledVol20)/10000
    
    
    gapPropN$formBro[i] <- length(fineLoamGaps[fineLoamGaps$dbh>=gapPropN$minSz[i] & fineLoamGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(fineLoamGaps$dbh)
    gapPropN$formMot[i] <- length(heavyClayGaps[heavyClayGaps$dbh>=gapPropN$minSz[i] & heavyClayGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(heavyClayGaps$dbh)
    gapPropN$formPal[i] <- length(swellingClayGaps[swellingClayGaps$dbh>=gapPropN$minSz[i] & swellingClayGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(swellingClayGaps$dbh)
    gapPropN$formRed[i] <- length(lightClayGaps[lightClayGaps$dbh>=gapPropN$minSz[i] & lightClayGaps$dbh<gapPropN$maxSz[i],"dbh"])/length(lightClayGaps$dbh)
    # gapPropN$formBro[i] <- length(fineLoamGaps[fineLoamGaps$dbh>=gapPropN$minSz[i] & fineLoamGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledBro18+areaSampledBro20)/10000
    # gapPropN$formMot[i] <- length(heavyClayGaps[heavyClayGaps$dbh>=gapPropN$minSz[i] & heavyClayGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledMot18+areaSampledMot20)/10000
    # gapPropN$formPal[i] <- length(swellingClayGaps[swellingClayGaps$dbh>=gapPropN$minSz[i] & swellingClayGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledPal18+areaSampledPal20)/10000
    # gapPropN$formRed[i] <- length(lightClayGaps[lightClayGaps$dbh>=gapPropN$minSz[i] & lightClayGaps$dbh<gapPropN$maxSz[i],"dbh"])/(areaSampledRed18+areaSampledRed20)/10000
  }
  
  # Make plots
    gapPropYear <- as.matrix(gapProp[,c("yr15to18","yr18to20")])
    gapPropAge <- as.matrix(gapProp[,c("oldGrowth","secondary")])
    gapPropParent <- as.matrix(gapProp[,c("parAnd","parBoh","parMar","parVol")])
    gapPropForm <- as.matrix(gapProp[,c("formBro","formMot","formPal","formRed")])
    
    gapPropYearN <- as.matrix(gapPropN[,c("yr15to18","yr18to20")])
    gapPropAgeN <- as.matrix(gapPropN[,c("oldGrowth","secondary")])
    gapPropParentN <- as.matrix(gapPropN[,c("parAnd","parBoh","parMar","parVol")])
    gapPropFormN <- as.matrix(gapPropN[,c("formBro","formMot","formPal","formRed")])
    
    yRange <- 0.02+max(apply(gapPropYear,2,"sum"),apply(gapPropAge,2,"sum"),apply(gapPropParent,2,"sum"),apply(gapPropForm,2,"sum"))
    gapCols <- viridis::viridis(nrow(gapProp))
    
    

    par(mar=c(3,4,1,1), las=1)
    barplot(gapPropYear,
            names = c("2015-2018","2018-2020"),
            ylim=c(0,yRange),
            col=gapCols,
            ylab="Proportion of area")
    
    legend(x=0.15,y=0.085,col=gapCols,pch=19,bty="n",horiz=F,xpd=T,ncol=5,
           paste0(gapBins[1:(length(gapBins)-1)],"-",gapBins[2:(length(gapBins))]),
           title="Gap size (m2)")

    par(mfrow=c(2,3), mar=c(0,1,1,0), oma=c(3,3,5,1))
    
    barplot(gapPropAgeN,
            ylim=c(0,yRange),
            col=gapCols,
            ylab="Proportion of disturbances",
            names=rep("",2),
            cex.axis=1.5)

    barplot(gapPropParentN,
            ylim=c(0,yRange),
            names=rep("",4),
            yaxt="n", ylab = NA,
            col=gapCols)
    
    legend(x=-4,y=1.45,col=gapCols,pch=19,bty="n",horiz=F,
           xpd=NA,
           ncol=5,
           c(paste0(gapBins[1:(length(gapBins)-2)],"-",gapBins[2:(length(gapBins)-1)]),"6400+"),
           cex = 1.5)
    
    barplot(gapPropFormN,
            ylim=c(0,yRange),
            names=rep("",4),
            yaxt="n", ylab = NA,
            col=gapCols)
    
    barplot(gapPropAge,
            names = c("Old growth","Secondary"),
            ylim=c(0,yRange),
            col=gapCols,
            ylab="Proportion of area",
            cex.axis=1.5)

    barplot(gapPropParent,
            names = c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
            ylim=c(0,yRange),
            yaxt="n", ylab = NA,
            col=gapCols)

    barplot(gapPropForm,
            names = c("Brown fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
            ylim=c(0,yRange),
            yaxt="n", ylab = NA,
            col=gapCols)
    
  