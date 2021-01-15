#### Read data and define forest age polygon ####

  # Read grid info
    gridInfo <- read.csv("gridInfo_QAQC.csv")
  
  # Read in forest age
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    ageUse <- age[!(age$TYPE=="Clearings"),]
  
  # Read canopy height change rasters
    dchm17to18 <- raster::raster("dCHM17to18.tif")
    dchm18to19 <- raster::raster("dCHM18to19.tif")
    dchm19to20 <- raster::raster("dCHM19to20.tif")
  
    allLayerList <- list(dchm17to18,dchm18to19,dchm19to20)
  
  # Read new gap rasters
    gaps17to18 <- raster::raster("newGaps17to18.tif")
    gaps18to19 <- raster::raster("newGaps18to19.tif")
    gaps19to20 <- raster::raster("newGaps19to20.tif")
      
    gapLayerList <- list(gaps17to18, gaps18to19, gaps19to20)
    
  # Read in new gap polygons
    gaps18sp <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")
    gaps19sp <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")
    gaps20sp <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")
    
  # Read polygon buffer 25 m inland from lake
    buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
    buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
    buffer <- raster::intersect(buffer, ageUse)
  
    # Make new age class in buffer polygon to combine all secondary forests
      buffer$AgeClass <- "Other"
      buffer$AgeClass[buffer$Mascaro_Co == "> 400"] <- "OldGrowth"
      buffer$AgeClass[buffer$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
      # Separate into two different polygon data frames  
        oldGrowth <- raster::aggregate(buffer[buffer$AgeClass=="OldGrowth",])
        secondary <- raster::aggregate(buffer[buffer$AgeClass=="Secondary",])

#### Find % area in new gaps by forest age for each interval ####
        
  # 2017-2018
        
    # Old growth    
      allOld18 <- raster::crop(dchm17to18,oldGrowth)
      gapsOld18 <- raster::crop(gaps17to18,oldGrowth)
      values_allOld18 <- raster::getValues(allOld18)
      values_gapsOld18 <- raster::getValues(gapsOld18)
      round(100*length(values_gapsOld18[!is.na(values_gapsOld18)])/length(values_allOld18[!is.na(values_allOld18)]),2)
      
    # Secondary    
      allSec18 <- raster::crop(dchm17to18,secondary)
      gapsSec18 <- raster::crop(gaps17to18,secondary)
      values_allSec18 <- raster::getValues(allSec18)
      values_gapsSec18 <- raster::getValues(gapsSec18)
      round(100*length(values_gapsSec18[!is.na(values_gapsSec18)])/length(values_allSec18[!is.na(values_allSec18)]),2)
      
  
  # 2018-2019
      
    # Old growth    
      allOld19 <- raster::crop(dchm18to19,oldGrowth)
      gapsOld19 <- raster::crop(gaps18to19,oldGrowth)
      values_allOld19 <- raster::getValues(allOld19)
      values_gapsOld19 <- raster::getValues(gapsOld19)
      round(100*length(values_gapsOld19[!is.na(values_gapsOld19)])/length(values_allOld19[!is.na(values_allOld19)]),2)
      
    # Secondary    
      allSec19 <- raster::crop(dchm18to19,secondary)
      gapsSec19 <- raster::crop(gaps18to19,secondary)
      values_allSec19 <- raster::getValues(allSec19)
      values_gapsSec19 <- raster::getValues(gapsSec19)
      round(100*length(values_gapsSec19[!is.na(values_gapsSec19)])/length(values_allSec19[!is.na(values_allSec19)]),2)
  
  # 2019-2020

    # Old growth    
      allOld20 <- raster::crop(dchm19to20,oldGrowth)
      gapsOld20 <- raster::crop(gaps19to20,oldGrowth)
      values_allOld20 <- raster::getValues(allOld20)
      values_gapsOld20 <- raster::getValues(gapsOld20)
      round(100*length(values_gapsOld20[!is.na(values_gapsOld20)])/length(values_allOld20[!is.na(values_allOld20)]),2)
      
    # Secondary    
      allSec20 <- raster::crop(dchm19to20,secondary)
      gapsSec20 <- raster::crop(gaps19to20,secondary)
      values_allSec20 <- raster::getValues(allSec20)
      values_gapsSec20 <- raster::getValues(gapsSec20)
      round(100*length(values_gapsSec20[!is.na(values_gapsSec20)])/length(values_allSec20[!is.na(values_allSec20)]),2)
  
  

  
#### Get gap size frequency distributions ####
      
  # Define a Metropolis-Hastings MCMC algorithm to esimates the posterior distribution of lambda    
  # Based on recommendations from: https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
  
  # Data: vector of individual gap sizes (in # of pixels)
  
  # Likelihood function: 
    dpl <- function(param, data){
      sum(log(data^-param/VGAM::zeta(param)))
    }
    
  # Define a prior for lambda
    prior <- function(param){
      lambdaPrior <- dunif(param, min=1.001, max=5, log=T)
      #lambdaPrior <- dgamma(param-1, 0.001, 0.001, log=T)
      return(lambdaPrior)
    }
  
  # Define the posterior, the sum of the likelihood function and the prior (because log-transformed)
    posterior <- function(data, param){
      return(dpl(param,data) + prior(param))
    }
  
  # Define proposal function (adjust to get acceptance rate ~ 20%)
    proposalFunction <- function(param){
      return(rnorm(1,  mean=param, sd=0.11))
      #return(rgamma(1, shape=2500, rate=2500/(param-1))+1)
    }
  
  # Define MCMC routine
    runMCMC <- function(startvalue, iterations, data){
      
      chain = array(dim = c(iterations+1,1))
      chain[1] = startvalue
      accept = array(dim= c(iterations,1))
      
      for (i in 1:iterations){
        proposal = proposalFunction(chain[i])
        
        # Proposal must be greater than 1
        while(proposal < 1){
          proposal = proposalFunction(chain[i])
        }
        
        q1 <- dgamma(chain[i]-1, shape=2500, rate=2500/(proposal-1), log=T)
        q2 <- dgamma((proposal-1), shape=2500, rate=2500/(chain[i]-1), log=T)
        
        probabNumerator = exp(posterior(data, param = proposal) - posterior(data, param = chain[i]))
        probabDenominator = exp(q2-q1)
        #probab <- probabNumerator/probabDenominator
        probab <- probabNumerator
        
        if (runif(1) < probab){
          chain[i+1,] = proposal
          accept[i] = 1
        }else{
          chain[i+1,] = chain[i,]
          accept[i] = 0
        }
      }
      
      return(list(chain,accept))
      
    }
    
  # Get gap size vectors in each forest type
    areaPerimThresh <- 0.6
    
    # 2017 - 2018
      gaps18spOld <- gaps18sp[!is.na(sp::over(gaps18sp, oldGrowth)),]
      gapSz18Old <- gaps18spOld$area[gaps18spOld$ratio > areaPerimThresh]
      gaps18spSec <- gaps18sp[!is.na(sp::over(gaps18sp, secondary)),]
      gapSz18Sec <- gaps18spSec$area[gaps18spSec$ratio > areaPerimThresh]
      
    # 2018 - 2019
      gaps19spOld <- gaps19sp[!is.na(sp::over(gaps19sp, oldGrowth)),]
      gapSz19Old <- gaps19spOld$area[gaps19spOld$ratio > areaPerimThresh]
      gaps19spSec <- gaps19sp[!is.na(sp::over(gaps19sp, secondary)),]
      gapSz19Sec <- gaps19spSec$area[gaps19spSec$ratio > areaPerimThresh]
      
      # 2019 - 2020
      gaps20spOld <- gaps20sp[!is.na(sp::over(gaps20sp, oldGrowth)),]
      gapSz20Old <- gaps20spOld$area[gaps20spOld$ratio > areaPerimThresh]
      gaps20spSec <- gaps20sp[!is.na(sp::over(gaps20sp, secondary)),]
      gapSz20Sec <- gaps20spSec$area[gaps20spSec$ratio > areaPerimThresh]
      
      
  # Run MCMC routine to estimate lambda  
    # Set seed so "random" results are repeatable
      set.seed(1)
    # Set number of simulations  
      sims <- 100000
    # Run MCMC function for each year
      # 2017 - 2018
        lambda18Old <- runMCMC(startvalue = optimize(dpl, 
                                                     data=gapSz18Old, 
                                                     lower = 1.0001, upper = 20, 
                                                     maximum = T)$maximum, 
                               iterations = sims, 
                               data = gapSz18Old)
        
        lambda18Sec <- runMCMC(startvalue = optimize(dpl, 
                                                     data=gapSz18Sec, 
                                                     lower = 1.0001, upper = 20, 
                                                     maximum = T)$maximum, 
                               iterations = sims, 
                               data = gapSz18Sec)
        
      # 2018 - 2019
        lambda19Old <- runMCMC(startvalue = optimize(dpl, 
                                                     data=gapSz19Old, 
                                                     lower = 1.0001, upper = 20, 
                                                     maximum = T)$maximum, 
                               iterations = sims, 
                               data = gapSz19Old)
        
        lambda19Sec <- runMCMC(startvalue = optimize(dpl, 
                                                     data=gapSz19Sec, 
                                                     lower = 1.0001, upper = 20, 
                                                     maximum = T)$maximum, 
                               iterations = sims, 
                               data = gapSz19Sec)
        
      # 2019 - 2020
        lambda20Old <- runMCMC(startvalue = optimize(dpl, 
                                                     data=gapSz20Old, 
                                                     lower = 1.0001, upper = 20, 
                                                     maximum = T)$maximum, 
                               iterations = sims, 
                               data = gapSz20Old)
        
        lambda20Sec <- runMCMC(startvalue = optimize(dpl, 
                                                     data=gapSz20Sec, 
                                                     lower = 1.0001, upper = 20, 
                                                     maximum = T)$maximum, 
                               iterations = sims, 
                               data = gapSz20Sec)
        
  # Summarize results for each year
    # Define parameters for sampling chain
      burnIn <- 5001
      skipN <- 25
        
    # 2017-2018  
      lambdaMin18Old <- round(quantile(lambda18Old[[1]][seq(burnIn,sims,skipN)], 0.025),3)
      lambdaMed18Old <- round(quantile(lambda18Old[[1]][seq(burnIn,sims,skipN)], 0.50),3)
      lambdaMax18Old <- round(quantile(lambda18Old[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      lambdaMin18Sec <- round(quantile(lambda18Sec[[1]][seq(burnIn,sims,skipN)], 0.025),3)
      lambdaMed18Sec <- round(quantile(lambda18Sec[[1]][seq(burnIn,sims,skipN)], 0.50),3)
      lambdaMax18Sec <- round(quantile(lambda18Sec[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      
    # 2018-2019  
      lambdaMin19Old <- round(quantile(lambda19Old[[1]][seq(burnIn,sims,skipN)], 0.025),3)
      lambdaMed19Old <- round(quantile(lambda19Old[[1]][seq(burnIn,sims,skipN)], 0.50),3)
      lambdaMax19Old <- round(quantile(lambda19Old[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      lambdaMin19Sec <- round(quantile(lambda19Sec[[1]][seq(burnIn,sims,skipN)], 0.025),3)
      lambdaMed19Sec <- round(quantile(lambda19Sec[[1]][seq(burnIn,sims,skipN)], 0.50),3)
      lambdaMax19Sec <- round(quantile(lambda19Sec[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      
    # 2019-2020  
      lambdaMin20Old <- round(quantile(lambda20Old[[1]][seq(burnIn,sims,skipN)], 0.025),3)
      lambdaMed20Old <- round(quantile(lambda20Old[[1]][seq(burnIn,sims,skipN)], 0.50),3)
      lambdaMax20Old <- round(quantile(lambda20Old[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      lambdaMin20Sec <- round(quantile(lambda20Sec[[1]][seq(burnIn,sims,skipN)], 0.025),3)
      lambdaMed20Sec <- round(quantile(lambda20Sec[[1]][seq(burnIn,sims,skipN)], 0.50),3)
      lambdaMax20Sec <- round(quantile(lambda20Sec[[1]][seq(burnIn,sims,skipN)], 0.975),3)