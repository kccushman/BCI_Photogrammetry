# Do preliminary gap analyses

#### READ RAW DATA ####

# Read in canopy height rasters:
  dsm17 <- raster::raster("DSM_2017_corrected.tif")
  dsm18 <- raster::raster("DSM_2018_corrected.tif")
  dsm19 <- raster::raster("DSM_2019_corrected.tif")
  dsm20 <- raster::raster("DSM_2020_corrected.tif")
  
# Read in forest age (used to exclude recent clearings)
  age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  ageUse <- age[!(age$TYPE=="Clearings"),]

# Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    
# Read in BCI DEM
  dem <- raster::raster("D:/BCI_Spatial/BCI_Topo/LidarDEM_BCI.tif")
  dem <- raster::crop(dem, raster::extent(ageUse))
  dem <- raster::resample(dem,dsm17)
  
# Read in QAQC data
  qaqc <- read.csv("gridInfo_QAQC.csv")
 
#### PROCESS DATA ####
  
# Remove raster areas outside BCI perimeter (exclude within 25 m of lake)
  dsm17 <- raster::mask(dsm17, buffer)  
  dsm18 <- raster::mask(dsm18, buffer)  
  dsm19 <- raster::mask(dsm19, buffer)  
  dsm20 <- raster::mask(dsm20, buffer)

# Remove raster areas in clearings
  dsm17 <- raster::mask(dsm17, ageUse)  
  dsm18 <- raster::mask(dsm18, ageUse)  
  dsm19 <- raster::mask(dsm19, ageUse)  
  dsm20 <- raster::mask(dsm20, ageUse)

# Crop to ensure each raster has same extent
  dsm17 <- raster::crop(dsm17, raster::extent(ageUse))  
  dsm18 <- raster::crop(dsm18, raster::extent(ageUse))  
  dsm19 <- raster::crop(dsm19, raster::extent(ageUse))  
  dsm20 <- raster::crop(dsm20, raster::extent(ageUse))    
  
# Subtract ground elevation
  chm17 <- dsm17-dem
  chm18 <- dsm18-dem
  chm19 <- dsm19-dem
  chm20 <- dsm20-dem
  
# Set areas that fail QAQC to "NA"
  for(i in 1:dim(qaqc)[1]){
    
    # make extent object for current tile
      x1 <- qaqc[i, "xmin"] 
      x2 <- qaqc[i, "xmax"]
      y1 <- qaqc[i, "ymin"] 
      y2 <- qaqc[i, "ymax"] 
      extent_i <- raster::extent(c(x1,x2,y1,y2))
      extent_i <- as(extent_i, 'SpatialPolygons')
    
    # set values within failed tiles to NA  
    if(qaqc$Use17[i]==F){
      chm17 <- raster::mask(chm17, extent_i, inverse = T)
    }
      
    if(qaqc$Use18[i]==F){
      chm18 <- raster::mask(chm18, extent_i, inverse = T)
    }
    
    if(qaqc$Use19[i]==F){
      chm19 <- raster::mask(chm19, extent_i, inverse = T)
    }
    
    if(qaqc$Use20[i]==F){
      chm20 <- raster::mask(chm20, extent_i, inverse = T)
    }

  }

# Make sure all years have the same extent

  chm18 <- raster::crop(chm18, raster::extent(chm17))
  chm19 <- raster::crop(chm19, raster::extent(chm17))
  chm20 <- raster::crop(chm20, raster::extent(chm17))
  
  
# Cloud masks for 2017-2020
  mask17 <- raster::raster("CloudMask_2017.tif")
  mask18 <- raster::raster("CloudMask_2018.tif")
  mask19 <- raster::raster("CloudMask_2019.tif")
  mask20 <- raster::raster("CloudMask_2020.tif")
  
  # Resample to extent and resolution of CHMs
  mask17 <- raster::resample(mask17, chm17)
  mask18 <- raster::resample(mask18, chm18)
  mask19 <- raster::resample(mask19, chm19)
  mask20 <- raster::resample(mask20, chm20)
  
  # Remove cloud pixels
  chm17[!(mask17==1)] <- NA
  chm18[!(mask18==1)] <- NA
  chm19[!(mask19==1)] <- NA
  chm20[!(mask20==1)] <- NA
    
# # Save
#   raster::writeRaster(chm17, "CHM_2017_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm18, "CHM_2018_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm19, "CHM_2019_QAQC.tif", overwrite=T)
#   raster::writeRaster(chm20, "CHM_2020_QAQC.tif", overwrite=T)
  

#### READ PROCESSED DATA ####  
  
  # Canopy height models for all years
    chm17 <- raster::raster("CHM_2017_QAQC.tif")
    chm18 <- raster::raster("CHM_2018_QAQC.tif")
    chm19 <- raster::raster("CHM_2019_QAQC.tif")
    chm20 <- raster::raster("CHM_2020_QAQC.tif")
  
  # Calculate the change in canopy height for each interval  
    d17to18 <- chm18-chm17
    d18to19 <- chm19-chm18
    d19to20 <- chm20-chm19
    
  # Mask out areas that are initially < 5 m in height and decrease between two years
    
    # 2017 - 2018
      short17 <- rep(0, length(raster::values(chm17)))
        short17[raster::values(chm17)<5 & !is.na(raster::values(chm17))] <- 1
      
      omit17to18 <- rep(0, length(raster::values(d17to18)))
        omit17to18[raster::values(d17to18) < 0 & !is.na(raster::values(d17to18))] <- 1
        omit17to18[short17==0] <- 0
        
      d17to18@data@values[omit17to18==1] <- NA  

    # 2018 - 2019
      short18 <- rep(0, length(raster::values(chm18)))
      short18[raster::values(chm18)<5 & !is.na(raster::values(chm18))] <- 1
      
      omit18to19 <- rep(0, length(raster::values(d18to19)))
      omit18to19[raster::values(d18to19) < 0 & !is.na(raster::values(d18to19))] <- 1
      omit18to19[short18==0] <- 0
      
      d18to19@data@values[omit18to19==1] <- NA  
    
    # 2018 - 2019
      short19 <- rep(0, length(raster::values(chm19)))
      short19[raster::values(chm19)<5 & !is.na(raster::values(chm19))] <- 1
      
      omit19to20 <- rep(0, length(raster::values(d19to20)))
      omit19to20[raster::values(d19to20) < 0 & !is.na(raster::values(d19to20))] <- 1
      omit19to20[short19==0] <- 0
      
      d19to20@data@values[omit19to20==1] <- NA  
      
#### BINARY NEW GAPS #### 
  
  # Use ForestGapR package to delineate new gaps
    
    # Define gap height threshold, min gap size, max gap size, and min area:perimeter ratio
      gapHtThresh <- -5
      gapSzMin <- 10
      gapSzMax <- 10^6
      areaPerimThresh <- 0.6
      
    # 2017 - 2018
        
      # Identify gaps  
        gaps17to18 <- ForestGapR::getForestGaps(d17to18,
                                                threshold = gapHtThresh ,
                                                size=c(gapSzMin,gapSzMax))
        
      # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps17to18sp <- ForestGapR::GapSPDF(gaps17to18)
      
      # Calculate the area and perimeter from each gap object
        gaps17to18sp@data$area <- NA
        gaps17to18sp@data$perimeter <- NA
        for(i in 1:length(gaps17to18sp)){
          gaps17to18sp[gaps17to18sp$gap_id==i,"area"] <- raster::area(gaps17to18sp[gaps17to18sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps17to18sp[gaps17to18sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps17to18sp[gaps17to18sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps17to18sp[gaps17to18sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps17to18sp[gaps17to18sp$gap_id==i,"perimeter"] <- sum(perims_j)
          
        }
        
      # Calculate the ratio of area to perimeter
        gaps17to18sp@data$ratio <- gaps17to18sp@data$area/gaps17to18sp@data$perimeter
      
      # Remove gaps that don't meet a minimum area/perimeter threshold
        for(i in 1:length(gaps17to18sp@data$ratio)){
          if(gaps17to18sp@data$ratio[i] < areaPerimThresh){
            gaps17to18[gaps17to18==gaps17to18sp@data$gap_id[i]] <- NA
          }
        }
      
    # 2018 - 2019
        
      # Identify gaps  
        gaps18to19 <- ForestGapR::getForestGaps(d18to19,
                                                threshold = gapHtThresh ,
                                                size=c(gapSzMin,gapSzMax))
      
      # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps18to19sp <- ForestGapR::GapSPDF(gaps18to19)
      
      # Calculate the area and perimeter from each gap object
        gaps18to19sp@data$area <- NA
        gaps18to19sp@data$perimeter <- NA
        for(i in 1:length(gaps18to19sp)){
          gaps18to19sp[gaps18to19sp$gap_id==i,"area"] <- raster::area(gaps18to19sp[gaps18to19sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps18to19sp[gaps18to19sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps18to19sp[gaps18to19sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps18to19sp[gaps18to19sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps18to19sp[gaps18to19sp$gap_id==i,"perimeter"] <- sum(perims_j)
          
        }
      
      # Calculate the ratio of area to perimeter
        gaps18to19sp@data$ratio <- gaps18to19sp@data$area/gaps18to19sp@data$perimeter
      
      # Remove gaps that don't meet a minimum area/perimeter threshold
        for(i in 1:length(gaps18to19sp@data$ratio)){
          if(gaps18to19sp@data$ratio[i] < areaPerimThresh){
            gaps18to19[gaps18to19==gaps18to19sp@data$gap_id[i]] <- NA
          }
        }
        
    # 2019 - 2020
        
      # Identify gaps  
        gaps19to20 <- ForestGapR::getForestGaps(d19to20,
                                                threshold = gapHtThresh ,
                                                size=c(gapSzMin,gapSzMax))
      
      # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps19to20sp <- ForestGapR::GapSPDF(gaps19to20)
      
      # Calculate the area and perimeter from each gap object
        gaps19to20sp@data$area <- NA
        gaps19to20sp@data$perimeter <- NA
        for(i in 1:length(gaps19to20sp)){
          gaps19to20sp[gaps19to20sp$gap_id==i,"area"] <- raster::area(gaps19to20sp[gaps19to20sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps19to20sp[gaps19to20sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps19to20sp[gaps19to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps19to20sp[gaps19to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps19to20sp[gaps19to20sp$gap_id==i,"perimeter"] <- sum(perims_j)
          
        }
      
      # Calculate the ratio of area to perimeter
        gaps19to20sp@data$ratio <- gaps19to20sp@data$area/gaps19to20sp@data$perimeter
      
      # Remove gaps that don't meet a minimum area/perimeter threshold
        for(i in 1:length(gaps19to20sp@data$ratio)){
          if(gaps19to20sp@data$ratio[i] < areaPerimThresh){
            gaps19to20[gaps19to20==gaps19to20sp@data$gap_id[i]] <- NA
          }
        }
  
  # Save gap shapefiles
      gaps17to18sp$use <- ifelse(gaps17to18sp$ratio<areaPerimThresh,
                                 F,T)
      rgdal::writeOGR(gaps17to18sp,
                      dsn = "gaps17to18_shapefile",
                      layer = "gaps17to18sp", 
                      driver = "ESRI Shapefile")
      
      gaps18to19sp$use <- ifelse(gaps18to19sp$ratio<areaPerimThresh,
                                 F,T)
      rgdal::writeOGR(gaps18to19sp,
                      dsn = "gaps18to19_shapefile",
                      layer = "gaps18to19sp", 
                      driver = "ESRI Shapefile")
      
      gaps19to20sp$use <- ifelse(gaps19to20sp$ratio<areaPerimThresh,
                                 F,T)
      rgdal::writeOGR(gaps19to20sp,
                      dsn = "gaps19to20_shapefile",
                      layer = "gaps19to20sp", 
                      driver = "ESRI Shapefile")
  
  # # Save rasters of new gap pixels
  #     raster::writeRaster(gaps17to18, file="newGaps17to18.tif", overwrite = T)
  #     raster::writeRaster(gaps18to19, file="newGaps18to19.tif", overwrite = T)
  #     raster::writeRaster(gaps19to20, file="newGaps19to20.tif", overwrite = T)
      
      
  # # Save rasters of canopy height change
  #     raster::writeRaster(d17to18, file="dCHM17to18.tif", overwrite = T)
  #     raster::writeRaster(d18to19, file="dCHM18to19.tif", overwrite = T)
  #     raster::writeRaster(d19to20, file="dCHM19to20.tif", overwrite = T)

#### GAPS NEIGHBORING NA CELLS ####
      
  # 2017 to 2018 
      
    gaps17to18sp <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")    
    d17to18 <- raster::raster("dCHM17to18.tif")    
    
    gaps17to18sp$borderNAs <- NA
    
    for(i in 1:length(gaps17to18sp)){
      # create a 1 m buffer (one cell) around gap polygon
        polyBuff <- raster::buffer(gaps17to18sp[i,], 1)
      # find raster cells within that polygon
        cells <- raster::cellFromPolygon(object = d17to18,
                                         p = polyBuff)[[1]]
      # count NA cells within buffered gap polygon
        gaps17to18sp$borderNAs[i] <- length(d17to18[cells][is.na(d17to18[cells])])
  
    }

  # 2018 to 2019 
    
    gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")    
    d18to19 <- raster::raster("dCHM18to19.tif")    
    
    gaps18to19sp$borderNAs <- NA
    
    for(i in 1:length(gaps18to19sp)){
      # create a 1 m buffer (one cell) around gap polygon
      polyBuff <- raster::buffer(gaps18to19sp[i,], 1)
      # find raster cells within that polygon
      cells <- raster::cellFromPolygon(object = d18to19,
                                       p = polyBuff)[[1]]
      # count NA cells within buffered gap polygon
      gaps18to19sp$borderNAs[i] <- length(d18to19[cells][is.na(d18to19[cells])])
      
    }  
    
  # 2019 to 2020 
    
    gaps19to20sp <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")    
    d19to20 <- raster::raster("dCHM19to20.tif")    
    
    gaps19to20sp$borderNAs <- NA
    
    for(i in 1:length(gaps19to20sp)){
      # create a 1 m buffer (one cell) around gap polygon
      polyBuff <- raster::buffer(gaps19to20sp[i,], 1)
      # find raster cells within that polygon
      cells <- raster::cellFromPolygon(object = d19to20,
                                       p = polyBuff)[[1]]
      # count NA cells within buffered gap polygon
      gaps19to20sp$borderNAs[i] <- length(d19to20[cells][is.na(d19to20[cells])])
      
    }  
    
#### SUMMARY GAP STATS PER YEAR ####

      # Make vectors of all sampled values in each interval, and of gap values
      # sampled in each interval
        # All values
          allVals18 <- raster::values(d17to18)
          allVals19 <- raster::values(d18to19)
          allVals20 <- raster::values(d19to20)
        # Gap values
          gapVals18 <- raster::values(gaps17to18)
          gapVals19 <- raster::values(gaps18to19)
          gapVals20 <- raster::values(gaps19to20)
      
      
      # Total area sampled in each interval (in ha)
        area18 <- length(allVals18[!is.na(allVals18)])/10000
        area19 <- length(allVals19[!is.na(allVals19)])/10000
        area20 <- length(allVals20[!is.na(allVals20)])/10000

      # Gap area (% yr-1)
      # Correct for measuring 13 months in 2019-2020  
        round(100*length(gapVals18[!is.na(gapVals18)])/length(allVals18[!is.na(allVals18)]),1)
        round(100*length(gapVals19[!is.na(gapVals19)])/length(allVals19[!is.na(allVals19)]),1)
        round(100*length(gapVals20[!is.na(gapVals20)])/length(allVals20[!is.na(allVals20)])/(13/12),1)

      
      # Gap number (gaps ha-1 yr-1)
        round(length(unique(gapVals18[!is.na(gapVals18)]))/(area18),1)
        round(length(unique(gapVals19[!is.na(gapVals19)]))/(area19),1)
        round(length(unique(gapVals20[!is.na(gapVals20)]))/(area20),1)

#### CALCULATE GAP SIZE FREQUENCY DISTRIBUTION ####
        
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
  
  # # Test MCMC routine on simulated data
  #   lambda <- 2.5
  #   gapValues <- c(1:10000)
  #   gapProbs <- (gapValues^-lambda)/VGAM::zeta(lambda)
  #   gapSims <- sample(size=1000, x=gapValues, prob=gapProbs, replace=T)
  #   
  #   fitdp <- optimize(dpl, data=gapSims, lower = 1.0001, upper = 20, maximum = T)
  # 
  #   fitSim <- runMCMC(startvalue = fitdp$maximum, iterations = 5000, data = gapSims)
  #   hist(fitSim[[1]][1001:5000]);abline(v=lambda,col="red")
  #   sum(fitSim[[2]]);mean(fitSim[[1]][1001:5000])-lambda
        
  # Get vector of gap sizes for each interval
    gapSz18 <- gaps17to18sp$area[gaps17to18sp$ratio > areaPerimThresh]
    gapSz19 <- gaps18to19sp$area[gaps18to19sp$ratio > areaPerimThresh]
    gapSz20 <- gaps19to20sp$area[gaps19to20sp$ratio > areaPerimThresh]
    
  # Run MCMC routine to estimate lambda  
    # Set seed so "random" results are repeatable
      set.seed(1)
    # Set number of simulations  
      sims <- 100000
    # Run MCMC function for each year
      lambda18MC <- runMCMC(startvalue = optimize(dpl, 
                                                  data=gapSz18, 
                                                  lower = 1.0001, upper = 20, 
                                                  maximum = T)$maximum, 
                            iterations = sims, 
                            data = gapSz18)
      
      lambda19MC <- runMCMC(startvalue = optimize(dpl, 
                                                  data=gapSz19, 
                                                  lower = 1.0001, upper = 20, 
                                                  maximum = T)$maximum, 
                            iterations = sims, 
                            data = gapSz19)
      
      lambda20MC <- runMCMC(startvalue = optimize(dpl, 
                                                  data=gapSz20, 
                                                  lower = 1.0001, upper = 20, 
                                                  maximum = T)$maximum, 
                            iterations = sims, 
                            data = gapSz20)
      
    # Summarize results for each year
      # Define parameters for sampling chain
      burnIn <- 5001
      skipN <- 25
      
      # 2017-2018  
        lambdaMin18 <- round(quantile(lambda18MC[[1]][seq(burnIn,sims,skipN)], 0.025),3)
        lambdaMed18 <- round(quantile(lambda18MC[[1]][seq(burnIn,sims,skipN)], 0.50),3)
        lambdaMax18 <- round(quantile(lambda18MC[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      # 2018-2019  
        lambdaMin19 <- round(quantile(lambda19MC[[1]][seq(burnIn,sims,skipN)], 0.025),3)
        lambdaMed19 <- round(quantile(lambda19MC[[1]][seq(burnIn,sims,skipN)], 0.50),3)
        lambdaMax19 <- round(quantile(lambda19MC[[1]][seq(burnIn,sims,skipN)], 0.975),3)
      # 2019-2020  
        lambdaMin20 <- round(quantile(lambda20MC[[1]][seq(burnIn,sims,skipN)], 0.025),3)
        lambdaMed20 <- round(quantile(lambda20MC[[1]][seq(burnIn,sims,skipN)], 0.50),3)
        lambdaMax20 <- round(quantile(lambda20MC[[1]][seq(burnIn,sims,skipN)], 0.975),3)
        