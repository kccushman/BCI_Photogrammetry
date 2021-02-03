#### Read data ####
  
  # Read validation gap shapefiles from Raquel's 50 ha plot analysis
    valGaps <- rgdal::readOGR("Raquel50haData/gaps_2014_2019_all/gaps_2014_2019_all.shp")
    
    # Convert date to date class
    valGaps$dateClass <- as.Date(valGaps$date)
    
    # Only keep gaps > 10 m2
    valGaps <- valGaps[valGaps$area_m2 > 10,]
    
    # Split into years that correspond to annual island flights
    date17 <- as.Date("2017-06-22")
    date18 <- as.Date("2018-06-20")
    date19 <- as.Date("2019-06-24")
    date20 <- as.Date("2020-07-31")
    
    valGaps17to18 <- valGaps[valGaps$dateClass>date17 & valGaps$dateClass<date18,]
    valGaps18to19 <- valGaps[valGaps$dateClass>date18 & valGaps$dateClass<date19,]
    valGaps19to20 <- valGaps[valGaps$dateClass>date19 & valGaps$dateClass<date20,]
    
  # Read canopy height change rasters
    dchm17to18 <- raster::raster("dCHM17to18.tif")
    dchm18to19 <- raster::raster("dCHM18to19.tif")
    dchm19to20 <- raster::raster("dCHM19to20.tif")
    
  # Read new gap shapefiles
    gaps17to18 <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")
    gaps18to19 <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")
    gaps19to20 <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")
    
    # Only keep gaps within the plot (ADD)
    
#### Make initial plots for each year ####
    
    # 2017 to 2018
    raster::plot(dchm17to18,
                 ext = raster::extent(valGaps))
    raster::plot(valGaps17to18, col = NA, 
                 border = "black", lwd = 1, add=T)
    raster::plot(gaps17to18[gaps17to18$use==T,], col = NA, 
                 border = "red", lty=2, add=T)
    
    raster::plot(dchm18to19,
                 ext = raster::extent(valGaps))
    raster::plot(valGaps18to19, col = NA, 
                 border = "black", lwd = 1, add=T)
    raster::plot(gaps18to19[gaps18to19$use==T,], col = NA, 
                 border = "red", lty=2, add=T)
    
    # DON'T USE LAST YEAR because Raquel's data don't extend this long.
    raster::plot(dchm19to20,
                 ext = raster::extent(valGaps))
    raster::plot(valGaps19to20, col = NA, 
                 border = "black", lwd = 1, add=T)
    raster::plot(gaps19to20[gaps19to20$use==T,], col = NA, 
                 border = "red", lty=2, add=T)
    
    
#### Calculate precision and recall ####
    
  # First, find which of Raquel's gaps were covered in annual data
    valGaps17to18$area_sampled <- raster::extract(dchm17to18, valGaps17to18, fun = mean, na.rm=T)[,1]
    valGaps18to19$area_sampled <- raster::extract(dchm18to19, valGaps18to19, fun = mean, na.rm=T)[,1]
  
    # Only keep gaps sampled
    valGaps17to18 <- valGaps17to18[!is.na(valGaps17to18$area_sampled),]
    valGaps18to19 <- valGaps18to19[!is.na(valGaps18to19$area_sampled),]
    
  # Manually make proj4 strings the same
    sp::proj4string(gaps17to18) <- sp::proj4string(valGaps17to18)
    sp::proj4string(gaps18to19) <- sp::proj4string(valGaps18to19)
    
  # Calculate precision: how many of my gaps are also in Raquel's data? (ADD)
    
  # Calculate recall: how many of Raquel's gaps are also in my data?
    
    valGaps17to18$observed <- NA
    for(i in 1:length(valGaps17to18)){
      valGaps17to18$observed[i] <- rgeos::gIntersects(valGaps17to18[i,], gaps17to18, prepared = F)
    }
    
    valGaps18to19$observed <- NA
    for(i in 1:length(valGaps18to19)){
      valGaps18to19$observed[i] <- rgeos::gIntersects(valGaps18to19[i,], gaps18to19, prepared = F)
    }
    
    # Proportion of gaps recalled
      recall17to18_n <- length(valGaps17to18[valGaps17to18$observed==T,])/length(valGaps17to18)
      recall18to19_n <- length(valGaps18to19[valGaps18to19$observed==T,])/length(valGaps18to19)
      
    # Proportion of gap area (based on T/F) recalled
      recall17to18_area <- sum(valGaps17to18$area_m2[valGaps17to18$observed==T])/sum(valGaps17to18$area_m2)
      recall18to19_area <- sum(valGaps18to19$area_m2[valGaps18to19$observed==T])/sum(valGaps18to19$area_m2)
      