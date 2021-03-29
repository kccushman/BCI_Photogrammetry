#### Read data ####
  
  # Read validation gap shapefiles from Raquel's 50 ha plot analysis
    valGaps <- rgdal::readOGR("Raquel50haData/gaps_2014_2019_all/gaps_2014_2019_all.shp")
    
    # Convert date to date class
    valGaps$dateClass <- as.Date(valGaps$date)
    
    # Only keep gaps > 25 m2
    valGaps <- valGaps[valGaps$area_m2 > 25,]
    
    # Split into years that correspond to annual island flights
    date15 <- as.Date("2015-07-1")
    date18 <- as.Date("2018-06-20")
    date20 <- as.Date("2020-07-31")
    
    valGaps15to18 <- valGaps[valGaps$dateClass>date15 & valGaps$dateClass<date18,]
    valGaps18to20 <- valGaps[valGaps$dateClass>date18 & valGaps$dateClass<date20,]
    
  # Read canopy height change rasters
    dchm15to18 <- raster::raster("dCHM15to18_tin.tif")
    dchm18to20 <- raster::raster("dCHM18to20_tin.tif")
    
  # Read new gap shapefiles
    gaps15to18 <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
    gaps18to20 <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")
    
    # Only keep gaps within the plot (ADD)
    gaps15to18 <- raster::crop(gaps15to18,raster::extent(valGaps))
    gaps18to20 <- raster::crop(gaps18to20,raster::extent(valGaps))
    
    # Only keep gaps that pass area:perimeter threshold
    gaps15to18 <- gaps15to18[gaps15to18$use==T,]
    gaps18to20 <- gaps18to20[gaps18to20$use==T,]
    
#### Make initial plots for each year ####
    
    # 2015 to 2018 (Don't use 2018-2020 period because Raquel's data don't extend that long)
    raster::plot(dchm15to18,
                 ext = raster::extent(valGaps))
    raster::plot(valGaps15to18, col = NA, 
                 border = "black", lwd = 1, add=T)
    raster::plot(gaps15to18, col = NA, 
                 border = "red", lty=2, add=T)
    
#### Calculate precision and recall ####
    
  # First, find which of Raquel's gaps were covered in annual data
    valGaps15to18$area_sampled <- raster::extract(dchm15to18, valGaps15to18, fun = mean, na.rm=T)[,1]

    # Only keep gaps sampled
    valGaps15to18 <- valGaps15to18[!is.na(valGaps15to18$area_sampled),]

    # Manually make proj4 strings the same
      sp::proj4string(gaps15to18) <- sp::proj4string(valGaps15to18)

  # Calculate precision: how many of my gaps are also in Raquel's data?
    gaps15to18$observed <- NA
    for(i in 1:length(gaps15to18)){
      gaps15to18$observed[i] <- rgeos::gIntersects(gaps15to18[i,], valGaps15to18)
    }
    
    # Proportion of gaps recalled
    precision15to18_n <- length(gaps15to18[gaps15to18$observed==T,])/length(gaps15to18)
    
    # Proportion of gap area (based on T/F) recalled
    precision15to18_area <- sum(gaps15to18$area[gaps15to18$observed==T])/sum(gaps15to18$area)
    
    
  # Calculate recall: how many of Raquel's gaps are also in my data?
    
    valGaps15to18$observed <- NA
    for(i in 1:length(valGaps15to18)){
      valGaps15to18$observed[i] <- rgeos::gIntersects(valGaps15to18[i,], gaps15to18)
    }
    
    # Proportion of gaps recalled
      recall15to18_n <- length(valGaps15to18[valGaps15to18$observed==T,])/length(valGaps15to18)

    # Proportion of gap area (based on T/F) recalled
      recall15to18_area <- sum(valGaps15to18$area_m2[valGaps15to18$observed==T])/sum(valGaps15to18$area_m2)
      
  
  # Save shapefiles to look in ArcGIS
      rgdal::writeOGR(valGaps15to18,
                      dsn = "Monthly_50haGaps_2015-2018",
                      layer = "MonthlyGaps", 
                      driver = "ESRI Shapefile")
      
      rgdal::writeOGR(gaps15to18,
                      dsn = "Annual_50haGaps_2015-2018",
                      layer = "AnnualGaps", 
                      driver = "ESRI Shapefile")
      
