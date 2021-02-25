#### Load packages ####

library("INLA")

#### Load data ####

  # Gap polygons
    gaps15to17sp <- rgdal::readOGR("gaps15to17_shapefile/gaps15to17sp.shp")
    gaps17to18sp <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")
    gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")
    gaps19to20sp <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")
  
  # Gap rasters
    gaps15to17 <- raster::raster("newGaps15to17.tif")
    gaps17to18 <- raster::raster("newGaps17to18.tif")
    gaps18to19 <- raster::raster("newGaps18to19.tif")
    gaps19to20 <- raster::raster("newGaps19to20.tif")
    
  # Canopy height change rasters where only possible gap values (> 5 m initially) are included  
    d15to17tall <- raster::raster("dCHM15to17tall.tif")
    d17to18tall <- raster::raster("dCHM17to18tall.tif")
    d18to19tall <- raster::raster("dCHM18to19tall.tif")
    d19to20tall <- raster::raster("dCHM19to20tall.tif")
    
  # Forest age polygon
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
      ageUse <- age[!(age$TYPE=="Clearings"),]
      ageUse$AgeClass <- "Other"
      ageUse$AgeClass[ageUse$Mascaro_Co == "> 400"] <- "OldGrowth"
      ageUse$AgeClass[ageUse$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
      
  # Soil type polygon  
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(d17to18tall))
    
  # Aspect raster
    aspectRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_21.tif")
    
  # Distance above drainage raster
    drainRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/distAboveStream_1000.tif")
    
#### Create SpatialPolygonsDataFrame across BCI with proper order for R-INLA ####
    
    # Define cell size (in m)
    cellSize <- 100
    
    # Count number of cells in each dimension
    nCellX <- ceiling((d17to18tall@extent@xmax-d17to18tall@extent@xmin)/cellSize)
    nCellY <- ceiling((d17to18tall@extent@ymax-d17to18tall@extent@ymin)/cellSize)
    
    # Define a SpatialGridDataFrame
    bci.pix <- spatstat::im(mat = 1:(nCellX*nCellY),
                            xcol = d17to18tall@extent@xmin + 0.5*cellSize + cellSize*(0:(nCellX-1)),
                            yrow = d17to18tall@extent@ymin + 0.5*cellSize + cellSize*(0:(nCellY-1)))
    
    bci.grid <- as(bci.pix,"SpatialGridDataFrame")
    
    # Convert to SpatialPolygons
    bci.poly <- as(bci.grid, "SpatialPolygons")
    
    # Convert gap polygons to points
    gapsPts15to17 <- sp::SpatialPoints(coords = gaps17to18sp@data[gaps15to17sp$use==T,c("X1","X2")])
    gapsPts17to18 <- sp::SpatialPoints(coords = gaps17to18sp@data[gaps17to18sp$use==T,c("X1","X2")])
    gapsPts18to19 <- sp::SpatialPoints(coords = gaps18to19sp@data[gaps18to19sp$use==T,c("X1","X2")])
    gapsPts19to20 <- sp::SpatialPoints(coords = gaps19to20sp@data[gaps19to20sp$use==T,c("X1","X2")])
    
    # 2015 - 2017
      #Number of gap observations per cell
      idx <- over(gapsPts15to17, bci.poly)
      tab.idx <- table(idx)
      
      #Add number of gaps
      d <- data.frame(Ngaps = rep(0, length(bci.poly)))
      row.names(d) <- paste0("g", 1:length(bci.poly))
      d$Ngaps[as.integer(names(tab.idx))] <- tab.idx
      
      #SpatialPolygonsDataFrame
      bci.gaps <- SpatialPolygonsDataFrame(bci.poly, d)
      
      ## Correct INLA mapping (assumes that data are sorted top-to-bottom by column)    
      idx.mapping <- as.vector(t(matrix(1:(nCellX*nCellY), nrow = nCellX, ncol = nCellY)))
      bci.gaps17 <- bci.gaps[idx.mapping, ]
      
    
    # 2017 - 2018
      #Number of gap observations per cell
      idx <- over(gapsPts17to18, bci.poly)
      tab.idx <- table(idx)
      
      #Add number of gaps
      d <- data.frame(Ngaps = rep(0, length(bci.poly)))
      row.names(d) <- paste0("g", 1:length(bci.poly))
      d$Ngaps[as.integer(names(tab.idx))] <- tab.idx
      
      #SpatialPolygonsDataFrame
      bci.gaps <- SpatialPolygonsDataFrame(bci.poly, d)
      
      ## Correct INLA mapping (assumes that data are sorted top-to-bottom by column)    
      idx.mapping <- as.vector(t(matrix(1:(nCellX*nCellY), nrow = nCellX, ncol = nCellY)))
      bci.gaps18 <- bci.gaps[idx.mapping, ]
    
    # 2018 - 2019
      #Number of gap observations per cell
      idx <- over(gapsPts18to19, bci.poly)
      tab.idx <- table(idx)
      
      #Add number of gaps
      d <- data.frame(Ngaps = rep(0, length(bci.poly)))
      row.names(d) <- paste0("g", 1:length(bci.poly))
      d$Ngaps[as.integer(names(tab.idx))] <- tab.idx
      
      #SpatialPolygonsDataFrame
      bci.gaps <- SpatialPolygonsDataFrame(bci.poly, d)
      
      ## Correct INLA mapping (assumes that data are sorted top-to-bottom by column)    
      idx.mapping <- as.vector(t(matrix(1:(nCellX*nCellY), nrow = nCellX, ncol = nCellY)))
      bci.gaps19 <- bci.gaps[idx.mapping, ]
    
    # 2019 - 2020
      #Number of gap observations per cell
      idx <- over(gapsPts19to20, bci.poly)
      tab.idx <- table(idx)
      
      #Add number of gaps
      d <- data.frame(Ngaps = rep(0, length(bci.poly)))
      row.names(d) <- paste0("g", 1:length(bci.poly))
      d$Ngaps[as.integer(names(tab.idx))] <- tab.idx
      
      #SpatialPolygonsDataFrame
      bci.gaps <- SpatialPolygonsDataFrame(bci.poly, d)
      
      ## Correct INLA mapping (assumes that data are sorted top-to-bottom by column)    
      idx.mapping <- as.vector(t(matrix(1:(nCellX*nCellY), nrow = nCellX, ncol = nCellY)))
      bci.gaps20 <- bci.gaps[idx.mapping, ]

#### Quantify area measured and proportion of gaps in each cell ####
    
  # Define a function to count the percent of gap area (%) and total area measured (ha) in each cell    
    propGap <- function(gapPoly, gapLayer, baseLayer){
      
      gapVals <- data.frame(sampHa = rep(NA, length(gapPoly)),
                            gapProp = rep(NA, length(gapPoly)))
      
      for(i in 1:length(gapPoly)){
        base_i <- raster::crop(baseLayer, raster::extent(gapPoly[i,]))
        gaps_i <- raster::crop(gapLayer, raster::extent(gapPoly[i,]))
        pixSampled <- length(raster::values(base_i)[!is.na(raster::values(base_i))])
        pixGaps <- length(raster::values(gaps_i)[!is.na(raster::values(gaps_i))])
        gapVals$sampHa[i] <- pixSampled/10000
        gapVals$gapProp[i] <- pixGaps/pixSampled*100
      }
      return(gapVals)
    }

    # 2015-2017

      # ~ 30 mins each for 1 ha pixels
      gapArea <- propGap(gapPoly = bci.gaps17,
                         gapLayer = gaps15to17,
                         baseLayer = d15to17tall)
      bci.gaps17$areaObs <- gapArea$sampHa
      bci.gaps17$gapProp <- gapArea$gapProp
      
    # 2017-2018
      
      # ~ 30 mins each for 1 ha pixels
      gapArea <- propGap(gapPoly = bci.gaps18,
                         gapLayer = gaps18to18,
                         baseLayer = d17to18tall)
      bci.gaps18$areaObs <- gapArea$sampHa
      bci.gaps18$gapProp <- gapArea$gapProp
      
    # 2018-2019

      # ~ 30 mins each for 1 ha pixels
      gapArea <- propGap(gapPoly = bci.gaps19,
                         gapLayer = gaps18to19,
                         baseLayer = d18to19tall)
      bci.gaps19$areaObs <- gapArea$sampHa
      bci.gaps19$gapProp <- gapArea$gapProp

    # 2019-2020
  
      # ~ 30 mins each for 1 ha pixels
      gapArea <- propGap(gapPoly = bci.gaps20,
                         gapLayer = gaps19to20,
                         baseLayer = d19to20tall)
      bci.gaps20$areaObs <- gapArea$sampHa
      bci.gaps20$gapProp <- gapArea$gapProp
      
#### Assign forest age and soil type to each cell ####

  # Define a function to assign a cell a categorical value from a polygon 
    assignCat <- function(gapPoly, polyLayer, layerName){
      
      polyVals <- rep(NA, length(gapPoly))
      
      for(i in 1:length(gapPoly)){
        cat_i <- raster::crop(polyLayer, raster::extent(gapPoly[i,]))
        
        # Find polygon with highest area
        if(length(cat_i) > 0){
          
          cat_areas <- rep(NA, length(cat_i))
          
          for(j in 1:length(cat_i)){
            cat_areas[j] <- cat_i@polygons[[j]]@area
          }
          
          polyVals[i] <- as.character(cat_i@data[which(cat_areas==max(cat_areas)),layerName])
        }
      }
      return(polyVals)
    }  
      
    # 2015-2017
      
      # < 1 min for 1 ha pixels
      bci.gaps17$soil <- assignCat(gapPoly = bci.gaps17,
                                   polyLayer = soil,
                                   layerName = "SOIL")
      
      # < 1 min for 1 ha pixels
      bci.gaps17$age <- assignCat(gapPoly = bci.gaps17,
                                  polyLayer = ageUse,
                                  layerName = "AgeClass")
      
    # Assume identical across all years
      
      # 2017 - 2018  
        bci.gaps18$age <- bci.gaps17$age
        bci.gaps18$soil <- bci.gaps17$soil
        
      # 2018-2019
        bci.gaps19$age <- bci.gaps17$age
        bci.gaps19$soil <- bci.gaps17$soil
        
      # 2019-2020
        bci.gaps20$age <- bci.gaps17$age
        bci.gaps20$soil <- bci.gaps17$soil
    
#### Quantify aspect and height above drainage in each cell ####
      
  # Define a function to take the mean and median topographic variable within each cell
  quantTopo <- function(gapPoly, baseLayer=NA, topoLayer){
    
    topoVals <- data.frame(meanVal = rep(NA, length(gapPoly)),
                           medVal = rep(NA, length(gapPoly)))
    
    for(i in 1:length(gapPoly)){
      topo_i <- raster::values(raster::crop(topoLayer, raster::extent(gapPoly[i,])))
      
      if(!is.na(baseLayer)){
        base_i <- raster::values(raster::crop(baseLayer, raster::extent(gapPoly[i,])))
        topoVals$meanVal[i] <- mean(topo_i[!is.na(base_i)])
        topoVals$medVal[i] <- median(topo_i[!is.na(base_i)])
      }
      
      if(is.na(baseLayer)){
        topoVals$meanVal[i] <- mean(topo_i, na.rm=T)
        topoVals$medVal[i] <- median(topo_i, na.rm=T)
      }
      
      
    }
    return(topoVals)
  }   
      
    # Assume identical across all years
      
      # ~ 30 mins each for 1 ha pixels
      aspectQuant <- quantTopo(gapPoly = bci.gaps17,
                               baseLayer = NA,
                               topoLayer = aspectRaster)
      
      bci.gaps17$aspectMean <- aspectQuant$meanVal
      bci.gaps17$aspectMed <- aspectQuant$medVal
      bci.gaps18$aspectMean <- aspectQuant$meanVal
      bci.gaps18$aspectMed <- aspectQuant$medVal
      bci.gaps19$aspectMean <- aspectQuant$meanVal
      bci.gaps19$aspectMed <- aspectQuant$medVal
      bci.gaps20$aspectMean <- aspectQuant$meanVal
      bci.gaps20$aspectMed <- aspectQuant$medVal
      
      # ~ 30 mins each for 1 ha pixels
      drainQuant <- quantTopo(gapPoly = bci.gaps17,
                              baseLayer = NA,
                              topoLayer = drainRaster)
      
      bci.gaps17$drainMean <- drainQuant$meanVal
      bci.gaps17$drainMed <- drainQuant$medVal
      bci.gaps18$drainMean <- drainQuant$meanVal
      bci.gaps18$drainMed <- drainQuant$medVal
      bci.gaps19$drainMean <- drainQuant$meanVal
      bci.gaps19$drainMed <- drainQuant$medVal
      bci.gaps20$drainMean <- drainQuant$meanVal
      bci.gaps20$drainMed <- drainQuant$medVal
      
#### Calculate expected gaps per cell from area measured ####
  # Calculate mean gaps/ha across the whole island
    allGapsSum <-   length(gapsPts15to17) + length(gapsPts17to18) + length(gapsPts18to19) + length(gapsPts19to20)
      allAreaSum <-   (length(raster::values(d15to17tall)[!is.na(raster::values(d15to17tall))])+
                         length(raster::values(d17to18tall)[!is.na(raster::values(d17to18tall))]) + 
                         length(raster::values(d18to19tall)[!is.na(raster::values(d18to19tall))]) + 
                         length(raster::values(d19to20tall)[!is.na(raster::values(d19to20tall))]))/10000
    avgRateAll <- allGapsSum/allAreaSum
    
  # Caclulate expect gaps for each cell based on area measured
    bci.gaps17$expectedGaps <- round(bci.gaps17$areaObs*avgRateAll,2)
    bci.gaps18$expectedGaps <- round(bci.gaps18$areaObs*avgRateAll,2)
    bci.gaps19$expectedGaps <- round(bci.gaps19$areaObs*avgRateAll,2)
    bci.gaps20$expectedGaps <- round(bci.gaps20$areaObs*avgRateAll,2)
    

#### Combine all years into a single data frame and save ####
    
    bci.gaps17$Year <- "2017"
    bci.gaps18$Year <- "2018"
    bci.gaps19$Year <- "2019"
    bci.gaps20$Year <- "2020"
    
    bci.gapsAll <- rbind(as.data.frame(bci.gaps17),
                         as.data.frame(bci.gaps18),
                         as.data.frame(bci.gaps19),
                         as.data.frame(bci.gaps20))
    
    save(bci.gaps17, bci.gaps18, bci.gaps19, bci.gaps20, bci.gapsAll, file="INLA_prelim.RData")
    
#### Very preliminary models ####
    
    load("INLA/INLA_prelim.RData")
    
    # Make any rows NA where the expected value is 0
    bci.gapsAll[bci.gapsAll$expectedGaps<0.5,"Ngaps"] <- NA
    
    # How many observations are left?
    dim(bci.gapsAll[!is.na(bci.gapsAll$Ngaps),]) # 3715 observations
    
    # Look at summary topo values
    hist(bci.gapsAll$curvMed)
    hist(bci.gapsAll$slopeMed)
    hist(log(bci.gapsAll$slopeMed))
    hist(bci.gapsAll$aspectMed)
    
    bci.gapsAll$logSlopeMed <- log(bci.gapsAll$slopeMed)
    
    #Log-Poisson regression
    m0 <- inla(Ngaps ~ curvMed + logSlopeMed + soil + Year, family = "poisson",
               data = bci.gapsAll,
               E = expectedGaps,
               control.compute = list(dic = TRUE))
    summary(m0)
    
    #Log-Poisson regression with random effects
    bci.gapsAll$ID <- rep(1:(nCellX*nCellY),3)
    m0.re <- inla(Ngaps ~ curvMed + logSlopeMed + soil + Year + f(ID),
                  family = "poisson",
                  data = bci.gapsAll,
                  E = expectedGaps,
                  control.compute = list(dic = TRUE))
    summary(m0.re)
    
    #RW2d
    m0.rw2d <- inla(Ngaps ~ curvMed + logSlopeMed + soil + Year + 
                      f(ID, model = "rw2d", nrow = nCellY, ncol = nCellX),
                    family = "poisson", 
                    data = bci.gapsAll,
                    E = expectedGaps,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE))
    
    summary(m0.rw2d)
    
    #Matern2D
    m0.m2d <- inla(Ngaps ~ curvMed + logSlopeMed + soil + Year +
                     f(ID, model = "matern2d", nrow = nCellY, ncol = nCellX),
                   family = "poisson", data = bci.gapsAll,
                   E = expectedGaps,
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE))
    summary(m0.m2d)
    
    m0$dic$dic
    m0.re$dic$dic
    m0.rw2d$dic$dic
    m0.m2d$dic$dic
    
