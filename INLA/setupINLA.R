#### Load packages ####

library("INLA")

#### Load data ####

  # Gap polygons
    gaps17to18sp <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")
    gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")
    gaps19to20sp <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")
  
  # Canopy height change rasters where only possible gap values (> 5 m initially) are included  
    d17to18tall <- raster::raster("dCHM17to18tall.tif")
    d18to19tall <- raster::raster("dCHM18to19tall.tif")
    d19to20tall <- raster::raster("dCHM19to20tall.tif")
    
  # Sample topography rasters
    slopeRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_15.tif")
    curvRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_21.tif")
    aspectRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_21.tif")
    
  # Soil type polygon  
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(d17to18tall))
    
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
    gapsPts17to18 <- sp::SpatialPoints(coords = gaps17to18sp@data[gaps17to18sp$use==T,c("X1","X2")])
    gapsPts18to19 <- sp::SpatialPoints(coords = gaps18to19sp@data[gaps18to19sp$use==T,c("X1","X2")])
    gapsPts19to20 <- sp::SpatialPoints(coords = gaps19to20sp@data[gaps19to20sp$use==T,c("X1","X2")])
    
    
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

#### Quantify topographic covariates in each cell ####
    
  # Define a function to count non-NA area in each cell (in ha)
    countPix <- function(gapPoly, baseLayer){
      
      sampHa <- rep(NA, length(gapPoly))
      
      for(i in 1:length(gapPoly)){
        base_i <- raster::crop(baseLayer, raster::extent(gapPoly[i,]))
        pixSampled <- length(raster::values(base_i)[!is.na(raster::values(base_i))])
        sampHa[i] <- pixSampled/10000
      }
      return(sampHa)
    }
    
  # Define a function to take the mean and median topographic variable within each cell
    quantTopo <- function(gapPoly, baseLayer, topoLayer){
      
      topoVals <- data.frame(meanVal = rep(NA, length(gapPoly)),
                             medVal = rep(NA, length(gapPoly)))
      
      for(i in 1:length(gapPoly)){
        base_i <- raster::values(raster::crop(baseLayer, raster::extent(gapPoly[i,])))
        topo_i <- raster::values(raster::crop(topoLayer, raster::extent(gapPoly[i,])))
        
        topoVals$meanVal[i] <- mean(topo_i[!is.na(base_i)])
        topoVals$medVal[i] <- median(topo_i[!is.na(base_i)])
        
      }
      return(topoVals)
    }  
    
  # Define a function to assign a cell a categorical value from a polygon 
    assignCat <- function(gapPoly, polyLayer){
      
      polyVals <- rep(NA, length(gapPoly))
      
      for(i in 1:length(gapPoly)){
        cat_i <- raster::crop(polyLayer, raster::extent(gapPoly[i,]))
        
        # Find polygon with highest area
        if(length(cat_i) > 0){
          
          cat_areas <- rep(NA, length(cat_i))
          
          for(j in 1:length(cat_i)){
            cat_areas[j] <- cat_i@polygons[[j]]@area
          }
        
          polyVals[i] <- as.character(cat_i@data$SOIL[which(cat_areas==max(cat_areas))])
        }
      }
      return(polyVals)
    }  
    
    # 2017-2018
      
      # < 1 min for 1 ha pixels
      bci.gaps18$soil <- assignCat(gapPoly = bci.gaps18,
                                   polyLayer = soil)
      # ~ 15 mins for 1 ha pixels
      bci.gaps18$areaObs <- countPix(gapPoly = bci.gaps18,
                                     baseLayer = d17to18tall)
      # ~ 27 mins each for 1 ha pixels
      curvQuant <- quantTopo(gapPoly = bci.gaps18,
                             baseLayer = d17to18tall,
                             topoLayer = curvRaster)
      slopeQuant <- quantTopo(gapPoly = bci.gaps18,
                              baseLayer = d17to18tall,
                              topoLayer = slopeRaster)
      aspectQuant <- quantTopo(gapPoly = bci.gaps18,
                               baseLayer = d17to18tall,
                               topoLayer = aspectRaster)
      
      bci.gaps18$curvMean <- curvQuant$meanVal
      bci.gaps18$curvMed <- curvQuant$medVal
      bci.gaps18$slopeMean <- slopeQuant$meanVal
      bci.gaps18$slopeMed <- slopeQuant$medVal
      bci.gaps18$aspectMean <- aspectQuant$meanVal
      bci.gaps18$aspectMed <- aspectQuant$medVal
      
    # 2018-2019
      
      # < 1 min for 1 ha pixels
      bci.gaps19$soil <- assignCat(gapPoly = bci.gaps19,
                                   polyLayer = soil)
      # ~ 15 mins for 1 ha pixels
      bci.gaps19$areaObs <- countPix(gapPoly = bci.gaps19,
                                     baseLayer = d18to19tall)
      # ~ 27 mins each for 1 ha pixels
      curvQuant <- quantTopo(gapPoly = bci.gaps19,
                             baseLayer = d18to19tall,
                             topoLayer = curvRaster)
      slopeQuant <- quantTopo(gapPoly = bci.gaps19,
                              baseLayer = d18to19tall,
                              topoLayer = slopeRaster)
      aspectQuant <- quantTopo(gapPoly = bci.gaps19,
                               baseLayer = d18to19tall,
                               topoLayer = aspectRaster)
      
      bci.gaps19$curvMean <- curvQuant$meanVal
      bci.gaps19$curvMed <- curvQuant$medVal
      bci.gaps19$slopeMean <- slopeQuant$meanVal
      bci.gaps19$slopeMed <- slopeQuant$medVal
      bci.gaps19$aspectMean <- aspectQuant$meanVal
      bci.gaps19$aspectMed <- aspectQuant$medVal
      
    # 2019-2020
      
      # < 1 min for 1 ha pixels
      bci.gaps20$soil <- assignCat(gapPoly = bci.gaps20,
                                   polyLayer = soil)
      # ~ 15 mins for 1 ha pixels
      bci.gaps20$areaObs <- countPix(gapPoly = bci.gaps20,
                                     baseLayer = d19to20tall)
      # ~ 27 mins each for 1 ha pixels
      curvQuant <- quantTopo(gapPoly = bci.gaps20,
                             baseLayer = d19to20tall,
                             topoLayer = curvRaster)
      slopeQuant <- quantTopo(gapPoly = bci.gaps20,
                              baseLayer = d19to20tall,
                              topoLayer = slopeRaster)
      aspectQuant <- quantTopo(gapPoly = bci.gaps20,
                               baseLayer = d19to20tall,
                               topoLayer = aspectRaster)
      
      bci.gaps20$curvMean <- curvQuant$meanVal
      bci.gaps20$curvMed <- curvQuant$medVal
      bci.gaps20$slopeMean <- slopeQuant$meanVal
      bci.gaps20$slopeMed <- slopeQuant$medVal
      bci.gaps20$aspectMean <- aspectQuant$meanVal
      bci.gaps20$aspectMed <- aspectQuant$medVal
      
#### Calculate expected gaps per cell from area measured ####
  # Calculate mean gaps/ha across the whole island
    allGapsSum <-   length(gapsPts17to18) + length(gapsPts18to19) + length(gapsPts19to20)
    allAreaSum <-   (length(raster::values(d17to18tall)[!is.na(raster::values(d17to18tall))]) + 
                       length(raster::values(d18to19tall)[!is.na(raster::values(d18to19tall))]) + 
                       length(raster::values(d19to20tall)[!is.na(raster::values(d19to20tall))]))/10000
    avgRateAll <- allGapsSum/allAreaSum
    
    # avgRate18 <- length(gapsPts17to18)/(length(raster::values(d17to18tall)[!is.na(raster::values(d17to18tall))])/10000)
    # avgRate19 <- length(gapsPts18to19)/(length(raster::values(d18to19tall)[!is.na(raster::values(d18to19tall))])/10000)
    # avgRate20 <- length(gapsPts19to20)/(length(raster::values(d19to20tall)[!is.na(raster::values(d19to20tall))])/10000)
      
  # Caclulate expect gaps for each cell based on area measured
    bci.gaps18$expectedGaps <- round(bci.gaps18$areaObs*avgRateAll,2)
    bci.gaps19$expectedGaps <- round(bci.gaps19$areaObs*avgRateAll,2)
    bci.gaps20$expectedGaps <- round(bci.gaps20$areaObs*avgRateAll,2)
    

#### Combine three years into a single data frame and save ####

    bci.gaps18$Year <- "2018"
    bci.gaps19$Year <- "2019"
    bci.gaps20$Year <- "2020"
    
    bci.gapsAll <- rbind(as.data.frame(bci.gaps18),
                         as.data.frame(bci.gaps19),
                         as.data.frame(bci.gaps20))
    
    save(bci.gaps18, bci.gaps19, bci.gaps20, bci.gapsAll, file="INLA_prelim.RData")
    
#### Very preliminary models ####
    
    # Make any rows NA where the expected value is 0
    bci.gapsAll[bci.gapsAll$expectedGaps<0.5,"Ngaps"] <- NA
    
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
                  control.compute = list(dic = TRUE))
    summary(m0.re)
    
    #RW2d
    m0.rw2d <- inla(Ngaps ~ curvMed + logSlopeMed + soil + Year + 
                      f(ID, model = "rw2d", nrow = nCellY, ncol = nCellX),
                    family = "poisson", 
                    data = bci.gapsAll,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE))
    
    summary(m0.rw2d)
    
    #Matern2D
    m0.m2d <- inla(Ngaps ~ curvMed + logSlopeMed + soil + Year +
                     f(ID, model = "matern2d", nrow = nCellY, ncol = nCellX),
                   family = "poisson", data = bci.gapsAll,
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE))
    summary(m0.m2d)
    
    m0$dic$dic
    m0.re$dic$dic
    m0.rw2d$dic$dic
    m0.m2d$dic$dic
    
