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
  
#### Create SpatialPolygonsDataFrame across BCI ####
    
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
    #Mapping
    idx.mapping <- as.vector(t(matrix(1:(nCellX*nCellY), nrow = nCellX, ncol = nCellY)))
    bci.gaps2 <- bci.gaps[idx.mapping, ]

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
    
    # ~ 15 mins for 1 ha pixels
    bci.gaps2$areaObs <- countPix(gapPoly = bci.gaps2,
                                  baseLayer = d19to20tall)
    # ~ 27 mins each for 1 ha pixels
    system.time(curvQuant <- quantTopo(gapPoly = bci.gaps2,
                                       baseLayer = d19to20tall,
                                       topoLayer = curvRaster))
    system.time(slopeQuant <- quantTopo(gapPoly = bci.gaps2,
                                       baseLayer = d19to20tall,
                                       topoLayer = slopeRaster))
    system.time(aspectQuant <- quantTopo(gapPoly = bci.gaps2,
                                       baseLayer = d19to20tall,
                                       topoLayer = aspectRaster))
    
    bci.gaps2$curvMean <- curvQuant$meanVal
    bci.gaps2$curvMed <- curvQuant$medVal
    bci.gaps2$slopeMean <- slopeQuant$meanVal
    bci.gaps2$slopeMed <- slopeQuant$medVal
    bci.gaps2$aspectMean <- aspectQuant$meanVal
    bci.gaps2$aspectMed <- aspectQuant$medVal
    
#### Calculate expected gaps per cell from area measured ####
  # Calculate mean gaps/ha across the whole island
    avgRate <- length(gapsPts19to20)/(length(raster::values(d19to20tall)[!is.na(raster::values(d19to20tall))])/10000)
    
  # Caclulate expect gaps for each cell based on area measured
    bci.gaps2$expectedGaps <- round(bci.gaps2$areaObs*avgRate)
    
  # Make any rows NA where the expected value is 0
    bci.gaps2[bci.gaps2$expectedGaps==0,"Ngaps"] <- NA

#### Very preliminary models ####
    
    # Look at summary topo values
    hist(bci.gaps2@data$curvMed)
    hist(bci.gaps2@data$slopeMed)
    hist(log(bci.gaps2@data$slopeMed))
    hist(bci.gaps2@data$aspectMed)
    
    bci.gaps2$logSlopeMed <- log(bci.gaps2@data$slopeMed)
    
    plot(Ngaps~curvMed, data=bci.gaps2@data, pch=19, col=adjustcolor("black",0.1))
    plot(Ngaps~slopeMed, data=bci.gaps2@data, pch=19, col=adjustcolor("black",0.1))
    plot(Ngaps~aspectMed, data=bci.gaps2@data, pch=19, col=adjustcolor("black",0.1))
    
    
    #Log-Poisson regression
    m0 <- inla(Ngaps ~ curvMed + logSlopeMed, family = "poisson",
               data = as.data.frame(bci.gaps2),
               E = expectedGaps,
               control.compute = list(dic = TRUE))
    summary(m0)
    
    #Log-Poisson regression with random effects
    bci.gaps2$ID <- 1:length(bci.gaps2)
    m0.re <- inla(Ngaps ~ curvMed + logSlopeMed + f(ID), family = "poisson",
                  data = as.data.frame(bci.gaps2),
                  control.compute = list(dic = TRUE))
    summary(m0.re)
    
    #RW2d
    m0.rw2d <- inla(Ngaps ~ curvMed + logSlopeMed +
                      f(ID, model = "rw2d", nrow = nCellY, ncol = nCellX),
                    family = "poisson", data = as.data.frame(bci.gaps2),
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE) )
    
    summary(m0.rw2d)
    
    #Matern2D
    m0.m2d <- inla(Ngaps ~ curvMed + logSlopeMed +
                     f(ID, model = "matern2d", nrow = nCellY, ncol = nCellX),
                   family = "poisson", data = as.data.frame(bci.gaps2),
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE))
    summary(m0.m2d)
    
    m0$dic$dic
    m0.re$dic$dic
    m0.rw2d$dic$dic
    m0.m2d$dic$dic
    
    save(bci.gaps2, file="INLA_prelim.RData")