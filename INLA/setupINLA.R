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
    # resample to same extent as gap rasters (adds NA area to edges)
    drainRaster <- raster::resample(drainRaster, gaps19to20)
    
#### Create SpatialPolygonsDataFrame across BCI with proper order for R-INLA ####
    
    # Define cell size (in m)
    cellSize <- 40
    
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
    gapsPts15to17 <- sp::SpatialPoints(coords = gaps15to17sp@data[gaps15to17sp$use==T,c("X1","X2")])
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
    
  # Define a function to count the percent of gap area (0 to 1) and total area measured (ha) in each cell    
      propGap <- function(gapPoly, gapLayer, baseLayer, cellSz, nX, nY){
        
        gapVals <- data.frame(sampHa = rep(NA, length(gapPoly)),
                              gapProp = rep(NA, length(gapPoly)))
        
        baseCrop <- raster::crop(baseLayer, raster::extent(gapPoly))
        baseCrop[!is.na(baseCrop)] <- 1
        gapCrop <- raster::crop(gapLayer, raster::extent(gapPoly))
        gapCrop[!is.na(gapCrop)] <- 1
        
        resampleBase <- raster::aggregate(baseCrop, cellSz, fun = sum)
        resampleGap <- raster::aggregate(gapCrop, cellSz, fun = sum)
        
        valuesBase <- raster::values(resampleBase)
        valuesGap <- raster::values(resampleGap)
        
        # reorder to proper order for INLA object
        gapVals$sampHa <- valuesBase[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]/10000
        gapVals$gapProp <- valuesGap[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]/valuesBase[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]
        
        return(gapVals)
      }

    # 2015-2017

      gapArea <- propGap(gapPoly = bci.gaps17,
                         gapLayer = gaps15to17,
                         baseLayer = d15to17tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
      bci.gaps17$areaObs <- gapArea$sampHa
      bci.gaps17$gapProp <- gapArea$gapProp
      
    # 2017-2018
      
      gapArea <- propGap(gapPoly = bci.gaps18,
                         gapLayer = gaps17to18,
                         baseLayer = d17to18tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
      bci.gaps18$areaObs <- gapArea$sampHa
      bci.gaps18$gapProp <- gapArea$gapProp
      
    # 2018-2019

      gapArea <- propGap(gapPoly = bci.gaps19,
                         gapLayer = gaps18to19,
                         baseLayer = d18to19tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
      bci.gaps19$areaObs <- gapArea$sampHa
      bci.gaps19$gapProp <- gapArea$gapProp

    # 2019-2020
  
      gapArea <- propGap(gapPoly = bci.gaps20,
                         gapLayer = gaps19to20,
                         baseLayer = d19to20tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
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
      
     soilVals <- assignCat(gapPoly = bci.gaps17,
                                   polyLayer = soil,
                                   layerName = "SOIL")
      
      ageVals <- assignCat(gapPoly = bci.gaps17,
                                  polyLayer = ageUse,
                                  layerName = "AgeClass")
      
    # Assume identical across all years
      # 2015 - 2017  
        bci.gaps17$age <- ageVals
        bci.gaps17$soil <- soilVals
      
      # 2017 - 2018  
        bci.gaps18$age <- ageVals
        bci.gaps18$soil <- soilVals
        
      # 2018-2019
        bci.gaps19$age <- ageVals
        bci.gaps19$soil <- soilVals
        
      # 2019-2020
        bci.gaps20$age <- ageVals
        bci.gaps20$soil <- soilVals
    
#### Quantify aspect and height above drainage in each cell ####
      
  # Define a function to take the mean and median topographic variable within each cell
    quantTopo <- function(gapPoly, topoLayer, cellSz, nX, nY){
      
      topoVals <- data.frame(meanVal = rep(NA, length(gapPoly)),
                             medVal = rep(NA, length(gapPoly)))
      
      rasterCrop <- raster::crop(topoLayer, raster::extent(gapPoly))
      resampleMean <- raster::aggregate(rasterCrop, cellSz)
      resampleMedian <- raster::aggregate(rasterCrop, cellSz, fun = median)
      
      valuesMean <- raster::values(resampleMean)
      valuesMedian <- raster::values(resampleMedian)
      
      # reorder to proper order for INLA object
      topoVals$meanVal <- valuesMean[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]
      topoVals$medVal <- valuesMedian[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]
      
      return(topoVals)
    }
      
    # Assume identical across all years
      
      aspectQuant <- quantTopo(gapPoly = bci.gaps17,
                               topoLayer = aspectRaster,
                               cellSz = cellSize,
                               nX = nCellX,
                               nY = nCellY)
      
      bci.gaps17$aspectMean <- aspectQuant$meanVal
      bci.gaps17$aspectMed <- aspectQuant$medVal
      bci.gaps18$aspectMean <- aspectQuant$meanVal
      bci.gaps18$aspectMed <- aspectQuant$medVal
      bci.gaps19$aspectMean <- aspectQuant$meanVal
      bci.gaps19$aspectMed <- aspectQuant$medVal
      bci.gaps20$aspectMean <- aspectQuant$meanVal
      bci.gaps20$aspectMed <- aspectQuant$medVal
      
      drainQuant <- quantTopo(gapPoly = bci.gaps17,
                              topoLayer = drainRaster,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
      
      bci.gaps17$drainMean <- drainQuant$meanVal
      bci.gaps17$drainMed <- drainQuant$medVal
      bci.gaps18$drainMean <- drainQuant$meanVal
      bci.gaps18$drainMed <- drainQuant$medVal
      bci.gaps19$drainMean <- drainQuant$meanVal
      bci.gaps19$drainMed <- drainQuant$medVal
      bci.gaps20$drainMean <- drainQuant$meanVal
      bci.gaps20$drainMed <- drainQuant$medVal

#### Calculate curvature and slope across a range of smoothing values ####

# Vector of different sigmas quantified    
  smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
    
# Curvature  
    
  for(i in 1:length(smoothScales)){
    curvQuant <- quantTopo(gapPoly = bci.gaps18,
                           topoLayer = raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",smoothScales[i],".tif")),
                           cellSz = cellSize,
                           nX = nCellX,
                           nY = nCellY)
    
    if(i==1){
      curvAll <- curvQuant
      names(curvAll)[names(curvAll)%in%c("meanVal","medVal")] <- paste0(c("curvMean","curvMed"),'_',smoothScales[i])
    }
    
    if(i>1){
      curvAll <- cbind(curvAll,curvQuant)
      names(curvAll)[names(curvAll)%in%c("meanVal","medVal")] <- paste0(c("curvMean","curvMed"),'_',smoothScales[i])
    }
    
  }
    
  # Assume identical across all years
    bci.gaps17@data <- cbind(bci.gaps17@data, curvAll)
    bci.gaps18@data <- cbind(bci.gaps18@data, curvAll)
    bci.gaps19@data <- cbind(bci.gaps19@data, curvAll)
    bci.gaps20@data <- cbind(bci.gaps20@data, curvAll)
    
# Slope  
    
    for(i in 1:length(smoothScales)){
      slopeQuant <- quantTopo(gapPoly = bci.gaps18,
                             topoLayer = raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",smoothScales[i],".tif")),
                             cellSz = cellSize,
                             nX = nCellX,
                             nY = nCellY)
      
      if(i==1){
        slopeAll <- slopeQuant
        names(slopeAll)[names(slopeAll)%in%c("meanVal","medVal")] <- paste0(c("slopeMean","slopeMed"),'_',smoothScales[i])
      }
      
      if(i>1){
        slopeAll <- cbind(slopeAll,slopeQuant)
        names(slopeAll)[names(slopeAll)%in%c("meanVal","medVal")] <- paste0(c("slopeMean","slopeMed"),'_',smoothScales[i])
      }
      
    }
    
    # Assume identical across all years
    bci.gaps17@data <- cbind(bci.gaps17@data, slopeAll)
    bci.gaps18@data <- cbind(bci.gaps18@data, slopeAll)
    bci.gaps19@data <- cbind(bci.gaps19@data, slopeAll)
    bci.gaps20@data <- cbind(bci.gaps20@data, slopeAll)
    
#### Combine all years into a single data frame and save ####
    
    bci.gaps17$Year <- "2017"
    bci.gaps18$Year <- "2018"
    bci.gaps19$Year <- "2019"
    bci.gaps20$Year <- "2020"
    
    bci.gapsAll <- rbind(as.data.frame(bci.gaps17),
                         as.data.frame(bci.gaps18),
                         as.data.frame(bci.gaps19),
                         as.data.frame(bci.gaps20))
    
#    save(bci.gaps17, bci.gaps18, bci.gaps19, bci.gaps20, bci.gapsAll, file="INLA/INLA_prelim_40m.RData")
    
#### Smoothing scale sensitivity analysis ####
    
    load("INLA/INLA_prelim_40m.RData")
    
    # Scale 2015-2017 and 2019-2020 values to account for different interval lengths
    bci.gapsAll[bci.gapsAll$Year=="2017","gapProp"] <- (1/2)*bci.gapsAll[bci.gapsAll$Year=="2017","gapProp"]
    bci.gapsAll[bci.gapsAll$Year=="2020","gapProp"] <- (12/13)*bci.gapsAll[bci.gapsAll$Year=="2020","gapProp"]
    
    # Make any rows NA where the area observed is less that 0.2 ha
    bci.gapsAll[!is.na(bci.gapsAll$areaObs) & bci.gapsAll$areaObs<0.12,"gapProp"] <- NA
    
    # Make all covariates NA where the proportion of gaps is NA
    bci.gapsAll[is.na(bci.gapsAll$gapProp),] <- NA
    
    # How many observations are left?
    nrow(bci.gapsAll[!is.na(bci.gapsAll$gapProp),]) # 13130 observations
    
    # Scale topographic covariates
    
    for(i in 8:57){
      bci.gapsAll$new <- NA
      bci.gapsAll[!is.na(bci.gapsAll$gapProp),"new"] <- scale(bci.gapsAll[!is.na(bci.gapsAll$gapProp),i])
      names(bci.gapsAll)[names(bci.gapsAll)=="new"] <- paste0("Sc_",names(bci.gapsAll)[i])
    }
    
    
    # Make a data frame to store results from sensitivity analysis
    # across all smoothing scales for 
    
    smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
    
    topoScaleResults <- data.frame(curvScale = rep(smoothScales,length(smoothScales)),
                                   slopeScale = rep(smoothScales, each = length(smoothScales)),
                                   margLik = NA,
                                   DIC = NA)
    
    # # for 100 m grid
    # nCellX <- 57
    # nCellY <- 53
    
    # for 40 m grid
    nCellX <- 142
    nCellY <- 131
    
    # Make an ID value for each cell
    bci.gapsAll$ID <- rep(1:(nCellX*nCellY),length(unique(bci.gapsAll[!is.na(bci.gapsAll$Year),"Year"])))
    
    # Add censoring value for beta distribution
    cens <- 0.002
    bci.gapsAll$gapPropCens <- bci.gapsAll$gapProp
    bci.gapsAll$gapPropCens[bci.gapsAll$gapPropCens <= cens] <- 0
    bci.gapsAll$gapPropCens[bci.gapsAll$gapPropCens >= 1-cens] <- 1
    
    # Remote NA rows
    bci.gapsAll_Sub <- bci.gapsAll[!is.na(bci.gapsAll$gapProp),]
    
    # Use mean values for topographic variables--with quadratic terms
    
    for(i in 1:nrow(topoScaleResults)){
      
      # For initial model comparison, use only fixed effects
      
      fixed_i <- paste0("Sc_curvMean_",topoScaleResults$curvScale[i]," + I(Sc_curvMean_",topoScaleResults$curvScale[i],"^2) + Sc_slopeMean_",topoScaleResults$slopeScale[i]," + I(Sc_slopeMean_",topoScaleResults$slopeScale[i],"^2) + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
      form_i <- formula(paste0("gapPropCens ~ ",fixed_i))
      #random_i <- "f(ID, model = \"matern2d\", nrow = nCellY, ncol = nCellX)"
      #form_i <- formula(paste0("gapPropCens ~ ",fixed_i," + ",random_i))
      
      model_i <- inla(form_i,
                      family = "beta",
                      data = bci.gapsAll_Sub,
                      control.compute = list(dic = TRUE),
                      control.family = list(beta.censor.value = cens)) 
      
      #summary(model_i)
      topoScaleResults$margLik[i] <- model_i$mlik[2]
      topoScaleResults$DIC[i] <- model_i$dic$dic
      print(i)
    }
    
    write.csv(topoScaleResults, "INLA/INLA_scaleTopo_MeanVals_Quad.csv", row.names = F)
  
    # Plot results
    
      # Read results
      resultsMeanQuad <- read.csv("INLA/INLA_scaleTopo_MeanVals_quad.csv")
    
      # Convert results to matrices--curvature scale along columns, slope scale along hows
      meanQuad_DIC <- matrix(data = resultsMeanQuad$margLik,
                           nrow = length(smoothScales),
                           ncol = length(smoothScales),
                           byrow = F,
                           dimnames = list(smoothScales,smoothScales))
      
      plotMin <- min(meanQuad_DIC)
      plotMax <- max(meanQuad_DIC)
      plotBreaks <- seq(plotMin,plotMax,0.5)
    
      plot(meanQuad_DIC, breaks = plotBreaks,
           ylab = "Curvature scale",
           xlab = "Slope scale")
    
      resultsMeanQuad[which(resultsMeanQuad$DIC==min(resultsMeanQuad$DIC)),]

#### Using best scale result, do variable selection ####
curvScale <- 45
slopeScale <- 63      
      
# STEP 1

  fixed_ref <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref)) 
  
  fixed_a <- paste0("I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a)) 
  
  fixed_b <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b)) 
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c)) 
  
  fixed_d <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d)) 
  
  fixed_e <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e)) 
  
  fixed_f <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + Sc_drainMean + I(Sc_drainMean^2) + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f)) 
  
  fixed_g <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + I(Sc_drainMean^2) + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g)) 
  
  fixed_h <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + Sc_drainMean + Year")
  form_h <- formula(paste0("gapPropCens ~ ",fixed_h)) 
  
  fixed_i <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soil + age + Sc_drainMean + I(Sc_drainMean^2)")
  form_i <- formula(paste0("gapPropCens ~ ",fixed_i)) 

  form_list1 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g, form_h, form_i)
  
  step1_results <- data.frame(model = c("ref","a","b","c","d","e","f","g","h","i"),
                               margLik = NA,
                               DIC = NA)
  
  for(i in 1:nrow(step1_results)){
    
    model_i <- inla(form_list1[[i]],
                    family = "beta",
                    data = bci.gapsAll_Sub,
                    control.compute = list(dic = TRUE),
                    control.family = list(beta.censor.value = cens)) 
    
    #summary(model_i)
    step1_results$margLik[i] <- model_i$mlik[2]
    step1_results$DIC[i] <- model_i$dic$dic
    print(i)
    
  }
  

#### OLD preliminary models ####
    
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
    

#### Look at spatial scales of rasters vs analysis ####
  smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
    
  # Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    

  
  curv <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",45,".tif")) 
  curv <- raster::crop(curv, raster::extent(c(626000,627000,1012500,1013500)))
  
  
  curvMed <- raster::aggregate(curv, fact = 40, fun = median)
  curvMean <- raster::aggregate(curv, fact = 40, fun = mean)
  
  
  raster::plot(curv)

  raster::plot(curvMed)

  raster::plot(curvMean)


  
  slope <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",smoothScales[i],".tif")) 
  
  slope <- raster::mask(slope, buffer)
  
  raster::plot(slope, ext = raster::extent(c(626000,627000,1012500,1013500)))
  raster::plot(bci.gaps, add=T)
#### OLD Calculate expected gaps per cell from area measured ####
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