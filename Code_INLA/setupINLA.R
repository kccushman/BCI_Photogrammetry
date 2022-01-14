#### Load data ####
library(INLA)

  # Gap polygons
    gaps15to18sp <- rgdal::readOGR("Data_GapShapefiles/gaps15to18_shapefile/gaps15to18sp.shp")
    gaps18to20sp <- rgdal::readOGR("Data_GapShapefiles/gaps18to20_shapefile/gaps18to20sp.shp")
  
  # Gap rasters
    gaps15to18 <- raster::raster("Data_HeightRasters/newGaps15to18.tif")
    gaps18to20 <- raster::raster("Data_HeightRasters/newGaps18to20.tif")
    
  # Canopy height change rasters where only likely gap values (> 10 m initially) are included  
    d15to18tall <- raster::raster("Data_HeightRasters/dCHM15to18tall.tif")
    d18to20tall <- raster::raster("Data_HeightRasters/dCHM18to20tall.tif")
   
  # Canopy height change raster omitting two very large gaps in second interval (7190 and 15124 m2; largest in 15to18 was 2742)
    # find gap id of big gaps
    bigIDs <- gaps18to20sp$gap_id[gaps18to20sp$area>3000]
    
    # remove these gaps from base raster  
    d18to20tall_alt <- raster::mask(d18to20tall, gaps18to20sp[gaps18to20sp$gap_id%in%bigIDs,], inverse=T)

  # Canopy height raster at the beginning of each interval
    chm15 <- raster::raster("Data_HeightRasters/CHM_2015_QAQC.tif")
    chm18 <- raster::raster("Data_HeightRasters/CHM_2018_QAQC.tif")
    
    
  # Forest age polygon
    age <- rgdal::readOGR("Data_Ancillary/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
      age$AgeClass <- "Other"
      age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
      age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
      ageUse <- age[!(age$AgeClass=="Other"),]
      
  # Soil type polygon  
    soil <- rgdal::readOGR("Data_Ancillary/BCI_Soils/BCI_Soils.shp")
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

  # Aspect raster
    aspectRaster <- raster::raster("Data_TopographyRasters/Aspect_smooth_21.tif")
    # resample to same extent as gap rasters (adds NA area to edges)
    aspectRaster <- raster::resample(aspectRaster, gaps18to20)
    
  # Distance above drainage raster
    drainRaster <- raster::raster("Data_TopographyRasters/distAboveStream_1000.tif")
    # resample to same extent as gap rasters (adds NA area to edges)
    drainRaster <- raster::resample(drainRaster, gaps18to20)
    drainRasterSq <- drainRaster^2

#### Create SpatialPolygonsDataFrame across BCI with proper order for R-INLA ####
    
    # Define cell size (in m)
    cellSize <- 40
    
    # Count number of cells in each dimension
    nCellX <- ceiling((d15to18tall@extent@xmax-d15to18tall@extent@xmin)/cellSize)
    nCellY <- ceiling((d15to18tall@extent@ymax-d15to18tall@extent@ymin)/cellSize)
    
    # Define a SpatialGridDataFrame
    bci.pix <- spatstat::im(mat = 1:(nCellX*nCellY),
                            xcol = d15to18tall@extent@xmin + 0.5*cellSize + cellSize*(0:(nCellX-1)),
                            yrow = d15to18tall@extent@ymin + 0.5*cellSize + cellSize*(0:(nCellY-1)))
    
    bci.grid <- as(bci.pix,"SpatialGridDataFrame")
    
    # Convert to SpatialPolygons
    bci.poly <- as(bci.grid, "SpatialPolygons")
    
    # Convert gap polygons to points
    gapsPts15to18 <- sp::SpatialPoints(coords = gaps15to18sp@data[,c("X1","X2")])
    gapsPts18to20 <- sp::SpatialPoints(coords = gaps18to20sp@data[,c("X1","X2")])
    
    # 2015 - 2018
      #Number of gap observations per cell
      idx <- over(gapsPts15to18, bci.poly)
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
      
    
    # 2018 - 2020
      #Number of gap observations per cell
      idx <- over(gapsPts18to20, bci.poly)
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
        baseCropNew <- raster::values(baseCrop)
        baseCropNew[!is.na(baseCropNew)] <- 1
        raster::values(baseCrop) <- baseCropNew
        
        gapCrop <- raster::crop(gapLayer, raster::extent(gapPoly))
        gapCropNew <- raster::values(gapCrop)
        gapCropNew[!is.na(gapCropNew)] <- 1
          # omit values that are NA in the base (tall) layer (i.e. don't count gaps that are initially < 10 m)
          gapCropNew[is.na(baseCropNew)] <- NA
        raster::values(gapCrop) <- gapCropNew
          
        
        resampleBase <- raster::aggregate(baseCrop, cellSz, fun = sum)
        resampleGap <- raster::aggregate(gapCrop, cellSz, fun = sum)
        
        valuesBase <- raster::values(resampleBase)
        valuesGap <- raster::values(resampleGap)
        valuesGap[is.na(valuesGap)] <- 0
        
        # reorder to proper order for INLA object
        gapVals$sampHa <- valuesBase[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]/10000
        gapVals$gapProp <- valuesGap[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]/valuesBase[as.vector(t(matrix(1:(nX*nY), nrow = nX, ncol = nY)))]
        
        return(gapVals)
      }

    # 2015-2018

      gapArea <- propGap(gapPoly = bci.gaps18,
                         gapLayer = gaps15to18,
                         baseLayer = d15to18tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
      bci.gaps18$areaObs <- gapArea$sampHa
      bci.gaps18$gapProp <- gapArea$gapProp

    # 2018-2020
  
      gapArea <- propGap(gapPoly = bci.gaps20,
                         gapLayer = gaps18to20,
                         baseLayer = d18to20tall,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
      bci.gaps20$areaObs <- gapArea$sampHa
      bci.gaps20$gapProp <- gapArea$gapProp
      
    # 2018-2020 (omitting two very large gaps)
      
      gapArea <- propGap(gapPoly = bci.gaps20,
                         gapLayer = gaps18to20,
                         baseLayer = d18to20tall_alt,
                         cellSz = cellSize,
                         nX = nCellX,
                         nY = nCellY)
      bci.gaps20$areaObs_alt <- gapArea$sampHa
      bci.gaps20$gapProp_alt <- gapArea$gapProp
      
      # make colums for first interval too
      bci.gaps18$areaObs_alt <- bci.gaps18$areaObs
      bci.gaps18$gapProp_alt <- bci.gaps18$gapProp
      
    # Normalize the proportion of gaps observed to per year
      nYr15to18 <- as.numeric(as.Date("2018-06-07") - as.Date("2015-06-26"))/365
      nYr18to20 <- as.numeric(as.Date("2020-07-31") - as.Date("2018-06-07"))/365
      
      bci.gaps18$gapProp <- bci.gaps18$gapProp/nYr15to18
      bci.gaps20$gapProp <- bci.gaps20$gapProp/nYr18to20
      bci.gaps18$gapProp_alt <- bci.gaps18$gapProp_alt/nYr15to18
      bci.gaps20$gapProp_alt <- bci.gaps20$gapProp_alt/nYr18to20
      
#### Assign forest age, soil parent material, and soil form to each cell ####

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
      
      soilParentVals <- assignCat(gapPoly = bci.gaps18,
                                  polyLayer = soil,
                                  layerName = "SoilParent")
      
      soilFormVals <- assignCat(gapPoly = bci.gaps18,
                                polyLayer = soil,
                                layerName = "SoilForm")
      
      ageVals <- assignCat(gapPoly = bci.gaps18,
                           polyLayer = ageUse,
                           layerName = "AgeClass")
      
    # Assume identical across all years
      # 2015 - 2018  
        bci.gaps18$age <- ageVals
        bci.gaps18$soilParent <- soilParentVals
        bci.gaps18$soilForm <- soilFormVals
        
      
      # 2019-2020
        bci.gaps20$age <- ageVals
        bci.gaps20$soilParent <- soilParentVals
        bci.gaps20$soilForm <- soilFormVals
        
      # Omit cells where age is "Other"
        bci.gaps18[bci.gaps18$age == "Other" & !is.na(bci.gaps18$age),c("age","gapProp")] <- NA
        bci.gaps20[bci.gaps20$age == "Other" & !is.na(bci.gaps20$age),c("age","gapProp")] <- NA
        
#### Quantify aspect, height above drainage, and initial canopy height in each cell ####
      
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
      
      aspectQuant <- quantTopo(gapPoly = bci.gaps18,
                               topoLayer = aspectRaster,
                               cellSz = cellSize,
                               nX = nCellX,
                               nY = nCellY)
      
      bci.gaps18$aspectMean <- aspectQuant$meanVal
      bci.gaps20$aspectMean <- aspectQuant$meanVal

      
    # Height above nearest drainage--linear
      drainQuant <- quantTopo(gapPoly = bci.gaps18,
                              topoLayer = drainRaster,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
      
      bci.gaps18$drainMean <- drainQuant$meanVal
      bci.gaps20$drainMean <- drainQuant$meanVal

    # Height above nearest drainage--squared
      drainQuantSq <- quantTopo(gapPoly = bci.gaps18,
                              topoLayer = drainRasterSq,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
      
      bci.gaps18$drainMean_Sq <- drainQuantSq$meanVal
      bci.gaps20$drainMean_Sq <- drainQuantSq$meanVal


    # Initial canopy height -- linear
      canopyHt15 <- quantTopo(gapPoly = bci.gaps18,
                              topoLayer = chm15,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
      
      bci.gaps18$initialHt <- canopyHt15$meanVal
      
      canopyHt18 <- quantTopo(gapPoly = bci.gaps18,
                              topoLayer = chm18,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
      
      bci.gaps20$initialHt <- canopyHt18$meanVal

#### Calculate curvature and slope across a range of smoothing values ####

# Vector of different sigmas quantified    
  smoothScales <- c(1,2,3,4,6,8,12,16,24,32,48,64)
    
# Curvature  
    
  for(i in 1:length(smoothScales)){
    curvQuant <- quantTopo(gapPoly = bci.gaps18,
                           topoLayer = raster::raster(paste0("Data_TopographyRasters/Curv_smooth_",smoothScales[i],".tif")),
                           cellSz = cellSize,
                           nX = nCellX,
                           nY = nCellY)
    
    curvQuant_Sq <- quantTopo(gapPoly = bci.gaps18,
                              topoLayer = (raster::raster(paste0("Data_TopographyRasters/Curv_smooth_",smoothScales[i],".tif")))^2,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
    
    if(i==1){
      curvAll <- curvQuant
      names(curvAll)[names(curvAll)%in%c("meanVal","medVal")] <- paste0(c("curvMean","curvMed"),'_',smoothScales[i])
      
      curvAll_Sq <- curvQuant_Sq
      names(curvAll_Sq)[names(curvAll_Sq)%in%c("meanVal","medVal")] <- paste0(c("curvMean","curvMed"),'_',smoothScales[i],"_Sq")
    }
    
    if(i>1){
      curvAll <- cbind(curvAll,curvQuant)
      names(curvAll)[names(curvAll)%in%c("meanVal","medVal")] <- paste0(c("curvMean","curvMed"),'_',smoothScales[i])
      
      curvAll_Sq <- cbind(curvAll_Sq,curvQuant_Sq)
      names(curvAll_Sq)[names(curvAll_Sq)%in%c("meanVal","medVal")] <- paste0(c("curvMean","curvMed"),'_',smoothScales[i],"_Sq")
    }
    
  }
    
  # Assume identical across all years
    bci.gaps18@data <- cbind(bci.gaps18@data, curvAll, curvAll_Sq)
    bci.gaps20@data <- cbind(bci.gaps20@data, curvAll, curvAll_Sq)
    
# Slope  
    
    for(i in 1:length(smoothScales)){
      slopeQuant <- quantTopo(gapPoly = bci.gaps18,
                             topoLayer = raster::raster(paste0("Data_TopographyRasters/Slope_smooth_",smoothScales[i],".tif")),
                             cellSz = cellSize,
                             nX = nCellX,
                             nY = nCellY)
      
      slopeQuant_Sq <- quantTopo(gapPoly = bci.gaps18,
                                topoLayer = (raster::raster(paste0("Data_TopographyRasters/Slope_smooth_",smoothScales[i],".tif")))^2,
                                cellSz = cellSize,
                                nX = nCellX,
                                nY = nCellY)
      
      if(i==1){
        slopeAll <- slopeQuant
        names(slopeAll)[names(slopeAll)%in%c("meanVal","medVal")] <- paste0(c("slopeMean","slopeMed"),'_',smoothScales[i])
        
        slopeAll_Sq <- slopeQuant_Sq
        names(slopeAll_Sq)[names(slopeAll_Sq)%in%c("meanVal","medVal")] <- paste0(c("slopeMean","slopeMed"),'_',smoothScales[i],"_Sq")
      }
      
      if(i>1){
        slopeAll <- cbind(slopeAll,slopeQuant)
        names(slopeAll)[names(slopeAll)%in%c("meanVal","medVal")] <- paste0(c("slopeMean","slopeMed"),'_',smoothScales[i])
        
        slopeAll_Sq <- cbind(slopeAll_Sq,slopeQuant_Sq)
        names(slopeAll_Sq)[names(slopeAll_Sq)%in%c("meanVal","medVal")] <- paste0(c("slopeMean","slopeMed"),'_',smoothScales[i],"_Sq")
      }
      
    }
    
    # Assume identical across all years
    bci.gaps18@data <- cbind(bci.gaps18@data, slopeAll, slopeAll_Sq)
    bci.gaps20@data <- cbind(bci.gaps20@data, slopeAll, slopeAll_Sq)
    
#### Combine all years into a single data frame and save ####
    
    # Add year
    bci.gaps18$Year <- "2018"
    bci.gaps20$Year <- "2020"
    
    # Add censoring value for beta distribution
    cens <- 0.0001
    
    bci.gaps18$gapPropCens <- bci.gaps18$gapProp
    bci.gaps18$gapPropCens[bci.gaps18$gapPropCens <= cens] <- 0
    bci.gaps18$gapPropCens[bci.gaps18$gapPropCens >= 1-cens] <- 1
    
    bci.gaps20$gapPropCens <- bci.gaps20$gapProp
    bci.gaps20$gapPropCens[bci.gaps20$gapPropCens <= cens] <- 0
    bci.gaps20$gapPropCens[bci.gaps20$gapPropCens >= 1-cens] <- 1
    
    bci.gaps18$gapPropCens_alt <- bci.gaps18$gapProp_alt
    bci.gaps18$gapPropCens_alt[bci.gaps18$gapPropCens_alt <= cens] <- 0
    bci.gaps18$gapPropCens_alt[bci.gaps18$gapPropCens_alt >= 1-cens] <- 1
    
    bci.gaps20$gapPropCens_alt <- bci.gaps20$gapProp_alt
    bci.gaps20$gapPropCens_alt[bci.gaps20$gapPropCens_alt <= cens] <- 0
    bci.gaps20$gapPropCens_alt[bci.gaps20$gapPropCens_alt >= 1-cens] <- 1
    
    # Make any rows NA where the area observed is less than half the pixel
    bci.gaps18[!is.na(bci.gaps18$areaObs) & bci.gaps18$areaObs<(0.5*cellSize*cellSize/10000),"gapPropCens"] <- NA
    bci.gaps20[!is.na(bci.gaps20$areaObs) & bci.gaps20$areaObs<(0.5*cellSize*cellSize/10000),"gapPropCens"] <- NA
    
    bci.gaps18[!is.na(bci.gaps18$areaObs_alt) & bci.gaps18$areaObs_alt<(0.5*cellSize*cellSize/10000),"gapPropCens_alt"] <- NA
    bci.gaps20[!is.na(bci.gaps20$areaObs_alt) & bci.gaps20$areaObs_alt<(0.5*cellSize*cellSize/10000),"gapPropCens_alt"] <- NA
    
    # Combine years
    bci.gapsAll <- rbind(as.data.frame(bci.gaps18),
                         as.data.frame(bci.gaps20))
    
    # How many observations are left?
    nrow(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),]) # 17059 observations
    nrow(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens_alt),]) # 17043 observations (omitting large gaps)
    
    # Remove median values from topographic covariates (not used, keep mean)
    medCols <- grepl(pattern="Med",names(bci.gapsAll))
    bci.gapsAll <- bci.gapsAll[,!medCols]
    
    # Scale topographic covariates
    for(i in c(10:60)){
      bci.gapsAll$new <- NA
      bci.gapsAll[!is.na(bci.gapsAll$age),"new"] <- scale(bci.gapsAll[!is.na(bci.gapsAll$age),i])
      names(bci.gapsAll)[names(bci.gapsAll)=="new"] <- paste0("Sc_",names(bci.gapsAll)[i])
    }


    save(bci.gaps18, bci.gaps20, bci.gapsAll, cellSize, nCellX, nCellY, cens, file="Data_INLA/INLA_40m.RData")
    
#### Smoothing scale sensitivity analysis ####
    
    # Run commented code if returning to script (not running from beginning to format data) 
    # library(INLA)
    # load("Data_INLA/INLA_40m.RData")
 
    
    # Make a data frame to store results from sensitivity analysis
    # across all smoothing scales for curvature and slope
    
    smoothScales <- c(1,2,3,4,6,8,12,16,24,32,48,64)
    
    topoScaleResults <- data.frame(curvScale = rep(smoothScales,length(smoothScales)),
                                   slopeScale = rep(smoothScales, each = length(smoothScales)),
                                   margLik = NA,
                                   DIC = NA)
    
    # Remove NA rows
    bci.gapsAll_Sub <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),]
    
    # Use mean values for topographic variables--with quadratic terms
    
    for(i in 1:nrow(topoScaleResults)){
      
      # For initial model comparison, use only fixed effects
      
      fixed_i <- paste0("Sc_curvMean_",topoScaleResults$curvScale[i]," + Sc_curvMean_",topoScaleResults$curvScale[i],"_Sq + Sc_slopeMean_",topoScaleResults$slopeScale[i]," + Sc_slopeMean_",topoScaleResults$slopeScale[i],"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
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
    }
    
    write.csv(topoScaleResults, "Data_INLA/INLA_SmoothingScaleResults.csv", row.names = F)
 
    # FIGURE WITH THESE RESULTS IN FILE MakeFigures.R

    # Best scale for curvature is 2, best scale for slope is 16
    resultsMeanQuad[which(resultsMeanQuad$DIC==min(resultsMeanQuad$DIC)),]

#### Using best scale result, do initial variable selection without Matern spatial autocorrelation term ####
  curvScale <- 2
  slopeScale <- 16
  
  fixed_ref <- paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref))
  
  fixed_a <-  paste0("Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a))
  
  fixed_b <-  paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b))
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c))
  
  fixed_d <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d))
  
  fixed_e <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e))
  
  fixed_f <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + soilParent + soilForm + age + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f))
  
  fixed_g <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilForm + age + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g))
  
  fixed_h <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + age + Year")
  form_h <- formula(paste0("gapPropCens ~ ",fixed_h))
  
  fixed_i <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + Year")
  form_i <- formula(paste0("gapPropCens ~ ",fixed_i))

  form_list1 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g, form_h, form_i)
  
  step1_results <- data.frame(model = c("ref","a","b","c","d","e","f","g","h","i"),
                               margLik = NA,
                               DIC = NA)
  
  for(i in 1:nrow(step1_results)){

    tryCatch(model_i <- inla(form_list1[[i]],
                    family = "beta",
                    data = bci.gapsAll,
                    control.compute = list(dic = TRUE),
                    control.family = list(beta.censor.value = cens),
                    verbose = F),
             error = function(e) {print(paste("iteration",i,"failed"))})
    
    #summary(model_i)
    step1_results$margLik[i] <- model_i$mlik[2]
    step1_results$DIC[i] <- model_i$dic$dic
    print(i)
    
  }
  
  step1_results$dDIC <- step1_results$DIC - min(step1_results$DIC,na.rm=T)
  step1_results <- step1_results[order(step1_results$dDIC),]
  

#### Using best scale result, do final variable selection with Matern spatial autocorrelation term ####
  curvScale <- 2
  slopeScale <- 16
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  # Define random effects matrix
  random_mat <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  
  # STEP 1

  fixed_ref <- paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref," + ",random_mat))
  
  fixed_a <-  paste0("Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a," + ",random_mat))
  
  fixed_b <-  paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b," + ",random_mat))
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c," + ",random_mat))
  
  fixed_d <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d," + ",random_mat))
  
  fixed_e <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e," + ",random_mat))
  
  fixed_f <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + soilParent + soilForm + age + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f," + ",random_mat))
  
  fixed_g <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilForm + age + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g," + ",random_mat))
  
  fixed_h <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + age + Year")
  form_h <- formula(paste0("gapPropCens ~ ",fixed_h," + ",random_mat))
  
  fixed_i <-  paste0("Sc_curvMean_",curvScale," + Sc_curvMean_",curvScale,"_Sq + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + Year")
  form_i <- formula(paste0("gapPropCens ~ ",fixed_i," + ",random_mat))
  
  form_list1 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g, form_h, form_i)
  
  step1b_results <- data.frame(model = c("ref","a","b","c","d","e","f","g","h","i"),
                              margLik = NA,
                              DIC = NA)
  
  for(i in 1:nrow(step1b_results)){
    
    tryCatch(model_i <- inla(form_list1[[i]],
                             family = "beta",
                             data = bci.gapsAll_Order,
                             control.compute = list(dic = TRUE),
                             control.family = list(beta.censor.value = cens),
                             verbose = F),
             error = function(e) {print(paste("iteration",i,"failed"))})
    
    #summary(model_i)
    step1b_results$margLik[i] <- model_i$mlik[2]
    step1b_results$DIC[i] <- model_i$dic$dic
    print(i)
    
  }
  
  step1b_results$dDIC <- step1b_results$DIC - min(step1b_results$DIC,na.rm=T)
  step1b_results <- step1b_results[order(step1b_results$dDIC),]
  
  # STEP 2
  
  fixed_ref <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref," + ",random_mat))
  
  fixed_a <- paste0("Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a," + ",random_mat))
  
  fixed_b <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b," + ",random_mat))
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c," + ",random_mat))
  
  fixed_d <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d," + ",random_mat))
  
  fixed_e <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + soilParent + soilForm + age + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e," + ",random_mat))
  
  fixed_f <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilForm + age + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f," + ",random_mat))
  
  fixed_g <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + age + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g," + ",random_mat))
  
  fixed_h <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + Year")
  form_h <- formula(paste0("gapPropCens ~ ",fixed_h," + ",random_mat))
  
  
  form_list2 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g, form_h)
  
  step2b_results <- data.frame(model = c("ref","a","b","c","d","e","f","g","h"),
                              margLik = NA,
                              DIC = NA)
  
  for(i in 1:nrow(step2b_results)){
    
    tryCatch(model_i <- inla(form_list2[[i]],
                             family = "beta",
                             data = bci.gapsAll_Order,
                             control.compute = list(dic = TRUE),
                             control.family = list(beta.censor.value = cens),
                             verbose = F),
             error = function(e) {print(paste("iteration",i,"failed"))})
    
    #summary(model_i)
    step2b_results$margLik[i] <- model_i$mlik[2]
    step2b_results$DIC[i] <- model_i$dic$dic
    print(i)
  }
  
  step2b_results$dDIC <- step2b_results$DIC - min(step2b_results$DIC,na.rm=T)
  step2b_results <- step2b_results[order(step2b_results$dDIC),]

  # # STEP 3
  # 
  # fixed_ref <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean + Sc_drainMean_sq + soilParent + age + Year")
  # form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref," + ",random_mat))
  # 
  # fixed_a <- paste0("Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean + Sc_drainMean_sq + soilParent + age + Year")
  # form_a <- formula(paste0("gapPropCens ~ ",fixed_a," + ",random_mat))
  # 
  # fixed_b <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean + Sc_drainMean_sq + soilParent + age + Year")
  # form_b <- formula(paste0("gapPropCens ~ ",fixed_b," + ",random_mat))
  # 
  # fixed_c <-  paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_drainMean + Sc_drainMean_sq + soilParent + age + Year")
  # form_c <- formula(paste0("gapPropCens ~ ",fixed_c," + ",random_mat))
  # 
  # fixed_d <-   paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean_sq + soilParent + age + Year")
  # form_d <- formula(paste0("gapPropCens ~ ",fixed_d," + ",random_mat))
  # 
  # fixed_e <-  paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean + soilParent + age + Year")
  # form_e <- formula(paste0("gapPropCens ~ ",fixed_e," + ",random_mat))
  # 
  # fixed_f <-   paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean + Sc_drainMean_sq + age + Year")
  # form_f <- formula(paste0("gapPropCens ~ ",fixed_f," + ",random_mat))
  # 
  # fixed_g <-  paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_sq + Sc_drainMean + Sc_drainMean_sq + soilParent + Year")
  # form_g <- formula(paste0("gapPropCens ~ ",fixed_g," + ",random_mat))
  # 
  # form_list3 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g)
  # 
  # step3b_results <- data.frame(model = c("ref","a","b","c","d","e","f","g"),
  #                              margLik = NA,
  #                              DIC = NA)
  # 
  # for(i in 1:nrow(step3b_results)){
  #   
  #   tryCatch(model_i <- inla(form_list3[[i]],
  #                            family = "beta",
  #                            data = bci.gapsAll_Order,
  #                            control.compute = list(dic = TRUE),
  #                            control.family = list(beta.censor.value = cens),
  #                            verbose = F),
  #            error = function(e) {print(paste("iteration",i,"failed"))})
  #   
  #   #summary(model_i)
  #   step3b_results$margLik[i] <- model_i$mlik[2]
  #   step3b_results$DIC[i] <- model_i$dic$dic
  #   print(i)
  # }
  # 
  # step3b_results$dDIC <- step3b_results$DIC - min(step3b_results$DIC,na.rm=T)
  # step3b_results <- step3b_results[order(step3b_results$dDIC),]    

#### Run full model with Matern spatial autocorrelation term ####
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  # Reorder factors so that "base" level is the group with the most data
  bci.gapsAll_Order$soilParent <- relevel(as.factor(bci.gapsAll_Order$soilParent), "Bohio")
  bci.gapsAll_Order$soilForm <- relevel(as.factor(bci.gapsAll_Order$soilForm), "BrownFineLoam")
  bci.gapsAll_Order$age <- relevel(as.factor(bci.gapsAll_Order$age), "OldGrowth")
  bci.gapsAll_Order$Year <- relevel(as.factor(bci.gapsAll_Order$Year), "2020")
  
  
  # Run the best (and next most simple) models with full spatial autocorrelation
  fixed_full <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_full <- formula(paste0("gapPropCens ~ ",fixed_full," + ",random_full))
  
  model_full <- inla(form_full,
                     family = "beta",
                     data = bci.gapsAll_Order,
                     control.compute = list(dic = TRUE),
                     control.family = list(beta.censor.value = cens))
  
  save(model_full, file = "Data_INLA/INLA_fullModelResult.RData")
  
#### Run alternate models separately for each year ####
  
  # Run commented code if returning to script (not running from beginning to format data) 
  # library(INLA)
  # load("Data_INLA/INLA_40m.RData")
  
  curvScale <- 2
  slopeScale <- 16
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  # Reorder factors so that "base" level is the group with the most data
  bci.gapsAll$soilParent <- relevel(as.factor(bci.gapsAll$soilParent), "Bohio")
  bci.gapsAll$soilForm <- relevel(as.factor(bci.gapsAll$soilForm), "BrownFineLoam")
  bci.gapsAll$age <- relevel(as.factor(bci.gapsAll$age), "OldGrowth")
  bci.gapsAll$Year <- relevel(as.factor(bci.gapsAll$Year), "2020")
  
  
  # Make new data frames for each year
  bci.gapsAll1 <- bci.gapsAll[bci.gapsAll$Year=="2018",]
  bci.gapsAll2 <- bci.gapsAll[bci.gapsAll$Year=="2020",]
  bci.gapsAll1$ID <- 1:nrow(bci.gapsAll1)
  bci.gapsAll2$ID <- 1:nrow(bci.gapsAll2)

  # Define model form, without year term
  fixed_full_sep <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age")
  random_full_sep <- "f(ID, model = \"matern2d\", nrow = nCellY, ncol = nCellX)"
  form_full_sep <- formula(paste0("gapPropCens ~ ",fixed_full_sep," + ",random_full_sep))
  
  model_full1 <- inla(form_full_sep,
                     family = "beta",
                     data = bci.gapsAll1,
                     control.compute = list(dic = TRUE),
                     control.family = list(beta.censor.value = cens))
  model_full2 <- inla(form_full_sep,
                      family = "beta",
                      data = bci.gapsAll2,
                      control.compute = list(dic = TRUE),
                      control.family = list(beta.censor.value = cens))
  
  save(model_full1, model_full2, file = "Data_INLA/INLA_fullModelResult_separate.RData")
  
#### Run alternate model for 2018-2020 without biggest gaps ####
  
  # Run commented code if returning to script (not running from beginning to format data) 
  # library(INLA)
  # load("Data_INLA/INLA_40m.RData")
  
  curvScale <- 2
  slopeScale <- 16
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  # Reorder factors so that "base" level is the group with the most data
  bci.gapsAll_Order$soilParent <- relevel(as.factor(bci.gapsAll_Order$soilParent), "Bohio")
  bci.gapsAll_Order$soilForm <- relevel(as.factor(bci.gapsAll_Order$soilForm), "BrownFineLoam")
  bci.gapsAll_Order$age <- relevel(as.factor(bci.gapsAll_Order$age), "OldGrowth")
  bci.gapsAll_Order$Year <- relevel(as.factor(bci.gapsAll_Order$Year), "2020")
  
  # Make new data frames for each year
  bci.gapsAll1 <- bci.gapsAll[bci.gapsAll$Year=="2018",]
  bci.gapsAll2 <- bci.gapsAll[bci.gapsAll$Year=="2020",]
  bci.gapsAll1$ID <- 1:nrow(bci.gapsAll1)
  bci.gapsAll2$ID <- 1:nrow(bci.gapsAll2)
  
  # Run the best (and next most simple) models with full spatial autocorrelation
  fixed_full <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_full_alt <- formula(paste0("gapPropCens_alt ~ ",fixed_full," + ",random_full))
  
  model_full_alt <- inla(form_full_alt,
                         family = "beta",
                         data = bci.gapsAll_Order,
                         control.compute = list(dic = TRUE),
                         control.family = list(beta.censor.value = cens))
  
  # Define model form, without year term
  fixed_full_sep <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age")
  random_full_sep <- "f(ID, model = \"matern2d\", nrow = nCellY, ncol = nCellX)"
  form_full_sepAlt <- formula(paste0("gapPropCens_alt ~ ",fixed_full_sep," + ",random_full_sep))
  
  model_full2_alt <- inla(form_full_sepAlt,
                      family = "beta",
                      data = bci.gapsAll2,
                      control.compute = list(dic = TRUE),
                      control.family = list(beta.censor.value = cens))
  
  save(model_full_alt, model_full2_alt, file = "Data_INLA/INLA_fullModelResult_noLargeGaps.RData")
  
#### Run alternate model with initial canopy height ####
  
 # Run commented code if returning to script (not running from beginning to format data) 
 # library(INLA)
 # load("Data_INLA/INLA_40m.RData")
  
  curvScale <- 2
  slopeScale <- 16
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  # Reorder factors so that "base" level is the group with the most data
  bci.gapsAll_Order$soilParent <- relevel(as.factor(bci.gapsAll_Order$soilParent), "Bohio")
  bci.gapsAll_Order$soilForm <- relevel(as.factor(bci.gapsAll_Order$soilForm), "BrownFineLoam")
  bci.gapsAll_Order$age <- relevel(as.factor(bci.gapsAll_Order$age), "OldGrowth")
  bci.gapsAll_Order$Year <- relevel(as.factor(bci.gapsAll_Order$Year), "2020")
  
  # Run the best models with full spatial autocorrelation AND add initial canopy height
  fixed_full_ht <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + soilParent + soilForm + age + Year + Sc_initialHt")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_full_ht <- formula(paste0("gapPropCens ~ ",fixed_full_ht," + ",random_full))
  
  model_full_ht <- inla(form_full_ht,
                     family = "beta",
                     data = bci.gapsAll_Order,
                     control.compute = list(dic = TRUE),
                     control.family = list(beta.censor.value = cens))
  
  save(model_full_ht, file = "Data_INLA/INLA_fullModelResult_initialHt.RData")
  
#### Run alternate models isolating soil, topography, and age terms ####
  
  # Run commented code if returning to script (not running from beginning to format data) 
  # library(INLA)
  # load("Data_INLA/INLA_40m.RData")
  
  curvScale <- 2
  slopeScale <- 16
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Interval 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    # Interval 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  # Reorder factors so that "base" level is the group with the most data
  bci.gapsAll_Order$soilParent <- relevel(as.factor(bci.gapsAll_Order$soilParent), "Bohio")
  bci.gapsAll_Order$soilForm <- relevel(as.factor(bci.gapsAll_Order$soilForm), "BrownFineLoam")
  bci.gapsAll_Order$age <- relevel(as.factor(bci.gapsAll_Order$age), "OldGrowth")
  bci.gapsAll_Order$Year <- relevel(as.factor(bci.gapsAll_Order$Year), "2020")
  
  # Model with only forest age and year
  fixed_ageOnly <- paste0("age + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_ageOnly <- formula(paste0("gapPropCens ~ ",fixed_ageOnly," + ",random_full))
  
  model_ageOnly <- inla(form_ageOnly,
                        family = "beta",
                        data = bci.gapsAll_Order,
                        control.compute = list(dic = TRUE),
                        control.family = list(beta.censor.value = cens))
  
  save(model_ageOnly, file = "Data_INLA/INLA_ModelResult_ageOnly.RData")
  

  
  # Model with only soils terms
  fixed_soilsOnly <- paste0("soilParent + soilForm + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_soilsOnly <- formula(paste0("gapPropCens ~ ",fixed_soilsOnly," + ",random_full))
  
  model_soilsOnly <- inla(form_soilsOnly,
                       family = "beta",
                       data = bci.gapsAll_Order,
                       control.compute = list(dic = TRUE),
                       control.family = list(beta.censor.value = cens))
  
  save(model_soilsOnly, file = "Data_INLA/INLA_ModelResult_soilsOnly.RData")
  
  # Model with only topography terms
  fixed_topoOnly <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_slopeMean_",slopeScale,"_Sq + Sc_drainMean + Sc_drainMean_Sq + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_topoOnly <- formula(paste0("gapPropCens ~ ",fixed_topoOnly," + ",random_full))
  
  model_topoOnly <- inla(form_topoOnly,
                          family = "beta",
                          data = bci.gapsAll_Order,
                          control.compute = list(dic = TRUE),
                          control.family = list(beta.censor.value = cens))
  
  save(model_topoOnly, file = "Data_INLA/INLA_ModelResult_topoOnly.RData")