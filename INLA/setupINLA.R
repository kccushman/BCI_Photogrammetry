#### Load data ####
library(INLA)

  # Gap polygons
    gaps15to18sp <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
    gaps18to20sp <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")
  
  # Gap rasters
    gaps15to18 <- raster::raster("newGaps15to18_tin.tif")
    gaps18to20 <- raster::raster("newGaps18to20_tin.tif")
    
  # Canopy height change rasters where only possible gap values (> 5 m initially) are included  
    d15to18tall <- raster::raster("dCHM15to18tall_tin.tif")
    d18to20tall <- raster::raster("dCHM18to20tall_tin.tif")
    
  # Forest age polygon
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
      ageUse <- age[!(age$TYPE=="Clearings"),]
      ageUse$AgeClass <- "Other"
      ageUse$AgeClass[ageUse$Mascaro_Co == "> 400"] <- "OldGrowth"
      ageUse$AgeClass[ageUse$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
      
  # Soil type polygon  
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
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
    aspectRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Aspect_smooth_21.tif")
    
  # Distance above drainage raster
    drainRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/distAboveStream_1000.tif")
    # resample to same extent as gap rasters (adds NA area to edges)
    drainRaster <- raster::resample(drainRaster, gaps18to20)
    
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
    gapsPts15to18 <- sp::SpatialPoints(coords = gaps15to18sp@data[gaps15to18sp$use==T,c("X1","X2")])
    gapsPts18to20 <- sp::SpatialPoints(coords = gaps18to20sp@data[gaps18to20sp$use==T,c("X1","X2")])
    
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
      
      aspectQuant <- quantTopo(gapPoly = bci.gaps18,
                               topoLayer = aspectRaster,
                               cellSz = cellSize,
                               nX = nCellX,
                               nY = nCellY)
      
      bci.gaps18$aspectMean <- aspectQuant$meanVal
      bci.gaps18$aspectMed <- aspectQuant$medVal
      bci.gaps20$aspectMean <- aspectQuant$meanVal
      bci.gaps20$aspectMed <- aspectQuant$medVal
      
      drainQuant <- quantTopo(gapPoly = bci.gaps18,
                              topoLayer = drainRaster,
                              cellSz = cellSize,
                              nX = nCellX,
                              nY = nCellY)
      
      bci.gaps18$drainMean <- drainQuant$meanVal
      bci.gaps18$drainMed <- drainQuant$medVal
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
    bci.gaps18@data <- cbind(bci.gaps18@data, curvAll)
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
    bci.gaps18@data <- cbind(bci.gaps18@data, slopeAll)
    bci.gaps20@data <- cbind(bci.gaps20@data, slopeAll)
    
#### Combine all years into a single data frame and save ####
    
    bci.gaps18$Year <- "2018"
    bci.gaps20$Year <- "2020"
    
    bci.gapsAll <- rbind(as.data.frame(bci.gaps18),
                         as.data.frame(bci.gaps20))
    
    # Scale 2015-2018 and 2018-2020 values to account for different interval lengths
    bci.gapsAll[bci.gapsAll$Year=="2018","gapProp"] <- (1/3)*bci.gapsAll[bci.gapsAll$Year=="2018","gapProp"]
    bci.gapsAll[bci.gapsAll$Year=="2020","gapProp"] <- (12/25)*bci.gapsAll[bci.gapsAll$Year=="2020","gapProp"]
    
    # Make any rows NA where the area observed is less than half the pixel
    bci.gapsAll$gapPropCens <- bci.gapsAll$gapProp
    bci.gapsAll[!is.na(bci.gapsAll$areaObs) & bci.gapsAll$areaObs<(0.5*cellSize*cellSize/10000),"gapPropCens"] <- NA
    
    # How many observations are left?
    nrow(bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),]) # 16276 observations
    
    # Scale topographic covariates
    for(i in 9:58){
      bci.gapsAll$new <- NA
      bci.gapsAll[!is.na(bci.gapsAll$age),"new"] <- scale(bci.gapsAll[!is.na(bci.gapsAll$age),i])
      names(bci.gapsAll)[names(bci.gapsAll)=="new"] <- paste0("Sc_",names(bci.gapsAll)[i])
    }
    
    # Add censoring value for beta distribution
    cens <- 0.0001
    bci.gapsAll$gapPropCens <- bci.gapsAll$gapProp
    bci.gapsAll$gapPropCens[bci.gapsAll$gapPropCens <= cens] <- 0
    bci.gapsAll$gapPropCens[bci.gapsAll$gapPropCens >= 1-cens] <- 1
    
    save(bci.gaps18, bci.gaps20, bci.gapsAll, cellSize, nCellX, nCellY, cens, file="INLA/INLA_prelim_40m_tin.RData")
    
#### Smoothing scale sensitivity analysis ####
    library(INLA)
    load("INLA/INLA_prelim_40m_tin.RData")
 
    
    # Make a data frame to store results from sensitivity analysis
    # across all smoothing scales for 
    
    smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
    
    topoScaleResults <- data.frame(curvScale = rep(smoothScales,length(smoothScales)),
                                   slopeScale = rep(smoothScales, each = length(smoothScales)),
                                   margLik = NA,
                                   DIC = NA)
    
    # Remove NA rows
    bci.gapsAll_Sub <- bci.gapsAll[!is.na(bci.gapsAll$gapProp),]
    
    # Use mean values for topographic variables--with quadratic terms
    
    for(i in 1:nrow(topoScaleResults)){
      
      # For initial model comparison, use only fixed effects
      
      fixed_i <- paste0("Sc_curvMean_",topoScaleResults$curvScale[i]," + I(Sc_curvMean_",topoScaleResults$curvScale[i],"^2) + Sc_slopeMean_",topoScaleResults$slopeScale[i]," + I(Sc_slopeMean_",topoScaleResults$slopeScale[i],"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
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
      meanQuad_DIC <- matrix(data = resultsMeanQuad$DIC,
                           nrow = length(smoothScales),
                           ncol = length(smoothScales),
                           byrow = F,
                           dimnames = list(smoothScales,smoothScales))
      meanQuad_DIC <- meanQuad_DIC-min(meanQuad_DIC)
      
      plotMin <- min(meanQuad_DIC)
      plotMax <- max(meanQuad_DIC)
      plotBreaks <- c(seq(0,30,1),plotMax)
    
      library(plot.matrix)
      par(mar=c(4,4,2,4))
      plot(meanQuad_DIC, breaks = plotBreaks,
           ylab = expression("Curvature scale ("~sigma~")"),
           xlab = expression("Slope scale ("~sigma~")"),
           main = expression(Delta~"DIC score"))
    
      resultsMeanQuad[which(resultsMeanQuad$DIC==min(resultsMeanQuad$DIC)),]

#### Using best scale result, do variable selection with Matern spatial autocorrelation term ####
curvScale <- 9
slopeScale <- 15

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

# STEP 1
  random_mat <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"

  fixed_ref <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref," + ",random_mat))

  fixed_a <- paste0("I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a," + ",random_mat))
  
  fixed_b <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b," + ",random_mat))
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c," + ",random_mat))
  
  fixed_d <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d," + ",random_mat)) 
  
  fixed_e <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e," + ",random_mat))
  
  fixed_f <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + soilParent + soilForm + age + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f," + ",random_mat))
  
  fixed_g <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilForm + age + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g," + ",random_mat))
  
  fixed_h <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + age + Year")
  form_h <- formula(paste0("gapPropCens ~ ",fixed_h," + ",random_mat))
  
  fixed_i <- paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + Year")
  form_i <- formula(paste0("gapPropCens ~ ",fixed_i," + ",random_mat))
  
  form_list1 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g, form_h, form_i)
  
  step1_results <- data.frame(model = c("ref","a","b","c","d","e","f","g","h","i"),
                               margLik = NA,
                               DIC = NA)
  
  for(i in 3:nrow(step1_results)){

    tryCatch(model_i <- inla(form_list1[[i]],
                    family = "beta",
                    data = bci.gapsAll_Order,
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
  
# STEP 2
  
  fixed_ref <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref," + ",random_mat))
  
  fixed_a <- paste0("Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a," + ",random_mat))
  
  fixed_b <- paste0("Sc_curvMean_",curvScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b," + ",random_mat))
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c," + ",random_mat))
  
  fixed_d <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d," + ",random_mat))
  
  fixed_e <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + soilParent + soilForm + age + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e," + ",random_mat))
  
  fixed_f <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilForm + age + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f," + ",random_mat))
  
  fixed_g <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + age + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g," + ",random_mat))
  
  fixed_h <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + Year")
  form_h <- formula(paste0("gapPropCens ~ ",fixed_h," + ",random_mat))
  
  form_list2 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g, form_h)
  
  step2_results <- data.frame(model = c("ref","a","b","c","d","e","f","g","h"),
                              margLik = NA,
                              DIC = NA)
  
  for(i in 1:nrow(step2_results)){
    
    tryCatch(model_i <- inla(form_list2[[i]],
                             family = "beta",
                             data = bci.gapsAll_Order,
                             control.compute = list(dic = TRUE),
                             control.family = list(beta.censor.value = cens),
                             verbose = F),
             error = function(e) {print(paste("iteration",i,"failed"))})
    
    #summary(model_i)
    step2_results$margLik[i] <- model_i$mlik[2]
    step2_results$DIC[i] <- model_i$dic$dic
    print(i)
  }
  
  step2_results$dDIC <- step2_results$DIC - min(step2_results$DIC,na.rm=T)
  step2_results <- step2_results[order(step2_results$dDIC),]
  
  # STEP 3
  
  fixed_ref <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_ref <- formula(paste0("gapPropCens ~ ",fixed_ref," + ",random_mat))
  
  fixed_a <- paste0("Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_a <- formula(paste0("gapPropCens ~ ",fixed_a," + ",random_mat))
  
  fixed_b <- paste0("Sc_curvMean_",curvScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_b <- formula(paste0("gapPropCens ~ ",fixed_b," + ",random_mat))
  
  fixed_c <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_c <- formula(paste0("gapPropCens ~ ",fixed_c," + ",random_mat))
  
  fixed_d <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + soilParent + soilForm + age + Year")
  form_d <- formula(paste0("gapPropCens ~ ",fixed_d," + ",random_mat))
  
  fixed_e <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilForm + age + Year")
  form_e <- formula(paste0("gapPropCens ~ ",fixed_e," + ",random_mat))
  
  fixed_f <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + age + Year")
  form_f <- formula(paste0("gapPropCens ~ ",fixed_f," + ",random_mat))
  
  fixed_g <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + I(Sc_slopeMean_",slopeScale,"^2) + I(Sc_drainMean^2) + soilParent + soilForm + Year")
  form_g <- formula(paste0("gapPropCens ~ ",fixed_g," + ",random_mat))
  
  form_list3 <- list(form_ref, form_a, form_b, form_c, form_d, form_e, form_f, form_g)
  
  step3_results <- data.frame(model = c("ref","a","b","c","d","e","f","g"),
                              margLik = NA,
                              DIC = NA)
  
  for(i in 1:nrow(step3_results)){
    
    tryCatch(model_i <- inla(form_list3[[i]],
                             family = "beta",
                             data = bci.gapsAll_Order,
                             control.compute = list(dic = TRUE),
                             control.family = list(beta.censor.value = cens),
                             verbose = F),
             error = function(e) {print(paste("iteration",i,"failed"))})
    
    #summary(model_i)
    step3_results$margLik[i] <- model_i$mlik[2]
    step3_results$DIC[i] <- model_i$dic$dic
    print(i)
  }
  
  step3_results$dDIC <- step3_results$DIC - min(step3_results$DIC,na.rm=T)
  step3_results <- step3_results[order(step3_results$dDIC),]
  
  
  # OLD
  
  for(i in 1:nrow(step2_results)){
    
    model_i <- inla(form_list2[[i]],
                    family = "beta",
                    data = bci.gapsAll_Sub,
                    control.compute = list(dic = TRUE),
                    control.family = list(beta.censor.value = cens)) 
    
    #summary(model_i)
    step2_results$margLik[i] <- model_i$mlik[2]
    step2_results$DIC[i] <- model_i$dic$dic
    print(i)
    
  }
  
  step2_results$dDIC <- step2_results$DIC - min(step2_results$DIC)
  step2_results <- step2_results[order(step2_results$dDIC),]
  
#### Run full model with Matern spatial autocorrelation term -- alternate autocorrelation ####
  
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
  
  # Run the best (and next most simple) models with full spatial autocorrelation
  fixed_full <-  paste0("Sc_curvMean_",curvScale," + I(Sc_curvMean_",curvScale,"^2) + Sc_slopeMean_",slopeScale," + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_full <- formula(paste0("gapPropCens ~ ",fixed_full," + ",random_full))
  
  fixed_full_simp <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  random_full <- "f(ID, model = \"matern2d\", nrow = nCellY*2, ncol = nCellX)"
  form_full_simp <- formula(paste0("gapPropCens ~ ",fixed_full_simp," + ",random_full))
  
  model_full <- inla(form_full,
                     family = "beta",
                     data = bci.gapsAll_Order,
                     control.compute = list(dic = TRUE),
                     control.family = list(beta.censor.value = cens))
  
  model_full_simp <- inla(form_full_simp,
                     family = "beta",
                     data = bci.gapsAll_Order,
                     control.compute = list(dic = TRUE),
                     control.family = list(beta.censor.value = cens))
  
  save(model_full, model_full_simp, file = "INLA/INLA_fullModelResult.RData")
  
  
  fixed_full_simp2 <- paste0("Sc_curvMean_",curvScale," + Sc_drainMean + I(Sc_drainMean^2) + soilParent + soilForm + age + Year")
  form_full_simp2 <- formula(paste0("gapPropCens ~ ",fixed_full_simp2," + ",random_full))
  
  model_full_simp2 <- inla(form_full_simp2,
                          family = "beta",
                          data = bci.gapsAll_Order,
                          control.compute = list(dic = TRUE),
                          control.family = list(beta.censor.value = cens))
  model_full_simp$dic$dic
  model_full_simp2$dic$dic
  
#### Run full model with RW2d spatial autocorrelation term -- alternate autocorrelation ####
  
  # Make an ID value for each cell
  bci.gapsAll$Order <- 1:nrow(bci.gapsAll)
  
  # Reorder so that INLA thinks there is one spatial pattern
  newOrder <- c()
  for(i in 1:nCellX){
    # Year 1
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY)):(i*nCellY))
    #Year 2
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + nCellX*nCellY):(i*nCellY + nCellX*nCellY))
    #Year 3
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + 2*nCellX*nCellY):(i*nCellY + 2*nCellX*nCellY))
    #Year 4
    newOrder <- c(newOrder,(1 + (i-1)*(nCellY) + 3*nCellX*nCellY):(i*nCellY + 3*nCellX*nCellY))
  }
  
  # Reorder and make a new ID column
  bci.gapsAll_Order <- bci.gapsAll[newOrder,]
  bci.gapsAll_Order$ID <- 1:nrow(bci.gapsAll_Order)
  
  fixed_full <- paste0("Sc_curvMean_",curvScale," + Sc_slopeMean_",slopeScale," + soil + age + Sc_drainMean + I(Sc_drainMean^2) + Year")
  random_full3 <- "f(ID, model = \"rw2d\", nrow = nCellY*4, ncol = nCellX)"
  form_full3 <- formula(paste0("gapPropCens ~ ",fixed_full," + ",random_full3))
  
  model_full3 <- inla(form_full3,
                      family = "beta",
                      data = bci.gapsAll_Order,
                      control.compute = list(dic = TRUE),
                      control.family = list(beta.censor.value = cens))
  
  save(model_full3, file = "INLA/INLA_fullModelResult3.RData")
  
#### Look at model residuals ####
  load("INLA/INLA_fullModelResult2.RData")
  
  # Get fitted values and reorder to match original data frame
  bci.gapsAll$pred <- model_full2$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  bci.gapsAll$resd <- bci.gapsAll$pred - bci.gapsAll$gapPropCens
  
  pred <- model_full2$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
  
  pred <- pred[!is.na(bci.gapsAll$gapPropCens)]
  obsv <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"gapPropCens"]
  resd <- obsv - pred
  
  
  plot(x = pred, y = obsv,
       col = adjustcolor("black", 0.1),
       ylab = "Observed",
       xlab = "Predicted",
       pch=19)
  abline(a=0,b=1,col="red")
  summary(lm(obsv~pred))
  
  
  par(mfrow=c(2,1), las=1, mar=c(4,4,0,2), oma=c(1,1,1,1))
  
  curvSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_21"], y = resd, df=50)
  plot(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"curvMean_21"], y = resd,
       pch=19,
       col = adjustcolor("black", 0.05),
       ylab = "Residual value",
       xlab = "Curvature value")
  abline(h=0,col="red")
  lines(x = curvSpline$x, y= curvSpline$y, col="orange",lty=2,lwd=2)
  
  curvRMSE <- data.frame(min = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[1:10],
                         max = seq(min(curvSpline$x),max(curvSpline$x),length.out=11)[2:11],
                         mid = NA,
                         RMSE = NA)
  curvRMSE$mid <- 0.5*(curvRMSE$min+curvRMSE$max)
  for(i in 1:nrow(curvRMSE)){
    res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                         & bci.gapsAll$curvMean_21 > curvRMSE$min[i]
                         & bci.gapsAll$curvMean_21 <= curvRMSE$max[i],"resd"]
    curvRMSE$RMSE[i] <- sqrt(mean(res_i^2))
  }
  
  plot(RMSE~mid, data=curvRMSE,
       ylab = "RMSE",
       xlab = "Curvature value",
       pch=19, ylim=c(0,0.05))
  
  
  
  slopeSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_9"], y = resd, df=50)
  plot(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"slopeMean_9"], y = resd,
       pch=19,
       col = adjustcolor("black", 0.05),
       ylab = "Residual value",
       xlab = "Slope (degree)")
  abline(h=0,col="red")
  lines(x = slopeSpline$x, y= slopeSpline$y, col="orange",lty=2,lwd=2)
  slopeRMSE <- data.frame(min = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[1:10],
                         max = seq(min(slopeSpline$x),max(slopeSpline$x),length.out=11)[2:11],
                         mid = NA,
                         RMSE = NA)
  slopeRMSE$mid <- 0.5*(slopeRMSE$min+slopeRMSE$max)
  for(i in 1:nrow(slopeRMSE)){
    res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                         & bci.gapsAll$slopeMean_9 > slopeRMSE$min[i]
                         & bci.gapsAll$slopeMean_9 <= slopeRMSE$max[i],"resd"]
    slopeRMSE$RMSE[i] <- sqrt(mean(res_i^2))
  }
  
  plot(RMSE~mid, data=slopeRMSE,
       ylab = "RMSE",
       xlab = "Slope (degree)",
       pch=19, ylim=c(0,0.05))
  
  
  
  aspectSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = resd, df=50)
  plot(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"aspectMean"], y = resd,
       pch=19,
       col = adjustcolor("black", 0.05),
       ylab = "Residual value",
       xlab = "Aspect (degree)")
  abline(h=0,col="red")
  lines(x = aspectSpline$x, y= aspectSpline$y, col="orange",lty=2,lwd=2)
  aspectRMSE <- data.frame(min = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[1:10],
                          max = seq(min(aspectSpline$x),max(aspectSpline$x),length.out=11)[2:11],
                          mid = NA,
                          RMSE = NA)
  aspectRMSE$mid <- 0.5*(aspectRMSE$min+aspectRMSE$max)
  for(i in 1:nrow(aspectRMSE)){
    res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                         & bci.gapsAll$aspectMean > aspectRMSE$min[i]
                         & bci.gapsAll$aspectMean <= aspectRMSE$max[i],"resd"]
    aspectRMSE$RMSE[i] <- sqrt(mean(res_i^2))
  }
  plot(RMSE~mid, data=aspectRMSE,
       ylab = "RMSE",
       xlab = "Aspect (degree)",
       pch=19, ylim=c(0,0.05))
  
  
  drainSpline <- smooth.spline(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = resd, df=50)
  plot(x = bci.gapsAll[!is.na(bci.gapsAll$gapPropCens),"drainMean"], y = resd,
       pch=19,
       col = adjustcolor("black", 0.05),
       ylab = "Residual value",
       xlab = " Height above drainage (m)")
  abline(h=0,col="red")
  lines(x = drainSpline$x, y= drainSpline$y, col="orange",lty=2,lwd=2)
  drainRMSE <- data.frame(min = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[1:10],
                           max = seq(min(drainSpline$x),max(drainSpline$x),length.out=11)[2:11],
                           mid = NA,
                           RMSE = NA)
  drainRMSE$mid <- 0.5*(drainRMSE$min+drainRMSE$max)
  for(i in 1:nrow(drainRMSE)){
    res_i <- bci.gapsAll[!is.na(bci.gapsAll$gapPropCens) 
                         & bci.gapsAll$drainMean > drainRMSE$min[i]
                         & bci.gapsAll$drainMean <= drainRMSE$max[i],"resd"]
    drainRMSE$RMSE[i] <- sqrt(mean(res_i^2))
  }
  plot(RMSE~mid, data=drainRMSE,
       ylab = "RMSE",
       xlab = "drain (degree)",
       pch=19, ylim=c(0,0.05))

#### Make plots of observed, average predicted, and SD predicted values for each year ####
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  
  
  allMeans <-  model_full2$summary.fitted.values$mean[order(bci.gapsAll_Order$Order)]
#  allMeans[is.na(bci.gapsAll$gapPropCens)] <- NA
  
  allSds <-  model_full2$summary.fitted.values$sd[order(bci.gapsAll_Order$Order)]
#  allSds[is.na(bci.gapsAll$gapPropCens)] <- NA
  
  bci.gaps17$predictedMean <- allMeans[1:(nCellX*nCellY)]
  bci.gaps18$predictedMean <- allMeans[(1 + nCellX*nCellY):(2*nCellX*nCellY)]
  bci.gaps19$predictedMean <- allMeans[(1 + 2*nCellX*nCellY):(3*nCellX*nCellY)]
  bci.gaps20$predictedMean <- allMeans[(1 + 3*nCellX*nCellY):(4*nCellX*nCellY)]
  
  bci.gaps17$predictedSD <- allSds[1:(nCellX*nCellY)]
  bci.gaps18$predictedSD <- allSds[(1 + nCellX*nCellY):(2*nCellX*nCellY)]
  bci.gaps19$predictedSD <- allSds[(1 + 2*nCellX*nCellY):(3*nCellX*nCellY)]
  bci.gaps20$predictedSD <- allSds[(1 + 3*nCellX*nCellY):(4*nCellX*nCellY)]

  # convert to a raster object
  
  # 2017
  predMeanRaster17 <- raster::raster(x = matrix(data = bci.gaps17$predictedMean,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps17)@xmin,
                                     xmx = raster::extent(bci.gaps17)@xmax,
                                     ymn = raster::extent(bci.gaps17)@ymin,
                                     ymx = raster::extent(bci.gaps17)@ymax)
  predMeanRaster17 <- raster::mask(predMeanRaster17, buffer)
  predMeanRaster17[predMeanRaster17<0] <- NA
  
  
  obsRaster17 <- raster::raster(x = matrix(data = bci.gaps17$gapProp,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps17)@xmin,
                                     xmx = raster::extent(bci.gaps17)@xmax,
                                     ymn = raster::extent(bci.gaps17)@ymin,
                                     ymx = raster::extent(bci.gaps17)@ymax)
  obsRaster17 <- raster::mask(obsRaster17, buffer)

  
  
  predSdRaster17 <- raster::raster(x = matrix(data = bci.gaps17$predictedSD,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps17)@xmin,
                                     xmx = raster::extent(bci.gaps17)@xmax,
                                     ymn = raster::extent(bci.gaps17)@ymin,
                                     ymx = raster::extent(bci.gaps17)@ymax)
  predSdRaster17[is.na(predMeanRaster17)] <- NA
  

  # 2018
  predMeanRaster18 <- raster::raster(x = matrix(data = bci.gaps18$predictedMean,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps18)@xmin,
                                     xmx = raster::extent(bci.gaps18)@xmax,
                                     ymn = raster::extent(bci.gaps18)@ymin,
                                     ymx = raster::extent(bci.gaps18)@ymax)
  predMeanRaster18 <- raster::mask(predMeanRaster18, buffer)
  predMeanRaster18[predMeanRaster18<0] <- NA
  
  
  obsRaster18 <- raster::raster(x = matrix(data = bci.gaps18$gapProp,
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps18)@xmin,
                                xmx = raster::extent(bci.gaps18)@xmax,
                                ymn = raster::extent(bci.gaps18)@ymin,
                                ymx = raster::extent(bci.gaps18)@ymax)
  obsRaster18 <- raster::mask(obsRaster18, buffer)
  
  
  
  predSdRaster18 <- raster::raster(x = matrix(data = bci.gaps18$predictedSD,
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps18)@xmin,
                                   xmx = raster::extent(bci.gaps18)@xmax,
                                   ymn = raster::extent(bci.gaps18)@ymin,
                                   ymx = raster::extent(bci.gaps18)@ymax)
  predSdRaster18[is.na(predMeanRaster18)] <- NA
  
  # 2019
  predMeanRaster19 <- raster::raster(x = matrix(data = bci.gaps19$predictedMean,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps19)@xmin,
                                     xmx = raster::extent(bci.gaps19)@xmax,
                                     ymn = raster::extent(bci.gaps19)@ymin,
                                     ymx = raster::extent(bci.gaps19)@ymax)
  predMeanRaster19 <- raster::mask(predMeanRaster19, buffer)
  predMeanRaster19[predMeanRaster19<0] <- NA
  
  
  obsRaster19 <- raster::raster(x = matrix(data = bci.gaps19$gapProp,
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps19)@xmin,
                                xmx = raster::extent(bci.gaps19)@xmax,
                                ymn = raster::extent(bci.gaps19)@ymin,
                                ymx = raster::extent(bci.gaps19)@ymax)
  obsRaster19 <- raster::mask(obsRaster19, buffer)
  
  
  
  predSdRaster19 <- raster::raster(x = matrix(data = bci.gaps19$predictedSD,
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps19)@xmin,
                                   xmx = raster::extent(bci.gaps19)@xmax,
                                   ymn = raster::extent(bci.gaps19)@ymin,
                                   ymx = raster::extent(bci.gaps19)@ymax)
  predSdRaster19[is.na(predMeanRaster19)] <- NA
  
  # 2020
  predMeanRaster20 <- raster::raster(x = matrix(data = bci.gaps20$predictedMean,
                                                nrow = nCellY,
                                                ncol = nCellX,
                                                byrow = F),
                                     xmn = raster::extent(bci.gaps20)@xmin,
                                     xmx = raster::extent(bci.gaps20)@xmax,
                                     ymn = raster::extent(bci.gaps20)@ymin,
                                     ymx = raster::extent(bci.gaps20)@ymax)
  predMeanRaster20 <- raster::mask(predMeanRaster20, buffer)
  predMeanRaster20[predMeanRaster20<0] <- NA
  
  
  obsRaster20 <- raster::raster(x = matrix(data = bci.gaps20$gapProp,
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps20)@xmin,
                                xmx = raster::extent(bci.gaps20)@xmax,
                                ymn = raster::extent(bci.gaps20)@ymin,
                                ymx = raster::extent(bci.gaps20)@ymax)
  obsRaster20 <- raster::mask(obsRaster20, buffer)
  
  
  
  predSdRaster20 <- raster::raster(x = matrix(data = bci.gaps20$predictedSD,
                                              nrow = nCellY,
                                              ncol = nCellX,
                                              byrow = F),
                                   xmn = raster::extent(bci.gaps20)@xmin,
                                   xmx = raster::extent(bci.gaps20)@xmax,
                                   ymn = raster::extent(bci.gaps20)@ymin,
                                   ymx = raster::extent(bci.gaps20)@ymax)
  predSdRaster20[is.na(predMeanRaster20)] <- NA
  
  # Plot rasters
  colBreaks <- seq(0,0.85,0.05)

  
  # 2015-2017
  jpeg(height=3000,width=2000,file = "Figure S10. Raster patterns.jpg")
  par(mfrow=c(4,3))
  raster::plot(obsRaster17,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Observed")
  raster::plot(buffer, add=T)
  
  raster::plot(predMeanRaster17,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Predicted (mean)")
  raster::plot(buffer, add=T)
  
  raster::plot(predSdRaster17,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(seq(0,0.1,0.01))),
               breaks = seq(0,0.1,0.01),
               main = "Predicted (SD)")
  raster::plot(buffer, add=T)
  mtext("2015 - 2017", outer=F, side = 3, line = -1.5)
  
  # 2017-2018
  raster::plot(obsRaster18,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Observed")
  raster::plot(buffer, add=T)
  
  raster::plot(predMeanRaster18,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Predicted (mean)")
  raster::plot(buffer, add=T)
  
  raster::plot(predSdRaster18,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(seq(0,0.1,0.01))),
               breaks = seq(0,0.1,0.01),
               main = "Predicted (SD)")
  raster::plot(buffer, add=T)
  mtext("2017 - 2018", outer=F, side = 3, line = -1.5)
  
  # 2018-2019
  raster::plot(obsRaster19,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Observed")
  raster::plot(buffer, add=T)
  
  raster::plot(predMeanRaster19,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Predicted (mean)")
  raster::plot(buffer, add=T)
  
  raster::plot(predSdRaster19,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(seq(0,0.1,0.01))),
               breaks = seq(0,0.1,0.01),
               main = "Predicted (SD)")
  raster::plot(buffer, add=T)
  mtext("2018 - 2019", outer=F, side = 3, line = -1.5)
  
  # 2019-2020
  raster::plot(obsRaster20,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Observed")
  raster::plot(buffer, add=T)
  
  raster::plot(predMeanRaster20,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(colBreaks)),
               breaks = colBreaks,
               main = "Predicted (mean)")
  raster::plot(buffer, add=T)
  
  raster::plot(predSdRaster20,
               bty = "n", box = F, xaxt="n", yaxt="n",
               col = viridis::viridis(length(seq(0,0.1,0.01))),
               breaks = seq(0,0.1,0.01),
               main = "Predicted (SD)")
  raster::plot(buffer, add=T)
  mtext("2019 - 2020", outer=F, side = 3, line = -1.5)
dev.off()

#### Calculate estimated value with only fixed effects ####

# Calculate total fixed effects for each pixel
  
  # Intercept
  bci.gapsAll$fix_int <- model_full2$summary.fixed$mean[model_full2$names.fixed=="(Intercept)"]
  
  # Curvature
  bci.gapsAll$fix_C <- bci.gapsAll$Sc_curvMean_21*(model_full2$summary.fixed$mean[model_full2$names.fixed=="Sc_curvMean_21"])
  
  # Slope
  bci.gapsAll$fix_S <- bci.gapsAll$Sc_slopeMean_9*(model_full2$summary.fixed$mean[model_full2$names.fixed=="Sc_slopeMean_9"])
  
  # Height above drainage
  bci.gapsAll$fix_H <- bci.gapsAll$Sc_drainMean*(model_full2$summary.fixed$mean[model_full2$names.fixed=="Sc_drainMean"])
  bci.gapsAll$fix_H2 <- (bci.gapsAll$Sc_drainMean^2)*(model_full2$summary.fixed$mean[model_full2$names.fixed=="I(Sc_drainMean^2)"])
  
  # Soil
  bci.gapsAll$fix_soil <- 0
  
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Barbour","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilBarbour"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Fairchild","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilFairchild"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Gross","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilGross"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Harvard","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilHarvard"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Hood","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilHood"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Lake","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilLake"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Lutz","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilLutz"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Marron","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilMarron"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Poacher","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilPoacher"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Standley","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilStandley"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Wetmore","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilWetmore"]
  bci.gapsAll[!is.na(bci.gapsAll$soil) & bci.gapsAll$soil=="Zetek","fix_soil"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="soilZetek"]
  
  # Age
  bci.gapsAll$fix_age <- 0
  
  bci.gapsAll[!is.na(bci.gapsAll$age) & bci.gapsAll$age=="Secondary","fix_age"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="ageSecondary"]
  
  # Year 
  bci.gapsAll$fix_year <- 0
  bci.gapsAll[bci.gapsAll$Year=="2018","fix_age"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="Year2018"]
  bci.gapsAll[bci.gapsAll$Year=="2019","fix_age"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="Year2019"]
  bci.gapsAll[bci.gapsAll$Year=="2020","fix_age"] <- model_full2$summary.fixed$mean[model_full2$names.fixed=="Year2020"]
  
  # All
  bci.gapsAll$fix_sum <- bci.gapsAll$fix_int + bci.gapsAll$fix_C + bci.gapsAll$fix_S + bci.gapsAll$fix_H + bci.gapsAll$fix_H2 + bci.gapsAll$fix_soil + bci.gapsAll$fix_age + bci.gapsAll$fix_year
  
  # Predicted value with just fixed effects
  bci.gapsAll$fix_pred <- exp(bci.gapsAll$fix_sum)/(1 + exp(bci.gapsAll$fix_sum))

# Plot 
  fixRaster17 <- raster::raster(x = matrix(data = bci.gapsAll[bci.gapsAll$Year=="2017","fix_pred"],
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps17)@xmin,
                                xmx = raster::extent(bci.gaps17)@xmax,
                                ymn = raster::extent(bci.gaps17)@ymin,
                                ymx = raster::extent(bci.gaps17)@ymax)
  fixRaster17 <- raster::mask(fixRaster17, buffer)

  fixRaster18 <- raster::raster(x = matrix(data = bci.gapsAll[bci.gapsAll$Year=="2018","fix_pred"],
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps18)@xmin,
                                xmx = raster::extent(bci.gaps18)@xmax,
                                ymn = raster::extent(bci.gaps18)@ymin,
                                ymx = raster::extent(bci.gaps18)@ymax)
  fixRaster18 <- raster::mask(fixRaster18, buffer)
  
  fixRaster19 <- raster::raster(x = matrix(data = bci.gapsAll[bci.gapsAll$Year=="2019","fix_pred"],
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps19)@xmin,
                                xmx = raster::extent(bci.gaps19)@xmax,
                                ymn = raster::extent(bci.gaps19)@ymin,
                                ymx = raster::extent(bci.gaps19)@ymax)
  fixRaster19 <- raster::mask(fixRaster19, buffer)

  fixRaster20 <- raster::raster(x = matrix(data = bci.gapsAll[bci.gapsAll$Year=="2020","fix_pred"],
                                           nrow = nCellY,
                                           ncol = nCellX,
                                           byrow = F),
                                xmn = raster::extent(bci.gaps20)@xmin,
                                xmx = raster::extent(bci.gaps20)@xmax,
                                ymn = raster::extent(bci.gaps20)@ymin,
                                ymx = raster::extent(bci.gaps20)@ymax)
  fixRaster20 <- raster::mask(fixRaster20, buffer)
  
  raster::plot(fixRaster20,
               bty="n", box=F, xaxt="n", yaxt="n",
               col = viridis::viridis(50))
  raster::plot(buffer,add=T)

#### Figure 5: Make plots of average spatial pattern across all years ####

  avgStack <- raster::stack(fixRaster17,fixRaster18,fixRaster19,fixRaster20)
  
  avgPredictedRaster <- raster::calc(avgStack, sum, na.rm=T)

  par(mfrow=c(1,1))
  avgBreaks <- c(seq(0,0.2,0.005),seq(0.25,0.5,0.05))
  raster::plot(avgPredictedRaster,
               col=gapCols(length(avgBreaks)),
               breaks = avgBreaks,
               main = "Average spatial pattern -- mean predicted")  
  
  avgBreaks <- c(seq(0,0.12,0.005))
  raster::plot(sdPredictedRaster,
               col=gapCols(length(avgBreaks)),
               breaks = avgBreaks,
               main = "Average spatial pattern -- SD predicted") 
  
#### Low canopy area and average spatial pettern ####
  
  blockData <- read.csv("bootstrapBlocks.csv")
  lo09 <- raster::raster("binaryLoCanopy.tif")

  blockData$propGap <- NA
  blockData$propLo <- NA
  for(i in 1:nrow(blockData)){
    gapSub <- raster::crop(fixRaster20, raster::extent(blockData$xmin[i],
                                                       blockData$xmax[i],
                                                       blockData$ymin[i],
                                                       blockData$ymax[i]))
    blockData$propGap[i] <- mean(gapSub@data@values, na.rm=T)
    
    loSub <- raster::crop(lo09, raster::extent(blockData$xmin[i],
                                                       blockData$xmax[i],
                                                       blockData$ymin[i],
                                                       blockData$ymax[i]))
    blockData$propLo[i] <- length(loSub@data@values[!is.na(loSub@data@values) & loSub@data@values==0])/length(loSub@data@values[!is.na(loSub@data@values) & loSub@data@values==1])
    print(i)
  }
  
  plot(propLo~propGap, data = blockData,
       pch = 19, 
       xlab = "Average predicted disturbance frequency",
       ylab = "Proportion of area with canopy < 10 m")
  summary(lm(propLo~propGap, data = blockData))
  
  
#### Figure 6: Fixed effects sizes ####  
  fixedResults <- model_full2$summary.fixed
  # Remove "Other" age class (with only NA pixels)
  fixedResults <- fixedResults[-16,]
  
  par(mfrow=c(1,1))
  plot(x =fixedResults[2:nrow(fixedResults),"mean"],
       y = nrow(fixedResults):2,
       xlim = range(fixedResults[2:nrow(fixedResults),c(3,5)]) + c(0,0.5),
       pch = 19, 
       cex = 1,
       xlab = "Marginal fixed effect",
       ylab = NA, yaxt = "n")
  arrows(x0 = fixedResults[2:nrow(fixedResults),"0.025quant"],
         x1 = fixedResults[2:nrow(fixedResults),"0.975quant"],
         y0 = nrow(fixedResults):2,
         y1 = nrow(fixedResults):2,
         angle = 90, code=3,
         length = 0.05)
  abline(v=0, lty=2)
  text(row.names(fixedResults)[2:nrow(fixedResults)],
        x = fixedResults[2:nrow(fixedResults),"0.975quant"]+0.3,
       y = nrow(fixedResults):2)
  
  
  
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
#### Fig S2. Spatial scales of rasters vs analysis ####
  smoothScales <- c(0,3,9,15,21,27,33,39,45,51,57,63)
  
  # Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
  
  dem1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",3,".tif")) 
  dem1 <- raster::crop(dem1, raster::extent(c(626000,627000,1012500,1013500)))
  dem2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",21,".tif")) 
  dem2 <- raster::crop(dem2, raster::extent(c(626000,627000,1012500,1013500)))
  dem3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",39,".tif")) 
  dem3 <- raster::crop(dem3, raster::extent(c(626000,627000,1012500,1013500)))
  dem4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",57,".tif")) 
  dem4 <- raster::crop(dem4, raster::extent(c(626000,627000,1012500,1013500)))
  
  curv1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",3,".tif")) 
  curv1 <- raster::crop(curv1, raster::extent(c(626000,627000,1012500,1013500)))
  curv2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",21,".tif")) 
  curv2 <- raster::crop(curv2, raster::extent(c(626000,627000,1012500,1013500)))
  curv3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",39,".tif")) 
  curv3 <- raster::crop(curv3, raster::extent(c(626000,627000,1012500,1013500)))
  curv4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",57,".tif")) 
  curv4 <- raster::crop(curv4, raster::extent(c(626000,627000,1012500,1013500)))
  
  slope1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",3,".tif")) 
  slope1 <- raster::crop(slope1, raster::extent(c(626000,627000,1012500,1013500)))
  slope2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",21,".tif")) 
  slope2 <- raster::crop(slope2, raster::extent(c(626000,627000,1012500,1013500)))
  slope3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",39,".tif")) 
  slope3 <- raster::crop(slope3, raster::extent(c(626000,627000,1012500,1013500)))
  slope4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",57,".tif")) 
  slope4 <- raster::crop(slope4, raster::extent(c(626000,627000,1012500,1013500)))
  
  jpeg(filename = "Figure S2. Smoothing scales example.jpg", width=1000,height = 1500)
  par(mfrow = c(4,3), mar=c(2,2,2,6), oma=c(1,1,1,1))
  
    raster::plot(dem1,
                 bty="n", box=F,
                 xaxt = "n", yaxt="n")
    mtext(expression(sigma~"= 3"),side=2,outer=F, cex=2)
    raster::plot(curv1,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::viridis(256))
    raster::plot(slope1,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = heat.colors(256))
    
    raster::plot(dem2,xaxt = "n", yaxt="n",
                 bty="n", box=F)
    mtext(expression(sigma~"= 21"),side=2,outer=F, cex=2)
    raster::plot(curv2,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::viridis(256))
    raster::plot(slope2,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = heat.colors(256))

    raster::plot(dem3,xaxt = "n", yaxt="n",
                 bty="n", box=F)
    mtext(expression(sigma~"= 39"),side=2,outer=F, cex=2)
    raster::plot(curv3,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::viridis(256))
    raster::plot(slope3,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = heat.colors(256))
    
    raster::plot(dem4,xaxt = "n", yaxt="n",
                 bty="n", box=F)
    mtext(expression(sigma~"= 57"),side=2,outer=F, cex=2)
    mtext(expression("Elevation (m)"),side=1,outer=F, cex=2)
    raster::plot(curv4,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::viridis(256))
    mtext(expression("Curvature (LaPlacian convexity)"),side=1,outer=F, cex=2)
    raster::plot(slope4,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = heat.colors(256))
    mtext(expression("Slope (degrees)"),side=1,outer=F, cex=2)
  
  dev.off()
  
  
  
  curvMed <- raster::aggregate(curv, fact = 40, fun = median)
  curvMean <- raster::aggregate(curv, fact = 40, fun = mean)
  
  
  raster::plot(curv)
  
  raster::plot(curvMed)
  
  raster::plot(curvMean)
  
  
  
  slope <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",smoothScales[i],".tif")) 
  
  slope <- raster::mask(slope, buffer)
  
  raster::plot(slope, ext = raster::extent(c(626000,627000,1012500,1013500)))
  raster::plot(bci.gaps, add=T)  