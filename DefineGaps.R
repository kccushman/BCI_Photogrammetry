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
 
#### PROCESS AND SAVE CANOPY HEIGHT DATA ####
  
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
  

#### MAKE AND SAVE CANOPY HEIGHT CHANGE RASTERS ####  
  
  # Canopy height models for all years
    chm17 <- raster::raster("CHM_2017_QAQC.tif")
    chm18 <- raster::raster("CHM_2018_QAQC.tif")
    chm19 <- raster::raster("CHM_2019_QAQC.tif")
    chm20 <- raster::raster("CHM_2020_QAQC.tif")
  
  # Calculate the change in canopy height for each interval  
    d17to18 <- chm18-chm17
    d18to19 <- chm19-chm18
    d19to20 <- chm20-chm19
    
  # Make additional rasters where values below and initial height threshold are
  # omitted 
    d17to18tall <- d17to18
    d18to19tall <- d18to19
    d19to20tall <- d19to20
    
  # Mask out areas that are initially < 5 m in height and decrease between two years
    
    # 2017 - 2018
      short17 <- rep(0, length(raster::values(chm17)))
        short17[raster::values(chm17)<5 & !is.na(raster::values(chm17))] <- 1
      
        d17to18tall@data@values[short17==1] <- NA

    # 2018 - 2019
        short18 <- rep(0, length(raster::values(chm18)))
        short18[raster::values(chm18)<5 & !is.na(raster::values(chm18))] <- 1
        
        d18to19tall@data@values[short18==1] <- NA
    
    # 2019 - 2020
        short19 <- rep(0, length(raster::values(chm19)))
        short19[raster::values(chm19)<5 & !is.na(raster::values(chm19))] <- 1
        
        d19to20tall@data@values[short19==1] <- NA
      
  # Save rasters of canopy height change, with and without 5 m initial height mask
      # raster::writeRaster(d17to18, file="dCHM17to18.tif", overwrite = T)
      # raster::writeRaster(d18to19, file="dCHM18to19.tif", overwrite = T)
      # raster::writeRaster(d19to20, file="dCHM19to20.tif", overwrite = T)
      # 
      # raster::writeRaster(d17to18tall, file="dCHM17to18tall.tif", overwrite = T)
      # raster::writeRaster(d18to19tall, file="dCHM18to19tall.tif", overwrite = T)
      # raster::writeRaster(d19to20tall, file="dCHM19to20tall.tif", overwrite = T)  
        
        
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
        
  # # Save rasters of new gap pixels
  #     raster::writeRaster(gaps17to18, file="newGaps17to18.tif", overwrite = T)
  #     raster::writeRaster(gaps18to19, file="newGaps18to19.tif", overwrite = T)
  #     raster::writeRaster(gaps19to20, file="newGaps19to20.tif", overwrite = T)
      
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
  
  # Calculate NA cells as a proportion of total gap perimeter
    gaps17to18sp$pctNAs <- gaps17to18sp$borderNAs/gaps17to18sp$perimeter
    gaps18to19sp$pctNAs <- gaps18to19sp$borderNAs/gaps18to19sp$perimeter
    gaps19to20sp$pctNAs <- gaps19to20sp$borderNAs/gaps19to20sp$perimeter
      
  # What proportion of gaps neighbor NA cells?
    pct17to18 <- 100*round(length(gaps17to18sp[gaps17to18sp$use==T & gaps17to18sp$borderNAs>0,])/length(gaps17to18sp[gaps17to18sp$use==T,]),3)
    pct18to19 <- 100*round(length(gaps18to19sp[gaps18to19sp$use==T & gaps18to19sp$borderNAs>0,])/length(gaps18to19sp[gaps18to19sp$use==T,]),3)    
    pct19to20 <- 100*round(length(gaps19to20sp[gaps19to20sp$use==T & gaps19to20sp$borderNAs>0,])/length(gaps19to20sp[gaps19to20sp$use==T,]),3)
  
  # What proportion of gaps have 10% or more of perimeter
    length(gaps17to18sp[gaps17to18sp$use==T & gaps17to18sp$pctNAs>0.1,])/length(gaps17to18sp[gaps17to18sp$use==T,])
    length(gaps18to19sp[gaps18to19sp$use==T & gaps18to19sp$pctNAs>0.1,])/length(gaps18to19sp[gaps18to19sp$use==T,])
    length(gaps19to20sp[gaps19to20sp$use==T & gaps19to20sp$pctNAs>0.1,])/length(gaps19to20sp[gaps19to20sp$use==T,])
      
    
    
  # Histograms of pct of perimeter in NA cells
    par(mfrow=c(1,3))
    hist(gaps17to18sp[gaps17to18sp$use==T & gaps17to18sp$pctNAs>0,]$pctNAs,
         main = paste0(pct17to18, "% in18to19"),
         xlab = NA,
         breaks=seq(0,1,.02),
         xlim=c(0.001,0.6),
         col="black",border="white")
    hist(gaps18to19sp[gaps18to19sp$use==T & gaps18to19sp$pctNAs>0,]$pctNAs,
         main = paste0(pct18to19, "% in 18to19"),
         xlab = NA,
         breaks=seq(0,1,.02),
         xlim=c(0.001,0.6),
         col="black",border="white")
    hist(gaps19to20sp[gaps19to20sp$use==T & gaps19to20sp$pctNAs>0,]$pctNAs,
         main = paste0(pct19to20, "% in 19to20"),
         xlab = NA,
         breaks=seq(0,2,.02),
         xlim=c(0.001,0.6),
         col="black",border="white")
    mtext("NA cells (proportion of perimeter)", side=1, outer=T, line=-1)
    

    # Relationship between gap size and NA cells
    par(mfrow=c(1,3))
    plot(x = gaps17to18sp[gaps17to18sp$use==T,]$area,
         y = gaps17to18sp[gaps17to18sp$use==T,]$borderNAs,
         xlim=c(10,1200),
         ylim=c(0,60),
         xlab = NA,
         ylab = NA,
         pch=19,
         col = adjustcolor("black",0.05))
    plot(x = gaps18to19sp[gaps18to19sp$use==T,]$area,
         y = gaps18to19sp[gaps18to19sp$use==T,]$borderNAs,
         xlim=c(10,1200),
         ylim=c(0,60),
         xlab = NA,
         ylab = NA,
         pch=19,
         col = adjustcolor("black",0.05))
    plot(x = gaps19to20sp[gaps19to20sp$use==T,]$area,
         y = gaps19to20sp[gaps19to20sp$use==T,]$borderNAs,
         xlim=c(10,1200),
         ylim=c(0,60),
         xlab = NA,
         ylab = NA,
         pch=19,
         col = adjustcolor("black",0.05))
    mtext("Gap size (m2)", side=1, outer=T, line=-2)
    mtext("Bordering NA pixels (#)", side=2, outer=T, line=-2)
    
    
#### ASSIGN GAP POLYGONS TO BOOTSTRAP GROUPS ####
  blocks <- read.csv("bootstrapBlocks.csv")
  
  gaps17to18sp$block <- NA 
  gaps18to19sp$block <- NA 
  gaps19to20sp$block <- NA
  for(i in 1:dim(blocks)[1]){
    gaps17to18sp@data[gaps17to18sp$X1>=blocks$xmin[i] & gaps17to18sp$X1 < blocks$xmax[i]
                      & gaps17to18sp$X2>=blocks$ymin[i] & gaps17to18sp$X2 < blocks$ymax[i],"block"] <- i
    gaps18to19sp@data[gaps18to19sp$X1>=blocks$xmin[i] & gaps18to19sp$X1 < blocks$xmax[i]
                      & gaps18to19sp$X2>=blocks$ymin[i] & gaps18to19sp$X2 < blocks$ymax[i],"block"] <- i
    gaps19to20sp@data[gaps19to20sp$X1>=blocks$xmin[i] & gaps19to20sp$X1 < blocks$xmax[i]
                      & gaps19to20sp$X2>=blocks$ymin[i] & gaps19to20sp$X2 < blocks$ymax[i],"block"] <- i
  }  
    
    
#### SAVE GAP POLYGONS ####
    
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
    
    # Load gap shapefiles
    gaps17to18sp <- rgdal::readOGR("gaps17to18_shapefile/gaps17to18sp.shp")
    gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile/gaps18to19sp.shp")
    gaps19to20sp <- rgdal::readOGR("gaps19to20_shapefile/gaps19to20sp.shp")
    