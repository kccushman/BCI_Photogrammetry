# Define new canopy disturbances

#### READ RAW DATA ####

# Read in canopy height rasters:
  dsm15 <- raster::raster("DSM_2015_corrected_tin.tif")
  dsm18 <- raster::raster("DSM_2018_corrected_tin2.tif")
  dsm20 <- raster::raster("DSM_2020_corrected_tin2.tif")
  
# Read in forest age (used to exclude recent clearings)
  age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
  ageUse <- age[!(age$TYPE=="Clearings"),]

# Read polygon buffer 25 m inland from lake
  buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
  buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    
# Read in BCI DEM
  dem <- raster::raster("D:/BCI_Spatial/BCI_Topo/LidarDEM_BCI.tif")
  dem <- raster::crop(dem, raster::extent(ageUse))
  dem <- raster::resample(dem,dsm15)
  
# Read in QAQC data
  qaqc <- read.csv("gridInfo_QAQC.csv")
 
#### PROCESS AND SAVE CANOPY HEIGHT DATA ####
  
# Remove raster areas outside BCI perimeter (exclude within 25 m of lake)
  dsm15 <- raster::mask(dsm15, buffer)  
  dsm18 <- raster::mask(dsm18, buffer)  
  dsm20 <- raster::mask(dsm20, buffer)

# Remove raster areas in clearings
  dsm15 <- raster::mask(dsm15, ageUse)  
  dsm18 <- raster::mask(dsm18, ageUse)  
  dsm20 <- raster::mask(dsm20, ageUse)

# Crop to ensure each raster has same extent
  dsm15 <- raster::crop(dsm15, raster::extent(ageUse))  
  dsm18 <- raster::crop(dsm18, raster::extent(ageUse))  
  dsm20 <- raster::crop(dsm20, raster::extent(ageUse))    
  
# Subtract ground elevation
  chm15 <- dsm15-dem
  chm18 <- dsm18-dem
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
    if(qaqc$Use15[i]==F){
      chm15 <- raster::mask(chm15, extent_i, inverse = T)
    }  
    
    if(qaqc$Use18[i]==F){
      chm18 <- raster::mask(chm18, extent_i, inverse = T)
    }
    
    if(qaqc$Use20[i]==F){
      chm20 <- raster::mask(chm20, extent_i, inverse = T)
    }

  }

# Make sure all years have the same extent
  chm18 <- raster::crop(chm18, raster::extent(chm15))
  chm20 <- raster::crop(chm20, raster::extent(chm15))
  
  
# Cloud masks for 2017-2020
  mask15 <- raster::raster("CloudMask_2015.tif")
  mask18 <- raster::raster("CloudMask_2018.tif")
  mask20 <- raster::raster("CloudMask_2020.tif")
  
  # Resample to extent and resolution of CHMs
  mask15 <- raster::resample(mask15, chm15)
  mask18 <- raster::resample(mask18, chm18)
  mask20 <- raster::resample(mask20, chm20)
  
  # Remove cloud pixels
  chm15[!(mask15>0.99)] <- NA
  chm18[!(mask18==1)] <- NA
  chm20[!(mask20==1)] <- NA
    
# # Save
    raster::writeRaster(chm15, "CHM_2015_QAQC_tin.tif", overwrite=T)
    raster::writeRaster(chm18, "CHM_2018_QAQC_tin.tif", overwrite=T)
    raster::writeRaster(chm20, "CHM_2020_QAQC_tin.tif", overwrite=T)
  

#### MAKE AND SAVE CANOPY HEIGHT CHANGE RASTERS ####  
  
  # Canopy height models for all years
    # chm15 <- raster::raster("CHM_2015_QAQC_tin.tif")
    # chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
    # chm19 <- raster::raster("CHM_2019_QAQC_tin.tif")
    # chm20 <- raster::raster("CHM_2020_QAQC_tin.tif")
  
  # Calculate the change in canopy height for each interval  
    d15to18 <- chm18-chm15
    d18to20 <- chm20-chm18
    
    
  # Make additional rasters where values below and initial height threshold are
  # omitted 
    d15to18tall <- d15to18
    d18to20tall <- d18to20
    
  # Mask out areas that are initially < 5 m in height and decrease between two years
    # 2015 - 2018
      short15 <- rep(0, length(raster::values(chm15)))
      short15[raster::values(chm15)<5 & !is.na(raster::values(chm15))] <- 1
      
      d15to18tall@data@values[short15==1] <- NA
        
    # 2018 - 2020
      short18 <- rep(0, length(raster::values(chm18)))
      short18[raster::values(chm18)<5 & !is.na(raster::values(chm18))] <- 1
        d18to20tall@data@values[short18==1] <- NA    
      
  # Save rasters of canopy height change, with and without 5 m initial height mask
     
     raster::writeRaster(d15to18, file="dCHM15to18_tin.tif", overwrite = T)
     raster::writeRaster(d18to20, file="dCHM18to20_tin.tif", overwrite = T)
     
     raster::writeRaster(d15to18tall, file="dCHM15to18tall_tin.tif", overwrite = T)
     raster::writeRaster(d18to20tall, file="dCHM18to20tall_tin.tif", overwrite = T)  
    
        
#### BINARY NEW GAPS #### 
  
  # Use ForestGapR package to delineate new gaps
    
    # Define gap height threshold, min gap size, max gap size, and min area:perimeter ratio
      gapHtThresh <- -5
      gapSzMin <- 25
      gapSzMax <- 10^6
      areaPerimThresh <- 0.6
      
      # 2015 - 2018
      
      # Identify gaps  
      gaps15to18 <- ForestGapR::getForestGaps(d15to18,
                                              threshold = gapHtThresh ,
                                              size=c(gapSzMin,gapSzMax))
      
      # Create a Spatial Polygon Data Frame object, where each polygon is a gap
      gaps15to18sp <- ForestGapR::GapSPDF(gaps15to18)
      
      # Calculate the area and perimeter from each gap object
      gaps15to18sp@data$area <- NA
      gaps15to18sp@data$perimeter <- NA
      for(i in 1:length(gaps15to18sp)){
        gaps15to18sp[gaps15to18sp$gap_id==i,"area"] <- raster::area(gaps15to18sp[gaps15to18sp$gap_id==i,])
        
        perims_j <- c()
        for(j in 1:length(gaps15to18sp[gaps15to18sp$gap_id==i,]@polygons[[1]]@Polygons)){
          
          coordsj <- gaps15to18sp[gaps15to18sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
          
          lengths_k <- c()
          for(k in 2:dim(coordsj)[1]){
            lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
            
          }
          perims_j[j] <- sum(lengths_k)
          
          if(gaps15to18sp[gaps15to18sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
            perims_j[j] <- 0
          }
          
        }
        gaps15to18sp[gaps15to18sp$gap_id==i,"perimeter"] <- sum(perims_j)
        
      }
      
      # Calculate the ratio of area to perimeter
      gaps15to18sp@data$ratio <- gaps15to18sp@data$area/gaps15to18sp@data$perimeter
      
      # Remove gaps that don't meet a minimum area/perimeter threshold
      for(i in 1:length(gaps15to18sp@data$ratio)){
        if(gaps15to18sp@data$ratio[i] < areaPerimThresh){
          gaps15to18[gaps15to18==gaps15to18sp@data$gap_id[i]] <- NA
        }
      }
      
  
    # 2018 - 2020
        
        # Identify gaps  
        gaps18to20 <- ForestGapR::getForestGaps(d18to20,
                                                threshold = gapHtThresh ,
                                                size=c(gapSzMin,gapSzMax))
        
        # Create a Spatial Polygon Data Frame object, where each polygon is a gap
        gaps18to20sp <- ForestGapR::GapSPDF(gaps18to20)
        
        # Calculate the area and perimeter from each gap object
        gaps18to20sp@data$area <- NA
        gaps18to20sp@data$perimeter <- NA
        for(i in 1:length(gaps18to20sp)){
          gaps18to20sp[gaps18to20sp$gap_id==i,"area"] <- raster::area(gaps18to20sp[gaps18to20sp$gap_id==i,])
          
          perims_j <- c()
          for(j in 1:length(gaps18to20sp[gaps18to20sp$gap_id==i,]@polygons[[1]]@Polygons)){
            
            coordsj <- gaps18to20sp[gaps18to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@coords
            
            lengths_k <- c()
            for(k in 2:dim(coordsj)[1]){
              lengths_k[k-1] <- sqrt((coordsj[k,1]-coordsj[k-1,1])^2 + (coordsj[k,2]-coordsj[k-1,2])^2)
              
            }
            perims_j[j] <- sum(lengths_k)
            
            if(gaps18to20sp[gaps18to20sp$gap_id==i,]@polygons[[1]]@Polygons[[j]]@hole==T){
              perims_j[j] <- 0
            }
            
          }
          gaps18to20sp[gaps18to20sp$gap_id==i,"perimeter"] <- sum(perims_j)
          
        }
        
        # Calculate the ratio of area to perimeter
        gaps18to20sp@data$ratio <- gaps18to20sp@data$area/gaps18to20sp@data$perimeter
        
        # Remove gaps that don't meet a minimum area/perimeter threshold
        for(i in 1:length(gaps18to20sp@data$ratio)){
          if(gaps18to20sp@data$ratio[i] < areaPerimThresh){
            gaps18to20[gaps18to20==gaps18to20sp@data$gap_id[i]] <- NA
          }
        }      
  # Create variable to store whether gaps pass area:perimeter threshold 
      gaps15to18sp$use <- ifelse(gaps15to18sp$ratio<areaPerimThresh,
                                 F,T)
 
      gaps18to20sp$use <- ifelse(gaps18to20sp$ratio<areaPerimThresh,
                                 F,T)  
        
  # Save rasters of new gap pixels
       raster::writeRaster(gaps15to18, file="newGaps15to18_tin.tif", overwrite = T)
       raster::writeRaster(gaps18to20, file="newGaps18to20_tin.tif", overwrite = T)
      
#### GAPS NEIGHBORING NA CELLS ####
        
  # 2015 to 2018 
    
    # gaps15to18sp <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")    
    # d15to18 <- raster::raster("dCHM15to18.tif")    
    
    gaps15to18sp$borderNAs <- NA
    
    for(i in 1:length(gaps15to18sp)){
      # create a 1 m buffer (one cell) around gap polygon
      polyBuff <- raster::buffer(gaps15to18sp[i,], 1)
      # find raster cells within that polygon
      cells <- raster::cellFromPolygon(object = d15to18,
                                       p = polyBuff)[[1]]
      # count NA cells within buffered gap polygon
      gaps15to18sp$borderNAs[i] <- length(d15to18[cells][is.na(d15to18[cells])])
      
    }
        

  # 2018 to 2020 
    
    # gaps18to20sp <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")    
    # d18to20 <- raster::raster("dCHM18to20.tif")    
    
    gaps18to20sp$borderNAs <- NA
    
    for(i in 1:length(gaps18to20sp)){
      # create a 1 m buffer (one cell) around gap polygon
      polyBuff <- raster::buffer(gaps18to20sp[i,], 1)
      # find raster cells within that polygon
      cells <- raster::cellFromPolygon(object = d18to20,
                                       p = polyBuff)[[1]]
      # count NA cells within buffered gap polygon
      gaps18to20sp$borderNAs[i] <- length(d18to20[cells][is.na(d18to20[cells])])
      
    } 
    
  # Calculate NA cells as a proportion of total gap perimeter
    gaps15to18sp$pctNAs <- gaps15to18sp$borderNAs/gaps15to18sp$perimeter
    gaps18to20sp$pctNAs <- gaps18to20sp$borderNAs/gaps18to20sp$perimeter
    
  # What proportion of gaps neighbor NA cells?
    pct15to18 <- 100*round(length(gaps15to18sp[gaps15to18sp$use==T & gaps15to18sp$borderNAs>0,])/length(gaps15to18sp[gaps15to18sp$use==T,]),3)
    pct18to20 <- 100*round(length(gaps18to20sp[gaps18to20sp$use==T & gaps18to20sp$borderNAs>0,])/length(gaps18to20sp[gaps18to20sp$use==T,]),3)
    
  # What proportion of gaps have 10% or more of perimeter
    length(gaps15to18sp[gaps15to18sp$use==T & gaps15to18sp$pctNAs>0.1,])/length(gaps15to18sp[gaps15to18sp$use==T,])
    length(gaps18to20sp[gaps18to20sp$use==T & gaps18to20sp$pctNAs>0.1,])/length(gaps18to20sp[gaps18to20sp$use==T,])
    
#### ASSIGN GAP POLYGONS TO BOOTSTRAP GROUPS ####
  blocks <- read.csv("bootstrapBlocks.csv")
    
  gaps15to18sp$block <- NA 
  gaps18to20sp$block <- NA
  
  for(i in 1:dim(blocks)[1]){
    gaps15to18sp@data[gaps15to18sp@data[,1]>=blocks$xmin[i] & gaps15to18sp@data[,1] < blocks$xmax[i]
                      & gaps15to18sp@data[,2]>=blocks$ymin[i] & gaps15to18sp@data[,2] < blocks$ymax[i],"block"] <- i
   
     gaps18to20sp@data[gaps18to20sp$X1>=blocks$xmin[i] & gaps18to20sp$X1 < blocks$xmax[i]
                      & gaps18to20sp$X2>=blocks$ymin[i] & gaps18to20sp$X2 < blocks$ymax[i],"block"] <- i
  }  
    
    
#### SAVE GAP POLYGONS ####
    
    # Save gap shapefiles
    rgdal::writeOGR(gaps15to18sp,
                    dsn = "gaps15to18_shapefile_tin",
                    layer = "gaps15to18sp", 
                    driver = "ESRI Shapefile")
    
    rgdal::writeOGR(gaps18to20sp,
                    dsn = "gaps18to20_shapefile_tin",
                    layer = "gaps18to20sp", 
                    driver = "ESRI Shapefile")
    
#### Figure 1: Example CHMs ####
    d18to19 <- raster::raster("dCHM18to19_tin.tif")   
    chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
    chm19 <- raster::raster("CHM_2019_QAQC_tin.tif")
    gaps18to19sp <- rgdal::readOGR("gaps18to19_shapefile_tin/gaps18to19sp.shp")
    
    sampleExt <- raster::extent(c(625356,625418,1011182,1011235))
    example18 <- raster::crop(chm18,sampleExt)
    example19 <- raster::crop(chm19,sampleExt)
    example18to19 <- raster::crop(d18to19, sampleExt)
    
    raster::plot(example18,
                 col=rev(terrain.colors(20)),
                 breaks=seq(0,40,2),
                 bty = "n", box = F,xaxt = "n")
    raster::plot(example19,
                 col=rev(terrain.colors(20)),
                 breaks=seq(0,40,2),
                 bty = "n", box = F,xaxt = "n")
    raster::plot(example18to19,
                 col=viridis::viridis(12),
                 breaks=seq(-35,20,5),
                 bty = "n", box = F,xaxt = "n")
    raster::plot(gaps18to19sp, add=T, border="black", lwd=3)
    
    
#### Figure 2: area sampled and gaps in each interval ####

# Load rasters    
d15to18 <- raster::raster("dCHM15to18_tin.tif")     
d18to20 <- raster::raster("dCHM18to20_tin.tif")  

gaps15to18 <- raster::raster("newGaps15to18_tin.tif")
gaps18to20 <- raster::raster("newGaps18to20_tin.tif")


buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  

    
  makePlotRaster <- function(dRaster, buffer){
    plotRaster <- dRaster
    plotRaster[is.na(plotRaster)] <- -1000
    plotRaster <- raster::mask(plotRaster, buffer)
    return(plotRaster)
  }  
   
  plot15to18 <- makePlotRaster(d15to18, buffer)  
  plot18to20 <- makePlotRaster(d18to20, buffer)  
  
  

    raster::plot(plot15to18,
                 bty="n", box=F,yaxt="n",
                 breaks = c(-1001,-999,1000),
                 col = c("lightgrey","white"),
                 main = "2015 - 2018")
    raster::plot(buffer,add=T)
    raster::plot(gaps15to18, col = "red",add=T)
    
 
    raster::plot(plot18to20,
                 bty="n", box=F,yaxt="n",
                 breaks = c(-1001,-999,1000),
                 col = c("lightgrey","white"),
                 main = "2018 - 2020")
    raster::plot(buffer,add=T)
    raster::plot(gaps18to20, col = "red",add=T)