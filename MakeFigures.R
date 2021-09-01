#### Define colors ####
  # Colors used in plots
  col18 <- "blue"
  col20 <- "#d95f02"
  colBro <- wesanderson::wes_palette("Rushmore1",5)[1]
  colMot <- wesanderson::wes_palette("Rushmore1",5)[4]
  colPal <- wesanderson::wes_palette("Rushmore1",5)[3]
  colRed <- wesanderson::wes_palette("Rushmore1",5)[5]
  colAnd <- wesanderson::wes_palette("Chevalier1",4)[2]
  colBoh <- wesanderson::wes_palette("Chevalier1",4)[4]
  colMar <- wesanderson::wes_palette("Chevalier1",4)[3]
  colVol <- wesanderson::wes_palette("Chevalier1",4)[1]
  colOld <- wesanderson::wes_palette("Moonrise2",4)[1]
  colSec <- wesanderson::wes_palette("Moonrise2",4)[2]
#### Figure 1: Example CHMs ####
d18to20 <- raster::raster("dCHM18to20_tin.tif")   
chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
chm20 <- raster::raster("CHM_2020_QAQC_tin.tif")
gaps18to20sp <- rgdal::readOGR("gaps18to20_shapefile_tin/gaps18to20sp.shp")

sampleExt <- raster::extent(c(625840,625896,1010870,1010909))
example18 <- raster::crop(chm18,sampleExt)
example20 <- raster::crop(chm20,sampleExt)
example18to20 <- raster::crop(d18to20, sampleExt)
examplesGaps <- raster::crop(gaps18to20sp, sampleExt)

par(mfrow=c(2,1), mar=c(1,1,0,2))
raster::plot(example18,
             col=rev(terrain.colors(20)),
             breaks=seq(0,40,2),
             bty = "n", box = F,xaxt = "n",yaxt="n")
raster::plot(example20,
             col=rev(terrain.colors(20)),
             breaks=seq(0,40,2),
             bty = "n", box = F,xaxt = "n",yaxt="n")
raster::plot(example18to20,
             col=viridis::viridis(12),
             breaks=seq(-35,20,5),
             bty = "n", box = F,xaxt = "n")
raster::plot(examplesGaps, add=T, border="black", lwd=3)


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


par(mfrow=c(1,2), mar=c(2,1,1,1), oma=c(1,1,1,3))
raster::plot(plot15to18,
             bty="n", box=F,yaxt="n",
             breaks = c(-1001,-999,1000),
             col = c("lightgrey","white"),
             main = NA,
             legend=F)
raster::plot(buffer,add=T)
raster::plot(gaps15to18, col = "red",add=T, legend=F)


raster::plot(plot18to20,
             bty="n", box=F,yaxt="n",
             breaks = c(-1001,-999,1000),
             col = c("lightgrey","white"),
             main = NA,
             legend=F)
raster::plot(buffer,add=T)
raster::plot(gaps18to20, col = "red",add=T, legend=F)

#### Figure S1: plot corrected and uncorrected dCHM #### 
  # Read polygon buffer 25 m inland from lake
    buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
    buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs") 

  # Make uncorrected height change rasters
    raw15 <- raster::raster("DSM_2015_raw_tin.tif")
    raw18 <- raster::raster("DSM_2018_raw_tin.tif")
    raw20 <- raster::raster("DSM_2020_raw_tin.tif")
    
    raw15to18 <- raw18-raw15
      raw15to18 <- raster::crop(raw15to18, buffer)
      raw15to18 <- raster::mask(raw15to18, buffer)
    raw18to20 <- raw20-raw18
      raw18to20 <- raster::crop(raw18to20, buffer)
      raw18to20 <- raster::mask(raw18to20, buffer)
    
  # Read corrected height change rasters
    cor15to18 <- raster::raster("dCHM15to18_tin.tif")
    cor18to20 <- raster::raster("dCHM18to20_tin.tif")
    
  # Set NA values within buffer to -10000 for plotting  
    makePlotRaster <- function(dRaster, buffer){
      plotRaster <- dRaster
      plotRaster[is.na(plotRaster)] <- -300
      plotRaster <- raster::mask(plotRaster, buffer)
      return(plotRaster)
    }  
    raw15to18 <- makePlotRaster(raw15to18, buffer)
    raw18to20 <- makePlotRaster(raw18to20, buffer)
    cor15to18 <- makePlotRaster(cor15to18, buffer)
    cor18to20 <- makePlotRaster(cor18to20, buffer)
    
    meanRaw15to18 <- mean(raster::values(raw15to18),na.rm=T)
    meanRaw18to20 <- mean(raster::values(raw18to20),na.rm=T)
    meanCor15to18 <- mean(raster::values(cor15to18),na.rm=T)
    meanCor18to20 <- mean(raster::values(cor18to20),na.rm=T)
    
    
  # Make plot of corrected and uncorrected height change
    colBrks <- c(-50,-25,-10,-5,5,10,25,50)
    colPal <- colorRampPalette(c("red","darksalmon","khaki2",
                                  "white",
                                  "skyblue1","skyblue3","darkblue"))
    
      par(mfrow=c(2,2), mar=c(0,1,1,1), oma=c(1,1,1,3))
      raster::plot(raw15to18 - meanRaw15to18,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA,
                   legend.width=1.5)
      raster::plot(buffer,add=T)
      
      raster::plot(cor15to18 ,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA,
                   legend = F)
      raster::plot(buffer,add=T)
      
      raster::plot(raw18to20- meanRaw18to20,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA,
                   legend = F)
      raster::plot(buffer,add=T)
      
      raster::plot(cor18to20,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA,
                   legend = F)
      raster::plot(buffer,add=T)


#### Figure S2 (2015 data correction): in "Validation_50HaPlot.R" ####
#### Figure S5: Results from disturbance validation in 50 ha plot ####
    
## Panel a: plot monthly disturbances
    valGaps15to18 <- rgdal::readOGR("toCheckRaquel/toCheckRaquel.shp")
    valGapsSzHt <- valGaps15to18[valGaps15to18$area2 >= 25 & valGaps15to18$htDrop <= -5 & !is.na(valGaps15to18$htDrop),]
    
    # make a jiggered date field to plot
    valGapsSzHt$datClss <- as.Date(valGapsSzHt$datClss)
    valGapsSzHt$datePlot <-  valGapsSzHt$datClss + rnorm(mean = 0, sd = 7, n = nrow(valGapsSzHt))
    
    par(mar=c(3,3,1,2), oma = c(1,2,2,0), mfrow=c(1,2), las=1)
    
    # observed with ht and area agreement
    plot(-htDrop~area2, data=valGapsSzHt[valGapsSzHt$vslChck %in% c(0,5),],
         xlim=c(20,1800),
         ylim=c(5,45),
         pch=20,
         col=adjustcolor("grey",0.99),
         cex = 0.75,
         log="xy",
         ylab=NA,
         xlab = NA)
    mtext(expression("Gap area " (m^2)),side=1,outer=T, las=0)
    mtext(expression("Mean height decrease (m)"),side=2,outer=F, las=0, line=2.5)
    mtext(expression("Disturbances detected in monthly"),side=3,outer=F, las=0, line=1)
    mtext(expression("high-resolution data"),side=3,outer=F, las=0, line=0)
    
    # Red: monthly data are correct but missed by me
    points(-htDrop~area2, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==1,],
           pch=20,
           cex = 0.75,
           col=adjustcolor("red",0.5))
    
    # Orange: false positive in monthly data
    points(-htDrop~area2, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==2,],
           pch=20,
           cex = 0.75,
           col=adjustcolor("orange",0.6))
    
    # Purple: unclear which data set is correct
    points(-htDrop~area2, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==3,],
           pch=20,
           cex = 0.75,
           col=adjustcolor("purple",0.6))
    
    # BLUE: likely size discrepancy (potential overestimation in monthly data)
      points(-htDrop~area2, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==4,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("blue",0.6))

      text("a",
           x=20,
           y=45)
      
    legend(x=90,
           y=52,
           bty="n",
           c("Observed in interannual data",
             "Missing in interannual data",
             "Potential size overestimate",
             "False positive",
             "Unclear"),
           col=adjustcolor(c("grey","red","blue","orange","purple"),c(0.6)),
          pch=20,
          pt.cex=1,
          cex=0.9)
      
    # # mismatch in space (group these with OK trees because not a detection problem)
    #   # Spatial alignment mismatch
    #   points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==5,],
    #          pch=20,
    #          col=adjustcolor("purple",0.6))
  
## Panel b: plot multiannual disturbances
    
    plotGaps <- rgdal::readOGR("toCheckKC_final/toCheckKC_final.shp")
    
    
    # Observed with ht and area agreement
    plot(ratiCrc~area2, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==T,],
         ylim=c(0.07, 0.85),
         xlim=c(20,1800),
         pch=20,
         col=adjustcolor("grey",0.99),
         cex = 0.75,
         log="x",
         ylab=NA,
         xlab = NA)
    mtext(expression("Gap circularity " (4*pi*"Area"/"Perim"^2)),side=2,outer=F, las=0, line=2.5)
    mtext(expression("Disturbances detected in multiannual"),side=3,outer=F, las=0, line=1)
    mtext(expression("lower-resolution data"),side=3,outer=F, las=0, line=0)    
    # mismatch in space/time [consider this OK because not an error in detection]
      # Spatial alignment mismatch
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck==5,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("grey",0.99))
      # Temporal mismatch (i.e. I see beginning of slowly dying tree and Raquel records the later fall)
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(8),],
             pch=20,
             cex = 0.75,
             col=adjustcolor("grey",0.99))
    
    # BLUE: observed by Raquel's data when constraints are not imposed
      # Light blue: Observed by Raquel's data are less than 25m2
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==F & plotGaps$obsrvdH==T & !(plotGaps$vslChck %in% c(5,8,9,10)),],
             pch=20,
             cex = 0.75,
             col=adjustcolor("blue",0.6))
      # Observed by Raquel's data have ht drop less than 5 m
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==F & !(plotGaps$vslChck %in% c(5,8,9,10)),],
             pch=20,
             cex = 0.75,
             col=adjustcolor("blue",0.6))
      # Observed by Raquel's data are less than 25m2 AND have ht drop less than 5m (NONE)
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==F & plotGaps$obsrvdH==F & !(plotGaps$vslChck %in% c(5,8,9,10)),],
             pch=20,
             cex = 0.75,
             col=adjustcolor("blue",0.5))
    
    # RED: false positive in my data
      # Total false positive
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$vslChck==2,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("red",0.6))
      # Overestimation of gap area
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$vslChck==10,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("orange",0.6))
    
    
    # PURPLE: my data are correct
      # Unclear why Raquel misses this gap
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$vslChck==1,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("darkorchid1",0.6))
      # Decaying tree: total miss
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$vslChck==7,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("magenta1",0.6))
      # Decaying tree: partial miss
      points(ratiCrc~area2, data=plotGaps@data[plotGaps$vslChck==9,],
             pch=20,
             cex = 0.75,
             col=adjustcolor("palevioletred1",0.6))
      
      text("b",
           x=20,
           y=0.85)
      
      legend(x=40,
             y=0.90,
             bty="n",
             c("Observed",
               "Observed with size/height discrepancy",
               "Missing in monthly data: total",
               "Missing in monthly data: partial",
               "Missing in monthly data: unclear",
               "False positive: no disturbance",
               "False positive: small disturbance"),
             col=adjustcolor(c("grey",
                               "blue",
                               "magenta1","palevioletred1","darkorchid",
                               "red","orange"),c(0.6)),
             pch=20,
             pt.cex=1,
             cex=0.9)

  ## Stacked barplots
    # Monthly data
    monthlySums <- as.table(matrix(c(sum(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(0,5)])/sum(valGapsSzHt$area2),
                                                sum(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(1)])/sum(valGapsSzHt$area2),
                                                sum(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(4)])/sum(valGapsSzHt$area2),
                                                sum(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(2)])/sum(valGapsSzHt$area2),
                                                sum(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(3)])/sum(valGapsSzHt$area2),
                                      length(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(0,5)])/length(valGapsSzHt$area2),
                                                  length(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(1)])/length(valGapsSzHt$area2),
                                                  length(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(4)])/length(valGapsSzHt$area2),
                                                  length(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(2)])/length(valGapsSzHt$area2),
                                                  length(valGapsSzHt$area2[valGapsSzHt$vslChck %in% c(3)])/length(valGapsSzHt$area2)),nrow=2,
                                   byrow=T))
    barplot(t(monthlySums),
            col=adjustcolor(c("grey","red","blue","orange","purple"),c(0.6)),
            names.arg = c("Gap area", "Gap number"))
    mtext("c", adj = 0.1, outer=T)
    par(las=0)
    mtext("Proportion of gap area or #", side = 2, outer=T)
    par(las=1)
    
    # Interannual data
    annualSums <- as.table(matrix(c(sum(plotGaps@data$area2[(plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==T & !(plotGaps$vslChck==9)) | (plotGaps$vslChck %in% c(5,8))])/sum(plotGaps$area2),
                                    sum(plotGaps@data$area2[(plotGaps$obsrvdAl==T & (plotGaps$obsrvdAr==F | plotGaps$obsrvdH==F)) & !(plotGaps$vslChck %in% c(5,8,9,10))])/sum(plotGaps$area2),
                                    sum(plotGaps@data$area2[plotGaps$vslChck==7])/sum(plotGaps$area2),
                                    sum(plotGaps@data$area2[plotGaps$vslChck==9])/sum(plotGaps$area2),
                                    sum(plotGaps@data$area2[plotGaps$vslChck==1])/sum(plotGaps$area2),
                                    sum(plotGaps@data$area2[plotGaps$vslChck==10])/sum(plotGaps$area2),
                                    sum(plotGaps@data$area2[plotGaps$vslChck==2])/sum(plotGaps$area2),
                                    length(plotGaps@data$area2[(plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==T & !(plotGaps$vslChck==9)) | (plotGaps$vslChck %in% c(5,8))])/length(plotGaps$area2),
                                    length(plotGaps@data$area2[(plotGaps$obsrvdAl==T & (plotGaps$obsrvdAr==F | plotGaps$obsrvdH==F)) & !(plotGaps$vslChck %in% c(5,8,9,10))])/length(plotGaps$area2),
                                    length(plotGaps@data$area2[plotGaps$vslChck==1])/length(plotGaps$area2),
                                    length(plotGaps@data$area2[plotGaps$vslChck==7])/length(plotGaps$area2),
                                    length(plotGaps@data$area2[plotGaps$vslChck==9])/length(plotGaps$area2),
                                    length(plotGaps@data$area2[plotGaps$vslChck==10])/length(plotGaps$area2),
                                    length(plotGaps@data$area2[plotGaps$vslChck==2])/length(plotGaps$area2)),nrow=2,
                                   byrow=T))
    barplot(t(annualSums),
            col=adjustcolor(c("grey",
                              "blue",
                              "magenta1","palevioletred1","darkorchid",
                              "red","orange"),c(0.6)),
            names.arg = c("Gap area", "Gap number"),
            yaxt="n")
    mtext("d", adj = 0.6, outer=T)
    

#### Figure S4. Smoothing scale example plots ####

    # Read polygon buffer 25 m inland from lake
    buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
    buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    
    dem1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",2,".tif")) 
    dem1 <- raster::crop(dem1, raster::extent(c(626000,627000,1012500,1013500)))
    dem2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",8,".tif")) 
    dem2 <- raster::crop(dem2, raster::extent(c(626000,627000,1012500,1013500)))
    dem3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",16,".tif")) 
    dem3 <- raster::crop(dem3, raster::extent(c(626000,627000,1012500,1013500)))
    dem4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",32,".tif")) 
    dem4 <- raster::crop(dem4, raster::extent(c(626000,627000,1012500,1013500)))
    
    curv1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",2,".tif")) 
    curv1 <- raster::crop(curv1, raster::extent(c(626000,627000,1012500,1013500)))
    curv2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",8,".tif")) 
    curv2 <- raster::crop(curv2, raster::extent(c(626000,627000,1012500,1013500)))
    curv3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",16,".tif")) 
    curv3 <- raster::crop(curv3, raster::extent(c(626000,627000,1012500,1013500)))
    curv4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",32,".tif")) 
    curv4 <- raster::crop(curv4, raster::extent(c(626000,627000,1012500,1013500)))
    
    slope1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",2,".tif")) 
    slope1 <- raster::crop(slope1, raster::extent(c(626000,627000,1012500,1013500)))
    slope2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",8,".tif")) 
    slope2 <- raster::crop(slope2, raster::extent(c(626000,627000,1012500,1013500)))
    slope3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",16,".tif")) 
    slope3 <- raster::crop(slope3, raster::extent(c(626000,627000,1012500,1013500)))
    slope4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",32,".tif")) 
    slope4 <- raster::crop(slope4, raster::extent(c(626000,627000,1012500,1013500)))
    
    legendWidth = 2
    legendCex = 2
    jpeg(filename = "Figure S4. Smoothing scales example.jpg", width=1000,height = 1100)
    par(mfrow = c(4,3), mar=c(2,2,2,5), oma=c(1,1,1,1))
    
    raster::plot(dem1,
                 bty="n", box=F,
                 xaxt = "n", yaxt="n",
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    mtext(expression(sigma~"= 2"),side=2,outer=F, cex=2)
    raster::plot(curv1,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::cividis(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    raster::plot(slope1,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::plasma(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    
    raster::plot(dem2,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    mtext(expression(sigma~"= 8"),side=2,outer=F, cex=2)
    raster::plot(curv2,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::cividis(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    raster::plot(slope2,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::plasma(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    
    raster::plot(dem3,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    mtext(expression(sigma~"= 16"),side=2,outer=F, cex=2)
    raster::plot(curv3,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::cividis(128),axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    raster::plot(slope3,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::plasma(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    
    raster::plot(dem4,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    mtext(expression(sigma~"= 32"),side=2,outer=F, cex=2)
    mtext(expression("Elevation (m)"),side=1,outer=F, cex=2, line=2)
    raster::plot(curv4,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::cividis(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    mtext(expression("Curvature (LaPlacian convexity)"),side=1,outer=F, cex=2, line=2)
    raster::plot(slope4,xaxt = "n", yaxt="n",
                 bty="n", box=F,
                 col = viridis::plasma(128),
                 axis.args=list(cex.axis=legendCex),
                 legend.width=legendWidth)
    mtext(expression("Slope (degrees)"),side=1,outer=F, cex=2, line=2)
    
    dev.off()
    
    
#### Figure S5. Spatial averaging at scale of INLA analysis ####
    curvRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_2.tif")
    slopeRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_16.tif")
    drainRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/distAboveStream_1000.tif")
    
    curvCrop <- raster::crop(curvRaster, raster::extent(c(626000,627000,1012500,1013500)))
    curvMean <- raster::aggregate(curvCrop, fact = 40, fun = mean)
    
    slopeCrop <- raster::crop(slopeRaster, raster::extent(c(626000,627000,1012500,1013500)))
    slopeMean <- raster::aggregate(slopeCrop, fact = 40, fun = mean)
    
    drainCrop <- raster::crop(drainRaster, raster::extent(c(626000,627000,1012500,1013500)))
    drainMean <- raster::aggregate(drainCrop, fact = 40, fun = mean)
    
    legendWidth = 2
    legendCex = 2
    
    jpeg(filename = "Figure S5. INLA scale example.jpg", width=1000,height = 600)
    par(mfrow = c(2,3), mar=c(2,2,2,5), oma=c(1,1,1,1))
      
      # 1 m scale
      raster::plot(curvCrop,col = viridis::cividis(128),
                   bty="n", box=F, yaxt="n", xaxt="n",
                   axis.args=list(cex.axis=legendCex),
                   legend.width=legendWidth)
      
      raster::plot(slopeCrop,col = viridis::plasma(128),
                   bty="n", box=F, yaxt="n", xaxt="n",
                   axis.args=list(cex.axis=legendCex),
                   legend.width=legendWidth)
      
      raster::plot(drainCrop,col = viridis::viridis(128),
                   bty="n", box=F, yaxt="n", xaxt="n",
                   axis.args=list(cex.axis=legendCex),
                   legend.width=legendWidth)
      
      # resampled
      raster::plot(curvMean, col = viridis::cividis(128),
                   bty="n", box=F, yaxt="n", xaxt="n",
                   axis.args=list(cex.axis=legendCex),
                   legend.width=legendWidth)
      mtext(expression("Curvature (LaPlacian convexity)"),side=1,outer=F, cex=1.8, line=2)
      
      raster::plot(slopeMean, col = viridis::plasma(128),
                   bty="n", box=F, yaxt="n", xaxt="n",
                   axis.args=list(cex.axis=legendCex),
                   legend.width=legendWidth)
      mtext(expression("Slope (degrees)"),side=1,outer=F, cex=1.8, line=2)
      
      raster::plot(drainMean, col = viridis::viridis(128),
                   bty="n", box=F, yaxt="n", xaxt="n",
                   axis.args=list(cex.axis=legendCex),
                   legend.width=legendWidth)
      mtext(expression("Height above drainage (m)"),side=1,outer=F, cex=1.8, line=2)
      
    dev.off()
    
    
#### Figure S6. Best smoothing scale for INLA analysis ####
    # Read results
    resultsMeanQuad <- read.csv("INLA/INLA_scaleTopo_MeanVals_quad.csv")
    smoothScales <- c(1,2,3,4,6,8,12,16,24,32,48,64)
    
    # Convert results to matrices--curvature scale along columns, slope scale along hows
    meanQuad_DIC <- matrix(data = resultsMeanQuad$DIC,
                           nrow = length(smoothScales),
                           ncol = length(smoothScales),
                           byrow = F,
                           dimnames = list(smoothScales,smoothScales))
    meanQuad_DIC <- meanQuad_DIC-min(meanQuad_DIC)
    
    plotMin <- min(meanQuad_DIC)
    plotMax <- ceiling(max(meanQuad_DIC))
    plotBreaks <- seq(0,plotMax+2,2)
    
    library(plot.matrix)
    par(mar=c(4,5,2,8), mfrow=c(1,1))
    plot(meanQuad_DIC, breaks = plotBreaks,
         col = rev(wesanderson::wes_palette("Zissou1", length(plotBreaks), type = "continuous")),
         ylab = expression("Curvature scale ("~sigma~")"),
         xlab = expression("Slope scale ("~sigma~")"),
         main = expression(Delta~"DIC score"))
    
#### Figure S4: Location of corrected 2015 island-wide data ####
    chm15_raw <- raster::raster("CHM_2015_QAQC_tin_wBias.tif")
    chm15_cor <- raster::raster("CHM_2015_QAQC_tin.tif")
    
    buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
    buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    
    values_raw <- raster::values(chm15_raw)
    values_cor <- raster::values(chm15_cor)
    values_dif <- which(!(values_cor==values_raw) & !is.na(values_cor))
    
    # What percent of values are changed?
    100*length(values_dif)/ length(values_cor[!is.na(values_cor)])
    
    chm15_cor_plot <- chm15_cor
    raster::values(chm15_cor_plot)[values_dif] <- -9999
    
    makePlotRaster <- function(dRaster, buffer){
      plotRaster <- dRaster
      plotRaster[is.na(plotRaster)] <- -1000
      plotRaster <- raster::mask(plotRaster, buffer)
      return(plotRaster)
    }  
    
    chm15_cor_plot <- makePlotRaster(chm15_cor_plot, buffer)  
    
    raster::plot(chm15_cor_plot,
                 breaks = c(-10000,-1001,-100,999),
                 col = c("orange","grey","lightgreen"))
    
    
#### Figure S3: Example of height change correction around gaps ####
    plotShp <- rgdal::readOGR("D:/BCI_Spatial/BCI50ha/BCI_50ha.shp")
    plotShp <- sp::spTransform(plotShp, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
    
    chm09 <- raster::raster("CHM_2009_QAQC.tif")
    chm15_raw <- raster::raster("CHM_2015_QAQC_tin_wBias.tif")
    chm15_cor <- raster::raster("CHM_2015_QAQC_tin.tif")
    chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
    
    # Crop rasters to plot outline
    chm09 <- raster::crop(chm09, plotShp)
    chm15_raw <- raster::crop(chm15_raw, plotShp)
    chm15_cor <- raster::crop(chm15_cor, plotShp)
    chm18 <- raster::crop(chm18, plotShp)
    
    # Find corrected and uncorrected gaps
      htChange_raw <- chm18 - chm15_raw
      htChange_cor <- chm18 - chm15_cor
    
      # Define new function based on original ForestGapR function that does not
      # include diagonal pixels in the same gap
      # Make new function
      getForestGaps <- function (chm_layer, threshold = 10, size = c(1, 10^4)) 
      {
        chm_layer[chm_layer > threshold] <- NA
        chm_layer[chm_layer <= threshold] <- 1
        gaps <- raster::clump(chm_layer, directions = 4, gap = FALSE)
        rcl <- raster::freq(gaps)
        rcl[, 2] <- rcl[, 2] * raster::res(chm_layer)[1]^2
        rcl <- cbind(rcl[, 1], rcl)
        z <- raster::reclassify(gaps, rcl = rcl, right = NA)
        z[is.na(gaps)] <- NA
        gaps[z > size[2]] <- NA
        gaps[z < size[1]] <- NA
        gaps <- raster::clump(gaps, directions = 4, gap = FALSE)
        names(gaps) <- "gaps"
        return(gaps)
      }
    
    
      # Define gap height threshold, min gap size, max gap size, and min area:perimeter ratio
      gapHtThresh <- -5
      gapSzMin <- 25
      gapSzMax <- 10^6
    
      # Identify gaps  
      gaps_raw <- getForestGaps(htChange_raw,
                                  threshold = gapHtThresh ,
                                  size=c(gapSzMin,gapSzMax))
      gaps_cor <- getForestGaps(htChange_cor,
                                threshold = gapHtThresh ,
                                size=c(gapSzMin,gapSzMax))
      gaps_rawsp <- ForestGapR::GapSPDF(gaps_raw)
      gaps_corsp <- ForestGapR::GapSPDF(gaps_cor)
      
      rgdal::writeOGR(gaps_rawsp,
                      dsn = "50haPlotGaps_Original",
                      layer = "50haPlotGaps_Original", 
                      driver = "ESRI Shapefile")
      rgdal::writeOGR(gaps_corsp,
                      dsn = "50haPlotGaps_Corrected",
                      layer = "50haPlotGaps_Corrected", 
                      driver = "ESRI Shapefile")
      
      
      # PLOT EXAMPLE 1
      ext1 <- raster::extent(gaps_rawsp[gaps_rawsp$gap_id==185,])
      ext1@xmin <- ext1@xmin - 10
      ext1@xmax <- ext1@xmax + 7
      ext1@ymin <- ext1@ymin - 1
      ext1@ymax <- ext1@ymax + 7
      
      # Plot 2009, 2015, and 2018 canopy heights
      par(mfrow=c(1,4), mar=c(1,1,1,1), oma=c(1,1,1,6))
      
      raster::plot(chm09,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   breaks=seq(0,45,3),
                   col = rev(terrain.colors(length(seq(0,45,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      raster::plot(chm15_raw,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   breaks=seq(0,45,3),
                   col = rev(terrain.colors(length(seq(0,45,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      raster::plot(chm18,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   breaks=seq(0,45,3),
                   col = rev(terrain.colors(length(seq(0,45,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      raster::plot(chm15_cor,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=T,
                   legend.width = 3,
                   breaks=seq(0,45,3),
                   col = rev(terrain.colors(length(seq(0,45,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      par(mfrow=c(1,3), mar=c(1,1,1,1), oma=c(1,1,1,6))
      raster::plot(chm15_raw-chm09,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   legend.width = 3,
                   breaks=c(-15,seq(-5,30,5)),
                   col = viridis::viridis(length(c(-15,seq(-5,30,5)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      raster::plot(chm18-chm15_raw,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   legend.width = 3,
                   breaks=c(-15,seq(-5,30,5)),
                   col = viridis::viridis(length(c(-15,seq(-5,30,5)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      raster::plot(chm18-chm15_cor,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=T,
                   legend.width = 3,
                   breaks=c(-15,seq(-5,30,5)),
                   col = viridis::viridis(length(c(-15,seq(-5,30,5)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==185,],add=T, lwd=2)
      
      # PLOT EXAMPLE 2
      ext1 <- raster::extent(gaps_rawsp[gaps_rawsp$gap_id==199,])
      ext1@xmin <- ext1@xmin - 15
      ext1@xmax <- ext1@xmax + 2
      ext1@ymin <- ext1@ymin - 5
      ext1@ymax <- ext1@ymax + 2
      
      # Plot 2009, 2015, and 2018 canopy heights
      par(mfrow=c(1,4), mar=c(1,1,1,1), oma=c(1,1,1,6))
      
      raster::plot(chm09,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   breaks=seq(0,36,3),
                   col = rev(terrain.colors(length(seq(0,36,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      
      raster::plot(chm15_raw,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   breaks=seq(0,36,3),
                   col = rev(terrain.colors(length(seq(0,36,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      
      raster::plot(chm18,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   breaks=seq(0,36,3),
                   col = rev(terrain.colors(length(seq(0,36,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      
      raster::plot(chm15_cor,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=T,
                   legend.width = 3,
                   breaks=seq(0,36,3),
                   col = rev(terrain.colors(length(seq(0,36,3)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      

      
      par(mfrow=c(1,3), mar=c(1,1,1,1), oma=c(1,1,1,6))
      raster::plot(chm15_raw-chm09,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   legend.width = 3,
                   breaks=c(-20,seq(-5,20,5)),
                   col = viridis::viridis(length(c(-20,seq(-5,20,5)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      
      raster::plot(chm18-chm15_raw,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=F,
                   legend.width = 3,
                   breaks=c(-20,seq(-5,20,5)),
                   col = viridis::viridis(length(c(-20,seq(-5,20,5)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      
      raster::plot(chm18-chm15_cor,
                   ext = ext1,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   legend=T,
                   legend.width = 3,
                   breaks=c(-20,seq(-5,20,5)),
                   col = viridis::viridis(length(c(-20,seq(-5,20,5)))))
      raster::plot(gaps_rawsp[gaps_rawsp$gap_id==199,],add=T, lwd=2)
      raster::plot(gaps_corsp[gaps_corsp$gap_id==162,],add=T, lwd=2, border="white", lty=2)
      

    
    
#### Figure S#: Correlations among topographic variables ####

  # Get values for whole island: load, resample, and mask rasters to island
      
      # Use gap polygon so that 40 m aggregation is the same used in the INLA analysis
      gaps15to18sp <- rgdal::readOGR("gaps15to18_shapefile_tin/gaps15to18sp.shp")
      
      # BCI outline
      buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
      buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
      
      # Curvature raster
      curvRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_2.tif")
      curvRaster <- raster::crop(curvRaster, raster::extent(gaps15to18sp))
      curvMean <- raster::aggregate(curvRaster, 40)
      curvMean <- raster::mask(curvMean,buffer)
      
      # Slope raster
      slopeRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_16.tif")
      slopeRaster <- raster::crop(slopeRaster, raster::extent(gaps15to18sp))
      slopeMean <- raster::aggregate(slopeRaster, 40)
      slopeMean <- raster::mask(slopeMean,buffer)
      
      # Distance above drainage raster
      drainRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/distAboveStream_1000.tif")
      # resample to same extent as other rasters (adds NA area to edges)
      drainRaster <- raster::resample(drainRaster, curvRaster)
      drainRaster <- raster::crop(drainRaster, raster::extent(gaps15to18sp))
      drainMean <- raster::aggregate(drainRaster, 40)
      drainMean <- raster::mask(drainMean,buffer)
      
  # Get values from rasters
      curvVals <- raster::values(curvMean); curvVals <- curvVals[!is.na(curvVals)]
      slopeVals <- raster::values(slopeMean); slopeVals <- slopeVals[!is.na(slopeVals)]
      drainVals <- raster::values(drainMean); drainVals <- drainVals[!is.na(drainVals)]
      
  # Plot correlations
      
    par(mfrow=c(2,2), mar=c(2,2,1,1), oma=c(2,2,0,0))
    plot(x = curvVals, y = slopeVals,
         pch = 19,
         col = adjustcolor("black",0.05),
         xaxt="n")
    mtext("Slope (degrees)", outer=F, line = 2.5, side=2)
    text(bquote(rho~"="~.(round(cor(x = curvVals, y = slopeVals),2))),
         x = 3.5, y = 30)
    
    plot(x = drainVals, y = slopeVals,
         pch = 19,
         col = adjustcolor("black",0.05),
         yaxt="n")
    mtext("HAND (m)", outer=F, line = 2.5, side=1)
    text(bquote(rho~"="~.(round(cor(x = drainVals, y = slopeVals),2))),
         x = 33, y = 30)
    
    plot(x = curvVals, y = drainVals,
         pch = 19,
         col = adjustcolor("black",0.05))
    mtext("Curvature (LaPlacian convexity)", outer=F, line = 2.5, side=1)
    mtext("HAND (m)", outer=F, line = 2.5, side=2)
    text(bquote(rho~"="~.(round(cor(x = curvVals, y = drainVals),2))),
         x = 3.5, y = 33)
      
      
#### Figure S#: Proportion of topography in each forest age and soil type ####
    
    # RUN CODE ABOVE FOR CORRELATIONS AMONG TOPOGRAPHIC VARIABLES #
    
  # Define forest age and soil type polygons
    
    # Forest age polygon
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    age$AgeClass <- "Other"
    age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
    age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    ageUse <- age[!(age$AgeClass=="Other"),]
    
    # Soil type polygon  
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(age))
    
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
    
  # Make curvature density plots
    # All island
    cDens <- density(raster::values(curvMean),na.rm=T)
    # By forest age
    cDensOld <- density(raster::values(raster::mask(curvMean,ageUse[ageUse$AgeClass=="OldGrowth",])),
                        na.rm=T)
    cDensSec <- density(raster::values(raster::mask(curvMean,ageUse[ageUse$AgeClass=="Secondary",])),
                        na.rm=T)
    # By soil parent material
    cDensAnd <- density(raster::values(raster::mask(curvMean,soil[soil$SoilParent=="Andesite",])),
                        na.rm=T)
    cDensBoh <- density(raster::values(raster::mask(curvMean,soil[soil$SoilParent=="Bohio",])),
                        na.rm=T)
    cDensVol <- density(raster::values(raster::mask(curvMean,soil[soil$SoilParent=="CaimitoVolcanic",])),
                        na.rm=T)
    cDensMar <- density(raster::values(raster::mask(curvMean,soil[soil$SoilParent=="CaimitoMarineSedimentary",])),
                        na.rm=T)
    # By soil form
    cDensPal <- density(raster::values(raster::mask(curvMean,soil[soil$SoilForm=="PaleSwellingClay",])),
                        na.rm=T)
    cDensRed <- density(raster::values(raster::mask(curvMean,soil[soil$SoilForm=="RedLightClay",])),
                        na.rm=T)
    cDensMot <- density(raster::values(raster::mask(curvMean,soil[soil$SoilForm=="MottledHeavyClay",])),
                        na.rm=T)
    cDensBro <- density(raster::values(raster::mask(curvMean,soil[soil$SoilForm=="BrownFineLoam",])),
                        na.rm=T)
    
  # Make slope density plots
    # All island
    sDens <- density(raster::values(slopeMean),na.rm=T)
    # By forest age
    sDensOld <- density(raster::values(raster::mask(slopeMean,ageUse[ageUse$AgeClass=="OldGrowth",])),
                        na.rm=T)
    sDensSec <- density(raster::values(raster::mask(slopeMean,ageUse[ageUse$AgeClass=="Secondary",])),
                        na.rm=T)
    # By soil parent material
    sDensAnd <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilParent=="Andesite",])),
                        na.rm=T)
    sDensBoh <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilParent=="Bohio",])),
                        na.rm=T)
    sDensVol <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilParent=="CaimitoVolcanic",])),
                        na.rm=T)
    sDensMar <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilParent=="CaimitoMarineSedimentary",])),
                        na.rm=T)
    # By soil form
    sDensPal <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilForm=="PaleSwellingClay",])),
                        na.rm=T)
    sDensRed <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilForm=="RedLightClay",])),
                        na.rm=T)
    sDensMot <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilForm=="MottledHeavyClay",])),
                        na.rm=T)
    sDensBro <- density(raster::values(raster::mask(slopeMean,soil[soil$SoilForm=="BrownFineLoam",])),
                        na.rm=T)  
    
  # Make HAND density plots
    # All island
    dDens <- density(raster::values(drainMean),na.rm=T)
    # By forest age
    dDensOld <- density(raster::values(raster::mask(drainMean,ageUse[ageUse$AgeClass=="OldGrowth",])),
                        na.rm=T)
    dDensSec <- density(raster::values(raster::mask(drainMean,ageUse[ageUse$AgeClass=="Secondary",])),
                        na.rm=T)
    # By soil parent material
    dDensAnd <- density(raster::values(raster::mask(drainMean,soil[soil$SoilParent=="Andesite",])),
                        na.rm=T)
    dDensBoh <- density(raster::values(raster::mask(drainMean,soil[soil$SoilParent=="Bohio",])),
                        na.rm=T)
    dDensVol <- density(raster::values(raster::mask(drainMean,soil[soil$SoilParent=="CaimitoVolcanic",])),
                        na.rm=T)
    dDensMar <- density(raster::values(raster::mask(drainMean,soil[soil$SoilParent=="CaimitoMarineSedimentary",])),
                        na.rm=T)
    # By soil form
    dDensPal <- density(raster::values(raster::mask(drainMean,soil[soil$SoilForm=="PaleSwellingClay",])),
                        na.rm=T)
    dDensRed <- density(raster::values(raster::mask(drainMean,soil[soil$SoilForm=="RedLightClay",])),
                        na.rm=T)
    dDensMot <- density(raster::values(raster::mask(drainMean,soil[soil$SoilForm=="MottledHeavyClay",])),
                        na.rm=T)
    dDensBro <- density(raster::values(raster::mask(drainMean,soil[soil$SoilForm=="BrownFineLoam",])),
                        na.rm=T)
    
    axisSz <- 1.4
  # Make plot for forest age
    par(mfrow=c(2,3), las=1, mar=c(3,3,1,1))
    
    # OLD GROWTH
    plot(cDens,
         main = NA,
         col  = adjustcolor("grey",1),
         ylim=range(c(cDens$y,cDensOld$y,cDensSec$y)),
         xlim=range(c(cDens$x,cDensOld$x,cDensSec$x)),
         cex.axis=axisSz,
         xaxt="n",
         lwd=2)
    lines(cDensOld,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],1),
          lwd=2)
    plot(sDens,
         main = NA,
         col  = adjustcolor("grey",1),
         ylim=range(c(sDens$y,sDensOld$y,sDensSec$y)),
         xlim=range(c(sDens$x,sDensOld$x,sDensSec$x)),
         xaxt="n",
         lwd=2,
         cex.axis=axisSz)
    lines(sDensOld,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],1),
          lwd=2)
    plot(dDens,
         main = NA,
         col  = adjustcolor("grey",1),
         ylim=range(c(dDens$y,dDensOld$y,dDensSec$y)),
         xlim=range(c(dDens$x,dDensOld$x,dDensSec$x)),
         xaxt="n",
         lwd=2,
         cex.axis=axisSz)
    lines(dDensOld,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[1],1),
          lwd=2)
    legend(x=10,y=0.12,
           c("All BCI","Old growth"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Moonrise2",4)[1]),1),
           bty="n",
           cex=axisSz,
           lwd=2)
    
    # SECONDARY
    plot(cDens,
         main = NA,
         col  = adjustcolor("grey",1),
         ylim=range(c(cDens$y,cDensOld$y,cDensSec$y)),
         xlim=range(c(cDens$x,cDensOld$x,cDensSec$x)),
         cex.axis=axisSz,
         lwd=2)
    mtext("Curvature (LaPlacian convexity)", side=1, outer=F, line=2.5)
    lines(cDensSec,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensOld$y,sDensSec$y)),
         xlim=range(c(cDens$x,cDensOld$x,sDensSec$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    mtext("Slope (degrees)", side=1, outer=F, line=2.5)
    lines(sDensSec,
          main = NA,
          ylim=range(c(sDens$y,sDensOld$y,sDensSec$y)),
          col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensOld$y,dDensSec$y)),
         xlim=range(c(cDens$x,cDensOld$x,dDensSec$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensSec,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Moonrise2",4)[2],1),
          lwd=2)
    mtext("HAND (m)", side=1, outer=F, line=2.5)
    legend(x=10,y=0.12,
           c("All BCI","Secondary"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Moonrise2",4)[2]),1),
           cex=axisSz,
           bty="n",
           lwd=2)
    par(las=0)
    mtext("Frequency",side=2,outer=T)
    par(las=1)
    
  # Make plot for soil parent material
    par(mfrow=c(4,3), mar=c(3,3,0,1),oma=c(2,1,1,1))
    
    # CaimitoVolcanic
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensVol$y,cDensAnd$y,cDensMar$y,cDensBoh$y)),
         xlim=range(c(cDens$x,cDensVol$x,cDensAnd$x,cDensMar$x,cDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensVol,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensVol$y,sDensAnd$y,sDensMar$y,sDensBoh$y)),
         xlim=range(c(sDens$x,sDensVol$x,sDensAnd$x,sDensMar$x,sDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensVol,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensVol$y,dDensAnd$y,dDensMar$y,dDensBoh$y)),
         xlim=range(c(dDens$x,dDensVol$x,dDensAnd$x,dDensMar$x,dDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensVol,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[1],1),
          lwd=2)
    legend(x=7,y=0.18,
           c("All BCI","Caimito volcanic"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Chevalier1",4)[1]),1),
           cex = axisSz-0.2,
           bty="n",
           lwd=2)
    
    # Andesite
    plot(cDens,
         main = NA,
         col  = adjustcolor("grey",1),
         ylim=range(c(cDens$y,cDensVol$y,cDensAnd$y,cDensMar$y,cDensBoh$y)),
         xlim=range(c(cDens$x,cDensVol$x,cDensAnd$x,cDensMar$x,cDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         lwd=2)
    lines(cDensAnd,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensVol$y,sDensAnd$y,sDensMar$y,sDensBoh$y)),
         xlim=range(c(sDens$x,sDensVol$x,sDensAnd$x,sDensMar$x,sDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensAnd,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensVol$y,dDensAnd$y,dDensMar$y,dDensBoh$y)),
         xlim=range(c(dDens$x,dDensVol$x,dDensAnd$x,dDensMar$x,dDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensAnd,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[2],1),
          lwd=2)
    legend(x=7,y=0.18,
           c("All BCI","Andesite"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Chevalier1",4)[2]),1),
           bty="n",
           cex = axisSz-0.2,
           lwd=2)  
    
    # CaimitoMarineSedimentary
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensVol$y,cDensAnd$y,cDensMar$y,cDensBoh$y)),
         xlim=range(c(cDens$x,cDensVol$x,cDensAnd$x,cDensMar$x,cDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensMar,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensVol$y,sDensAnd$y,sDensMar$y,sDensBoh$y)),
         xlim=range(c(sDens$x,sDensVol$x,sDensAnd$x,sDensMar$x,sDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensMar,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensVol$y,dDensAnd$y,dDensMar$y,dDensBoh$y)),
         xlim=range(c(dDens$x,dDensVol$x,dDensAnd$x,dDensMar$x,dDensBoh$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensMar,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[3],1),
          lwd=2)
    legend(x=7,y=0.18,
           c("All BCI","Caimito marine"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Chevalier1",4)[3]),1),
           bty="n",
           cex = axisSz-0.2,
           lwd=2)
    
    # Bohio
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensVol$y,cDensAnd$y,cDensMar$y,cDensBoh$y)),
         xlim=range(c(cDens$x,cDensVol$x,cDensAnd$x,cDensMar$x,cDensBoh$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensBoh,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],1),
          lwd=2)
    mtext("Curvature (LaPlacian convexity)", side=1, outer=F, line=2.5)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensVol$y,sDensAnd$y,sDensMar$y,sDensBoh$y)),
         xlim=range(c(sDens$x,sDensVol$x,sDensAnd$x,sDensMar$x,sDensBoh$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensBoh,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],1),
          lwd=2)
    mtext("Slope (degrees)", side=1, outer=F, line=2.5)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensVol$y,dDensAnd$y,dDensMar$y,dDensBoh$y)),
         xlim=range(c(dDens$x,dDensVol$x,dDensAnd$x,dDensMar$x,dDensBoh$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensBoh,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Chevalier1",4)[4],1),
          lwd=2)
    mtext("HAND (m)", side=1, outer=F, line=2.5)
    legend(x=7,y=0.18,
           c("All BCI","Bohio"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Chevalier1",4)[4]),1),
           bty="n",
           cex = axisSz-0.2,
           lwd=2)  
    par(las=0)
    mtext("Frequency",side=2,outer=T)
    par(las=1)
    
  # Make plot for soil form
    par(mfrow=c(4,3), mar=c(3,3,0,1),oma=c(2,1,1,1))
    
    # BrownFineLoam
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensBro$y,cDensPal$y,cDensMot$y,cDensRed$y)),
         xlim=range(c(cDens$x,cDensBro$x,cDensPal$x,cDensMot$x,cDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensBro,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensBro$y,sDensPal$y,sDensMot$y,sDensRed$y)),
         xlim=range(c(sDens$x,sDensBro$x,sDensPal$x,sDensMot$x,sDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensBro,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensBro$y,dDensPal$y,dDensMot$y,dDensRed$y)),
         xlim=range(c(dDens$x,dDensBro$x,dDensPal$x,dDensMot$x,dDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensBro,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[1],1),
          lwd=2)
    legend(x=3,y=0.18,
           cex = axisSz - 0.2,
           c("All BCI","Brown fine loam"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Rushmore1",5)[1]),1),
           bty="n",
           lwd=2)
    
    # PaleSwellingClay
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensBro$y,cDensPal$y,cDensMot$y,cDensRed$y)),
         xlim=range(c(cDens$x,cDensBro$x,cDensPal$x,cDensMot$x,cDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensPal,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensBro$y,sDensPal$y,sDensMot$y,sDensRed$y)),
         xlim=range(c(sDens$x,sDensBro$x,sDensPal$x,sDensMot$x,sDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensPal,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensBro$y,dDensPal$y,dDensMot$y,dDensRed$y)),
         xlim=range(c(dDens$x,dDensBro$x,dDensPal$x,dDensMot$x,dDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensPal,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[3],1),
          lwd=2)
    legend(x=3,y=0.18,
           cex = axisSz - 0.2,
           c("All BCI","Pale swelling clay"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Rushmore1",5)[3]),1),
           bty="n",
           lwd=2)  
    
    # MottledHeavyClay
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensBro$y,cDensPal$y,cDensMot$y,cDensRed$y)),
         xlim=range(c(cDens$x,cDensBro$x,cDensPal$x,cDensMot$x,cDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensMot,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],1),
          lwd=2)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensBro$y,sDensPal$y,sDensMot$y,sDensRed$y)),
         xlim=range(c(sDens$x,sDensBro$x,sDensPal$x,sDensMot$x,sDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensMot,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],1),
          lwd=2)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensBro$y,dDensPal$y,dDensMot$y,dDensRed$y)),
         xlim=range(c(dDens$x,dDensBro$x,dDensPal$x,dDensMot$x,dDensRed$x)),
         cex.axis=axisSz,
         xaxt="n",
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensMot,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[4],1),
          lwd=2)
    legend(x=3,y=0.18,
           cex = axisSz - 0.2,
           c("All BCI","Mottled heavy clay"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Rushmore1",5)[4]),1),
           bty="n",
           lwd=2)
    
    # RedLightClay
    plot(cDens,
         main = NA,
         ylim=range(c(cDens$y,cDensBro$y,cDensPal$y,cDensMot$y,cDensRed$y)),
         xlim=range(c(cDens$x,cDensBro$x,cDensPal$x,cDensMot$x,cDensRed$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(cDensRed,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],1),
          lwd=2)
    mtext("Curvature (LaPlacian convexity)", side=1, outer=F, line=2.5)
    plot(sDens,
         main = NA,
         ylim=range(c(sDens$y,sDensBro$y,sDensPal$y,sDensMot$y,sDensRed$y)),
         xlim=range(c(sDens$x,sDensBro$x,sDensPal$x,sDensMot$x,sDensRed$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(sDensRed,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],1),
          lwd=2)
    mtext("Slope (degrees)", side=1, outer=F, line=2.5)
    plot(dDens,
         main = NA,
         ylim=range(c(dDens$y,dDensBro$y,dDensPal$y,dDensMot$y,dDensRed$y)),
         xlim=range(c(dDens$x,dDensBro$x,dDensPal$x,dDensMot$x,dDensRed$x)),
         cex.axis=axisSz,
         col  = adjustcolor("grey",1),
         lwd=2)
    lines(dDensRed,
          main = NA,
          col = adjustcolor(wesanderson::wes_palette("Rushmore1",5)[5],1),
          lwd=2)
    mtext("HAND (m)", side=1, outer=F, line=2.5)
    legend(x=3,y=0.18,
           cex = axisSz - 0.2,
           c("All BCI","Red light clay"),
           col  = adjustcolor(c("grey",wesanderson::wes_palette("Rushmore1",5)[5]),1),
           bty="n",
           lwd=2)    
    
    par(las=0)
    mtext("Frequency",side=2,outer=T)
    par(las=1)
#### Figure S6: Proportion of area with height decrease >=5 vs initial canopy height ####

# Look at initial canopy height and transitions per height class

# Load rasters    
gaps15to18 <- raster::raster("newGaps15to18_tin.tif")
gaps18to20 <- raster::raster("newGaps18to20_tin.tif")

chm09 <- raster::raster("CHM_2009_QAQC.tif")
chm15 <- raster::raster("CHM_2015_QAQC_tin.tif")
chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
chm20 <- raster::raster("CHM_2020_QAQC_tin.tif")

# Make sure all have same extents (2009 slightly shorter)
chm15 <- raster::crop(chm15, raster::extent(chm09))
chm18 <- raster::crop(chm18, raster::extent(chm09))
chm20 <- raster::crop(chm20, raster::extent(chm09))

chm09_vals <- raster::values(chm09)
chm15_vals <- raster::values(chm15)
chm18_vals <- raster::values(chm18)
chm20_vals <- raster::values(chm20)

gaps15to18_vals <- raster::values(gaps15to18)
gaps18to20_vals <- raster::values(gaps18to20)

chm15_vals[is.na(chm18_vals)] <- NA
chm18_vals[is.na(chm20_vals)] <- NA


allStart <- floor(max(c(chm15_vals,chm18_vals),na.rm=T))

propGap <- data.frame(start = 5:allStart,
                      n15 = NA,
                      propGap15 = NA,
                      n18 = NA,
                      propGap18 = NA,
                      n20=NA,
                      n09=NA)

for(i in 1:nrow(propGap)){
  
  vals15 <- which(!is.na(chm15_vals) & chm15_vals>=propGap$start[i] & chm15_vals<(propGap$start[i]+1))
  gaps15 <- gaps15to18_vals[vals15]
  propGap$n15[i] <- length(vals15)/length(chm15_vals[!is.na(chm15_vals)])
  propGap$propGap15[i] <- length(gaps15[!is.na(gaps15)])/length(vals15)

  vals18 <- which(!is.na(chm18_vals) & chm18_vals>=propGap$start[i] & chm18_vals<(propGap$start[i]+1))
  gaps18 <- gaps18to20_vals[vals18]
  propGap$n18[i] <- length(vals18)/length(chm18_vals[!is.na(chm18_vals)])
  propGap$propGap18[i] <- length(gaps18[!is.na(gaps18)])/length(vals18)
  
  vals09 <- which(!is.na(chm09_vals) & chm09_vals>=propGap$start[i] & chm09_vals<(propGap$start[i]+1))
  propGap$n09[i] <- length(vals09)/length(chm09_vals[!is.na(chm09_vals)])
  
  vals20 <- which(!is.na(chm20_vals) & chm20_vals>=propGap$start[i] & chm20_vals<(propGap$start[i]+1))
  propGap$n20[i] <- length(vals20)/length(chm20_vals[!is.na(chm20_vals)])
  
}


# Normalize the proportion of gaps observed to per year
nYr15to18 <- as.numeric(as.Date("2018-06-07") - as.Date("2015-06-26"))/365
nYr18to20 <- as.numeric(as.Date("2020-07-31") - as.Date("2018-06-07"))/365

col18 <- "blue"
col20 <- "#d95f02"

# Plot canopy height distributions
par(mfrow=c(2,1), mar=c(0,4,0,1),oma=c(4,1,2,1))
plot(n15~start,
     xlim = c(5,50),
     data = propGap,
     type="l",
     ylim=c(0,0.06),
     xlab = NA,
     xaxt="n",
     bty= "n",
     ylab = "Proportion of total area",
     col = adjustcolor("black",0.6),
     lty=1,
     lwd=3)
abline(h=0)
# lines(n09~start,
#       data = propGap,
#       col=adjustcolor("orange",0.6), lwd=3)

lines(n18~start,
      data = propGap,
      col=adjustcolor(col18,0.6), lwd=3)

lines(n20~start,
      data = propGap,
      col=adjustcolor(col20,0.6), lwd=3)

legend(x=35,y=0.06,
       bty="n",
       c("2015","2018","2020"),
       lty=c(1,1,1),
       col=adjustcolor(c("black",col18,col20),0.6), lwd=3)



# Plot probability of becoming a gap for area in each height bin


plot(propGap15/nYr15to18~start,
     data = propGap,
     col = adjustcolor(col18,0.6),
     type="l",
     ylim=c(0,0.045),
     xlim=c(5,50),
     xlab = "Initial canopy height (m)",
     ylab = "Proportion of height area decreasing > 5m",
     lwd=3)

lines(propGap18/nYr18to20~start, data=propGap,
      col=adjustcolor(col20,0.6), lwd=3)

legend(x=5,y=0.04,
       bty="n",
       c("2015-2018","2018-2020"),
       col=adjustcolor(c(col18,col20),0.6), lwd=3, lty=c(1,1))








#### Figure S#: CDF of gap size and distribution of 2009 canopy height per forest type ####

# Define forest age and soil type polygons

  # Forest age polygon
    age <- rgdal::readOGR("D:/BCI_Spatial/Enders_Forest_Age_1935/Ender_Forest_Age_1935.shp")
    age$AgeClass <- "Other"
    age$AgeClass[age$Mascaro_Co == "> 400"] <- "OldGrowth"
    age$AgeClass[age$Mascaro_Co %in% c("80-110", "120-130")] <- "Secondary"
    ageUse <- age[!(age$AgeClass=="Other"),]
  
  # Soil type polygon  
    soil <- rgdal::readOGR("D:/BCI_Spatial/BCI_Soils/BCI_Soils.shp")
    soil <- sp::spTransform(soil,raster::crs(age))
  
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

# Read 2009 CHM
  chm09 <- raster::raster("CHM_2009_QAQC.tif")

# Make density plots
  
  # Whole island
    chm_all <- density(raster::values(chm09),
                       n = 512, from = 0, to = 70, na.rm=T)

  # By forest age
    chm_Old <- density(raster::values(raster::mask(chm09, ageUse[ageUse$AgeClass=="OldGrowth",])), 
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Sec <- density(raster::values(raster::mask(chm09, ageUse[ageUse$AgeClass=="Secondary",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    
  # By soil parent material
    chm_And <- density(raster::values(raster::mask(chm09, soil[soil$SoilParent=="Andesite",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Boh <- density(raster::values(raster::mask(chm09, soil[soil$SoilParent=="Bohio",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Mar <- density(raster::values(raster::mask(chm09, soil[soil$SoilParent=="CaimitoMarineSedimentary",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Vol <- density(raster::values(raster::mask(chm09, soil[soil$SoilParent=="CaimitoVolcanic",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    
  # By soil form
    chm_Bro <- density(raster::values(raster::mask(chm09, soil[soil$SoilForm=="BrownFineLoam",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Mot <- density(raster::values(raster::mask(chm09, soil[soil$SoilForm=="MottledHeavyClay",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Pal <- density(raster::values(raster::mask(chm09, soil[soil$SoilForm=="PaleSwellingClay",])),
                       n = 512, from = 0, to = 70, na.rm=T)
    chm_Red <- density(raster::values(raster::mask(chm09, soil[soil$SoilForm=="RedLightClay",])),
                       n = 512, from = 0, to = 70, na.rm=T)

  # Make canopy height CDFs
    htRange <- seq(0,70,length.out = 512)
    
    htCDFs <- data.frame(ht = htRange,
                          Age_OldGrowth = NA,
                          Age_Secondary = NA,
                          Parent_Bohio = NA,
                          Parent_CaimitoVolcanic=NA,
                          Parent_CaimitoMarineSedimentary=NA,
                          Parent_Andesite=NA,
                          Form_RedLightClay=NA,
                          Form_BrownFineLoam=NA,
                          Form_PaleSwellingClay=NA,
                          Form_MottledHeavyClay=NA)
    for(i in 1:length(htRange)){
      htCDFs$Age_OldGrowth[i] <- sum(chm_Old$y[chm_Old$x <= htCDFs$ht[i]])/sum(chm_Old$y)
      htCDFs$Age_Secondary[i] <- sum(chm_Sec$y[chm_Sec$x <= htCDFs$ht[i]])/sum(chm_Sec$y)
      
      htCDFs$Parent_Bohio[i] <- sum(chm_Boh$y[chm_Boh$x <= htCDFs$ht[i]])/sum(chm_Boh$y)
      htCDFs$Parent_CaimitoVolcanic[i] <- sum(chm_Vol$y[chm_Vol$x <= htCDFs$ht[i]])/sum(chm_Vol$y)
      htCDFs$Parent_CaimitoMarineSedimentary[i] <- sum(chm_Mar$y[chm_Mar$x <= htCDFs$ht[i]])/sum(chm_Mar$y)
      htCDFs$Parent_Andesite[i] <- sum(chm_And$y[chm_And$x <= htCDFs$ht[i]])/sum(chm_And$y)
      
      htCDFs$Form_RedLightClay[i] <- sum(chm_Red$y[chm_Red$x <= htCDFs$ht[i]])/sum(chm_Red$y)
      htCDFs$Form_BrownFineLoam[i] <- sum(chm_Bro$y[chm_Bro$x <= htCDFs$ht[i]])/sum(chm_Bro$y)
      htCDFs$Form_PaleSwellingClay[i] <- sum(chm_Pal$y[chm_Pal$x <= htCDFs$ht[i]])/sum(chm_Pal$y)
      htCDFs$Form_MottledHeavyClay[i] <- sum(chm_Mot$y[chm_Mot$x <= htCDFs$ht[i]])/sum(chm_Mot$y)
    }
       
  # Calculate cumulative proportion of area in gaps with initial canopy height

    # Load rasters    
      gaps15to18 <- raster::raster("newGaps15to18_tin.tif")
      gaps18to20 <- raster::raster("newGaps18to20_tin.tif")
      gaps15to18_vals <- raster::values(gaps15to18)
      gaps18to20_vals <- raster::values(gaps18to20)
      
      chm15 <- raster::raster("CHM_2015_QAQC_tin.tif")
      chm18 <- raster::raster("CHM_2018_QAQC_tin.tif")
      
    # Normalize the proportion of gaps observed to per year
      nYr15to18 <- as.numeric(as.Date("2018-06-07") - as.Date("2015-06-26"))/365
      nYr18to20 <- as.numeric(as.Date("2020-07-31") - as.Date("2018-06-07"))/365  

      
    # Make separate rasters separately for forest age, parent material, and soil form
      chm15_OldGrowth <- raster::values(raster::mask(chm15, ageUse[ageUse$AgeClass=="OldGrowth",])); chm15_OldGrowth[chm15_OldGrowth<10] <- NA
      chm15_Secondary <- raster::values(raster::mask(chm15, ageUse[ageUse$AgeClass=="Secondary",])); chm15_Secondary[chm15_Secondary<10] <- NA
      chm18_OldGrowth <- raster::values(raster::mask(chm18, ageUse[ageUse$AgeClass=="OldGrowth",])); chm18_OldGrowth[chm18_OldGrowth<10] <- NA
      chm18_Secondary <- raster::values(raster::mask(chm18, ageUse[ageUse$AgeClass=="Secondary",])); chm18_Secondary[chm18_Secondary<10] <- NA
      
      chm15_Bohio <- raster::values(raster::mask(chm15, soil[soil$SoilParent=="Bohio",])); chm15_Bohio[chm15_Bohio<10] <- NA
      chm15_CaimitoVolcanic <- raster::values(raster::mask(chm15, soil[soil$SoilParent=="CaimitoVolcanic",])); chm15_CaimitoVolcanic[chm15_CaimitoVolcanic<10] <- NA
      chm15_CaimitoMarineSedimentary <- raster::values(raster::mask(chm15, soil[soil$SoilParent=="CaimitoMarineSedimentary",])); chm15_CaimitoMarineSedimentary[chm15_CaimitoMarineSedimentary<10] <- NA
      chm15_Andesite <- raster::values(raster::mask(chm15, soil[soil$SoilParent=="Andesite",])); chm15_Andesite[chm15_Andesite<10] <- NA
      chm18_Bohio <- raster::values(raster::mask(chm18, soil[soil$SoilParent=="Bohio",])); chm18_Bohio[chm18_Bohio<10] <- NA
      chm18_CaimitoVolcanic <- raster::values(raster::mask(chm18, soil[soil$SoilParent=="CaimitoVolcanic",])); chm18_CaimitoVolcanic[chm18_CaimitoVolcanic<10] <- NA
      chm18_CaimitoMarineSedimentary <- raster::values(raster::mask(chm18, soil[soil$SoilParent=="CaimitoMarineSedimentary",])); chm18_CaimitoMarineSedimentary[chm18_CaimitoMarineSedimentary<10] <- NA
      chm18_Andesite <- raster::values(raster::mask(chm18, soil[soil$SoilParent=="Andesite",])); chm18_Andesite[chm18_Andesite<10] <- NA
      
      chm15_RedLightClay <- raster::values(raster::mask(chm15, soil[soil$SoilForm=="RedLightClay",])); chm15_RedLightClay[chm15_RedLightClay<10] <- NA
      chm15_BrownFineLoam <- raster::values(raster::mask(chm15, soil[soil$SoilForm=="BrownFineLoam",])); chm15_BrownFineLoam[chm15_BrownFineLoam<10] <- NA
      chm15_PaleSwellingClay <- raster::values(raster::mask(chm15, soil[soil$SoilForm=="PaleSwellingClay",])); chm15_PaleSwellingClay[chm15_PaleSwellingClay<10] <- NA
      chm15_MottledHeavyClay <- raster::values(raster::mask(chm15, soil[soil$SoilForm=="MottledHeavyClay",])); chm15_MottledHeavyClay[chm15_MottledHeavyClay<10] <- NA
      chm18_RedLightClay <- raster::values(raster::mask(chm18, soil[soil$SoilForm=="RedLightClay",])); chm18_RedLightClay[chm18_RedLightClay<10] <- NA
      chm18_BrownFineLoam <- raster::values(raster::mask(chm18, soil[soil$SoilForm=="BrownFineLoam",])); chm18_BrownFineLoam[chm18_BrownFineLoam<10] <- NA
      chm18_PaleSwellingClay <- raster::values(raster::mask(chm18, soil[soil$SoilForm=="PaleSwellingClay",])); chm18_PaleSwellingClay[chm18_PaleSwellingClay<10] <- NA
      chm18_MottledHeavyClay <- raster::values(raster::mask(chm18, soil[soil$SoilForm=="MottledHeavyClay",])); chm18_MottledHeavyClay[chm18_MottledHeavyClay<10] <- NA
      
      htRange <- seq(10,50,length.out = 128)
      
      gapCDFs <- data.frame(ht = htRange,
                            Age_OldGrowth = NA,
                            Age_Secondary = NA,
                            Parent_Bohio = NA,
                            Parent_CaimitoVolcanic=NA,
                            Parent_CaimitoMarineSedimentary=NA,
                            Parent_Andesite=NA,
                            Form_RedLightClay=NA,
                            Form_BrownFineLoam=NA,
                            Form_PaleSwellingClay=NA,
                            Form_MottledHeavyClay=NA)
      
      for(i in 1:nrow(gapCDFs)){
        gapCDFs$Age_OldGrowth[i] <- (length(chm15_OldGrowth[!is.na(chm15_OldGrowth) & chm15_OldGrowth <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_OldGrowth[!is.na(chm18_OldGrowth) & chm18_OldGrowth <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_OldGrowth[!is.na(chm15_OldGrowth)]) + length(chm18_OldGrowth[!is.na(chm18_OldGrowth)]))
        gapCDFs$Age_Secondary[i] <- (length(chm15_Secondary[!is.na(chm15_Secondary) & chm15_Secondary <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_Secondary[!is.na(chm18_Secondary) & chm18_Secondary <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_Secondary[!is.na(chm15_Secondary)]) + length(chm18_Secondary[!is.na(chm18_Secondary)]))
        
        gapCDFs$Parent_Bohio[i] <- (length(chm15_Bohio[!is.na(chm15_Bohio) & chm15_Bohio <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_Bohio[!is.na(chm18_Bohio) & chm18_Bohio <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_Bohio[!is.na(chm15_Bohio)]) + length(chm18_Bohio[!is.na(chm18_Bohio)]))
        gapCDFs$Parent_CaimitoVolcanic[i] <- (length(chm15_CaimitoVolcanic[!is.na(chm15_CaimitoVolcanic) & chm15_CaimitoVolcanic <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_CaimitoVolcanic[!is.na(chm18_CaimitoVolcanic) & chm18_CaimitoVolcanic <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_CaimitoVolcanic[!is.na(chm15_CaimitoVolcanic)]) + length(chm18_CaimitoVolcanic[!is.na(chm18_CaimitoVolcanic)]))
        gapCDFs$Parent_CaimitoMarineSedimentary[i] <- (length(chm15_CaimitoMarineSedimentary[!is.na(chm15_CaimitoMarineSedimentary) & chm15_CaimitoMarineSedimentary <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_CaimitoMarineSedimentary[!is.na(chm18_CaimitoMarineSedimentary) & chm18_CaimitoMarineSedimentary <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_CaimitoMarineSedimentary[!is.na(chm15_CaimitoMarineSedimentary)]) + length(chm18_CaimitoMarineSedimentary[!is.na(chm18_CaimitoMarineSedimentary)]))
        gapCDFs$Parent_Andesite[i] <- (length(chm15_Andesite[!is.na(chm15_Andesite) & chm15_Andesite <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_Andesite[!is.na(chm18_Andesite) & chm18_Andesite <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_Andesite[!is.na(chm15_Andesite)]) + length(chm18_Andesite[!is.na(chm18_Andesite)]))
        
        gapCDFs$Form_RedLightClay[i] <- (length(chm15_RedLightClay[!is.na(chm15_RedLightClay) & chm15_RedLightClay <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_RedLightClay[!is.na(chm18_RedLightClay) & chm18_RedLightClay <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_RedLightClay[!is.na(chm15_RedLightClay)]) + length(chm18_RedLightClay[!is.na(chm18_RedLightClay)]))
        gapCDFs$Form_BrownFineLoam[i] <- (length(chm15_BrownFineLoam[!is.na(chm15_BrownFineLoam) & chm15_BrownFineLoam <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_BrownFineLoam[!is.na(chm18_BrownFineLoam) & chm18_BrownFineLoam <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_BrownFineLoam[!is.na(chm15_BrownFineLoam)]) + length(chm18_BrownFineLoam[!is.na(chm18_BrownFineLoam)]))
        gapCDFs$Form_PaleSwellingClay[i] <- (length(chm15_PaleSwellingClay[!is.na(chm15_PaleSwellingClay) & chm15_PaleSwellingClay <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_PaleSwellingClay[!is.na(chm18_PaleSwellingClay) & chm18_PaleSwellingClay <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_PaleSwellingClay[!is.na(chm15_PaleSwellingClay)]) + length(chm18_PaleSwellingClay[!is.na(chm18_PaleSwellingClay)]))
        gapCDFs$Form_MottledHeavyClay[i] <- (length(chm15_MottledHeavyClay[!is.na(chm15_MottledHeavyClay) & chm15_MottledHeavyClay <= gapCDFs$ht[i] & !is.na(gaps15to18_vals)])/nYr15to18 + length(chm18_MottledHeavyClay[!is.na(chm18_MottledHeavyClay) & chm18_MottledHeavyClay <= gapCDFs$ht[i] & !is.na(gaps18to20_vals)])/nYr18to20)/(
          length(chm15_MottledHeavyClay[!is.na(chm15_MottledHeavyClay)]) + length(chm18_MottledHeavyClay[!is.na(chm18_MottledHeavyClay)]))
        
        print(i)
      }
    
    yLimVal_a <- range(gapCDFs[,-1])
    yLimVal_b <- range(c(chm_all$y, chm_Old$y, chm_Sec$y,
                         chm_And$y, chm_Boh$y, chm_Mar$y, chm_Vol$y,
                         chm_Bro$y, chm_Mot$y, chm_Pal$y, chm_Red$y))  
    cxAxis = 1.2
    
    par(mfrow=c(2,3), mar=c(1,1,0,1), oma=c(3,5,2,1), las=1)
  
  # GAP CUMULATIVE AREA
    
    # age
      plot(Age_OldGrowth~ht, data = gapCDFs,
           type = "l",
           xaxt="n",
           xlim = c(0,50),
           ylim = yLimVal_a,
           lwd = 2,
           col = adjustcolor(colOld, 0.8),
           cex.axis = cxAxis)
      lines(Age_Secondary~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colSec, 0.8))
      mtext("Forest age", side=3, outer=F, line=0.5)
      legend(x=0,y=0.028,
             c("Old growth","Secondary"),
             col=adjustcolor(c(colOld,colSec),0.8),
             lwd=2,
             bty="n")
      mtext("Cumulative disturbance area (% yr-1)", side=2, outer=F, line=4, las=0, cex=0.8)
      
      
      # parent material 
      plot(Parent_Bohio~ht, data = gapCDFs,
           type = "l",
           xaxt="n",yaxt="n",
           xlim = c(0,50),
           ylim = yLimVal_a,
           lwd = 2,
           col = adjustcolor(colBoh, 0.8),
           cex.axis = cxAxis)
      lines(Parent_CaimitoVolcanic~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colVol, 0.8))
      lines(Parent_CaimitoMarineSedimentary~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colMar, 0.8))
      lines(Parent_Andesite~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colAnd, 0.8))
      mtext("Soil parent material", side=3, outer=F, line=0.5)
      legend(x=0,y=0.028,
             c("Andesite","Bohio","Caimito marine","Caimito volcanic"),
             col=adjustcolor(c(colAnd,colBoh,colMar,colVol),0.8),
             lwd=2,
             bty="n")
      
      # soil form
      plot(Form_RedLightClay~ht, data = gapCDFs,
           type = "l",
           xaxt="n",yaxt="n",
           xlim = c(0,50),
           ylim = yLimVal_a,
           lwd = 2,
           col = adjustcolor(colRed, 0.8),
           cex.axis = cxAxis)
      lines(Form_BrownFineLoam~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colBro, 0.8))
      lines(Form_PaleSwellingClay~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colPal, 0.8))
      lines(Form_MottledHeavyClay~ht, data = gapCDFs,
            lwd = 2,
            col = adjustcolor(colMot, 0.8))
      mtext("Soil form", side=3, outer=F, line=0.5)
      legend(x=0,y=0.028,
             c("Fine loam","Mottled heavy clay","Pale swelling clay","Red light clay"),
             col=adjustcolor(c(colBro,colMot,colPal,colRed),0.8),
             lwd=2,
             bty="n")
      
    # CANOPY HEIGHT DISTRIBUTIONS
      # age
      plot(chm_Old,
           type = "l",
           main = NA,
           xlim = c(0,50),
           ylim=yLimVal_b,
           lwd = 2,
           col = adjustcolor(colOld, 0.8),
           cex.axis = cxAxis)
      lines(chm_Sec,
            lwd = 2,
            col = adjustcolor(colSec, 0.8))
      mtext("Canopy height distribution", side=2, outer=F, line=4, las=0, cex=0.8)
      
    
      # parent material
      plot(chm_Boh,
           type = "l",
           yaxt = "n",
           main = NA,
           xlim = c(0,50),
           ylim=yLimVal_b,
           lwd = 2,
           col = adjustcolor(colBoh, 0.8),
           cex.axis = cxAxis)
      lines(chm_Vol,
            lwd = 2,
            col = adjustcolor(colVol, 0.8))
      lines(chm_Mar,
            lwd = 2,
            col = adjustcolor(colMar, 0.8))
      lines(chm_And,
            lwd = 2,
            col = adjustcolor(colAnd, 0.8))
      
      # soil form    
        plot(chm_Red,
             type = "l",
             yaxt = "n",
             main = NA,
             xlim = c(0,50),
             ylim=yLimVal_b,
             lwd = 2,
             col = adjustcolor(colRed, 0.8),
             cex.axis = cxAxis)
        lines(chm_Bro,
              lwd = 2,
              col = adjustcolor(colBro, 0.8))
        lines(chm_Pal,
              lwd = 2,
              col = adjustcolor(colPal, 0.8))
        lines(chm_Mot,
              lwd = 2,
              col = adjustcolor(colMot, 0.8))
        
        mtext("Canopy height (m)", side=1, outer=T, line=1)
    
    
  # OLD PLOTS   
  # # Make plot per class
  #   # Find y-axis limit
  #   yLimVal <- range(c(chm_all$y, chm_Old$y, chm_Sec$y,
  #                      chm_And$y, chm_Boh$y, chm_Mar$y, chm_Vol$y,
  #                      chm_Bro$y, chm_Mot$y, chm_Pal$y, chm_Red$y))
  #   
  #   axSize <- 1.2
  #   
  #   par(mfcol=c(4,3), mar=c(1,1,0,1), oma=c(3,5,2,1))
  #   
  #   # Plot forest age
  #     plot(chm_all,
  #          xaxt="n",
  #          ylim=yLimVal,
  #          col = adjustcolor("grey",0.6),
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Old, 
  #           col = adjustcolor(colOld,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Old growth"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colOld),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     mtext("Forest age", side=3, line = 0.5)
  #     
  #     
  #     plot(chm_all,
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Sec, 
  #           col = adjustcolor(colSec,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Secondary"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colSec),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #       # Add empty plots
  #     plot(0,type='n',axes=FALSE,ann=FALSE)
  #     plot(0,type='n',axes=FALSE,ann=FALSE)
  #   
  #   # Plot soil parent material
  #     plot(chm_all,
  #          xaxt="n",
  #          yaxt="n",
  #          ylim=yLimVal,
  #          col = adjustcolor("grey",0.6),
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Vol, 
  #           col = adjustcolor(colVol,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Caimito volcanic"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colVol),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     mtext("Soil parent material", side=3, line = 0.5)
  #     
  #     
  #     plot(chm_all,
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          xaxt="n",
  #          yaxt="n",
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_And, 
  #           col = adjustcolor(colAnd,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Andesite"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colAnd),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     
  #     plot(chm_all,
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          xaxt="n",
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Mar, 
  #           col = adjustcolor(colMar,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Caimito marine"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colMar),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     
  #     plot(chm_all,
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Boh, 
  #           col = adjustcolor(colBoh,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Bohio"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colBoh),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     
  #   # Plot soil form
  #     plot(chm_all,
  #          xaxt="n",
  #          yaxt="n",
  #          ylim=yLimVal,
  #          col = adjustcolor("grey",0.6),
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Bro, 
  #           col = adjustcolor(colBro,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Brown fine loam"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colBro),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     mtext("Soil form", side=3, line = 0.5)
  #     
  #     
  #     plot(chm_all,
  #          yaxt="n",
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          xaxt="n",
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Pal, 
  #           col = adjustcolor(colPal,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Pale swelling clay"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colPal),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     
  #     plot(chm_all,
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          xaxt="n",
  #          yaxt="n",
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Mot, 
  #           col = adjustcolor(colMot,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Mottled heavy clay"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colMot),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     
  #     plot(chm_all,
  #          yaxt="n",
  #          col = adjustcolor("grey",0.6),
  #          ylim=yLimVal,
  #          lwd=2,
  #          main = NA,
  #          cex.axis = axSize)
  #     lines(chm_Red, 
  #           col = adjustcolor(colRed,0.6),
  #           lwd=2)
  #     legend(x = 30, y = 0.05,
  #            c("All BCI", "Red light clay"),
  #            bty = "n",
  #            col = adjustcolor(c("grey",colRed),0.6),
  #            lwd = 2,
  #            cex = axSize)
  #     mtext("Canopy height (m)", outer=T, side = 1, line = 1.5)
  #     par(las=0)
  #     mtext("Frequency", outer=T, side = 2, line = 2.5)
  #     par(las=1)
      