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
    
  # Make plot of corrected and uncorrected height change
    colBrks <- c(-300,-200,-30,-20,-10,-3,3,10,20,30,200)
    colPal <- colorRampPalette(c("grey","red","orangered","darksalmon","khaki2",
                                  "white",
                                  "skyblue1","skyblue3","cornflowerblue","darkblue"))
    
    pdf("FigureS1_output.pdf",width=8,height=7)
      par(mfrow=c(2,2), mar=c(0,1,1,1), oma=c(1,1,1,3))
      raster::plot(raw15to18,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA)
      raster::plot(buffer,add=T)
      
      raster::plot(cor15to18,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA)
      raster::plot(buffer,add=T)
      
      raster::plot(raw18to20,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA)
      raster::plot(buffer,add=T)
      
      raster::plot(cor18to20,
                   bty="n", box=F,yaxt="n",xaxt="n",
                   col = colPal(length(colBrks)-1),
                   breaks = colBrks,
                   main = NA)
      raster::plot(buffer,add=T)
    dev.off()
      


#### Figure S2 (2015 data correction): in "Validation_50HaPlot.R" ####
#### Figure S3: Results from disturbance validation in 50 ha plot ####
    
## Panel a: plot monthly disturbances
    valGaps15to18 <- rgdal::readOGR("toCheckRaquel/toCheckRaquel.shp")
    valGapsSzHt <- valGaps15to18[valGaps15to18$area2 >= 25 & valGaps15to18$htDrop <= -5 & !is.na(valGaps15to18$htDrop),]
    
    # make a jiggered date field to plot
    valGapsSzHt$datClss <- as.Date(valGapsSzHt$datClss)
    valGapsSzHt$datePlot <-  valGapsSzHt$datClss + rnorm(mean = 0, sd = 7, n = nrow(valGapsSzHt))
    
    # Only keep good gaps (remove false positives) before plotting/calculating summary stats
    valGapsUse <- valGapsSzHt[!(valGapsSzHt$vslChck==2),]
    
    par(mar=c(3,3,1,1), oma = c(1,2,1,0), mfrow=c(1,2), las=1)
    
    # observed with ht and area agreement
    plot(area2~datePlot, data=valGapsUse[valGapsUse$vslChck %in% c(0,3:5),],
         ylim=c(20,1400),
         pch=20,
         col=adjustcolor("grey",0.6),
         log="y",
         ylab=NA,
         xlab = NA)
    mtext(expression("Gap area " (m^2)),side=2,outer=T, las=0)
    mtext(expression("Disturbance date"),side=1,outer=F, las=0, line=2.5)
    mtext(expression("Monthly disturbances"),side=3,outer=F, las=0, line=0)
    
    # red: monthly data are correct but missed by me
    points(area2~datePlot, data=valGapsUse[valGapsUse$observd==F & valGapsUse$vslChck==1,],
           pch=20,
           col=adjustcolor("red",0.6))
    
    legend(x=as.Date("2017-01-01"),
           y=1800,
           bty="n",
           c("Observed",
             "Not observed"),
           col=adjustcolor(c("grey","red"),c(0.6)),
           pch=20,
           pt.cex=1.5,
           cex=0.8)
    text("a",
         x=as.Date("2015-07-01"),
         y=1400)
    
    # Don't plot these details for monthly data
    # # BLUE: likely size discrepancy (potential overestimation in monthly data)
    #   points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==c(3,4),],
    #          pch=20,
    #          col=adjustcolor("blue",0.6))
    # 
    # # RED: false positive in monthly data
    #   # Total false positive
    #   points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==2,],
    #          pch=20,
    #          col=adjustcolor("red",0.6))
    #   
    # # PURPLE: mismatch in space
    #   # Spatial alignment mismatch
    #   points(area2~datePlot, data=valGapsSzHt[valGapsSzHt$observd==F & valGapsSzHt$vslChck==5,],
    #          pch=20,
    #          col=adjustcolor("purple",0.6))
    
    
    
    # legend(x=as.Date("2017-01-01"),
    #        y=650,
    #        bty="n",
    #        c("Correctly observd in annual data",
    #          "Likely size discrepancy (monthly overestimate?)",
    #          "False positive",
    #          "Not observd in annual data and correct",
    #          "Spatial misalignment with annual data"),
    #        col=adjustcolor(c("grey","lightblue","red","green","purple"),c(0.6)),
    #       pch=20,
    #       pt.cex=1.5,
    #       cex=0.8)
    
    
    sum(valGapsUse$area2[valGapsUse$vslChck %in% c(0,3:5)])/sum(valGapsUse$area2)
    sum(valGapsUse$area2[valGapsUse$observd==F])/sum(valGapsUse$area2)
    sum(valGapsUse$area2[valGapsUse$observd==F  & valGapsUse$vslChck %in% c(3,4)])/sum(valGapsUse$area2)
    sum(valGapsUse$area2[valGapsUse$observd==F  & valGapsUse$vslChck %in% c(1)])/sum(valGapsUse$area2)
    sum(valGapsUse$area2[valGapsUse$observd==F  & valGapsUse$vslChck %in% c(5)])/sum(valGapsUse$area2)
    
    length(valGapsUse$area2[valGapsUse$vslChck %in% c(0,2:5)])/length(valGapsUse$area2)
    length(valGapsUse$area2[valGapsUse$observd==F])/length(valGapsUse$area2)
    length(valGapsUse$area2[valGapsUse$observd==F  & valGapsUse$vslChck %in% c(3,4)])/length(valGapsUse$area2)
    length(valGapsUse$area2[valGapsUse$observd==F  & valGapsUse$vslChck %in% c(1)])/length(valGapsUse$area2)
    length(valGapsUse$area2[valGapsUse$observd==F  & valGapsUse$vslChck %in% c(5)])/length(valGapsUse$area2)
    
    
## Panel b: plot multiannual disturbances
    
    plotGaps <- rgdal::readOGR("toCheckKC_final/toCheckKC_final.shp")
    
    
    # Observed with ht and area agreement
    plot(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==T,],
         ylim=c(20,1400),
         xlim=c(0.05,0.65),
         pch=20,
         col=adjustcolor("grey",0.6),
         log="y",
         ylab=NA,
         xlab = NA)
    mtext(expression("Gap circularity " (4*pi*"Area"/"Perim"^2)),side=1,outer=F, las=0, line=2.5)
    mtext(expression("Multiannual disturbances"),side=3,outer=F, las=0, line=0)
    
    # mismatch in space/time [consider this OK because not an error in either]
    # Spatial alignment mismatch
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck==5,],
           pch=20,
           col=adjustcolor("grey",0.6))
    # Temporal mismatch (i.e. I see beginning of slowly dying tree and Raquel records the later fall)
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(6,8),],
           pch=20,
           col=adjustcolor("grey",0.6))
    
    # BLUE: observed by Raquel's data when constraints are not imposed
    # Observed by Raquel's data are less than 25m2
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==F & plotGaps$obsrvdH==T,],
           pch=20,
           col=adjustcolor("grey",0.6))
    # Observed by Raquel's data have ht drop less than 5 m
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==F & !(plotGaps$vslChck==9),],
           pch=20,
           col=adjustcolor("grey",0.6))
    # Observed by Raquel's data are less than 25m2 AND have ht drop less than 5m (NONE)
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==F & plotGaps$obsrvdH==F & !(plotGaps$vslChck==9),],
           pch=20,
           col=adjustcolor("grey",0.5))
    
    # RED: false positive in my data
    # Total false positive
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck==2,],
           pch=20,
           col=adjustcolor("red",0.6))
    # Overestimation of gap area
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck==10,],
           pch=20,
           col=adjustcolor("red",0.6))
    
    
    # PURPLE: my data are correct
    # Unclear why Raquel misses this gap
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck==1,],
           pch=20,
           col=adjustcolor("purple",0.6))
    # Decaying tree: total miss
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$obsrvdAl==F & plotGaps$vslChck==7,],
           pch=20,
           col=adjustcolor("purple",0.6))
    # Decaying tree: partial miss
    points(area2~ratiCrc, data=plotGaps@data[plotGaps$vslChck==9,],
           pch=20,
           col=adjustcolor("purple",0.6))
    
    legend(x=0.2,
           y=1800,
           bty="n",
           c("True: observed in monthly data",
             "True: visually confirmed",
             "False positive"),
           col=adjustcolor(c("grey","purple","red"),c(0.6)),
           pch=20,
           pt.cex=1.5,
           cex=0.8)
    text("b",
         x=0.05,
         y=1400)
    
    # By area  
    # Total area in true gaps-- from comparison with monthly data OR confirmed visually
    sum(plotGaps$area2[plotGaps$obsrvdAl==T | plotGaps$vslChck %in% c(1,4:9)])/sum(plotGaps$area2)
    
    sum(plotGaps$area2[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==T])/sum(plotGaps$area2)
    sum(plotGaps$area2[plotGaps$obsrvdAl==T])/sum(plotGaps$area2)
    
    # False positive area
    sum(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(2,10)])/sum(plotGaps$area2)
    
    # Area missing from monthly data
    sum(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(1)])/sum(plotGaps$area2)
    sum(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(7,9)])/sum(plotGaps$area2)
    # Spatial or temporal mismatch
    sum(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(5,6,8)])/sum(plotGaps$area2)
    
    # By number of gaps
    # Total area in true gaps-- from comparison with monthly data OR confirmed visually
    length(plotGaps$area2[plotGaps$obsrvdAl==T | plotGaps$vslChck %in% c(1,4:9)])/length(plotGaps$area2)
    
    # Total area not observed
    length(plotGaps$area2[plotGaps$obsrvdAl==F])/length(plotGaps$area2)
    
    length(plotGaps$area2[plotGaps$obsrvdAl==T & plotGaps$obsrvdAr==T & plotGaps$obsrvdH==T])/length(plotGaps$area2)
    length(plotGaps$area2[plotGaps$obsrvdAl==T])/length(plotGaps$area2)
    
    # False positive area
    length(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(2,10)])/length(plotGaps$area2)
    # Area missing from monthly data
    length(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(1,7,9)])/length(plotGaps$area2)
    # Spatial or temporal mismatch
    length(plotGaps$area2[plotGaps$obsrvdAl==F & plotGaps$vslChck %in% c(5,6,8)])/length(plotGaps$area2)
    
#### Figure S4. Smoothing scale example plots ####

    # Read polygon buffer 25 m inland from lake
    buffer <- rgdal::readOGR("D:/BCI_Spatial/BCI_Outline_Minus25.shp")
    buffer <- sp::spTransform(buffer,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")  
    
    dem1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",2,".tif")) 
    dem1 <- raster::crop(dem1, raster::extent(c(626000,627000,1012500,1013500)))
    dem2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",8,".tif")) 
    dem2 <- raster::crop(dem2, raster::extent(c(626000,627000,1012500,1013500)))
    dem3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",24,".tif")) 
    dem3 <- raster::crop(dem3, raster::extent(c(626000,627000,1012500,1013500)))
    dem4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/DEM_smooth_",48,".tif")) 
    dem4 <- raster::crop(dem4, raster::extent(c(626000,627000,1012500,1013500)))
    
    curv1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",2,".tif")) 
    curv1 <- raster::crop(curv1, raster::extent(c(626000,627000,1012500,1013500)))
    curv2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",8,".tif")) 
    curv2 <- raster::crop(curv2, raster::extent(c(626000,627000,1012500,1013500)))
    curv3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",24,".tif")) 
    curv3 <- raster::crop(curv3, raster::extent(c(626000,627000,1012500,1013500)))
    curv4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Curv_smooth_",48,".tif")) 
    curv4 <- raster::crop(curv4, raster::extent(c(626000,627000,1012500,1013500)))
    
    slope1 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",2,".tif")) 
    slope1 <- raster::crop(slope1, raster::extent(c(626000,627000,1012500,1013500)))
    slope2 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",8,".tif")) 
    slope2 <- raster::crop(slope2, raster::extent(c(626000,627000,1012500,1013500)))
    slope3 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",24,".tif")) 
    slope3 <- raster::crop(slope3, raster::extent(c(626000,627000,1012500,1013500)))
    slope4 <- raster::raster(paste0("D:/BCI_Spatial/BCI_Topo/Slope_smooth_",48,".tif")) 
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
    mtext(expression(sigma~"= 24"),side=2,outer=F, cex=2)
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
    mtext(expression(sigma~"= 48"),side=2,outer=F, cex=2)
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
    curvRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Curv_smooth_8.tif")
    slopeRaster <- raster::raster("D:/BCI_Spatial/BCI_Topo/Slope_smooth_24.tif")
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
    par(mar=c(4,5,2,8))
    plot(meanQuad_DIC, breaks = plotBreaks,
         col = rev(wesanderson::wes_palette("Zissou1", length(plotBreaks), type = "continuous")),
         ylab = expression("Curvature scale ("~sigma~")"),
         xlab = expression("Slope scale ("~sigma~")"),
         main = expression(Delta~"DIC score"))
    
#### Figure S#: ####

# Look at initial canopy height and transitions per height class

# Load rasters    
d15to18 <- raster::raster("dCHM15to18_tin.tif")     
d18to20 <- raster::raster("dCHM18to20_tin.tif")

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

d15to18_vals <- raster::values(d15to18)
d18to20_vals <- raster::values(d18to20)

chm09_vals[is.na(d15to18_vals) | is.na(d18to20_vals)] <- NA
chm15_vals[is.na(d15to18_vals) | is.na(d18to20_vals)] <- NA
chm18_vals[is.na(d15to18_vals) | is.na(d18to20_vals)] <- NA
chm20_vals[is.na(d15to18_vals) | is.na(d18to20_vals)] <- NA

toChange <- which(!is.na(chm15_vals) & (chm15_vals-chm09_vals) > 3 & (chm18_vals-chm15_vals) < -1)
newVals <- chm09_vals[toChange]
#newVals <- 0.5*(chm09_vals[toChange]+chm18_vals[toChange])
chm15c <- chm15
raster::values(chm15c)[toChange] <- newVals
chm15_c_vals <-  raster::values(chm15c)
chm15_c_vals[is.na(d15to18_vals) | is.na(d18to20_vals)] <- NA
length(toChange)/length(chm15c_vals[!is.na(chm15c_vals)])



allStart <- floor(max(c(chm15_vals,chm18_vals),na.rm=T))

propGap <- data.frame(start = 5:allStart,
                      n15 = NA,
                      propGap15 = NA,
                      propNA15 = NA,
                      n18 = NA,
                      propGap18 = NA,
                      propNA18 = NA,
                      n20=NA,
                      n15_c = NA,
                      propGap15_c = NA,
                      propNA15_c = NA,
                      n09=NA)

for(i in 1:nrow(propGap)){
  
  vals15 <- which(!is.na(chm15_vals) & chm15_vals>=propGap$start[i] & chm15_vals<(propGap$start[i]+1))
  d15 <- chm18_vals[vals15] - chm15_vals[vals15]
  propGap$n15[i] <- length(vals15)/length(chm15_vals[!is.na(chm15_vals)])
  propGap$propGap15[i] <- length(d15[!is.na(d15) & d15<=-5])/length(vals15)
  propGap$propNA15[i] <- length(d15[is.na(d15)])/length(vals15)
  
  vals18 <- which(!is.na(chm18_vals) & chm18_vals>=propGap$start[i] & chm18_vals<(propGap$start[i]+1))
  d18 <- chm20_vals[vals18] - chm18_vals[vals18]
  propGap$n18[i] <- length(vals18)/length(chm18_vals[!is.na(chm18_vals)])
  propGap$propGap18[i] <- length(d18[!is.na(d18) & d18<=-5])/length(vals18)
  propGap$propNA18[i] <- length(d18[is.na(d18)])/length(vals18)
  
  vals09 <- which(!is.na(chm09_vals) & chm09_vals>=propGap$start[i] & chm09_vals<(propGap$start[i]+1))
  propGap$n09[i] <- length(vals09)/length(chm09_vals[!is.na(chm09_vals)])
  
  vals20 <- which(!is.na(chm20_vals) & chm20_vals>=propGap$start[i] & chm20_vals<(propGap$start[i]+1))
  propGap$n20[i] <- length(vals20)/length(chm20_vals[!is.na(chm20_vals)])
  
  vals15_c <- which(!is.na(chm15_c_vals) & chm15_c_vals>=propGap$start[i] & chm15_c_vals<(propGap$start[i]+1))
  d15_c <- chm18_vals[vals15_c] - chm15_c_vals[vals15_c]
  propGap$n15_c[i] <- length(vals15_c)/length(chm15_c_vals[!is.na(chm15_c_vals)])
  propGap$propGap15_c[i] <- length(d15_c[!is.na(d15_c) & d15_c<=-5])/length(vals15_c)
  propGap$propNA15_c[i] <- length(d15_c[is.na(d15_c)])/length(vals15_c)
}


# Normalize the proportion of gaps observed to per year
nYr15to18 <- as.numeric(as.Date("2018-06-07") - as.Date("2015-06-26"))/365
nYr18to20 <- as.numeric(as.Date("2020-07-31") - as.Date("2018-06-07"))/365

# Plot canopy height distributions
plot(n15~start,
     data = propGap,
     type="l",
     ylim=c(0,0.06),
     xlim=c(5,50),
     xlab = "Initial canopy height (m)",
     ylab = "Proportion of area",
     main = "Proportion of total area at height",
     col = adjustcolor("black",0.6),
     lty=3,
     lwd=3)

lines(n15_c~start,
      col = adjustcolor("black",0.6),
      data=propGap,lty=1,lwd=3)  

lines(n09~start,
      data = propGap,
      col=adjustcolor("orange",0.6), lwd=3)

lines(n18~start,
      data = propGap,
      col=adjustcolor("red",0.6), lwd=3)

lines(n20~start,
      data = propGap,
      col=adjustcolor("blue",0.6), lwd=3)

legend(x=35,y=0.06,
       bty="n",
       c("2009","2015 (original)","2015 (corrected)","2018","2020"),
       lty=c(1,3,1,1,1),
       col=adjustcolor(c("orange","black","black","red","blue"),0.6), lwd=3)



# Plot probability of becoming a gap for area in each height bin
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(propGap15/nYr15to18~start,
     data = propGap,
     type="l",
     ylim=c(0,0.06),
     xlim=c(5,50),
     xlab = "Initial canopy height (m)",
     ylab = "Proportion of area",
     main = "Proportion of area at height that decreases > 5 m",
     lwd=2)

lines(propGap18/nYr18to20~start, data=propGap,
      col="red", lwd=2)
lines(propGap15_c/nYr15to18~start, data=propGap,
      col="black",lty=2, lwd=2)

legend(x=25,y=0.06,
       bty="n",
       c("2015-2018 (original)","2015-2018 (corrected)","2018-2020"),
       col=c("black","black","red"), lwd=2, lty=c(1,2,1))







