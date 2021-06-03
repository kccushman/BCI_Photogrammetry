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







