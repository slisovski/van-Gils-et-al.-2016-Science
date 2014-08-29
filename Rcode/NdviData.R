##########################################################################################################
#### NDVI data ###########################################################################################
##########################################################################################################

library(raster)
library(rgdal)
library(maptools)
  data(wrld_simpl)
library(maps)
library(splancs)


# set working directory to folder containing the ndvi .hdf files downloaded
# from: ftp://ftp.orbit.nesdis.noaa.gov/pub/corp/scsb/wguo/GVIx/GVIx_VH_16km/NVI/ 

# requires rgdal library with hdf4 driver

# Step 1: get file paths and date index

fls <- list.files(pattern="*.hdf")
fls.year <- as.numeric(sapply(fls,substring,22,25))
fls.week <- as.numeric(sapply(fls,substring,28,29))

# Step 2: index raster with extent of ndvi data

tmp <- suppressWarnings(readGDAL(fls[25]))
r   <- raster(tmp)
extent(r) <- extent(-180, 180, -56.1, 75.324)  


# Step 3: Breeding Area Polygon ()

load("~/area.longlat.RData")

# Step 4: Create matrix with nrow = length of ndvi files and ncol = number of cells within breeding area

ind <- unlist(extract(r, area.longlat))
ndvi <- matrix(nrow = length(fls), ncol = length(ind))

# Step 5: Extract values for breeding area and fill ndvi matrix

for(i in 1547:length(fls)) {
  tmp <- suppressWarnings(readGDAL(fls[i]))
  r   <- raster(tmp)
  extent(r) <- extent(-180, 180, -56.1, 75.324)  
  ndvi[i,] <- unlist(extract(r, area.longlat))
}

ndvi <- ifelse(ndvi<(-50), NA, ndvi) # remove values below -50 (water and snow)


### Calculate phenological values as described in vanGils et al. XXX

ind1 <- paste(fls.year, fls.week, sep="_") 

tab <- matrix(c(rep(c(1983:2013), each = 52), rep(c(1:52), length(c(1983:2013)))), ncol = 2)
tab <- cbind(tab,ndvi[match(paste(tab[,1], tab[,2], sep="_"), ind1),])

out <- data.frame(year = 1983:2013, onset = NA, end = NA, duration = NA, max = NA, area = NA)


# opar <- par(mfrow = c(7, 5), mar = c(0.5,0.5,0.5,0.5), oma = c(6,6,4,1))
for(i in 1:length(1983:2013)){
  
  tmp  <- tab[tab[,1]==c(1983:2013)[i],]  
  y    <- matrix(c(rep(tmp[,1], ncol(tmp)-3), rep(tmp[,2], ncol(tmp)-3), rep(tmp[,3], ncol(tmp)-3), as.vector(tmp[,-c(1:3)])), ncol = 4) 
  
  if(sum(!is.na(y[y[,2]%in%c(20:40),4]))>100) {
    # plot(y[,2], y[,4], pch=16, col="grey70", ylim = c(-80, 500), xaxt="n", yaxt = "n")
    # axis(1, tcl = 0.2, labels = ifelse(i%in%c(26, 28), TRUE, FALSE), lwd = 1.1)
    # axis(2, tcl = 0.2, labels = ifelse(i%in%c(1, 11, 21), TRUE, FALSE), lwd = 1.1, las = 2)
    # box(lwd = 1.1)
    
    ls <- loess(y[,4]~y[,2], span = 0.3)
    pr <- predict(ls, newdata = c(1:52))
    
    # lines(1:52, pr, col = "orange", lwd = 5)
    
    
    tmp2 <- cbind(c(1:52), pr)
    nr <- 2
    ct1 <- FALSE
    repeat{
      if(!is.na(tmp2[nr,2]) & !is.na(tmp2[(nr-1),2])){
        if(tmp2[nr,2]<tmp2[(nr-1),2] & tmp2[nr,2]<tmp2[(nr+1),2] &
             mean(tmp2[(nr-(ifelse(nr>5, 5, nr))):nr,2], na.rm = T)<mean(tmp2[nr:(nr+5),2], na.rm = T)) {
          ct1 <- nr
          break
        }
      }
      nr <- nr + 1
    }
    
    # if(is.numeric(ct1)) abline(v=ct1, lty = 2, lwd = 2, col = "grey20")
    
    nr <- 51
    ct2 <- FALSE
    repeat{
      if(is.numeric(ct1) & nr == ct1) break
      if(!is.na(tmp2[nr,2]) & !is.na(tmp2[(nr+1),2])){
        if(tmp2[nr,2]<tmp2[(nr-1),2] & tmp2[nr,2]<tmp2[(nr+1),2] &
             length(which(tmp2[nr:52,2]>c(NA,tmp2[nr:52,2][-length(tmp2[nr:52,2])]) &
                            tmp2[nr:52,2]>c(tmp2[nr:52,2][-1],NA)))==0) {
          ct2 <- nr
          break
        }
      }
      nr <- nr - 1
    }    
    
    # if(is.numeric(ct2)) abline(v=ct2, lty = 2, lwd = 2, col = "grey20")    
    
    ndvi2 <- y[,4]
    if(is.numeric(ct1)) ndvi2[which(y[,2]<ct1)] <- NA
    if(is.numeric(ct2)) ndvi2[which(y[,2]>ct2)] <- NA    
    
    range <- c(ifelse(is.numeric(ct1), ct1, min(tmp[!is.na(tmp[,4]),2])),
               ifelse(is.numeric(ct2), ct2, max(tmp[!is.na(tmp[,4]),2])))
    
    
    # points(y[,2], ndvi2, pch = 16, col = "orange")
    
    pr2 <- predict(ls, newdata = range[1]:range[2])
    # lines(range[1]:range[2], pr2, lwd=4)
    
    ## Data
    ind3 <- which.max(pr2)
    out$max[i] <- pr2[ind3]
    
    ## Onset
    onc <- approxfun(x = pr2[1:ind3], y = range[1]:(range[1] + (ind3-1)))
    out$onset[i] <- onc(0) 
    # points(onc(0), 0, pch = 16, cex = 3, col = "firebrick")
    
    ## end
    end <- approxfun(x = pr2[ind3:length(pr2)], y = c(range[1]:range[2])[ind3:length(pr2)])
    out$end[i] <- end(100)
    # abline(h = c(0, 100), lty = 3, col = "grey40", lwd = 2)

    
    if(any(is.na(c(onc(0), end(100))))==FALSE){
      ## duration
      out$duration[i] <- out$end[i] - out$onset[i]
      
      
      ## area
      curve <- cbind(c(seq(onc(0), end(100), length = 30), onc(0)), predict(ls, newdata = c(seq(onc(0), end(100), length = 30), onc(0))))
      # polygon(curve[,1], curve[,2], border = "transparent", col = rgb(0.4, 0.4, 0.2, alpha = 0.4))
      out$area[i] <- areapl(curve)
      # points(end(100), 100, pch = 16, cex = 3, col = "cornflowerblue")
      # points(onc(0), 0, pch = 16, cex = 3, col = "firebrick") 
    }
    # mtext(c(1983:2011)[i], 3, line = -1.5, cex =1, at = 9)
  } else  {
    # plot(NA, ylim = c(-75, 500), xaxt="n", yaxt = "n")
    # mtext(c(1982:2011)[i], 3, line = -1.5, cex =1, at = 9)
  }
}
# plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# legend("topleft", c("Raw data", "Smooth line", "Smooth sensible \nperiod", "Borders sensible \nperiod", "Phenology \nthresholds", "Start", "End"),
#        pch = c(16, NA, NA, NA, NA, 16, 16), lty = c(NA,1,1,2,3,NA,NA), col = c("grey50", "orange", "black", "black", "black", "firebrick", "cornflowerblue"), 
#        bty = "n", cex = 0.6, y.intersp  = 1.4)
# mtext("Date (week of the year)", side = 1, outer = T, line = 3.5, cex = 1.3)
# mtext("Normalized Vegetation Index (NDVI)", side = 2, outer = T, line = 3.5, cex = 1.3)
# mtext("NDVI Data Manipulation (smoothing and phenology extraction)", 3, line = 1.5, outer = T, cex = 1.5)
# par(opar)

output <- out