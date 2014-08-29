library(raster)
library(rgdal)
library(maptools)
  data(wrld_simpl)
library(maps)
library(splancs)

##########################################################################################################
#### Breeding area polygon ###############################################################################
##########################################################################################################

## Equal area projection for snow cover data
load("area.laea.RData")

## lonlat projection
load("area.longlat.RData")



##########################################################################################################
#### Temperature data ####################################################################################
##########################################################################################################

load('climate.Taimyr.RData')

year.vector <- unique(daily$YEAR)
  daily$julian.date  <- as.POSIXlt(as.POSIXct(daily$DATE))$yday + 1
  daily.selected.range <- subset(daily, RANGE=="calcan")

D<- c()
par(mfrow = c(2, 2), mar = c(2,2,2,2), oma = c(5,5,0,0), las = 1, cex.axis = 1.4)
for (i in 1:length(year.vector)) {
  
  selected.dataset  =  daily.selected.range[daily.selected.range$YEAR==year.vector[i],]
  q2.model=lm(TEMPERATUR~julian.date+I(julian.date^2), selected.dataset)
  
  fnc <- approxfun(x = predict(q2.model)[1:which.max(predict(q2.model))], 
                   y = selected.dataset$julian.date[1:which.max(predict(q2.model))])
  D[i]=fnc(0)
if(i%in%c(30, 38, 40)){

with(selected.dataset, plot(julian.date, TEMPERATUR, pch = 16, col = "grey80", cex = 1.5, ylim = c(-22, 15), xlab = "", ylab = ""))
lines(selected.dataset$julian.date, predict(q2.model),
 	  col = "orange", lwd = 3)
points(fnc(0), 0, pch = 17, col = "cornflowerblue", cex = 3) 	  
text(220, -17, year.vector[i], cex = 1.4)
}
}


mtext("Temperature (degrees celsius)", 2, outer = T, las = 3, line = 2, cex = 2)
mtext("Day of the year", 1, outer = T, las = 1, line = 2, cex = 2)

plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("Raw data", "Model fit", "Start"),
       pch = c(16, NA, 17), col = c("grey80", "orange", "cornflowerblue"), lty = c(NA, 1, NA),
       bty = "n", cex = 2)




temp <- as.data.frame(cbind(year.vector, D))
names(temp) <- c('year','temp.start')
save(temp, file = "temp.RData")
load("temp.RData")

##########################################################################################################
#### NDVI data ###########################################################################################
##########################################################################################################
## data downloaded from: ftp://ftp.orbit.nesdis.noaa.gov/pub/corp/scsb/wguo/GVIx/GVIx_VH_16km/NVI/       #
##########################################################################################################

setwd("/Volumes/SLISOVSK/Data/NDVI")

fls <- list.files(pattern="*.hdf")
fls.year <- as.numeric(sapply(fls,substring,22,25))
fls.week <- as.numeric(sapply(fls,substring,28,29))
fls <- fls[!duplicated(paste(fls.year,fls.week,sep="_"))]
fls.year <- as.numeric(sapply(fls,substring,22,25))
fls.week <- as.numeric(sapply(fls,substring,28,29))

tmp <- suppressWarnings(readGDAL(fls[25]))

r   <- raster(tmp)
extent(r) <- extent(-180, 180, -56.1, 75.324)  

ind <- unlist(extract(r, area.longlat))
ndvi <- matrix(nrow = length(fls), ncol = length(ind))

for(i in 1547:length(fls)) {
  tmp <- suppressWarnings(readGDAL(fls[i]))
  r   <- raster(tmp)
  extent(r) <- extent(-180, 180, -56.1, 75.324)  
  ndvi[i,] <- unlist(extract(r, area.longlat))
}
setwd("~/Dropbox/Science/KnotNDVI/RS_KnotNDVI")


tmp <- suppressWarnings(readGDAL("~/Desktop/newNDVI/GVIX_NN_G16_C07_NVI_Y2013_P52.hdf"))


ndvi <- ifelse(ndvi<(-50), NA, ndvi)
  ind1 <- paste(fls.year, fls.week, sep="_") 

tab <- matrix(c(rep(c(1983:2013), each = 52), rep(c(1:52), length(c(1983:2013)))), ncol = 2)
tab <- cbind(tab,ndvi[match(paste(tab[,1], tab[,2], sep="_"), ind1),])

out <- data.frame(year = 1983:2013, onset = NA, end = NA, duration = NA, max = NA, area = NA)

opar <- par(mfrow = c(7, 5), mar = c(0.5,0.5,0.5,0.5), oma = c(6,6,4,1))
for(i in 1:length(1983:2013)){
  
  tmp  <- tab[tab[,1]==c(1983:2013)[i],]  
  y    <- matrix(c(rep(tmp[,1], ncol(tmp)-3), rep(tmp[,2], ncol(tmp)-3), rep(tmp[,3], ncol(tmp)-3), as.vector(tmp[,-c(1:3)])), ncol = 4) 
  
  if(sum(!is.na(y[y[,2]%in%c(20:40),4]))>100) {
    plot(y[,2], y[,4], pch=16, col="grey70", ylim = c(-80, 500), xaxt="n", yaxt = "n")
    axis(1, tcl = 0.2, labels = ifelse(i%in%c(26, 28), TRUE, FALSE), lwd = 1.1)
    axis(2, tcl = 0.2, labels = ifelse(i%in%c(1, 11, 21), TRUE, FALSE), lwd = 1.1, las = 2)
    box(lwd = 1.1)
    ls <- loess(y[,4]~y[,2], span = 0.3)
    pr <- predict(ls, newdata = c(1:52))
    lines(1:52, pr, col = "orange", lwd = 5)
    
    
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
    if(is.numeric(ct1)) abline(v=ct1, lty = 2, lwd = 2, col = "grey20")
    
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
    if(is.numeric(ct2)) abline(v=ct2, lty = 2, lwd = 2, col = "grey20")    
    
    ndvi2 <- y[,4]
    if(is.numeric(ct1)) ndvi2[which(y[,2]<ct1)] <- NA
    if(is.numeric(ct2)) ndvi2[which(y[,2]>ct2)] <- NA    
    
    range <- c(ifelse(is.numeric(ct1), ct1, min(tmp[!is.na(tmp[,4]),2])),
               ifelse(is.numeric(ct2), ct2, max(tmp[!is.na(tmp[,4]),2])))
    
    
    # points(y[,2], ndvi2, pch = 16, col = "orange")
    
    pr2 <- predict(ls, newdata = range[1]:range[2])
    lines(range[1]:range[2], pr2, lwd=4)
    
    ## Data
    ind3 <- which.max(pr2)
    out$max[i] <- pr2[ind3]
    
    ## Onset
    onc <- approxfun(x = pr2[1:ind3], y = range[1]:(range[1] + (ind3-1)))
    out$onset[i] <- onc(0) 
    points(onc(0), 0, pch = 16, cex = 3, col = "firebrick")
    
    ## end
    end <- approxfun(x = pr2[ind3:length(pr2)], y = c(range[1]:range[2])[ind3:length(pr2)])
    out$end[i] <- end(100)
    abline(h = c(0, 100), lty = 3, col = "grey40", lwd = 2)

    
    if(any(is.na(c(onc(0), end(100))))==FALSE){
      ## duration
      out$duration[i] <- out$end[i] - out$onset[i]
      
      
      ## area
      curve <- cbind(c(seq(onc(0), end(100), length = 30), onc(0)), predict(ls, newdata = c(seq(onc(0), end(100), length = 30), onc(0))))
      polygon(curve[,1], curve[,2], border = "transparent", col = rgb(0.4, 0.4, 0.2, alpha = 0.4))
      out$area[i] <- areapl(curve)
      points(end(100), 100, pch = 16, cex = 3, col = "cornflowerblue")
      points(onc(0), 0, pch = 16, cex = 3, col = "firebrick") 
    }
    mtext(c(1983:2011)[i], 3, line = -1.5, cex =1, at = 9)
  } else  {
    plot(NA, ylim = c(-75, 500), xaxt="n", yaxt = "n")
    mtext(c(1982:2011)[i], 3, line = -1.5, cex =1, at = 9)
  }
}
plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("Raw data", "Smooth line", "Smooth sensible \nperiod", "Borders sensible \nperiod", "Phenology \nthresholds", "Start", "End"),
       pch = c(16, NA, NA, NA, NA, 16, 16), lty = c(NA,1,1,2,3,NA,NA), col = c("grey50", "orange", "black", "black", "black", "firebrick", "cornflowerblue"), 
       bty = "n", cex = 0.6, y.intersp  = 1.4)
mtext("Date (week of the year)", side = 1, outer = T, line = 3.5, cex = 1.3)
mtext("Normalized Vegetation Index (NDVI)", side = 2, outer = T, line = 3.5, cex = 1.3)
mtext("NDVI Data Manipulation (smoothing and phenology extraction)", 3, line = 1.5, outer = T, cex = 1.5)
par(opar)

ndvi <- out
save(ndvi, file = "ndvi.RData")


##########################################################################################################
#### Snow Data ###########################################################################################
##########################################################################################################
## data downloaded from: ftp://sidads.colorado.edu/pub/DATASETS/nsidc0046_weekly_snow_seaice/            #
##########################################################################################################

### maximum likelihood function
gaussMLE <- function(week, size, prob) {
  require(bbmle)
  
  tab <- data.frame(week = week, size = size, p = prob)
  tab <- merge(data.frame(week = 1:52), tab, all.x = T)
  
  gauss.curve <- function(parms, intv = 0.1) {
    t <- seq(1, 52, intv)
    parms <- as.list(parms)
    fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
    fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
    c(fit1, fit2[-1])
  }
  
  gauss.loglik <- function(a1, a2, a3, a4, a5) {
    fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), 1)  
    fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
    # cat(paste(c(a1, a2, a3, a4, a5), sep = "  "), "\r")
    -sum(dbinom(x = round(tab[,3]*tab[,2],0), size = rep(100, length(fit)), prob = fit, log=TRUE), na.rm=T)
  }
  
  mle <- suppressWarnings(mle2(gauss.loglik, method="L-BFGS-B",
                               start=list(a1 = 33, a2 = 4,  a3 = 9,  a4 = 4, a5 = 9),
                               lower=list(a1 = 2,  a2 = 1.5,  a3 = 0.5,  a4 = 1.5,  a5 = 0.5),
                               upper=list(a1 = 50, a2 =  Inf,  a3 =  Inf, a4 =  Inf, a5 =  Inf),
  ))
  
  t <- seq(1, 52, 0.01)
  fit <- gauss.curve(coef(mle), intv = 0.01)  
  
  start <- approxfun(x = fit[1:which(t==round(coef(mle)[1],0))], y = t[1:which(t == round(coef(mle)[1]))])
  end   <- approxfun(x = fit[which(t==round(coef(mle)[1],0)):length(t)], y = t[which(t == round(coef(mle)[1],0)):length(t)])
  
  list(result = data.frame(week = t, fit = fit), start = start, end = end)
} 
## end maximum likelihood function

lon <-  raster("snow/longitude.tif")
lat <-  raster("snow/latitude.tif")
snow <- brick("snow/snow_cover_extent.tif")

  start <- as.POSIXct("1966-10-04", "GMT")
  end <- as.POSIXct("2014-02-03", "GMT")
  date <- seq(start, end, by = "7 days")
  week <- as.numeric(format(date, "%U"))
  year <- as.numeric(format(date, "%Y"))
  
  # Projection
  projSnow <- "+proj=stere +lat_0=90 +lon_0=10"
  xy <- project(cbind(values(lon), values(lat)), projSnow)
  area <- spTransform(area.longlat, CRS(projSnow))

  # Rasterize ([[1]])
  r <- raster(xmn = min(xy[,1]), xmx = max(xy[,1]), ymn = min(xy[,2]), ymx = max(xy[,2]), nrow=88, ncol=88)
  tmp <- rasterize(xy, y = tmp1, field = values(snow[[17]]))
  tmp2 <- resample(tmp, 
                   raster(raster(xmn = min(xy[,1]), xmx = max(xy[,1]), ymn = min(xy[,2]), ymx = max(xy[,2]), nrow=88*3, ncol=88*3)),
                   method = "ngb")


out <- matrix(ncol = 53, nrow = length(unique(year)))

for(i in which(year>1980)){

  cat(round(i/length(names(snow)),4), "%   ", "\r")
  tmp1 <- rasterize(xy, y = r, field = values(snow[[i]]))
  tmp2 <- crop(tmp1, area) 
  tmp3 <- resample(tmp2, 
                   raster(extent(tmp2), ncol = 6, nrow = 12),
                   method = "ngb")
  out[which(unique(year)==year[i]), which(c(0:52)==week[i])] <- extract(tmp3, area, mean, na.rm=T)

}

save(out, file="out.RData")
load("out.RData")

snowP <- data.frame(year = rep(min(unique(year)):max(unique(year)), 53), 
                    week = rep(0:52, each = length(unique(year))),
                    p =    as.vector(out))

phen <- matrix(nrow = length(1983:2013), ncol = 2) 

opar <- par(mfrow = c(7, 5), mar = c(0.5,0.5,0.5,0.5), oma = c(6,6,4,1))
for(i in 1983:2013) {
  tmp1 <- subset(snowP, year == i & week>0)
  plot(tmp1$week, tmp1$p, type = "o", pch = 16, col = "grey50", lwd = 2, cex = 0.5, ylim = c(-0.1, 1.1), xlim = c(-1,52),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(1, tcl = 0.2, labels = ifelse(i%in%c(2008, 2011), TRUE, FALSE), lwd = 1.1)
  axis(2, tcl = 0.2, labels = ifelse(i%in%c(1983, 1993, 2003), TRUE, FALSE), lwd = 1.1, las = 2)
  box(lwd = 1.1)
  # points(tmp1$week, tmp1$p2, pch = 16, type = "b", cex = 0.5)
  
  mod <- gaussMLE(week = tmp1$week, size = rep(100, nrow(tmp1)), prob = tmp1$p)
  lines(mod$result, lwd = 2, col = "orange")
  
  points(mod$start(0.66), 0.66, pch = 18, cex = 2, col = "firebrick")
  points(mod$end(0.33), 0.33, pch = 18, cex = 2, col = "cornflowerblue")
  text(8, 0.05, i, font = 2, cex = 1.4)
  phen[which(c(1983:2013)==i), 1] <- mod$start(0.66)
  phen[which(c(1983:2013)==i), 2] <- mod$end(0.33)
}
plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("Raw data", "Filtered data", "Model fit", "Start interp.", "End. interp."),
       pch = c(16, 16, NA, 18, 18), lty = c(1,1,1,NA,NA), col = c("grey50", "black", "orange", "firebrick", "cornflowerblue"), 
       bty = "n", cex = 1.3)
mtext("Date (week of the year)", side = 1, outer = T, line = 3.5, cex = 1.3)
mtext("Percentage of snow cover at breeding area", side = 2, outer = T, line = 3.5, cex = 1.3)
mtext("Snow Cover Data Manipulation (Curve fitting)", 3, line = 1.5, outer = T, cex = 1.5)
par(opar)


day <- function(year, week){
  tm <- seq(as.POSIXct(paste(year, "-01-01", sep = ""), "GMT"), as.POSIXct(paste(year, "-12-31", sep = ""), "GMT"),
            by = "day")
  wk <- seq(1, max(as.numeric(format(tm, "%U")))+1, length = length(tm))
  f <- approxfun(x = wk, y = tm)
  return(f(week))
}

snow <- data.frame(year = c(1983:2013))
snow$start <- round(as.POSIXct(apply(cbind(as.numeric(1983:2013), phen[,1]), 1, function(x) day(x[1],x[2])), origin = "1970-01-01"), "day")
snow$end <- round(as.POSIXct(apply(cbind(as.numeric(1983:2013), phen[,2]), 1, function(x) day(x[1],x[2])), origin = "1970-01-01"), "day") 
snow$yday.start <- as.POSIXlt(snow$start, "GMT")$yday + 1
snow$yday.end   <- as.POSIXlt(snow$end, "GMT")$yday + 1

save(snow, file = "snow.RData")
load("snow.RData")

##########################################################################################################
#### Red Knot Juvenile Bill Length #######################################################################
##########################################################################################################
load("data.Poland.RData")

  ## Bill data
  data.Poland$bill.all[is.na(data.Poland$BILL)]=data.Poland$BILL1[is.na(data.Poland$BILL)]
  data.Poland$bill.all[is.na(data.Poland$BILL)==F]=data.Poland$BILL[is.na(data.Poland$BILL)==F]
  data.Poland <- subset(data.Poland, AGE=="JUV")

n <- aggregate(data.frame(year = data.Poland$YEAR, n = rep(1, nrow(data.Poland))), by = list(data.Poland$YEAR), sum)[,c(1,3)]
b <- aggregate(data.frame(year = data.Poland$YEAR, b = data.Poland$bill.all), by = list(data.Poland$YEAR), mean, na.rm = T)[,c(1,3)]

bill <- merge(n, b, by = "Group.1", all.x = T, all.y = T)
names(bill) <- c("year", "n", "bill")

# only data with sample size >= 10
bill$bill <- ifelse(bill$n<10, NA, bill$bill)

## all bill lengths
bill.all <- data.Poland
bill.all <- data.frame(year = data.Poland$YEAR, bill=data.Poland$bill.all)

save(bill.all, file = "bill.all.RData")
save(bill, file = "bill.RData")


##########################################################################################################
#### Merge Data ##########################################################################################
##########################################################################################################
load("temp.RData")
  names(temp) <- c("year", "temp.start")
load("ndvi.RData")
  names(ndvi) <- c("year", "ndvi.start", "ndvi.end", "ndvi.duration", "ndvi.max", "ndvi.area")
load("snow.RData")
  names(snow) <- c("year", "snow.start.date", "snow.end.date", "snow.start.jD", "snow.end.jD")
load("bill.RData")

load("bill.all.RData")

knot.data <- merge(bill, temp, by = "year", all.x = T)
  knot.data <- merge(knot.data, snow, by = "year", all.x = T)
  knot.data <- merge(knot.data, ndvi, by = "year", all.x = T)

save(knot.data, file = "knot.data.RData")


knot.data.all <- merge(bill.all, temp, by = "year", all.x = T)
knot.data.all <- merge(knot.data.all, snow, by = "year", all.x = )
knot.data.all <- merge(knot.data.all, ndvi, by = "year", all.x = T)

save(knot.data.all, file = "knot.data.all.RData")
