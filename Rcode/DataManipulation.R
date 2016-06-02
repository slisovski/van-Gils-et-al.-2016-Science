library(raster)
library(rgeos)
library(gdalUtils)
library(rgdal)
library(RNetCDF)
library(maptools)
  data(wrld_simpl)
library(maps)
library(splancs)
library(bbmle)


## wd
setwd("~/Dropbox/Science/Projects/RedKnot_BillLength/Analysis")


##########################################################################################################
#### Breeding area polygon ###############################################################################
##########################################################################################################

## lonlat projection
load("Data/Basic/area.longlat.RData")

## stereographic projection
areaSP <- spTransform(area.longlat, CRSobj = CRS("+proj=stere +lat_0=90, +lon_0=0"))

## azimuth projection
areaAP <- spTransform(area.longlat, CRSobj = CRS("+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84"))



##########################################################################################################
#### Temperature data ####################################################################################
##########################################################################################################

load('Data/Basic/temp.RData')

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
# save(temp, file = "Data//Results/temp.phen.RData")
load("Data/Results/temp.phen.RData")

##########################################################################################################
#### NDVI data ###########################################################################################
##########################################################################################################
## data downloaded from: ftp://ftp.orbit.nesdis.noaa.gov/pub/corp/scsb/wguo/data/VHP_16km/VH/            #
##########################################################################################################

setwd("~/Desktop/VH")

fls <- list.files(pattern="*ND.hdf")
fls.year <- as.numeric(sapply(fls,substring,17,20))
fls.week <- as.numeric(sapply(fls,substring,21,23))
fls <- fls[!duplicated(paste(fls.year,fls.week,sep="_"))]
fls.year <- as.numeric(sapply(fls,substring,17,20))
fls.week <- as.numeric(sapply(fls,substring,21,23))

r   <- gdal_translate(fls[10], "tmp.tiff", sd_index = 2, output_Raster = TRUE)[[1]]
extent(r) <- extent(-180, 180, -56.1, 75.324)  

ind <- unlist(extract(r, area.longlat))
ndvi <- matrix(nrow = length(fls), ncol = length(ind))

for(i in 1:length(fls)) {
  if(i/100 == floor(i/100))  cat(i, " of ", length(fls), "\n")
  r <- gdal_translate(fls[i], "tmp.tiff", sd_index = 2, output_Raster = TRUE)[[1]]
  extent(r) <- extent(-180, 180, -56.1, 75.324)  
  ndvi[i,] <- unlist(extract(r, area.longlat))
}

setwd("~/Dropbox/Science/Projects/RedKnot_BillLength/Analysis")
# save(ndvi, file = "Data/Basic/ndvi_raw.RData")
load("Data/Basic/ndvi_raw.RData")


ndvi <- ifelse(ndvi<(-0.15), NA, ndvi)
ind1 <- paste(fls.year, fls.week, sep="_") 

tab <- matrix(c(rep(c(1983:2015), each = 52), rep(c(1:52), length(c(1983:2015)))), ncol = 2)
tab <- cbind(tab,ndvi[match(paste(tab[,1], tab[,2], sep="_"), ind1),])

out <- data.frame(year = 1983:2015, onset = NA, end = NA, duration = NA, max = NA, area = NA)

opar <- par(mfrow = c(7, 5), mar = c(0.5,0.5,0.5,0.5), oma = c(6,6,4,1))
for(i in 1:length(1983:2015)){
  
  tmp  <- tab[tab[,1]==c(1983:2015)[i],]  
  y    <- matrix(c(rep(tmp[,1], ncol(tmp)-3), rep(tmp[,2], ncol(tmp)-3), rep(tmp[,3], ncol(tmp)-3), as.vector(tmp[,-c(1:3)])), ncol = 4) 
  plot(y[,2], y[,4], pch=16, col="grey70", ylim = c(-0.1, 0.75), xaxt="n", yaxt = "n")
  
  # if((sum(!is.na(y[y[,2]%in%c(20:40),4]))/ncol(ndvi))>19) {
  axis(1, tcl = 0.2, labels = ifelse(i%in%c(31, 33), TRUE, FALSE), lwd = 1.1)
  axis(2, tcl = 0.2, labels = ifelse(i%in%c(1, 11, 21, 31), TRUE, FALSE), lwd = 1.1, las = 2)
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
  abline(h = out$max[i])
  abline(h=old$Amp[old$Year==c(1983:2015)[i]], col = "blue")
  ## Onset
  onc <- approxfun(x = pr2[1:ind3], y = range[1]:(range[1] + (ind3-1)))
  # if(c(1983:2015)[i]==2009) {
  
  out$onset[i] <- onc(0)
  
  ## end
  end <- approxfun(x = pr2[ind3:length(pr2)], y = c(range[1]:range[2])[ind3:length(pr2)])
  if(c(1983:2015)[i]==2009) out$end[i] <- end(0.13) else out$end[i] <- end(0.1)
  abline(h = c(0, 0.1), lty = 3, col = "grey40", lwd = 2)
  
  if(!any(is.na(c(out$onset[i], out$end[i])))){
    ## duration
    out$duration[i] <- out$end[i] - out$onset[i]
    
    ## area
    curve <- cbind(c(seq(out$onset[i], out$end[i], length = 30), out$onset[i]), predict(ls, newdata = c(seq(out$onset[i], out$end[i], length = 30), onc(0))))
    polygon(curve[,1], curve[,2], border = "transparent", col = rgb(0.4, 0.4, 0.2, alpha = 0.4))
    out$area[i] <- areapl(curve)
    points(end(0.1), 0.1, pch = 16, cex = 3, col = "cornflowerblue")
    points(onc(0), 0, pch = 16, cex = 3, col = "firebrick") 
  }
  mtext(c(1983:2015)[i], 3, line = -1.5, cex =1, at = 9)
  
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
# save(ndvi, file = "Data/Results/ndvi.phen.RData")
load("Data/Results/ndvi.phen.RData")



##########################################################################################################
#### Snow Data ###########################################################################################
##########################################################################################################
## data downloaded from: http://www.ncdc.noaa.gov/thredds/catalog/cdr/snowcover/                         #
#                        catalog.html?dataset=cdr/snowcover/nhsce_v01r01_19661004_20160104.nc            #
##########################################################################################################

### start: maximum likelihood function=
gaussMLE <- function(day, size, prob) {
  
  tab <- data.frame(day = day-6, size = size, p = prob)
  # tab <- merge(data.frame(week = 1:52), tab, all.x = T)
  
  gauss.curve <- function(parms, intv = 1) {
    t <- seq(1, 366, intv)
    parms <- as.list(parms)
    fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
    fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
    c(fit1, fit2[-1])
  }
  
  gauss.loglik <- function(a1, a2, a3, a4, a5) {
    fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), 1)  
    fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
    # cat(paste(c(a1, a2, a3, a4, a5), sep = "  "), "\r")
    -sum(dbinom(x = round(tab[,3]*tab[,2],0), size = rep(100, length(fit)), prob = fit[day], log=TRUE), na.rm=T)
  }
  
  mle <- suppressWarnings(mle2(gauss.loglik, method="L-BFGS-B",
                               start=list(a1 = 225, a2 = 40,  a3 = 9,  a4 = 40, a5 = 9),
                               lower=list(a1 = 50,  a2 = 5,   a3 = 0.5,a4=5,  a5 = 0.5),
                               upper=list(a1 = 225, a2 =  Inf,  a3 =  Inf, a4 =  Inf, a5 =  Inf),
  ))
  
  t <- seq(1, 366, 1)
  fit <- gauss.curve(coef(mle), intv = 1)  
  
  start <- approxfun(x = fit[1:which(t==round(coef(mle)[1],0))], y = t[1:which(t == round(coef(mle)[1]))])
  end   <- approxfun(x = fit[which(t==round(coef(mle)[1],0)):length(t)], y = t[which(t == round(coef(mle)[1],0)):length(t)])
  
  list(result = data.frame(week = t, fit = fit), start = start, end = end)
} 
## end maximum likelihood function




### Read netCDF file
fid <- open.nc("Data/Basic/nhsce_v01r01_19661004_20160104.nc")
print.nc(fid)

dat<-read.nc(fid)
## "days since 1966-10-03"
t<-dat$time
z<-dat$snow_cover_extent
ylat<-dat$longitude
xlon<-dat$latitude
close.nc(fid)

### Projection and raster extent
lon <- raster(ylat)
lat <- raster(xlon)
snow <- brick(z)

start <- as.POSIXct("1966-10-03", "GMT")
  date <- start + ((t-6)*24*60*60)

# Projection
projSnow <- "+proj=stere +lat_0=90 +lon_0=10"
xy <- project(cbind(values(lon), values(lat)), projSnow)
area <- spTransform(area.longlat, CRS(projSnow))

# Rasterize
r <- raster(xmn = min(xy[,1]), xmx = max(xy[,1]), ymn = min(xy[,2]), ymx = max(xy[,2]), nrow=88, ncol=88)
tmp <- rasterize(xy, y = r, field = values(snow[[1]]))

# plot(tmp1)
# wrld <- spTransform(subset(wrld_simpl, NAME != "Antarctica"), CRS(projSnow))
# plot(as(wrld, "SpatialLines"), add = TRUE, col = "green")


out <- data.frame(Date = date, Z = NA)

for(i in 1982:2015) {
  
  cat(round(i/length(names(snow)),4), "%   ", "\r")
  
  for(w in 1:length(date[as.numeric(format(date, "%Y"))==i])) {
    
    ind <- which(as.numeric(format(date, "%Y"))==i)[w]
    
  tmp1 <- rasterize(xy, y = tmp, field = values(snow[[ind]]))
  tmp2 <- crop(tmp1, area) 
  tmp3 <- resample(tmp2, 
                   raster(extent(tmp2), ncol = 6, nrow = 12),
                   method = "ngb")
  out[ind,2] <- extract(tmp3, area, mean, na.rm=T)
  
  }
}

# save(out, file = "Data/Basic/snow.RData")
load("Data/Basic/snow.RData")


phen <- matrix(nrow = length(1983:2015), ncol = 2) 

opar <- par(mfrow = c(7, 5), mar = c(0.5,0.5,0.5,0.5), oma = c(6,6,4,1))
for(i in 1983:2015) {
  tmp1 <- out[which(as.numeric(format(out[,1], "%Y"))==i),]
  
  p2 <- tmp1[,2]
  for(j in 2:(length(p2-1))) {
    if(!is.na(tmp1[j+1,2]) & tmp1[j-1,2]==tmp1[j+1,2] & tmp1[j,2]!=tmp1[j-1,2] & abs(tmp1[j-1,2] - tmp1[j,2])>0.25) p2[j] <- tmp1[j-1,2]
  } 
  
  days <- as.numeric(format(tmp1[,1], "%j"))
    
  plot(days, tmp1[,2], type = "o", pch = 16, col = "grey50", lwd = 3, cex = 1.1, ylim = c(-0.1, 1.1),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(1, tcl = 0.2, labels = ifelse(i%in%c(2008, 2011), TRUE, FALSE), lwd = 1.1)
  axis(2, tcl = 0.2, labels = ifelse(i%in%c(1983, 1993, 2003), TRUE, FALSE), lwd = 1.1, las = 2)
  box(lwd = 1.1)
  points(days, p2, pch = 16, type = "o", cex = 0.5, col = 1)
  
  mod <- gaussMLE(day = days, size = rep(100, nrow(tmp1)), prob = p2)
  lines(mod$result, lwd = 2, col = "orange")
  
  points(mod$start(0.75), 0.75, pch = 18, cex = 2, col = "firebrick")
  points(mod$end(0.25), 0.25, pch = 18, cex = 2, col = "cornflowerblue")
  text(50, 0.05, i, font = 2, cex = 1.4)
  
  phen[which(c(1983:2015)==i), 1] <- mod$start(0.75)
  phen[which(c(1983:2015)==i), 2] <- mod$end(0.25)
}
plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("Raw data", "Model fit", "Start interp.", "End. interp."),
pch = c(16, NA, 18, 18), lty = c(1,1,NA,NA), col = c("grey50", "orange", "firebrick", "cornflowerblue"), 
bty = "n", cex = 1.3)
mtext("Date (day of the year)", side = 1, outer = T, line = 3.5, cex = 1.3)
mtext("Percentage of snow cover at breeding area", side = 2, outer = T, line = 3.5, cex = 1.3)
mtext("Snow Cover Data Manipulation (Curve fitting)", 3, line = 1.5, outer = T, cex = 1.5)
par(opar)

snow.phen       <- data.frame(year = c(1983:2015))
snow.phen$start <- format(as.POSIXct(paste0(snow.phen$year, "-01-01"))+(phen[,1]-1)*24*60*60, "%b-%d")
snow.phen$end   <- format(as.POSIXct(paste0(snow.phen$year, "-01-01"))+(phen[,2]-1)*24*60*60, "%b-%d") 
snow.phen$yday.start <- phen[,1]
snow.phen$yday.end   <- phen[,2]

# save(snow.phen, file = "Data/Results/snow.phen.RData")
load("Data/Results/snow.phen.RData")

