##########################################################################################################
#### Snow Data ###########################################################################################
##########################################################################################################

# Snow cover data (netCDF) was downloaded from: http://data.ncdc.noaa.gov/thredds/catalog/cdr/snowcover/catalog.html
# geotif: longitude.tif, latitude.tif and snow_cover_extent.tif were produced using GDAL library (gdal_translate HDF5) in command line tools


library(bbmle)
library(raster)
library(rgdal)

load("~/area.longlat.RData")


### start: maximum likelihood function=
gaussMLE <- function(week, size, prob) {
  
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


snowP <- data.frame(year = rep(min(unique(year)):max(unique(year)), 53), 
                    week = rep(0:52, each = length(unique(year))),
                    p =    as.vector(out))

phen <- matrix(nrow = length(1983:2013), ncol = 2) 

# opar <- par(mfrow = c(7, 5), mar = c(0.5,0.5,0.5,0.5), oma = c(6,6,4,1))
for(i in 1983:2013) {
  tmp1 <- subset(snowP, year == i & week>0)
  # plot(tmp1$week, tmp1$p, type = "o", pch = 16, col = "grey50", lwd = 2, cex = 0.5, ylim = c(-0.1, 1.1), xlim = c(-1,52),
  #     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  # axis(1, tcl = 0.2, labels = ifelse(i%in%c(2008, 2011), TRUE, FALSE), lwd = 1.1)
  # axis(2, tcl = 0.2, labels = ifelse(i%in%c(1983, 1993, 2003), TRUE, FALSE), lwd = 1.1, las = 2)
  # box(lwd = 1.1)
  # points(tmp1$week, tmp1$p2, pch = 16, type = "b", cex = 0.5)
  
  mod <- gaussMLE(week = tmp1$week, size = rep(100, nrow(tmp1)), prob = tmp1$p)
  # lines(mod$result, lwd = 2, col = "orange")
  
  # points(mod$start(0.66), 0.66, pch = 18, cex = 2, col = "firebrick")
  # points(mod$end(0.33), 0.33, pch = 18, cex = 2, col = "cornflowerblue")
  # text(8, 0.05, i, font = 2, cex = 1.4)
  phen[which(c(1983:2013)==i), 1] <- mod$start(0.66)
  phen[which(c(1983:2013)==i), 2] <- mod$end(0.33)
}
# plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# legend("topleft", c("Raw data", "Filtered data", "Model fit", "Start interp.", "End. interp."),
       # pch = c(16, 16, NA, 18, 18), lty = c(1,1,1,NA,NA), col = c("grey50", "black", "orange", "firebrick", "cornflowerblue"), 
       # bty = "n", cex = 1.3)
# mtext("Date (week of the year)", side = 1, outer = T, line = 3.5, cex = 1.3)
# mtext("Percentage of snow cover at breeding area", side = 2, outer = T, line = 3.5, cex = 1.3)
# mtext("Snow Cover Data Manipulation (Curve fitting)", 3, line = 1.5, outer = T, cex = 1.5)
# par(opar)

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

retuen <- snow


