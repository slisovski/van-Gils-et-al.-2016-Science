library(bbmle)
library(maptools)
  data(wrld_simpl)
library(maps)
library(splancs)

### Figures supplementary materials

setwd("C:/Users/Simeon/Dropbox/Science/Projects/RedKnot_BillLength/Analysis")

### Temperature

load('Data/Basic/temp.RData')

year.vector <- unique(daily$YEAR)
daily$julian.date  <- as.POSIXlt(as.POSIXct(daily$DATE))$yday + 1
daily.selected.range <- subset(daily, RANGE=="calcan")

pdf("Figures/Temp_SM.pdf", height = 10, width = 10)
par(mfrow = c(2, 2), mar = c(3,2,3,2), oma = c(5,5,0,0), las = 1, cex.axis = 1.7)
for (i in 1:length(year.vector)) {
  
  selected.dataset  =  daily.selected.range[daily.selected.range$YEAR==year.vector[i],]
  q2.model=lm(TEMPERATUR~julian.date+I(julian.date^2), selected.dataset)
  
  fnc <- approxfun(x = predict(q2.model)[1:which.max(predict(q2.model))], 
                   y = selected.dataset$julian.date[1:which.max(predict(q2.model))])
  if(year.vector[i]%in%c(1995, 2005, 2015)){
    
    with(selected.dataset, plot(julian.date, TEMPERATUR, pch = 16, col = "grey80", cex = 1.5, 
                                xaxt = "n", ylim = c(-22, 15), xlab = "", ylab = ""))
    lines(selected.dataset$julian.date, predict(q2.model),
          col = "orange", lwd = 5)
    points(fnc(0), 0, pch = 24, col = 1, bg = "cornflowerblue", cex = 3, lwd = 2) 	  
    text(230, -20, year.vector[i], cex = 2.1)
    
    tm.ax <- seq(as.POSIXct(paste0(unique(selected.dataset$YEAR), "-05-01")), as.POSIXct(paste0(unique(selected.dataset$YEAR), "-11-01")), by = "month")
    
    axis(1, at = as.numeric(format(tm.ax, "%j")),
            labels = paste0("1 ", format(tm.ax, "%b")), las = 2, cex.axis = 1.7) 
  }
}

mtext(expression(paste("Temperature (",degree,"C)")), 2, outer = T, las = 3, line = 2, cex = 2)
mtext("Date", 1, outer = T, las = 1, line = 2.7, cex = 2)

plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("raw data", "model fit", expression(italic("D")["T0"])),
       pch = c(16, NA, 24), col = c("grey80", "orange", 1), pt.bg = c(NA, NA, "cornflowerblue"), lty = c(NA, 1, NA),
       bty = "n", cex = 2.2, lwd = c(1, 5, 2))
dev.off()



### NDVI

fls <- list.files("~/Desktop/VH", pattern="*ND.hdf")
fls.year <- as.numeric(sapply(fls,substring,17,20))
fls.week <- as.numeric(sapply(fls,substring,21,23))
fls <- fls[!duplicated(paste(fls.year,fls.week,sep="_"))]
fls.year <- as.numeric(sapply(fls,substring,17,20))
fls.week <- as.numeric(sapply(fls,substring,21,23))

load("Data/Basic/ndvi_raw.RData")

ndvi <- ifelse(ndvi<(-0.15), NA, ndvi)
ind1 <- paste(fls.year, fls.week, sep="_") 

tab <- matrix(c(rep(c(1983:2015), each = 52), rep(c(1:52), length(c(1983:2015)))), ncol = 2)
tab <- cbind(tab,ndvi[match(paste(tab[,1], tab[,2], sep="_"), ind1),])

pdf("Figures/NDVI_SM.pdf", height = 10, width = 10)
opar <- par(mfrow = c(2, 2), mar = c(2,2,2,2), oma = c(5,5,0,0), las = 1, cex.axis = 1.4)

out <- data.frame(year = c(1995, 2005, 2015), onset = NA, end = NA, duration = NA, max = NA, area = NA)

for(i in c(1995, 2005, 2015)){
  
  tmp  <- tab[tab[,1]==i,]  
  y    <- matrix(c(rep(tmp[,1], ncol(tmp)-3), rep(tmp[,2], ncol(tmp)-3), rep(tmp[,3], ncol(tmp)-3), as.vector(tmp[,-c(1:3)])), ncol = 4) 
  
    plot(y[,2], y[,4], pch=16, col=rgb(0.7,0.7,0.7, alpha = 0.1), ylim = c(-0.1, 0.5), xaxt="n")
    
    tm.ax <- seq(as.POSIXct(paste0(i, "-02-01")), as.POSIXct(paste0(i, "-12-01")), by = "month")
    
    axis(1, at = as.numeric(format(tm.ax, "%U")), labels = format(tm.ax, "%b"), las = 2)
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
    
    
    pr2 <- predict(ls, newdata = range[1]:range[2])
    lines(range[1]:range[2], pr2, lwd=4)
    
    ## Data
    ind3 <- which.max(pr2)
    out$max[which(c(2000, 2008, 2010)==i)] <- pr2[ind3]
    
    ## Onset
    onc <- approxfun(x = pr2[1:ind3], y = range[1]:(range[1] + (ind3-1)))
      out$onset[which(c(2000, 2008, 2010)==i)] <- onc(0.005) 

    
    ## end
    end <- approxfun(x = pr2[ind3:length(pr2)], y = c(range[1]:range[2])[ind3:length(pr2)])
    out$end[which(c(2000, 2008, 2010)==i)] <- end(0.1)
    abline(h = c(0, 0.1), lty = 3, col = "grey40", lwd = 2)
    
      ## duration
      out$duration[which(c(2000, 2008, 2010)==i)] <- out$end[which(c(2000, 2008, 2010)==i)] - out$onset[which(c(2000, 2008, 2010)==i)]
      
      ## area
      curve <- cbind(c(seq(onc(0), end(0.1), length = 30), onc(0)), predict(ls, newdata = c(seq(onc(0), end(0.1), length = 30), onc(0))))
      polygon(curve[,1], curve[,2], border = "transparent", col = rgb(0.4, 0.4, 0.2, alpha = 0.4))
      out$area[which(c(2000, 2008, 2010)==i)] <- areapl(curve)
      points(end(0.1), 0.1, pch = 21, cex = 3, lwd = 2, col = 1, bg = "cornflowerblue")
      points(onc(0.005), 0, pch = 21, lwd = 2, cex = 3, col = 1, bg = "firebrick") 
      mtext(i, 3, line = -1.5, cex =1, at = 9)
}
plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("raw NDVI data", "smooth line (span=0.3)", "smooth line (excl. cut-off)", "cut-off dates (albedo)", "phenology thresholds", expression("D1"["DNVI"]), expression("D2"["NDVI"])),
       pch = c(16, NA, NA, NA, NA, 21, 21), lty = c(NA,1,1,2,3,NA,NA), col = c("grey50", "orange", "black", "black", "black", 1, 1),
       pt.bg = c(rep(NA, 5), "firebrick", "cornflowerblue"),
       bty = "n", cex = 1.3, y.intersp  = 1.4)
mtext("date", side = 1, outer = T, line = 3.5, cex = 1.3)
mtext("normalized difference vegetation index (NDVI)", side = 2, outer = T, line = 3.5, cex = 1.3, las = 3)
par(opar)
dev.off()


### snow

### start: maximum likelihood function=
gaussMLE <- function(day, size, prob) {
  
  tab <- data.frame(day = day, size = size, p = prob)
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


load("Data/Basic/snow.RData")


phen <- matrix(nrow = length(1983:2015), ncol = 2) 

pdf("Figures/Snow_SM.pdf", height = 10, width = 10)

opar <- par(mfrow = c(2, 2), mar = c(3,2,2,2), oma = c(5,5,0,0), las = 1, cex.axis = 1.7)
for(i in c(1995, 2005, 2015)) {

  tmp1 <- out[which(as.numeric(format(out[,1], "%Y"))==i),]
  
  p2 <- tmp1[,2]
  for(j in 2:(length(p2-1))) {
    if(!is.na(tmp1[j+1,2]) & tmp1[j-1,2]==tmp1[j+1,2] & tmp1[j,2]!=tmp1[j-1,2] & abs(tmp1[j-1,2] - tmp1[j,2])>0.25) p2[j] <- tmp1[j-1,2]
  } 
  
  days <- as.numeric(format(tmp1[,1], "%j"))
  
  plot(days, tmp1[,2]*100, type = "n", pch = 16, col = "grey50", lwd = 2, cex = 0.5, ylim = c(-5, 105),
       xaxt = "n",xlab = "", ylab = "")
  tm.ax <- seq(as.POSIXct(paste0(i, "-01-01")), as.POSIXct(paste0(i, "-12-01")), by = "month")
  axis(1, at = as.numeric(format(tm.ax, "%j")), labels = paste0("1 ", format(tm.ax, "%b")), las = 2)
    
  mod <- gaussMLE(day = days, size = rep(100, length(tmp1)), prob = p2)
  lines(mod$result[,1], mod$result[,2]*100, lwd = 7, col = "orange")
  lines(days, tmp1[,2]*100, type = "o", pch = 16, cex = 0.9, lwd = 4, col = rgb(0.6, 0.6, 0.6, alpha = 0.7))
  lines(days, p2*100, type = "o", pch = 16, cex = 0.5, lwd = 0.6)
  
  points(mod$start(0.66), 66, pch = 21, cex = 3, col = 1, bg = "firebrick")
  points(mod$end(0.33), 33, pch = 21, cex = 3, col = 1, bg = "cornflowerblue")
  text(50, 0.05, i, font = 2, cex = 1.7)
}
plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topleft", c("raw data", "filtered data", "model fit", expression(italic("D")["SM"]), expression(italic("D")["SF"])),
       pch = c(16, 16, NA, 21, 21), lty = c(1,1,1,NA,NA), col = c("grey50", "black", "orange", "firebrick", "cornflowerblue"), 
       bty = "n", cex = 2.2, lwd = c(2,2,4,2,2))
mtext("Date", side = 1, outer = T, line = 3.2, cex = 1.8)
mtext("Snow Cover (%)", side = 2, outer = T, line = 3.2, cex = 1.8, las = 3)
par(opar)

dev.off()