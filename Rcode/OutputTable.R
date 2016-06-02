### OutputTable

## wd
setwd("~/Dropbox/Science/Projects/RedKnot_BillLength/Analysis")

## temp
load("Data/Results/temp.phen.RData")
  T0  <- data.frame(Year = temp$year, T0.yday = temp$temp.start, 
                   T0.date = format(as.POSIXct(paste0(temp$year, "-01-01"), tz = "GMT")+(temp$temp.start-1)*24*60*60, "%b %d"))

## ndvi

load("Data/Results/ndvi.phen.RData")
  NDVI  <- data.frame(Year = ndvi$year) 
  tmp01 <- t(apply(ndvi, 1, function(w) {
                          f <- approxfun(x = as.numeric(format(seq(as.POSIXct(paste0(w[1], "-01-01")), 
                                                                   as.POSIXct(paste0(w[1], "-12-31")), by = "day"), "%U")),
                                         y = seq(1:length(seq(as.POSIXct(paste0(w[1], "-01-01")), 
                                                              as.POSIXct(paste0(w[1], "-12-31")), by = "day"))))
                          cbind(f(w[2]), f(w[3]), w[4], w[5])
                          }))
  NDVI$D1_yday <- tmp01[,1] 
  NDVI$D1_date <- format(as.POSIXct(paste0(NDVI$Year, "-01-01"))+(tmp01[,1]-1)*24*60*60, "%b %d")

  NDVI$D2_yday <- tmp01[,2] 
  NDVI$D2_date <- format(as.POSIXct(paste0(NDVI$Year, "-01-01"))+(tmp01[,2]-1)*24*60*60, "%b %d")

  NDVI$Max  <- tmp01[,3]
  NDVI$Area <- tmp01[,4]*100

## snow
load("Data/Results/snow.phen1.RData")
load("Data/Results/snow.phen2.RData")


#phen1

  Snow1 <- data.frame(Year = snow.phen1$year, Dsm_yday1 = snow.phen1$yday.start, Dsm_date1 = format(as.POSIXct(snow.phen1$start), "%b %d"),
                      Dsf_yday1 = snow.phen1$yday.end, Dsm_date1 = format(as.POSIXct(snow.phen1$end), "%b %d"))

  Snow2 <- data.frame(Year = snow.phen2$year, Dsm_yday2 = snow.phen2$ydya.start, 
                      Dsm_date2 = format(as.POSIXct(paste0(snow.phen2$year, "-01-01"))+(snow.phen2$ydya.start*40*60*60)-1, "%b %d"),
                      Dsf_yday2 = snow.phen2$yday.end, 
                      Dsf_date2 = format(as.POSIXct(paste0(snow.phen2$year, "-01-01"))+(snow.phen2$yday.end*40*60*60)-1, "%b %d"))


out <- merge(T0, NDVI, by = "Year", all.x = T)
       out <- merge(out, Snow1, all.x = T)
       out <- merge(out, Snow2, all.x = T)

write.csv(out, "Data/Results/Summary.csv", row.names = F)


#####

load("Data/Basic/knot.data.all.RData")

out$Bill <- knot.data$bill[match(out$Year, knot.data$year)]
out$N    <- knot.data$n[match(out$Year, knot.data$year)]  

plot(ifelse(!is.na(out$Dsm_yday1), out$Dsm_yday1, out$Dsm_yday2), out$Bill, type = "n")
text(ifelse(!is.na(out$Dsm_yday1), out$Dsm_yday1, out$Dsm_yday2), out$Bill, substring(out$Year, 3,4))

