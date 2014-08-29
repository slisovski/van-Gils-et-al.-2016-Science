##########################################################################################################
#### Temperature data ####################################################################################
##########################################################################################################

# Temp = data.frame containing year (year), julian date (jDate) and temperature (temp)

D <- c()

for (i in unique(Temp$year)) {
	
	selected.dataset  =  Temp[Temp$year==i, ]
	
	q2.model=lm(temp~jDate+I(jDate^2), selected.dataset)
	
	fnc <- approxfun(x = predict(q2.model)[1:which.max(predict(q2.model))], 
                     y = selected.dataset$jDate[1:which.max(predict(q2.model))])
                     
    D <- c(D, fnc(0))

}                 

output <- data.frame(year = unique(Temp$year), Dt0 = D)