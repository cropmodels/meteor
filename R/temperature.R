# Author: Robert J. Hijmans
# License GPL3
# Version 0.2  January 2009, June 2016


hourlyFromDailyTemp <- function(tmin, tmax, doy, latitude) {
		d <- cbind(tmin, tmax, doy, latitude)
    	.Call('_meteor_hourlyFromDailyTemp', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4])
}


hourlyFromDailyRelh <- function(relh, tmin, tmax, doy, latitude) {
		d <- cbind(relh, tmin, tmax, doy, latitude)
    	.Call('_meteor_hourlyFromDailyRH', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4], d[,5])
}



dayTemp <- function(tmin, tmax, doy, latitude) {
	d <- cbind(tmin, tmax, doy, latitude)
    .Call('_meteor_daytimeTemperature', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4])
}


...diurnalTemp <- function(lat, date, tmin, tmax) {
	TC <- 4.0
    P <- 1.5
	dayl <- .daylength(lat, doyFromDate(date))
	nigthl <- 24 - dayl
    sunris <- 12 - 0.5 * dayl
    sunset <- 12 + 0.5 * dayl
	hrtemp <- vector(length=24)
	for (hr in 1:24) {
#    period a: dhour between midnight and sunrise;
		if ( hr < sunris)  {
			tsunst <- tmin+(tmax-tmin)*sin(pi*(dayl/(dayl+2*P)))
			hrtemp[hr] <- (tmin-tsunst*exp(-nigthl/TC)+(tsunst-tmin)*exp(-(hr+24-sunset)/TC))/(1-exp(-nigthl/TC))
		} else if ( hr < (12+P) ) {
#  period b: dhour between sunrise and normal time that tmax is reached (after noon)
			hrtemp[hr] <- tmin+(tmax-tmin)*sin(pi*(hr-sunris)/(dayl+2*P))
		} else if (hr < sunset) {
#  period c: dhour between time of tmax and sunset;
			hrtemp[hr] <- tmin+(tmax-tmin)*sin(pi*(hr-sunris)/(dayl+2*P))
		} else {
#  period d: dhour between sunset and midnight;
			tsunst <- tmin+(tmax-tmin)*sin(pi*(dayl/(dayl+2*P)))
			hrtemp[hr] <- (tmin-tsunst*exp(-nigthl/TC)+(tsunst-tmin)*exp(-(hr-sunset)/TC))/(1-exp(-nigthl/TC))
		}
	}
	return(hrtemp)
}
