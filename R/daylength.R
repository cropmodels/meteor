# Author: Robert J. Hijmans
# License GPL3


photoperiod <- function(doy, latitude) {
	d <- cbind(doy, latitude)
    .Call('_meteor_Photoperiod', PACKAGE = 'meteor', d[,1], d[,2])
}


.daylength <- function(lat, doy) {
	if (class(doy) == 'Date' | class(doy) == 'character') {
		doy <- doyFromDate(doy)
	}
	lat[lat > 90 | lat < -90] <- NA
	doy <- doy %% 365

#Ecological Modeling_, volume 80 (1995) pp. 87-95, called "A Model
#Comparison for Daylength as a Function of Latitude and Day of the Year."
	P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860*(doy-186)))))
	a <-  (sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(P)) / (cos(lat * pi/180) * cos(P));
	a <- pmin(pmax(a, -1), 1)
	DL <- 24 - (24/pi) * acos(a)
	return(DL)
}

.daylength2 <- function(lat, doy) {
	if (class(doy) == 'Date' | class(doy) == 'character') {
		doy <- doyFromDate(doy)
	}
	lat[lat > 90 | lat < -90] <- NA
	doy[doy==366] <- 365
	doy[doy < 1] <- 365 + doy[doy < 1]
	doy[doy > 365] <- doy[doy > 365] - 365
	if (isTRUE(any(doy<1))  | isTRUE(any(doy>365))) {
		stop('cannot understand value for doy')
	}

# after Goudriaan and Van Laar
	RAD <- pi/180
#  Sine and cosine of latitude (LAT)
    SINLAT <- sin(RAD * lat);
    COSLAT <- cos(RAD * lat);
# Maximal sine of declination;}
    SINDCM <- sin(RAD * 23.45)
#{Sine and cosine of declination (Eqns 3.4, 3.5);}
    SINDEC <- -SINDCM * cos(2*pi*(doy+10)/365)
    COSDEC <- sqrt(1-SINDEC*SINDEC);
#The terms A and B according to Eqn 3.3;}
    A <- SINLAT*SINDEC;
    B <- COSLAT*COSDEC;
    C <- A/B;
#Daylength according to Eqn 3.6; arcsin(c) = arctan(c/sqrt(c*c+1))}
    DAYL <- 12* (1+(2/pi)* atan(C/sqrt(C*C+1)))
	return(DAYL)
}
