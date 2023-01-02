# Author: Robert J. Hijmans
# License GPL3

if (!isGeneric("photoperiod")) {setGeneric("photoperiod", function(x, ...) standardGeneric("photoperiod"))}


setMethod("photoperiod", signature(x="Date"),
	function(x, latitude) {
		x <- fromDate(x, "doy")
		x <- cbind(x, latitude)
		.Call('_meteor_Photoperiod', PACKAGE = 'meteor', x[,1], x[,2])
	}
)

setMethod("photoperiod", signature(x="numeric"),
	function(x, latitude) {
		x <- cbind(x, latitude)
		.Call('_meteor_Photoperiod', PACKAGE = 'meteor', x[,1], x[,2])
	}
)

setMethod("photoperiod", signature(x="data.frame"),
	function(x) {
		if (!all(c("date", "latitude") %in% names(x))) {
			stop("x must have variables 'date', and 'latitude'")		
		}
		doy <- fromDate(x$date, "doy")
		.Call('_meteor_Photoperiod', PACKAGE = 'meteor', x$date, x$latitude)
	}
)

setMethod("photoperiod", signature(x="SpatRaster"),
	function(x, filename="", overwrite=FALSE, ...) {
		d <- terra::time(x)
		if (all(is.na(d))) {
			stop("the layers in x have no time stamps")
		}
		r <- terra::rast(x)[[1]]
		if (!terra::is.lonlat(x)) {
			lat <- terra::as.points(r, values=FALSE, na.rm=FALSE)
			lat <- terra::project(lat, crs="+proj=longlat")
			lat <- terra::crds(lat)[,2]
			dd <- data.frame(latitude=rep(lat, length(d)), date=rep(d, each=length(lat)))
			r <- terra::rast(x)
			terra::values(r) <- photoperiod(dd)
			if (filename != "") {
				r <- terra::writeRaster(r, filename=filename, overwrite=overwrite, ...)
			}
			r
		} else {
			r <- r[,1,drop=FALSE]
			terra::ext(r) <- terra::ext(x)
			lat <- terra::yFromRow(r, 1:nrow(r))
			dd <- data.frame(latitude=rep(lat, length(d)), date=rep(d, each=length(lat)))
			r <- terra::rast(x)
			terra::values(r) <- photoperiod(dd)
			terra::disagg(r, c(1, ncol(x)), filename=filename, overwrite=overwrite, ...)
		}		
	}
)



.daylength <- function(lat, doy) {

	if (inherits(doy, 'Date') | inherits(doy, 'character')) {
		doy <- fromDate(doy, "doy")
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
	if (inherits(doy, 'Date') | inherits(doy, 'character')) {
		doy <- doyFromDate(doy)
	}
	lat[lat > 90 | lat < -90] <- NA
	doy <- doy %% 365

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
