# Author: Robert J. Hijmans
# License GPL3

if (!isGeneric("Tg")) {setGeneric("Tg", function(x, ...) standardGeneric("Tg"))}
if (!isGeneric("Tnwb")) {setGeneric("Tnwb", function(x, ...) standardGeneric("Tnwb"))}
if (!isGeneric("WBGT")) {setGeneric("WBGT", function(x, ...) standardGeneric("WBGT"))}


setMethod("Tg", signature(x="data.frame"),
	function(x, latitude) {
		if (!all(c("temp", "rhum", "wind", "srad", "date") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', 'srad', and 'date'")		
		}
		year <- fromDate(x$date, "year")
		doy <- fromDate(x$date, "doy")
		.Tg1(x$temp, x$rhum, x$wind, x$srad, year, doy, latitude)
	}
)

setMethod("Tg", signature(x="SpatRasterDataset"),
	function(x, filename="", overwrite=FALSE, ...) {
		if (!all(c("temp", "rhum", "wind", "srad") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', 'srad'")		
		}
		date <- terra::time(x$temp)
		if (any(is.na(date))) {
			stop("x$temp does not have valid time stamps")
		}
		year <- fromDate(date, "year")
		doy <- fromDate(date, "doy")
		r <- terra::rast(x$temp)
		nc <- ncol(r)
		out <- terra::rast(r)
		terra::readStart(x)
		b <- terra::writeStart(out, filename, overwrite, wopt=list(...), n=5, sources=terra::sources(x))
		for (i in 1:b$n) {
			temp <- as.vector(t(terra::readValues(x$temp, b$row[i], b$nrows[i], 1, nc, mat=TRUE)))
			rhum <- as.vector(t(terra::readValues(x$rhum, b$row[i], b$nrows[i], 1, nc, mat=TRUE)))
			srad <- as.vector(t(terra::readValues(x$srad, b$row[i], b$nrows[i], 1, nc, mat=TRUE)))
			wind <- as.vector(t(terra::readValues(x$wind, b$row[i], b$nrows[i], 1, nc, mat=TRUE)))
			c1 <- (b$row[i]-1) * nc + 1 
			c2 <- (b$row[i]+b$nrows[i]-1) * nc + 1 
			lat <- terra::yFromCell(r, c1:c2)
			tg <- .Tg2(temp, rhum, wind, srad, year, doy, lat)
			terra::writeValues(out, tg, b$row[i], b$nrows[i])
		}
		terra::readStop(x)
		terra::writeStop(out)
	}
)


setMethod("Tnwb", signature(x="data.frame"),
	function(x, latitude) {
		if (!all(c("temp", "rhum", "wind", "srad", "date") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', 'srad', and 'date'")		
		}
		year <- fromDate(x$date, "year")
		doy <- fromDate(x$date, "doy")
		.Tnwb1(x$temp, x$rhum, x$wind, x$srad, year, doy, latitude)
	}
)

setMethod("WBGT", signature(x="data.frame"),
	function(x, latitude) {
		if (!all(c("temp", "rhum", "wind", "srad", "date") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', 'srad', and 'date'")		
		}
		year <- fromDate(x$date, "year")
		doy <- fromDate(x$date, "doy")
		.Tnwb1(x$temp, x$rhum, x$wind, x$srad, year, doy, latitude, FALSE)
	}
)