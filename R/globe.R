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
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, overwrite, wopt=list(...), n=5, sources=terra::sources(x))
		for (i in 1:b$n) {
			temp <- terra::readValues(x$temp, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			rhum <- terra::readValues(x$rhum, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			srad <- terra::readValues(x$srad, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			wind <- terra::readValues(x$wind, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			c1 <- (b$row[i]-1) * nc + 1 
			c2 <- (b$row[i]+b$nrows[i]-1) * nc 
			lat <- terra::yFromCell(r, c1:c2)
			tg <- .Tg2(temp, rhum, wind, srad, year, doy, lat)
			terra::writeValues(out, tg, b$row[i], b$nrows[i])
		}
		terra::writeStop(out)
	}
)


setMethod("Tnwb", signature(x="data.frame"),
	function(x, latitude, kelvin=FALSE) {
		if (!all(c("temp", "rhum", "wind", "srad", "date") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', 'srad', and 'date'")		
		}
		year <- fromDate(x$date, "year")
		doy <- fromDate(x$date, "doy")
		.Tnwb1(x$temp, x$rhum, x$wind, x$srad, year, doy, latitude, kelvin, TRUE)
	}
)


setMethod("Tnwb", signature(x="SpatRasterDataset"),
	function(x, kelvin=FALSE, filename="", overwrite=FALSE, ...) {
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
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, overwrite, wopt=list(...), n=5, sources=terra::sources(x))
		for (i in 1:b$n) {
			temp <- terra::readValues(x$temp, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			rhum <- terra::readValues(x$rhum, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			srad <- terra::readValues(x$srad, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			wind <- terra::readValues(x$wind, b$row[i], b$nrows[i], 1, nc, mat=TRUE)
			c1 <- (b$row[i]-1) * nc + 1 
			c2 <- (b$row[i]+b$nrows[i]-1) * nc 
			lat <- terra::yFromCell(r, c1:c2)
			tnwb <- .Tnwb2(temp, rhum, wind, srad, year, doy, lat, kelvin, TRUE)
			tnwb <- as.vector(matrix(tnwb, ncol=ncol(temp), nrow=nrow(temp), byrow=TRUE))
			terra::writeValues(out, tnwb, b$row[i], b$nrows[i])
		}
		terra::writeStop(out)
	}
)


setMethod("WBGT", signature(x="data.frame"),
	function(x, latitude, kelvin=FALSE) {
		if (!all(c("temp", "rhum", "wind", "srad", "date") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', 'srad', and 'date'")		
		}
		year <- fromDate(x$date, "year")
		doy <- fromDate(x$date, "doy")
		.Tnwb1(x$temp, x$rhum, x$wind, x$srad, year, doy, latitude, kelvin, FALSE)
	}
)


setMethod("WBGT", signature(x="SpatRasterDataset"),
	function(x, kelvin=FALSE, mask=NULL, filename="", overwrite=FALSE, ...) {
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
		stopifnot(terra::is.lonlat(r))
		out <- terra::rast(r)
		wopt <- list(...)
		if (is.null(wopt$names)) {
			wopt$names <- paste0("wbgt_", 1:nlyr(out))
		}
		if (!is.null(mask)) {
			if (nlyr(mask) > 1) {
				stop("mask should have a single layer")			
			}
			if (!compareGeom(out, mask)) {
				stop("the geometry of mask does not match the geometry of x")
			}
			terra::readStart(mask)
		}
		nc <- ncol(r)
		terra::readStart(x)
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, overwrite, wopt=wopt, n=5, sources=terra::sources(x))
		for (i in 1:b$n) {
			temp <- terra::readValues(x$temp, b$row[i], b$nrows[i], mat=TRUE)
			rhum <- terra::readValues(x$rhum, b$row[i], b$nrows[i], mat=TRUE)
			srad <- terra::readValues(x$srad, b$row[i], b$nrows[i], mat=TRUE)
			wind <- terra::readValues(x$wind, b$row[i], b$nrows[i], mat=TRUE)
			c1 <- (b$row[i]-1) * nc + 1 
			c2 <- (b$row[i]+b$nrows[i]-1) * nc 
			lat <- terra::yFromCell(r, c1:c2)
			if (!is.null(mask)) {
				m <- terra::readValues(mask, b$row[i], b$nrows[i], mat=TRUE)
				temp[is.na(m), ] <- NA 
			}
			tnwb <- .Tnwb2(temp, rhum, wind, srad, year, doy, lat, kelvin, FALSE)
			tnwb <- as.vector(matrix(tnwb, ncol=ncol(temp), nrow=nrow(temp), byrow=TRUE))
			terra::writeValues(out, tnwb, b$row[i], b$nrows[i])
		}
		if (!is.null(mask)) {
			terra::readStop(mask)
		}
		terra::writeStop(out)
	}
)

