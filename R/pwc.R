# Author: Robert J. Hijmans
# License GPL3

if (!isGeneric("pwc")) {setGeneric("pwc", function(x, ...) standardGeneric("pwc"))}
if (!isGeneric("utci")) {setGeneric("utci", function(x, ...) standardGeneric("utci"))}


setMethod("pwc", signature(x="numeric"),
	function(x, input="wbgt") {
		input <- match.arg(tolower(input), c("wbgt", "utci"))		
		if (input == "wbgt") {
			.pwc_wbgt(x)
		} else {
			.pwc_utci(x)		
		}
	}
)

setMethod("pwc", signature(x="SpatRaster"),
	function(x, input="wbgt", filename="", overwrite=FALSE, ...) {
		input <- match.arg(tolower(input), c("wbgt", "utci"))		
		out <- terra::rast(x)
		terra::readStart(x)
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, overwrite, wopt=list(...), n=5, sources=terra::sources(x))
		if (input == "wbgt") {
			for (i in 1:b$n) {
				v <- terra::readValues(x, b$row[i], b$nrows[i], mat=TRUE)
				v <- .pwc_wbgt(v)			
				terra::writeValues(out, v, b$row[i], b$nrows[i])
			}
		} else {
			for (i in 1:b$n) {
				v <- terra::readValues(x, b$row[i], b$nrows[i], mat=TRUE)
				v <- .pwc_ucti(v)			
				terra::writeValues(out, v, b$row[i], b$nrows[i])
			}
		}
		terra::writeStop(out)	
	}
)


setMethod("utci", signature(x="data.frame"),
	function(x) {
		if (!all(c("temp", "rhum", "tglb", "wind") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', and 'tglb'")		
		}
		.utci(x$temp, x$tglb, x$wind, x$rhum)
	}
)

