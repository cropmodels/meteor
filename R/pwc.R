# Author: Robert J. Hijmans
# License GPL3

if (!isGeneric("pwc")) {setGeneric("pwc", function(x, ...) standardGeneric("pwc"))}
if (!isGeneric("utci")) {setGeneric("utci", function(x, ...) standardGeneric("utci"))}


setMethod("pwc", signature(x="numeric"),
	function(x, input="wbgt") {
		input = match.arg(tolower(input), c("wbgt", "utci"))
		if (input == "wbgt") {
			.pwc_wbgt(x)
		} else {
			.pwc_utci(x)		
		}
	}
)

setMethod("pwc", signature(x="SpatRaster"),
	function(x, input="wbgt", filename="", overwrite=FALSE, ...) {
		input = match.arg(tolower(input), c("wbgt", "utci"))
		if (input == "wbgt") {
			terra::app(x, .pwc_wbgt, filename=filename, overwrite=overwrite, wopt=list(...))
		} else {
			terra::app(x, .pwc_utci, filename=filename, overwrite=overwrite, wopt=list(...))
		}
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

