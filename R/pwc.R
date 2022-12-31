# Author: Robert J. Hijmans
# License GPL3

if (!isGeneric("pwc")) {setGeneric("pwc", function(x, ...) standardGeneric("pwc"))}
if (!isGeneric("ucti")) {setGeneric("ucti", function(x, ...) standardGeneric("ucti"))}


setMethod("pwc", signature(x="numeric"),
	function(x, input="wbgt") {
		input = match.arg(tolower(input), c("wbgt", "ucti"))
		if (input == "wbgt") {
			.pwc_wbgt(x)
		} else {
			.pwc_ucti(x)		
		}
	}
)

setMethod("pwc", signature(x="SpatRaster"),
	function(x, input="wbgt", filename="", overwrite=FALSE, ...) {
		input = match.arg(tolower(input), c("wbgt", "ucti"))
		if (input == "wbgt") {
			app(x, .pwc_wbgt, filename=filename, overwrite=overwrite, wopt=list(...))
		} else {
			app(x, .pwc_ucti, filename=filename, overwrite=overwrite, wopt=list(...))
		}
	}
)


setMethod("ucti", signature(x="data.frame"),
	function(x, input="wbgt") {
		if (!all(c("temp", "rhum", "tglb", "wind", "rhum") %in% names(x))) {
			stop("x must have variables 'temp', 'rhum', 'wind', and 'tglb'")		
		}
		.ucti(x$temp, x$tglb, x$wind, x$rhum)
	}
)

