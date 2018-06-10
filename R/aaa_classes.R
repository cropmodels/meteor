
if (!isGeneric("crop<-")) { setGeneric("crop<-", function(x, value) standardGeneric("crop<-")) }	
if (!isGeneric("soil<-")) { setGeneric("soil<-", function(x, value) standardGeneric("soil<-")) }	
if (!isGeneric("control<-")) { setGeneric("control<-", function(x, value) standardGeneric("control<-")) }	
if (!isGeneric("weather<-")) { setGeneric("weather<-", function(x, value) standardGeneric("weather<-")) }	
if (!isGeneric("run")) { setGeneric("run", function(x, ...) standardGeneric("run")) }	


setClass("Weather",
	slots = c(data = "data.frame", 
	
			longitude = "numeric",
            latitude = "numeric",
            elevation = "numeric",
	        name = "character",
			country = "character",
			ID = "character"),
			
	
	prototype = list(
		data = data.frame(date=as.Date("2000-01-01"), tmin = 0, tmax=0, srad=0, prec=0, wind=0, svap=0),
		longitude = 0,
        latitude = 0,
        elevation = as.numeric(NA)
		),

	validity = function(object) {
		if (abs(latitude) > 90 )
			return("impossible latitude")
		return(TRUE)
	}
)


setMethod("$", "Weather", 
	function(x, name) {
		x@data[[name]]
	}
)


setReplaceMethod("$", "Weather", 
	function(x, name, value) { 
		x@data[[name]] = value 
		x 
	}
)


setMethod("[[", "Weather", 
	function(x, name) {
		x@data[[name]]
	}
)


setReplaceMethod("[[", "Weather", 
	function(x, name, value) { 
		x@data[[name]] = value 
		x 
	}
)


setMethod ('show' , 'Weather', 
	function(object) {
		#cat('class       :' , class(object), '\n')
		cat('longitude:' , object@longitude, ' latitude:' , object@latitude, ' elevation:', object@elevation, '\n')
		cat('\n')
		print(head(object@data, 3))
		cat('...\n')
		print(tail(object@data, 3))
	}
)	
	

