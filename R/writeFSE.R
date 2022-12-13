# Author: Robert J. Hijmans
# License GPL3


readFSEwth <- function(f) {
  d <- .trim(readLines(f))
  i <- substr(d, 1, 1) == '*'
  d <- d[!i ]
  h <- d[1]
  h <- unlist(strsplit(h, ' '))
  h <- as.numeric(h[h != ""])
  location <- h[1:3]
  d <- d[-1]
  x <- sapply(d, function(i)strsplit(i, ' '), USE.NAMES = FALSE)
  x <- lapply(x, function(y) as.numeric(y[y!=''] ))
  x <- do.call(rbind, x)
  date <- dateFromDoy(x[,3], x[,2])
  df <- data.frame(date, x[,4:9])
  colnames(df) <- c("date", "srad", "tmin", "tmax", "vapr", "wind", "prec")
  attr(df, 'location') <- location
  df 
}



writeFSEwth <- function(w, country="AAA", station=1, lon=0, lat=0, elev=0,  path=".") {
	
	nms <- c("date", "srad", "tmin", "tmax", "wind", "prec", "vapr")
	if (!all(nms %in% colnames(w))) {
		stop(paste("w does not have all names (", paste(nms, collapse=", "), ")"))
	}
	if (!dir.exists(path)) {
		error("path does not exist")
	}
	
	w$year <- yearFromDate(w$date)
	years <- unique(w$year)

	fnames <- file.path(path, paste0(country, station, '.', substr(years, 2, 4)))

	for (i in seq_along(years)) {

		yr <- years[i]

		thefile <- file(fnames[i], "w")
		
		cat("*-----------------------------------------------------------", "\n", file = thefile)
		cat("*  Created by the R package 'weather'\n", file = thefile)
		cat("*", "\n", file = thefile)
		cat("*  Column    Daily Value\n", file = thefile)
		cat("*     1      Station number\n", file = thefile)
		cat("*     2      Year\n", file = thefile)
		cat("*     3      Day\n", file = thefile)
		cat("*     4      irradiance         KJ m-2 d-1\n", file = thefile)
		cat("*     5      min temperature            oC\n", file = thefile)
		cat("*     6      max temperature            oC\n", file = thefile)
		cat("*     7      vapor pressure            kPa\n", file = thefile)
		cat("*     8      mean wind speed         m s-1\n", file = thefile)
		cat("*     9      precipitation          mm d-1\n", file = thefile)
		cat("*\n", file = thefile)
		cat("** WCCDESCRIPTION=meteor\n", file = thefile)
		cat("** WCCFORMAT=2\n", file = thefile)
		cat("** WCCYEARNR=", yr, "\n", file = thefile)
		cat("*-----------------------------------------------------------", "\n", file = thefile)

		cat(lon, lat, elev, '  0.00  0.00 \n', file = thefile)

		yw <- w[w$year==yr, ]
		yw[is.na(yw)] <- -9999
		
		expected <- diff(as.Date(paste0(yr, c("-01-01", "-12-31")))) + 1
		
		if (expected != nrow(yw)) {
			warning(paste("the data for year", yr, "is incomplete"))
		}
		
		for (d in 1:nrow(yw)) {
			cat("1  ", sprintf("%6.0f", yr), sprintf("%5.0f", d), sprintf("%10.0f", yw$srad[d]), sprintf("%8.1f", yw$tmin[d]), sprintf("%8.1f", yw$tmax[d]), sprintf("%8.1f", yw$vapr[d]), sprintf("%8.1f", yw$wind[d]), sprintf("%8.1f", yw$prec[d]), "\n", file=thefile)
		}
		close(thefile)
    }
	return(invisible(fnames))
}		

#wth$srad = wth$srad * 1000
#writeFSEwth(wth, 'NLD', 3, 5.67,51.97, 7 )


