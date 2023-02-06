# Author: Robert J. Hijmans
# License GPL3
# Version 1.0  December 2012


.trim <- function(x) { gsub("^\\s+|\\s+$", "", x) }

.downloadNASAPower <- function(lon, lat) {
	filename <- tempfile()
	vars <- c("toa_dwn", "swv_dwn", "lwv_dwn", "T2M", "T2MN", "T2MX", "RH2M", "DFP2M", "RAIN", "WS10M")
	d <- as.Date(Sys.time())
	eday <- dayFromDate(d)
	emon <- monthFromDate(d)
	eyr <- yearFromDate(d)
	part1 <- "http://earth-www.larc.nasa.gov/cgi-bin/cgiwrap/solar/agro.cgi?email=&step=1&lat="
	part2 <- paste(lat, "&lon=", lon, "&sitelev=&ms=1&ds=1&ys=1983&me=", emon, "&de=", eday, "&ye=", eyr, sep="")
	part3 <- paste("&p=", vars, sep='', collapse='')
	part3 <- paste(part3, "&submit=Submit", sep="")
	theurl <- paste(part1, part2, part3, sep="")
	utils::download.file(url=theurl, destfile=filename, method="auto", quiet=TRUE, mode="wb", cacheOK=TRUE)
	return(filename)
}

.processPower <- function(inf, outf) {

	lns <- readLines(inf)
	hdr <- lns[1:30]
	end <- which(hdr=="-END HEADER-")
	if (length(end) < 1) {
		stop('file not good')
	}
	
	lonlat <- strsplit(hdr[4], ' ')[[1]][c(7,3)]
	elevation <- strsplit(hdr[6], '= ')[[1]][2]
	d <- lns[(end-1):length(lns)]
	d <- d[-2]
	
	d <- strsplit ( gsub("[[:space:]]+", " ", gsub("[[:space:]]+$", "", d))  , " ")
	v <- do.call(rbind, d)
	v[v == '-'] <- NA
	cn <- v[1,]
	v <- v[-1,]
	d <- matrix(as.numeric(v), ncol=ncol(v))

	cns <- c('YEAR', 'DOY', 'toa_dwn', 'swv_dwn', 'lwv_dwn', 'T2M', 'T2MN', 'T2MX', 'RH2M', 'DFP2M', 'RAIN', 'WS10M')
	stopifnot(all(cn == cns))
	#colnames(d) <- cn
	colnames(d) <- c("year", "doy", "erad", "srad", "lwav", "tavg", "tmin", "tmax", "rhum", "tdew", "prec", "wind")	
	d <- data.frame(d)
	#rhnx <- rhMinMax2(d)
	#d <- cbind(d, rhnx)
	
	d$date <- as.Date(d$doy, origin=paste(d$year-1, "-12-31", sep=''))
	d$vapr <- d$rhum * .SVP(d$tavg) / 1000
	
	d <- d[, c("date", "srad", "lwav", "tavg", "tmin", "tmax", "tdew", "vapr", "rhum", "rhmn", "rhmx", "prec", "wind")]
	attr(d, 'location') <- as.numeric(c(lonlat, elevation))
	saveRDS(d, file=outf)
}


.cellFromLL <- function(lon, lat, res=1) {
	res <- 1 / res
	nrows <- 180 * res
    ncols <- 360 * res
	
	row <- floor((90 - lat) * res) 
	row[row < 0] <- 0
	row[row > (nrows-1)] <- nrows - 1

	col <- floor((lon + 180) * res)
	col[col < 0] <- 0
	col[col > (ncols-1)] <- ncols - 1

	(row) * ncols + col + 1
}

.llFromCell <- function(cell, res=1) {
	nrows <- 180 / res
    ncols <- 360 / res
	cell <- cell -1
	
    col = cell %% ncols
    row = trunc(cell / ncols);

    lon <- (col + 0.5) * res - 180
    lat <- 90 - (row + 0.5) * res
	
	cbind(lon, lat)
}



.wthPower <- function(lon, lat, folder=file.path(getwd(), 'power'), ...) {
	cell <- .cellFromLL(lon, lat)
	if (is.na(cell)) {	stop("invalid coordinates") }
	xy <- .llFromCell(cell)
	lon <- xy[1]
	lat <- xy[2]
	tile <- .cellFromLL(lon, lat, 30)
	
	folder <- file.path(folder, paste0("tile_", tile))
	if (!file.exists(folder)) {
		dir.create(folder, recursive=TRUE, showWarnings=FALSE)
	}
	fname <- file.path(folder, paste0(cell, ".rds"))
	if (! (file.exists(fname)) ) {
			theurl <- paste0("http://biogeo.ucdavis.edu/data/climate/daily/nasatiles/", tile, ".zip")
			tfilename <- paste0(folder, "/tile_", tile, ".zip", sep="")
			#if (! file.exists(tfilename)) {
			message(paste('downloading to', folder))
			utils::download.file(url=theurl, destfile=tfilename, method="auto", quiet=TRUE, mode = "wb")
			#}
			utils::unzip(tfilename, exdir=folder)
	}
	if (file.exists(fname)) {
        readRDS(fname)
    } else {
		stop('could not do it')
	}
}



.wthPowerOne <- function(lon, lat, folder=file.path(getwd(), 'power'), overwrite=FALSE, source='Davis', ...) {
	cell <- .cellFromLL(lon, lat)
	if (is.na(cell)) {	stop("invalid coordinates") }
	xy <- .llFromCell(cell)
	lon <- xy[1]
	lat <- xy[2]
	
	fname <- file.path(folder, paste0("nasa-power_", cell, ".rds"))
	
	if (! (file.exists(fname) | overwrite ) ) {
		dir.create(folder, FALSE, TRUE)
	
		source <- tolower(source)
		if (source == 'davis') {
			theurl <- paste0("http://biogeo.ucdavis.edu/data/climate/daily/nasapower/", cell, ".rds")
			utils::download.file(url=theurl, destfile=fname, method="auto", quiet=TRUE, mode="wb", cacheOK=TRUE)
		} else {
			tmpfile <- NULL
			tmpfile <- .downloadNASAPower(lon, lat)
			if (is.null(tmpfile)) {
				stop('it did not work out')
			}
			.processPower(tmpfile, fname)
		}
	}
	
	if (file.exists(fname)) {
        readRDS(fname)
    } else {
		stop('could not do it')
	}
}



..ICASAstyle <- function(lns) {

	lns <- lns[10:(length(lns)-1)]
	lns <- strsplit ( gsub("[[:space:]]+", " ", gsub("[[:space:]]+$", "", lns))  , " ")
	lns <- data.frame(matrix(as.numeric(unlist(lns)), ncol=length(lns[[1]]), byrow=T))
	#colnames(lns) <- h2[[1]]
	lns <- lns[,-1]
	colnames(lns) <- c("year", "doy", "srad", "tmax", "tmin", "prec", "wind", "tdew", "tmp", "rhum")

	#rhnx <- rhMinMax(lns[,'rhum'], lns[,'tmin'], lns[,'tmax'], lns[,'tmp']) 
	#vapr <- lns[,'rhum'] * .SVP(lns[,'tmp']) / 1000     # 100 for % and 10 to go from hPa to kPa
	#lns <- cbind(lns, rhnx, vapr)
	
	date <- dateFromDoy(lns[,'doy'], lns[,'year'])
	#lns <- cbind(as.data.frame(date), lns)
	lns[,-c(1:2)]
}

..getWthFileNASA <- function(filename) {

	lns <- readLines(filename)
	hdr <- lns[1:30]
	end <- which(hdr=="@ WEYR WEDAY  SRAD   TMAX   TMIN   RAIN   WIND   TDEW    T2M   RH2M")

	if (length(end) > 0) {
		h <- which(hdr == "@ INSI   WTHLAT   WTHLONG  WELEV   TAV   AMP  REFHT  WNDHT") + 1
		h1 <- hdr[h]
		y <- as.numeric(substr(h1, 9, 15))
		x <- as.numeric(substr(h1, 18, 25))
		alt <- as.numeric(substr(h1, 27, 32))

		h2 <- lns[end]
		h2 <- strsplit ( gsub("[[:space:]]+", " ", gsub("[[:space:]]+$", "", h2))  , " ")
		end2 <- which(lns=="</PRE></BODY></HTML>") - 1
		
		lns <- .trim( lns[ c((end+1):end2) ] )
	
		lns <- strsplit ( gsub("[[:space:]]+", " ", gsub("[[:space:]]+$", "", lns))  , " ")
		v <- unlist(lns)
		
		v <- as.numeric(v)
		v[v == -99] <- NA
		v <- data.frame(matrix(v, ncol=length(lns[[1]]), byrow=T))

		colnames(v) <- c("year", "doy", "srad", "tmax", "tmin", "prec", "wind", "tdew", "tmp", "rhum")	
		
	} else {

		end <- which(hdr=="-END HEADER-")
		if (end != 18) { warning('strange file') }
		
		hdr <- hdr[1:end]
		loc <- which(substr(hdr, 1, 9)=='Location:')
		loc <- substr(hdr[loc], 11, 100)
		loc <- unlist(strsplit(loc, "   "))
		loc <- unlist(strsplit(loc, " "))
		x <- as.numeric(loc[2])
		y <- as.numeric(loc[4])
		alt <- hdr[which(substr(hdr, 1, 9)=='Elevation')]
		alt <- as.numeric(unlist(strsplit(alt, "="))[2])
		
		lns <- lns[-c(1:end)]
		lns <- strsplit ( gsub("[[:space:]]+", " ", gsub("[[:space:]]+$", "", lns))  , " ")
		v <- unlist(lns)
		v[v=="-"] <- NA
		v <- data.frame(matrix(as.numeric(v), ncol=length(lns[[1]]), byrow=T))

		nicevars <- c("year", "doy", "srad", "tmp", "tmin", "tmax", "rhum", "prec")
		colnames(v) <- nicevars
	}
	
#	rhnx <- rhMinMax2(v) 
	v$vapr <- v$rhum * .SVP(v$tmp) / 1000     # 100 for % and 10 to go from hPa to kPa
#	rh <- cbind(v$rhum, rhnx, vapr)
#	colnames(rh) <- c("rh", "rhmin", "rhmax", "vapr")
#	i <- which(colnames(v) == 'rhum')
#	v <- cbind(v[, -i], rh)

	return(v)
	
#date <- dateFromDoy(v[,'doy'], v[,'year'])

#	r <- raster()
#	cell <- cellFromXY(r, c(x, y))
#	v <- v[, -c(1:2)]
#	makeWeather(cell, x, y, alt, date, v)
}

