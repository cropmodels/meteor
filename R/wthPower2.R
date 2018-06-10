# Author: Robert J. Hijmans
# License GPL3
# Version 1.0  December 2012


.trim <- function(x) { gsub("^\\s+|\\s+$", "", x) }



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

