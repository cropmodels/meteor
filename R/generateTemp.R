# R.J. Hijmans
# Licence GPL v3

.monthToYearV <- function(x) {
	y <- rep(NA, 365)
	y[1:31] <- rep(x[1], 31)
	y[32:59] <- rep(x[2], 28)
	y[60:90] <- rep(x[3], 31)
	y[91:120] <- rep(x[4], 30)
	y[121:151] <- rep(x[5], 31)
	y[152:181] <- rep(x[6], 30)
	y[182:212] <- rep(x[7], 31)
	y[213:243] <- rep(x[8], 31)
	y[244:273] <- rep(x[9], 30)
	y[274:304] <- rep(x[10], 31)
	y[305:334] <- rep(x[11], 30)
	y[335:365] <- rep(x[12], 31)
	y
}


.monthToYearM <- function(x) {
	y <- matrix(NA, nrow=365, ncol=ncol(x))
	colnames(y) <- colnames(x)
	y[1:31,] <- rep(x[1,], each=31)
	y[32:59,] <- rep(x[2,], each=28)
	y[60:90,] <- rep(x[3,], each=31)
	y[91:120,] <- rep(x[4,], each=30)
	y[121:151,] <- rep(x[5,], each=31)
	y[152:181,] <- rep(x[6,], each=30)
	y[182:212,] <- rep(x[7,], each=31)
	y[213:243,] <- rep(x[8,], each=31)
	y[244:273,] <- rep(x[9,], each=30)
	y[274:304,] <- rep(x[10,], each=31)
	y[305:334,] <- rep(x[11,], each=30)
	y[335:365,] <- rep(x[12,], each=31)
	y
}


.genTmp <- function(tmin, tmax, reps=1, std=.1, autocor=5) {
	stopifnot(length(tmin) == length(tmax))
	n <- length(tmin)
	nyrs <- n / 12
	if (nyrs %% 1 != 0) {
		stop('monthly climate data must be complete for a year')
	}
	reps <- round(reps)
	stopifnot(reps > 0)
	std <- std[1]
	
	nd <- c(31,28,31,30,31,30,31,31,30,31,30,31)
	nd <- rep(nd, nyrs)
	ds <- cumsum(nd) - nd/2
	ds <- c(ds[1]-31, ds, ds[length(ds)]+31)
	m <- 365*nyrs
	
	res <- list()
	for (r in 1:reps) {
		x <- stats::filter(stats::rnorm(m, sd=std[1]), filter=rep(1,autocor), circular=TRUE)
		y <- stats::filter(stats::runif(m, -0.5*std[1], 0.5*std[1]), filter=rep(1,autocor), circular=TRUE)
		tn <- c(tmin[length(tmin)], tmin, tmin[1])
		tx <- c(tmax[length(tmax)], tmax, tmax[1])
		tmin <- stats::spline(ds, tn, xout=1:m)$y + x - y
		tmax <- stats::spline(ds, tx, xout=1:m)$y + x + y
		res[[r]] <- cbind(tmin,tmax)
	}
	res
}


#g <- .genTmp(tmn, tmx)
#plot(g[[1]])
 
	