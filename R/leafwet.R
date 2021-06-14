# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3

.eLW <- function(rhmin, rhmax, tmin) {
# empirical leaf wetness estimation according to Hijmans, Forbes and Walker, 2001
    ewhr <- exp(-8.093137318+0.11636662*rhmax-0.03715678*rhmin+0.000358713*rhmin*rhmin)
    if (rhmin < 52) {
      ewhr52 <- exp(-8.093137318+0.11636662*rhmax-0.03715678*52+0.000358713*52*52);
      ewhr <- ewhr52 - (ewhr - ewhr52);
	}
    ewhr <- max(0, min(ewhr, 24))
    if (tmin < 0) {
		ewhr <- 0
	}
	return(ewhr)
}


.leafWet <- function(rhmn, rhmx, tmin, tmax, lat, date, simple=TRUE) {
	rh <- ( rhmn + rhmx) / 2
	rh <- .diurnalRH(rh, tmin, tmax, lat, date)
	if (simple) {
		lw <- length(rh[rh>=90])
	} else {
		w <- rh
		x <- (rh - 80) / (95 - 80)
		w[rh > 95] <- 1
		w[rh < 95] <- x[rh < 95]
		w[rh < 80] <- 0
		lw <- sum(w)
	}
	return(lw)
}


.leafWetWithRain <- function(rhmn, rhmx, prec, simple=TRUE) {
	lw <- .leafWet(rhmn, rhmx, simple=simple)
	prec[is.na(prec)] <- 0
	prhrs <- pmin(12, prec / 5)
	return(lw + (1 - lw/24) * prhrs)
}
