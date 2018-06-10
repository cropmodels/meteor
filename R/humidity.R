# Author: Robert J. Hijmans
# License GPL3


.saturatedVaporPressure <- function(tmp) {
   .611 * 10^(7.5 * tmp / (237.7 + tmp))  #kpa
}



.vaporPressureDeficit <- function(tmp, rh) {
    svp <- .saturatedVaporPressure(tmp)
    (1-(rh/100)) * svp
}


.rhMinMax <- function(rh, tmin, tmax) {

	tmin <- pmax(tmin, -5)
	tmax <- pmax(tmax, -5)
	tmp <- (tmin + tmax) / 2

	es <- .saturatedVaporPressure(tmp)
	vp <- rh / 100 * es

	es <- .saturatedVaporPressure(tmax)
	rhmn <- 100 * vp / es;
	rhmn <- pmax(0, pmin(100, rhmn))

	es <- .saturatedVaporPressure(tmin)
	rhmx <- 100*vp/es;
	rhmx <- pmax(0, pmin(100, rhmx))
	cbind(rhmn, rhmx)
}


.rhMinMax2 <- function(tmin, tmax, rhum) {

	tmin <- pmax(tmin, -5)
	tmax <- pmax(tmax, -5)
	tmp <- (tmin + tmax) / 2

	es <- .saturatedVaporPressure(tmp)
	vp <- rhum / 100 * es

	es <- .saturatedVaporPressure(tmax)
	rhmn <- 100 * vp / es;
	rhmn <- pmax(0, pmin(100, rhmn))

	es <- .saturatedVaporPressure(tmin)
	rhmx <- 100*vp/es;
	rhmx <- pmax(0, pmin(100, rhmx))
	cbind(rhmn, rhmx)
}

.diurnalRH <- function(rh, tmin, tmax, lat, date) {
	tmin <- pmax(tmin, -5)
	tmax <- pmax(tmax, -5)
	tmp <- (tmin + tmax) / 2
	vp <- .saturatedVaporPressure(tmp) * rh / 100

	hrtemp <- ...diurnalTemp(lat, date, tmin, tmax)
	hr <- 1:24
	es <- .saturatedVaporPressure(hrtemp[hr])
	rh <- 100*vp/es
	rh <- pmin(100, pmax(0, rh))
	return(rh)
}




.tDew <- function(temp, rh) {
	temp - (100 - rh)/5
}


.FtoC <- function(x) {(5/9)*(x-32) }
.CtoF <- function(x) { x*9/5 + 32 }

.atmp <- function(alt) {
  101.325 * (1 - 2.25577 * 10^-5 * alt) ^ 5.25588   # kPa
}


.rel2abshum <- function(rh, t) {
	es <- .saturatedVaporPressure(t)
	ea <- rh * es / 100
	M <- 18.02 # g/mol
	R <- 8.314472 # Pa?m?/(mol?K)
	T <- t + 273.15  # C to K
	hum <- ea*M/(T*R)
	return(hum)
}


.abs2rhumum <- function(hum, t) {
	M <- 18.02 # g/mol
	R <- 8.314472 # Pa?m?/(mol?K)
	T <- t + 273.15  # C to K
	ea <- hum / (M/(T*R))
	es <- .saturatedVaporPressure(t)
	rh <- 100 * ea / es
	rh  <- pmin(rh, 100)
	return(rh)
}



.rel2spechum <- function(rh, t, alt) {
	es <- .saturatedVaporPressure(t)
	ea <- es * (rh / 100)
	p <- .atmp(0)
	0.62198*ea / (p - ea)
}

.spec2rhumum <- function(spec, t, alt) {
	es <- .saturatedVaporPressure(t)
	100 * (spec * .atmp(alt)) / ((0.62198 + spec) * es)
}
