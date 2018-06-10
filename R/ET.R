# Author: Robert J. Hijmans
# License GPL3

ET0_PenmanMonteith <- function(temp, relh, atmp, Rn, G, ra, rs) {
		d <- cbind(temp, relh, atmp, Rn, G, ra, rs)
    .Call('_meteor_ET0_PenmanMonteith', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4], d[,5], d[,6], d[,7])
}

ET0_PriestleyTaylor <- function(temp, relh, atmp, Rn, G) {
		d <- cbind(temp, relh, atmp, Rn, G)
    .Call('_meteor_ET0_PriestleyTaylor', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4], d[,5])
}

ET0_Makkink <- function(temp, relh, atmp, Rs) {
		d <- cbind(temp, relh, atmp, Rs)
    .Call('_meteor_ET0_Makkink', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4])
}


ET0_ThornthwaiteWilmott <- function(temp, doy, latitude) {
		d <- cbind(temp, doy, latitude)
		.Call('_meteor_ET0_ThornthwaiteWilmott', PACKAGE = 'meteor', d[,1], d[,2], d[,3])
}

ET0_ThornthwaiteWilmottCamargo <- function(tmin, tmax, doy, latitude, Pereira=FALSE) {
		d <- cbind(tmin, tmax, doy, latitude)
		.Call('_meteor_ET0_ThornthwaiteWilmottCamargo', PACKAGE = 'meteor', d[,1], d[,2], d[,3], d[,4], Pereira)
}
