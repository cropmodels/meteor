
SVP <- function(temp) {
    .Call(`_meteor_SVP`, temp)
}

VP <- function(temp, relh) {
	x <- cbind(temp, relh)
    .Call(`_meteor_VP`, x[,1], x[,2])
}

VPD <- function(temp, relh) {
	x <- cbind(temp, relh)
    .Call(`_meteor_VPD`, x[,1], x[,2])
}