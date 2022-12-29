
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



#svp = .611 * 10^(7.5 * tmp / (237.7 + tmp))
# vp = .611 * 10^(7.5 * tdw / (237.7 + tdw))
# vp / .611 =  10^(7.5 * tdw / (237.7 + tdw))
# log10(vp / .611) = 7.5 * tdw / (237.7 + tdw)
# log10(vp /.611) / 7.5 = tdw / (237.7 + tdew)
## y = x / (x + a)  =>   x = ya / (1-y)

tDew <- function(temp, relh) {
	relh <- pmin(relh, 100)
	relh <- pmax(relh, 0)
	svp <- .saturatedVaporPressure(temp)
	vp <- svp * relh / 100
	y <- log10(vp /.611) / 7.5
	(y * 237.7) / (1 - y)
}
