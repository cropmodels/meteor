# Robert Hijmans and Jorrel Khalil S. Aunario
# Date :  8 May 2012
# Version 0.2
# Licence GPL v3


ExtraTerrestrialRadiation <- function(doy, latitude, sc=1367.7, FAO=FALSE) {
	d <- cbind(doy, latitude)
    .Call('_meteor_ExtraTerrestrialRadiation', PACKAGE = 'meteor', d[,1], d[,2], sc[1], FAO[1])
}


# TODO 
.sunhoursToSRad <- function(sunhour, doy, lat, anga, angb){
    deg2rad <- pi / 180
    dec <- -1 * asin(sin(23.45*deg2rad)*cos(2*pi*(doy+10)/365))
    sc <- 1370*(1+0.033*cos(2*pi*doy/365))
    
    sinld <- sin(deg2rad*lat)*sin(dec)
    cosld <- cos(deg2rad*lat)*cos(dec)
    aob <- sinld/cosld  
    
    dayl <- 12*(1+2*asin(aob)/pi)
    dsinb <- 3600*(dayl*sinld+24*cosld*sqrt(1-(aob^2))/pi)
    angot <- sc*dsinb
    
    return(angot*(anga+angb*((sunhour)/dayl))/1000000) 
}


.hargreavesSRad <- function(doy, tmax, tmin, latdd){
    DR <- 1+0.033*cos(pi*2/366*doy)    
    LatRad <- pi/180*latdd
    SRdec <- 0.409*sin(2*pi/366*doy-1.98)
    SHA <- acos(-tan(LatRad)*tan(SRdec))
    R_a <- 24*60/pi*0.082*DR*(SHA*sin(LatRad)*sin(SRdec)+cos(LatRad)*cos(SRdec)*sin(SHA))
    return(R_a*0.16*sqrt(tmax-tmin))       
}
