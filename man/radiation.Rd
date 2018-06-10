\name{ExtraTerrestrialRadiation}

\alias{ExtraTerrestrialRadiation}


\title{Extra-terrestrial Radiation}

\description{
Compute incoming radiation (J day-1 m-2 ) at the top of the atmosphere and photoperiod (daylength, sunshine duration).
}

\usage{
ExtraTerrestrialRadiation(doy, latitude, sc=1367.7, FAO=FALSE)
}


\arguments{
  \item{doy}{integer. Day of the year}
  \item{latitude}{numeric. Latitude}
  \item{sc}{numeric. The solar constant}  
  \item{FAO}{logical. If \code{TRUE} the algorithm described by Allen et al (1998) is used. If \code{FALSE} the approach by Goudriaan and Van Laar (1995) is used}
}

\value{
matrix with incoming radiation (J/day) and 
}


\references{
Goudriaan and Van Laar, 1995.

R.G. Allen, L.S. Pereira, D. Raes and M. Smith (1998). Crop Evaporation - Guidelines for computing crop water requirements. Irrigation and drainage paper 56. FAO, Rome, Italy. \url{http://www.fao.org/docrep/x0490e/x0490e07.htm}
}

\author{
Robert Hijmans, based on Python meteolib by Maarten J. Waterloo and J. Delsman \url{http://python.hydrology-amsterdam.nl/}
}

\examples{
ExtraTerrestrialRadiation(50, 60)
ExtraTerrestrialRadiation(50, 60, FAO=TRUE)
}

