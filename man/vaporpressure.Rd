\name{vaporpressure}

\alias{SVP}
\alias{VP}
\alias{VPD}


\title{Vapor pressure}

\description{
Functions to compute the saturated vapor pressure (SVP), actual vapor pressure (VP), and vapor pressure deficit (VPD) in Pascal.

For temperature < 0C the saturation vapour pressure equation for ice is used accoring to Goff and Gratch (1946), whereas for temperature >=0C that of Goff (1957) is used.
}

\usage{
SVP(temp)
VP(temp, relh)
VPD(temp, relh)
}


\arguments{
  \item{temp}{numeric. Temperature in degrees C}
  \item{relh}{relative humidity (percent)}
}

\value{
numeric vector (Pascal).
}


\references{
Goff, J.A.,and S. Gratch, 1946. Low-pressure properties of water from -160 to 212 F. Transactions of the American society of heating and ventilating engineers, p. 95-122, presented at the 52nd annual meeting of the American society of heating and ventilating engineers, New York, USA.

Goff, J. A. 1957. Saturation pressure of water on the new Kelvin temperature scale, Transactions of the American society of heating and ventilating engineers, pp 347-354, presented at the semi-annual meeting of the American society of heating and ventilating engineers, Murray Bay, Quebec, Canada.
}

\author{
Robert Hijmans, based on Python meteolib by Maarten J. Waterloo and J. Delsman \url{http://python.hydrology-amsterdam.nl/}
}

\examples{
temperature <- seq(-10,30,10)
SVP(temperature)
VP(temperature, 60)
VPD(temperature, 60)
}
