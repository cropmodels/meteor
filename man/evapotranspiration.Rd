\name{evapotranspiration}

\alias{ET0_PenmanMonteith}
\alias{ET0_PriestleyTaylor}
\alias{ET0_Makkink}
\alias{ET0_ThornthwaiteWilmott}
\alias{ET0_ThornthwaiteWilmottCamargo}


\title{Reference evapo-transpiration }

\description{
Functions to compute the reference evapotranspiration (ET0) from meteorological data. ET0 is a representation of the atmospheric water demand. The equations estimate the evapotranspiration rate of a short green crop (grass), completely shading the ground, of uniform height and with adequate water status in the soil profile. Actual evapotranspiration is equal reference evapotranspiration when there is ample water, but taller crops could have an evapotranspiration rate that is higher than ET0.
}

\usage{
ET0_PenmanMonteith(temp, relh, atmp, Rn, G, ra, rs)
ET0_PriestleyTaylor(temp, relh, atmp, Rn, G)
ET0_Makkink(temp, relh, atmp, Rs)
ET0_ThornthwaiteWilmott(temp, doy, latitude)
ET0_ThornthwaiteWilmottCamargo(tmin, tmax, doy, latitude, Pereira=FALSE)
}


\arguments{
  \item{temp}{temperature (degrees C)}
  \item{relh}{relative humidity (percent)}
  \item{atmp}{air pressure (hPa}
  \item{Rn}{net radiation (J m-2 day-1)}
  \item{Rs}{incoming solar radiation (J m-2 day-1)}
  \item{G}{soil heat flux (J m-2 day-1)}
  \item{ra}{aerodynamic resistance (s m-1)}
  \item{rs}{surface resistance (s m-1)}
  \item{doy}{integer. Day of the year (between 1--365)}
  \item{latitude}{numeric. Latitude}
  \item{tmin}{numeric. minimum temperature (C)}
  \item{tmax}{numeric. maximum temperature (C)}
  \item{Pereira}{logical. If \code{TRUE}, the Pereira adjustment for photoperiod is used}
}

\value{
vector with evaporation values (mm)
}

\references{
Allen, R.G., L.S. Pereira, D. Raes and M. Smith, 1998. Crop evapotranspiration. Guidelines for computing crop water requirements. FAO Irrigation and drainage paper 56. FAO - Food and Agriculture Organization of the United Nations, Rome, 1998. (http://www.fao.org/docrep/x0490e/x0490e07.htm)

Thornthwaite, C.W., 1948. An approach toward a rational classification of climate. Geogr. Rev. 38:55-94.

Willmott, C.J., Rowe, C.M. and Mintz, Y., 1985. Climatology of the terrestrial seasonal water cycle. J. Climatol. 5:589-606.

Camargo, A.P., Marin, F.R., Sentelhas, P.C. and Picini, A.G., 1999. Adjust of the Thornthwaite's method to estimate the potential evapotranspiration for arid and superhumid climates, based on daily temperature amplitude. Rev. Bras. Agrometeorol. 7(2):251-257

Pereira, A.R. and W.O. Pruitt, 2004. Adaptation of the Thornthwaite scheme for estimating daily reference evapotranspiration. Agricultural Water Management 66: 251-257
}

\author{
Robert Hijmans, partly based on Python evapolib by Maarten J. Waterloo \url{http://python.hydrology-amsterdam.nl/}
}

\examples{
ET0_PenmanMonteith(21.67, 67, 1013, 14100000, 500000, 104, 70)
ET0_PriestleyTaylor(21.65, 67, 1013, 18200000, 600000)
ET0_Makkink(21.65, 67, 1013, 24200000)
}
