\name{dayTemp}

\alias{dayTemp}


\title{Estimate the temperature during the day}

\description{
Estimate the mean temperature during the day (between sunrise and sunset) from daily minimum and maximum temperature.

The function requires the day of the year and latitude to compute the photoperiod.
}


\usage{
dayTemp(tmin, tmax, doy, latitude)
}


\arguments{
  \item{tmin}{numeric. minimum temperature (any unit)}
  \item{tmax}{numeric. maximum temperature (any unit)}
  \item{doy}{integer. Day of the year (between 1 and 365)}
  \item{latitude}{numeric. Latitude}  
}

\value{
numeric
}

\examples{
dayTemp(c(20,22), c(28,34), c(150,151), 52)
}
