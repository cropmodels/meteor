\name{dailyToHourly}

\alias{hourlyFromDailyTemp}
\alias{hourlyFromDailyRelh}


\title{Estimate hourly values from daily values }

\description{
Estimate hourly temperature from daily minimum and maximum temperature, or hourly relative humidity from average relative humidity and minimum and maximum temperature.

The functions require the day of the year and latitude to compute the photoperiod.
}


\usage{
hourlyFromDailyTemp(tmin, tmax, doy, latitude)
hourlyFromDailyRelh(relh, tmin, tmax, doy, latitude)
}


\arguments{
  \item{tmin}{numeric. minimum temperature (must be in C for hourlyFromDailyRelh)}
  \item{tmax}{numeric. maximum temperature (must be in C for hourlyFromDailyRelh)}
  \item{relh}{relative humidity (percent)}
  \item{doy}{integer. Day of the year (between 1 and 365)}
  \item{latitude}{numeric. Latitude}
}

\value{
matrix
}

\examples{
hourlyFromDailyTemp(c(20,22), c(28,34), c(150,151), 52)
hourlyFromDailyRelh(80, c(20,22), c(28,34), c(150,151), 52)
}
