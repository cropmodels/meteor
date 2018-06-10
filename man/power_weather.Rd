\name{power_weather}

\alias{power_weather}

\title{
Global weather data estimated from satellite data and models. 
}

\description{
This functions returns a data.frame with weather data from the NASA POWER database. It has the date, incoming solar radiation (srad, kJ m-2 day-1), minimum temperature (tmin, degrees C) and maximum temperature (tmax, degrees C), vapor pressure (vapr, Pa), precipitation (prec, mm), and windspeed (wind, m/s)

The data are from 1983-01-01 to 2016-12-31

Missing values for radiation (Jan to June 1983 and ...) and a few inbetween were replaced by the long term averages.

There are no precipitation values before 1997-01-01. Missing values for precipitaion after that date were estimated as the long term average (i.e., not a particularly good method). 

The data are at 1 degree spatial resolution. That is, they are the average for a large grid cell. 

These are estimates. They can give a good general impression, but they are not ground observations. 

Data are downloaded as-needed by tile. By default to folder called "power" in your working directory. 
}

\usage{
power_weather(lon, lat, folder=file.path(getwd(), 'power'), tiles=FALSE, ...)
}

\arguments{
\item{lon}{numeric}
\item{lat}{numeric}
\item{folder}{character}
\item{tiles}{logical. Download by tile?}
\item{...}{additional arguments}
}

\value{
data.frame
}


\examples{
\dontrun{
w <- power_weather(5, 50)
w$srad <- w$srad * 1000
wth <- subset(w, date > as.Date('2012-01-01'))
head(wth)
}
}
