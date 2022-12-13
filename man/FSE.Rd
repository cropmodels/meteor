\name{FSE weather}

\alias{readFSEwth}
\alias{writeFSEwth}
\alias{example_weather} 


\title{
Read FSE formatted weather data
}

\description{
Read or write FSE formatted weather data
}

\usage{
readFSEwth(f)
writeFSEwth(w, country='AAA', station=1, lon=0, lat=0, elev=0,  path=".", ...)

example_weather()
}

\arguments{
\item{f}{character. filename}
\item{w}{data.frame with daily weather data. Must include the following variables: "date", "srad", "tmin", "tmax", "wind", "prec", "vapr". The data must be sorted by date in ascending order. "date" must be a \code{\link{Date}}, the other variables must be numeric}
\item{country}{character code for a country (up to three letters)}
\item{station}{positive integer. Station number for the country}
\item{lon}{numeric. Longitude of the weather station (not used by the models)}
\item{lat}{numeric. Latitude of the weather station}
\item{elev}{numeric. Elevation of the weather station}
\item{path}{character. Folder where you want to write the files. It must exist}
}

\value{
\code{readFSEwth}: data.frame

\code{writeFSEwth}: character (invisibly) with the filenames
}
