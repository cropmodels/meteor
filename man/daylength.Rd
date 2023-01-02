\name{photoperiod}

\alias{photoperiod}
\alias{photoperiod,numeric-method}
\alias{photoperiod,Date-method}
\alias{photoperiod,data.frame-method}
\alias{photoperiod,SpatRaster-method}

\title{ photoperiod}

\description{
Compute photoperiod (daylength, sunshine duration) at a given latitude and day of the year.
}

\usage{
\S4method{photoperiod}{Date}(x, latitude)

\S4method{photoperiod}{data.frame}(x)

\S4method{photoperiod}{SpatRaster}(x, filename="", overwrite=FALSE, ...)
}

\arguments{
  \item{x}{Date, integer (day of the year), or data.frame (with variables "date" and "latitude", or SpatRaster}
  \item{latitude}{numeric. Latitude}
  \item{filename}{character. Output filename}
  \item{overwrite}{logical. If \code{TRUE}, \code{filename} is overwritten}  
  \item{...}{additional arguments for writing files as in \code{\link[terra]{writeRaster}}}     
}

\value{
double. Photoperiod in hours
}

\references{
Forsythe, W.C., E.J. Rykiel Jr., R.S. Stahl, H. Wu, R.M. Schoolfield, 1995. A model comparison for photoperiod as a function of latitude and day of the year. Ecological Modeling 80: 87-95.
}


\examples{
photoperiod(50, 52)
photoperiod(50, 5)
photoperiod(180, 55)

p <- photoperiod(1:365, 52)
d <- dateFromDoy(1:365, 2001)
plot(d, p)
}
