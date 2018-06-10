\name{photoperiod}

\alias{photoperiod}

\title{ photoperiod}

\description{
Compute photoperiod (daylength, sunshine duration) at a given latitude and day of the year.
}

\usage{
photoperiod(doy, latitude)
}

\arguments{
  \item{doy}{integer. Day of the year. Integer between 1 and 365}
  \item{latitude}{numeric. Latitude, in degrees, should be between -90 and 90 }
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
