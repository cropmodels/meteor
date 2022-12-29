\name{from}

\alias{fromDate}
\alias{fromYear}

\title{ Date manipulation }

\description{
Helper functions for manipulation of dates, including conversion between (day of year) (DOY) to date and back, and extraction of parts of a date.
}

\usage{
fromDate(x, v)
fromYear(y, v)
}

\arguments{
  \item{x}{Date or POSIX}
  \item{y}{integer (year)}

  \item{v}{character. Output variable. With \code{fromDate} it can be one of: "year", "month", "day", "doy". With \code{fromYear} it can be one of "leap", "ndays", "days"}
}

\value{
integer
}

\examples{
d <- as.Date("1999-12-30") + 1:2
d
fromDate(d, "month")
fromDate(d, "doy")
}
