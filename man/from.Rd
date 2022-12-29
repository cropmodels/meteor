\name{from}

\alias{fromDate}
\alias{fromYear}
\alias{fromDoy}

\title{ Date manipulation }

\description{
Helper functions for extracting information from dates. Or to create a Date from a day number.
}

\usage{
fromDate(x, v)
fromYear(y, v)
fromDoy(doy, y)
}

\arguments{
  \item{x}{Date or POSIX}
  \item{y}{integer (year)}
  \item{v}{character. Output variable. With \code{fromDate} it can be one of: "year", "month", "week", "day", "doy". With \code{fromYear} it can be one of "leap", "ndays", "days"}
  \item{doy}{integer (day of the year)}

}

\value{
integer or Date
}

\examples{
d <- as.Date("1999-12-30") + 1:2
d
fromDate(d, "month")
fromDate(d, "doy")

fromDoy(10, 2000)
}
