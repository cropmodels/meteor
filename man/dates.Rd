\name{dates}

\alias{dateFromDoy}
\alias{doyFromDate}
\alias{dayFromDate}
\alias{monthFromDate}
\alias{yearFromDate}
\alias{isLeapYear}
\alias{daysInYear}
\alias{daysOfYear}

\title{ date manipulation }

\description{
Helper functions for manipulation of dates, includding conversion between (day of year) (DOY) to date and back, and extraction of parts of a date.
}

\usage{
dateFromDoy(doy, year)
doyFromDate(date)
dayFromDate(date)
monthFromDate(date)
yearFromDate(date)
isLeapYear(year)
daysInYear(year)
}

\arguments{
  \item{doy}{integer. Day of the year (1..365) or (1..366) for leap years }
  \item{year}{integer. Year, e.g. 1982 }
  \item{date}{Date object or character formatted 'yyyy-mm-dd', e.g. '1982-11-23' }  
}

\value{
integer or Date
}

\examples{
doy <- 88
year <- 1970
date <- dateFromDoy(doy, year)
date
dateFromDoy(-15, 2000)
doyFromDate(date)
isLeapYear(2000)
daysInYear(2000)
daysInYear(1999)
}
