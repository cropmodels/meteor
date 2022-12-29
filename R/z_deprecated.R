# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1  January 2009

yearFromDate <- function(date) {
# date is a string like "2007-7-10"    YYYY-M-D
# to avoid date shifts because of your local time zone if date is a POSIX. 
	warning("this function will be removed. Use 'fromDate' instead")
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%Y"))
}


monthFromDate <- function(date) {
	warning("this function will be removed. Use 'fromDate' instead")
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%m"))
}

dayFromDate <- function(date) {
	warning("this function will be removed. Use 'fromDate' instead")
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%d"))
}

doyFromDate <- function(date) {
	warning("this function will be removed. Use 'fromDate' instead")
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%j"))
}

isLeapYear <- function(year) {
	warning("this function will be removed. Use 'fromYear' instead")
	year <- round(year)
    return( ((year %% 100 != 0) & (year %%4 ==0)) | (year %% 400==0) )
}

daysInYear <- function(year) {
	warning("this function will be removed. Use 'fromYear' instead")
	ifelse(isLeapYear(year), 366, 365)
}

daysOfYear <- function(year) {
	warning("this function will be removed. Use 'fromYear' instead")
	if (length(year) > 1) {
		stop('this function only accepts a single year as an argument')
	}
	firstday <- as.Date(paste(year, "-1-1", sep=""))
	lastday <- as.Date(paste(year, "-12-31", sep=""))
	d <- seq(firstday, to=lastday, by=1)
	return(d)
}	

dateFromDoy <- function(doy, year) {
	warning("this function will be removed. Use 'fromDoy' instead")
	year <- round(year)
	doy <- round(doy)
	return(as.Date(doy, origin=paste(year-1, "-12-31", sep='')))
}
