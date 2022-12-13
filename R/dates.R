# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1  January 2009

yearFromDate <- function(date) {
# date is a string like "2007-7-10"    YYYY-M-D
# to avoid date shifts because of your local time zone if date is a POSIX. 
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%Y"))
}

weekFromDate <- function(date) {
	date <- as.character(date)
	strftime(date, format = "%V")
}

monthFromDate <- function(date) {
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%m"))
}

dayFromDate <- function(date) {
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%d"))
}

doyFromDate <- function(date) {
	date <- as.character(date)
	as.numeric(format(as.Date(date), "%j"))
}

isLeapYear <- function(year) {
	year <- round(year)
    return( ((year %% 100 != 0) & (year %%4 ==0)) | (year %% 400==0) )
}

daysInYear <- function(year) {
	ifelse(isLeapYear(year), 366, 365)
}

daysOfYear <- function(year) {
	if (length(year) > 1) {
		stop('this function only accepts a single year as an argument')
	}
	firstday <- as.Date(paste(year, "-1-1", sep=""))
	lastday <- as.Date(paste(year, "-12-31", sep=""))
	d <- seq(firstday, to=lastday, by=1)
	return(d)
}	

dateFromDoy <- function(doy, year) {
	year <- round(year)
	doy <- round(doy)
	return(as.Date(doy, origin=paste(year-1, "-12-31", sep='')))
}
