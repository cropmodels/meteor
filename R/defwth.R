# Author: Robert J. Hijmans
# License GPL3
# Version 1.0  December 2012

example_weather <- function() {
	f <- system.file("extdata/18419.rds", package="meteor")
	readRDS(f)
}
