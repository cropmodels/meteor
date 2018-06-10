# R.J. Hijmans
# Licence GPL v3


.simRain <- function(rain, rainydays, years, markov=0.75, seed) {

	stopifnot(length(rain)==12)
	stopifnot(all(rain >= 0))
	stopifnot(length(rainydays)==12)
	stopifnot(all(rainydays >= 0))

	years <- as.integer(max(1, years))
	markov <- max(0.5, min(markov, 1))

	if (missing(seed)) {
		seed <- sample(.Machine$integer.max, 1)
	} else {
		seed <- as.integer(max(0, seed))
	}

	.markov_rain(rain, as.integer(rainydays), years, markov, seed = seed)

}
