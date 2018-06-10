/*
C++ code by Robert Hijmans, 2016

This function generates years of daily rainfall data on the basis of the given long term average monthly rainfall
(rain) and number of rainy days (raindays) using the method by Geng et al. (1986).

The Markov chain uses the transitional probabilities of a wet day after a dry day. The natural clustering of rainy days can be
described by a Markov parameter value of 0.75. In the absence of clustering Markov=1.00.
The parameters alpha and beta of the gamma distribution function are derived from the mean rainfall per wet day.

Derived from FORTRAN code in WOFOST model
C++ implementation and modifications by Robert Hijmans, June 2016
modifications
- simulate n years
- use cpp random number generators
- adjust daily rainfall to assure smooth changes in long term averages between days
*/


using namespace std;
#include <random>

//int main() { return 0; }

std::vector<double> simmeteo_rain(std::vector<double> rain, std::vector<double> raindays, int years, double markov, unsigned seed) {

      std::default_random_engine generator(seed);
      std::uniform_real_distribution<> runiform(0, 1);

      double A = 2.16;
      double B = 1.83;
//  double MARKOV = max(0.5, min(markov, 1.));

      int startday, endday, ndays;
      double rainWD, rand, alpha, beta, Pwd, Pww;

//     last day of preceding month
      double monthdays[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
      double Mlimit = 0.999 * A / (0.999 * B - 1.0);

      std::vector<double> simRain(365 * years, 0.);


      int prevm, nextm, midday;
      double adjust[365];
      double raindif;
      double nd = 30.;

	  int raind = 0;
	  
      for (int m=0; m<12; m++) {

          if (raindays[m] <= 0) {
              raind = 0;
          } else  {
              startday = monthdays[m];
              endday = monthdays[m+1];
              ndays = endday - startday;
             // amount of rain per wet day
              rainWD = rain[m] / raindays[m];
              if (rainWD >= Mlimit)  {
              // regression equation is valid
                  beta = B * rainWD - A;
                  alpha = rainWD / beta;
              } else {
              // adaptation to small beta
                  alpha = 0.999;
                  beta = rainWD / alpha;
              }
              // conditional rainfall probabilities
              Pwd = markov * raindays[m] / ndays;
              Pww = (1.00 - markov) + Pwd;

              gamma_distribution<double> rgamma(alpha, beta);

// Robert Hijmans, June 2016
// adjustment to assure that the change in average daily rainfall between days is smooth
// without this, each day of a given month has the same long term average, with sharp breaks between months

              prevm = m-1;
              nextm = m+1;
              if (m == 0) {
                prevm = 11;
              } else if (m == 11) {
                nextm = 0;
              } else if (m == 1) {
                nd = 29.;
              }
              midday = startday + 0.5 * ndays;
              raindif = (rain[m] - rain[prevm]) / nd;
              for (int d = startday; d < midday; d++) {
                adjust[d] = (rain[m] - (midday-d) * raindif) / rain[m];
              }
              raindif = (rain[nextm] - rain[m]) / nd;
              for (int d=midday; d < endday; d++) {
                adjust[d] = (rain[m] + (d-midday) * raindif) / rain[m];
              }

// end adjustment code

              for (int y=0; y<years; y++) {
                  int yy = y*365;
                  for (int d = startday; d < endday; d++) {
                      rand = runiform(generator);
                      if ((raind == 0  &&  rand <= Pwd) || (raind == 1  &&  rand <= Pww)) {
                          raind = 1;
                          simRain[d+yy] = rgamma(generator) * adjust[d];
                      } else {  // no rain
                          raind = 0;
                      }
                  } // generation
              } // years
          } // raindays > 0
      } // months
      return(simRain);
}

