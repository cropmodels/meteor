#include <Rcpp.h>


double sd(const std::vector<double> &x) {
  double m = accumulate(x.begin(), x.end(), 0.0) / 12; 
  double s = 0;
  for (size_t i=0; i<x.size(); i++) {
        double d = (x[i] - m);
	      s += (d * d);
  }
  return sqrt(s / 11);
}

double sdm(const std::vector<double> &x, const double& m) {
  double s = 0;
  for (size_t i=0; i<x.size(); i++) {
    double d = (x[i] - m);
    s += (d * d);
  }
  return sqrt(s / 11);
}


double raincv(std::vector<double> x) {
//# coefficient of variation (expressed as a percentage)
// the "+1" is to avoid strange CVs for areas where mean rainfaill is < 1)
	for (size_t i=0; i<12; i++) x[i]++;
	double m = accumulate(x.begin(), x.end(), 0.0) / 12; 
	return(100 * sdm(x, m) / m);
}

std::vector<double> window_sum(std::vector<double> x)  { 
	x.insert(x.end(), x.begin(), x.begin()+3);
	for (size_t i=0; i<12; i++) {
		for (size_t j=1; j<3; j++) {
			x[i] += x[i+j];
		}
	}
	x.resize(12);
	return x;
}

std::vector<double> window_mean(std::vector<double> x)  { 
  x.insert(x.end(), x.begin(), x.begin()+3);
  for (size_t i=0; i<12; i++) {
    for (size_t j=1; j<3; j++) {
      x[i] += x[i+j];
    }
    x[i] /= 3;
  }
  x.resize(12);
  return x;
}


// [[Rcpp::export(name = ".bcppvars")]]
std::vector<double> bcppvars(std::vector<double> prec, std::vector<double> tmin, std::vector<double> tmax) {

	std::vector<double> p(19);

  std::vector<double> tavg(12), trng(12);
	for (size_t i=0; i<12; i++) {
		tavg[i] = (tmin[i] + tmax[i]) / 2;
		trng[i] = tmax[i] - tmin[i];
	}
	
// P1. Annual Mean Temperature 
	p[0] = accumulate(tavg.begin(), tavg.end(), 0.0) / 12; 
// P2. Mean Diurnal Range(Mean(period max-min)) 
	p[1] = accumulate(trng.begin(), trng.end(), 0.0) / 12; 
// P4. Temperature Seasonality (standard deviation) 
	p[3] = 100 * sd(tavg);
// P5. Max Temperature of Warmest Period 
	p[4] = *max_element(tmax.begin(), tmax.end());
// P6. Min Temperature of Coldest Period 
	p[5] = *min_element(tmin.begin(), tmin.end());
// P7. Temperature Annual Range (P5-P6) 
	p[6] = p[4] - p[5];
// P3. Isothermality (P2 / P7) 
	p[2] = 100 * p[1] / p[6];
// P12. Annual Precipitation 
	p[11] = accumulate(prec.begin(), prec.end(), 0.0);
// P13. Precipitation of Wettest Period 
  p[12] = *max_element(prec.begin(), prec.end());
// P14. Precipitation of Driest Period 
	p[13] = *min_element(prec.begin(), prec.end());
	 // P15. Precipitation Seasonality(Coefficient of Variation) 
	p[14] = raincv(prec);

// precip by quarter (3 months)		
  std::vector<double> rain = window_sum(prec);
	std::vector<double>::iterator maxrain, minrain;
	
// P16. Precipitation of Wettest Quarter 
	maxrain = max_element(rain.begin(), rain.end());
	p[15] = *maxrain;
// P17. Precipitation of Driest Quarter 
	minrain = min_element(rain.begin(), rain.end());
	p[16] = *minrain;

// P8. Mean Temperature of Wettest Quarter 
	size_t wetqrt = std::distance(rain.begin(), maxrain);
	wetqrt = wetqrt == 11 ? 0 : wetqrt+1;
	p[7] = tavg[wetqrt];
// P9. Mean Temperature of Driest Quarter 
	size_t dryqrt = std::distance(rain.begin(), minrain);
	dryqrt = dryqrt == 11 ? 0 : dryqrt+1;
		p[8] = tavg[dryqrt];

	// P10 Mean Temperature of Warmest Quarter 
	std::vector<double> tmp = window_mean(tavg);
	std::vector<double>::iterator maxtmp, mintmp;
	maxtmp = max_element(tmp.begin(), tmp.end());
	p[9] = *maxtmp;

// P11 Mean Temperature of Coldest Quarter
	mintmp = min_element(tmp.begin(), tmp.end());
	p[10] = *mintmp;

// P18. Precipitation of Warmest Quarter 
	size_t warmqrt = std::distance(tmp.begin(), maxtmp);
	Rcpp::Rcout << warmqrt << std::endl;
	p[17] = rain[warmqrt];

// P19. Precipitation of Coldest Quarter 
	size_t coldqrt = std::distance(tmp.begin(), mintmp);
	p[18] = rain[coldqrt];

	return(p);
}



/*** R
tmin <- c(10,12,14,16,18,20,22,21,19,17,15,12)
tmax <- tmin + 5
tavg <- (tmin + tmax) / 2
prec <- c(3,2,10,30,80,160,80,20,40,60,20,10)
round(dismo::biovars(prec, tmin, tmax), 1)
round(bcppvars(prec, tmin, tmax), 1)
round(raster::movingFun(prec, 3, sum, circular=T), 1)
*/
