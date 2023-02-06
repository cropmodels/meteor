#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include "weather_lib.h"


#include "simmeteo.h"
#include <random>
#include <chrono>

//Rcpp::compileAttributes('e:/bitbucket/models/meteor/meteor')
//Rcpp::compileAttributes('/home/rhijmans/bitbucket/models/meteor/meteor')



// [[Rcpp::export(name = ".markov_rain")]]
NumericMatrix markov_rain(NumericVector rain, NumericVector rainydays, int years, double markov, unsigned seed) {

	years = max(1, years);
	std::vector<double> out = simmeteo_rain(Rcpp::as<std::vector<double> >(rain), Rcpp::as<std::vector<double> >(rainydays), years, markov, seed);

	NumericMatrix r(365, years);
	int size = 365 * years;

	for (int i = 0; i < size; i++) {
		r[i] = out[i];
	}

  return(r);
}



// [[Rcpp::export(name = ".hourlyFromDailyTemp")]]
NumericMatrix hourlyFromDailyTemp(NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector  latitude) {
	NumericMatrix out(tmin.size(), 24);
  std::vector<double>  d(24);
	for (int i=0; i < tmin.size(); i++) {
			d = dailyToHourlyTemperature(tmin[i], tmax[i], doy[i], latitude[i]);
			for (int j = 0; j<24; j++) {
    		out(i,j) = d[j];
			}
	}
	return(out);
}

// [[Rcpp::export(name =".hourlyFromDailyRH")]]
NumericMatrix hourlyFromDailyRH(NumericVector relh, NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector latitude) {
	NumericMatrix out(tmin.size(), 24);
	std::vector<double>  d(24);
	for (int i=0; i < tmin.size(); i++) {
			d = dailyToHourlyRelhum(relh[i], tmin[i], tmax[i], doy[i], latitude[i]);
			for (int j = 0; j<24; j++) {
    		out(i,j) = d[j];
			}
	}
	return(out);
}

// [[Rcpp::export(name =".daytimeTemperature")]]
NumericVector daytimeTemperature(NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector latitude) {
	NumericVector out(tmin.size());
	for (int i=0; i < tmin.size(); i++) {
		out[i] = dayTemperature(tmin[i], tmax[i], doy[i], latitude[i]);
	}
	return(out);
}



// [[Rcpp::export(name =".Photoperiod")]]
NumericVector Photoperiod(NumericVector doy, NumericVector latitude) {
	NumericVector out(doy.size());
	int d;
	for (int i=0; i < out.size(); i++) {
		if (std::isnan(doy[i]) || std::isnan(latitude[i])) {
			out[i] = NAN;
		} else {
			d = int(doy[i]) % 365;
			out[i] = photoperiod(d, latitude[i]);
		}
	}
	return(out);
}



//// ET /////
// [[Rcpp::export(name =".ET0_ThornthwaiteWilmott")]]
NumericVector ET0_ThornthwaiteWilmott(NumericVector temp, NumericVector doy, NumericVector latitude) {
	NumericVector out(temp.size());
	for (int i=0; i < out.size(); i++) {
		out[i] = EThornthwaiteWilmott(temp[i], doy[i], latitude[i]);
	}
	return(out);
}


// [[Rcpp::export(name =".ET0_ThornthwaiteWilmottCamargo")]]
NumericVector ET0_ThornthwaiteWilmottCamargo(NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector latitude, bool Pereira) {
		NumericVector out(tmin.size());
		for (int i=0; i < out.size(); i++) {
			out[i] = EThornthwaiteWilmottCamargo(tmin[i], tmax[i], doy[i], latitude[i], Pereira);
		}
		return(out);
}



// [[Rcpp::export(name =".ET0_PenmanMonteith")]]
NumericVector ET0_PenmanMonteith(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rn, NumericVector G, NumericVector ra, NumericVector rs) {
		NumericVector out(temp.size());
		for (int i=0; i < out.size(); i++) {
			out[i] = Epm(temp[i], relh[i], atmp[i], Rn[i], G[i], ra[i], rs[i]);
		}
		return(out);
}


// [[Rcpp::export(name =".ET0_PriestleyTaylor")]]
NumericVector ET0_PriestleyTaylor(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rn, NumericVector G) {
	NumericVector out(temp.size());
	for (int i=0; i < out.size(); i++) {
		out[i] = Ept(temp[i], relh[i], atmp[i], Rn[i], G[i]);
	}
	return(out);
}


// [[Rcpp::export(name =".ET0_Makkink")]]
NumericVector ET0_Makkink(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rs){
	NumericVector out(temp.size());
	for (int i=0; i < out.size(); i++) {
		out[i] = Em(temp[i], relh[i], atmp[i], Rs[i]);
	}
	return(out);
}

// [[Rcpp::export(name =".E_Penman")]]
NumericVector E_Penman(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rs, NumericVector Rext, NumericVector u, NumericVector alpha, NumericVector Z) {
	NumericVector out(temp.size());
	for (int i=0; i < out.size(); i++) {
		out[i] = Penman_E0(temp[i], relh[i], atmp[i], Rs[i], Rext[i], u[i], alpha[i], Z[i]);
	}
	return(out);
}


/////// vapor pressure

// [[Rcpp::export(name =".SVP")]]
NumericVector SVP(NumericVector temp) {
//saturated vapor pressure (Pa)
	NumericVector out(temp.size());
	for (int i=0; i < out.size(); i++) {
		out[i] = ES(temp[i]) ;
	}
	return(out);
}

// [[Rcpp::export(name =".VP")]]
NumericVector VP(NumericVector temp, NumericVector relh) {
// actual vapor pressure (Pa)
	NumericVector out(temp.size());
	double svp;
	for (int i=0; i < out.size(); i++) {
		svp = ES(temp[i]);
		out[i] = svp * relh[i] / 100;
	}
	return(out);
}

// [[Rcpp::export(name =".VPD")]]
NumericVector VPD(NumericVector temp, NumericVector relh){
// vapor pressure deficit (Pa)
	NumericVector out(temp.size());
	double svp;
	for (int i=0; i < out.size(); i++) {
		svp = ES(temp[i]);
		out[i] = svp - ( svp * relh[i] / 100 );
	}
	return(out);
}


///////////////
// radiation
// [[Rcpp::export(name =".ExtraTerrestrialRadiation")]]
NumericMatrix ExtraTerrestrialRadiation(NumericVector doy, NumericVector latitude, double sc=1367.7, bool FAO=false){
	NumericMatrix out(doy.size(), 2);
	colnames(out) = CharacterVector::create("Radiation", "Photoperiod");
	std::vector<double>  x;
	if (FAO) {
		for (int i=0; i < doy.size(); i++) {
			x = sun_NR(doy[i], latitude[i], sc);
			out(i,0) = x[0];
			out(i,1) = x[1];
		}
	} else {
		for (int i=0; i < doy.size(); i++) {
			x = potrad_dl(doy[i], latitude[i], sc);
			out(i,0) = x[0];
			out(i,1) = x[1];
		}
	}
	return(out);
}
