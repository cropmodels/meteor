/*
#include <vector>
#include <string>
#include "Rcpp.h"


// new fwbgTemp function that reads in Tg rather than calculating it
void fwbgTemp(std::vector<double> &out,
              const Rcpp::NumericVector tas, const Rcpp::NumericVector hurs,
              const Rcpp::NumericVector wind, const Rcpp::NumericVector srad, const Rcpp::NumericVector Tg, const Rcpp::NumericVector lat_deg, 
              const Rcpp::NumericVector year, const Rcpp::NumericVector doy, const double& tolerance, const int& optim=2, std::string output="pwc_wbgt_out") {
  size_t n = tas.size();
  
  // potential values of output
	bool pwc_wbgt_out = output == "pwc_wbgt_out";
	bool pwc_utci_out = output == "pwc_utci_out";
	bool utci_out = output == "utci_out";
	bool wbgt_out = output == "wbgt_out";
	bool fwbgTnwb_out = output == "fwbgTnwb_out";
  
	for (size_t i=0; i<n; i++) {
		if (std::isnan(tas[i]) || std::isnan(hurs[i]) || std::isnan(wind[i]) || std::isnan(srad[i])) { 
		  out.push_back(NAN);
		  continue; 
		}
		double radiation = srad[i];
		double zenith_rad = calZenith(doy[i], year[i], lat_deg[i]); 
		double relh = hurs[i] * 0.01;
		double Tair = tas[i] + kVal;
		double eair = relh * esat(Tair);	
		double emis_atm_out = emis_atm(Tair, relh);
		double density = Pair * 100./(Tair * r_air); // 
		double visc = viscosity(Tair);
		double Tnwb;
		
		Tnwb = optim_fTnwb_2steps(Tair, hurs[i], wind[i], radiation, zenith_rad, visc, emis_atm_out, eair, density, tolerance);
		  
		if (pwc_wbgt_out) {
		  // if else below sets pwc = 100 when temp is below lower value of range defined in Foster, et al 2021, Table 3
			if (tas[i] < 12.) {
				out.push_back(100.); 
			} else {
				double wbgt = 0.7 * Tnwb + 0.2 * Tg[i] + 0.1 * tas[i];
				double u = pwc_wbgt(wbgt); 
				out.push_back(u);
			}
		} else if (pwc_utci_out) { 
			if (tas[i] < 12.) {
				out.push_back(100.); 
			// Foster, et al 2021, Table 3
			} else {
				double utci = utci(tas[i], Tg[i], wind[i], hurs[i]);
				double u = pwc_utci(utci);
				out.push_back(u); }
			}
		} else if (wbgt_out) { 
			double wbgt = 0.7 * Tnwb + 0.2 * Tg[i] + 0.1 * tas[i];
			double u = wbgt;
			out.push_back(u);
		} else if (utci_out) { 
			double utci = utci(tas[i], Tg[i], wind[i], hurs[i]);
			out.push_back(utci);
		} else if (fwbgTnwb_out) { 
			out.push_back(Tnwb);
		} else {
			Rcpp::Rcout << "invalid output " << output <<  std::endl;
		}
  }
}

//function called from R
std::vector<double> pwc_lapp(Rcpp::NumericMatrix tas, Rcpp::NumericMatrix hurs,
                             Rcpp::NumericMatrix wind, Rcpp::NumericMatrix srad, Rcpp::NumericMatrix Tg, Rcpp::NumericMatrix lat,
                             Rcpp::NumericMatrix year, Rcpp::NumericMatrix doy,
                             double tolerance, int optim=2, std::string output="pwc") {
  
  std::vector<double> out;
  // size_t ds = dates.size();
  size_t n = tas.ncol();
  //  size_t m = tas.nrow();
  //  size_t l = lat.size();
  out.reserve(tas.size());
  
  // now we pass "out" by reference, and have it filled by fwbgTemp  
  for (size_t i=0; i<n; i++) { 
    fwbgTemp(out, tas(Rcpp::_, i), hurs(Rcpp::_, i), wind(Rcpp::_, i), srad(Rcpp::_, i), Tg(Rcpp::_, i), lat(Rcpp::_, i), year(Rcpp::_, i), doy(Rcpp::_, i), tolerance, optim, output);   
    
  }
  return out;
}

*/
