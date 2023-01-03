// authors: Gerald Nelson and Robert Hijmans
// based on code in the HeatStress package by Ana Casanueva

// RH in % and relh is a fraction

#include <vector>
#include <string>
#include <cmath>
#include "Rcpp.h"

// Physical constants
const double stefanb  =  0.000000056696;
const double cp  =  1003.5; // heat capacity of dry air at constant pressure 
const double m_air  =  28.97; 
const double m_h2o  =  18.015;
const double r_gas  =  8314.34;
const double r_air  =  r_gas / m_air;
const double  ratio  =  cp * m_air/ m_h2o;
const double  Pr  =  cp / (cp + (1.25 * r_air));

// Wick constants
const double emis_wick = 0.95; // emissivity;
const double  alb_wick = 0.4; // albedo;
const double diam_wick = 0.007; // diameter (in m)
const double len_wick = 0.0254; // length (in m)

// Globe constants
const double emis_globe = 0.95; // emissivity
const double   alb_globe = 0.05; // albedo
const double   diam_globe = 0.0508; // diameter (in m)

// Surface constants
const double   emis_sfc = 0.999;
const double SurfAlbedo = 0.4;
//const double   alb_sfc = SurfAlbedo;
const double propDirect = 0.8; //Assume a proportion of direct radiation = direct/(diffuse + direct)
const double Pair = 1010; // Atmospheric pressure in hPa

const double pi = 3.14159265358979323846;
//const double pi2 = 6.28318530717958647692;
const double min_speed = 0.1;
const double kVal = 273.15;

//diffusivity constants
const double pcrit13 = pow((36.4 * 218.), (1. / 3.));
const double  Tcrit512 = pow((132. * 647.3),  (5. / 12.));
const double  Tcrit12 = pow((132. * 647.3), 0.5);
const double Mmix = pow((1. / 28.97 + 1. / 18.015), 0.5);

const double toRad = pi/180.; 

inline double diffusivity(const double& Tk) { // saturation vapor pressure
	return( 0.000364 * pow((Tk / Tcrit12), 2.334) * pcrit13 * Tcrit512 * Mmix / (Pair / 1013.25) * 0.0001 );
}

inline double esat(const double& Tk) { // saturation vapor pressure 
  // units are hPa
  // Ta units C, esat return in C
  // Tk value of air temperature in Kelvin.
  // return saturation vapor pressure (hPa).
  //  return  exp(18.956 - 4030.18/(ta + 235.0));
  //saturation vapor pressure (hPa) over water.  
  // Reference: Buck's (1981) approximation (eqn 3) of Wexler's (1976) formula over liquid water.
	double satvp  = 1.004 * 6.1121 * exp(17.502 * (Tk - kVal) / (Tk - 32.18));
	return satvp;
}

inline double emis_atm(const double& Tk, const double& relh) { 
	double e = relh * esat(Tk);
	return 0.575 * pow(e, 0.143);
}


inline double Tdew(const double& ktas, const double& RH) {
  // Source: Lawrence, M. G. (2005). The Relationship between Relative Humidity and the Dewpoint Temperature in Moist Air: A Simple Conversion and Applications. Bulletin of the American Meteorological Society, 86(2), 225–234. doi: 10.1175/BAMS-86-2-225
  // constants for new dewpoint temperature
	double a1 = 17.625;
	double b1 = 243.04;
	double relh = RH/100;
	double logrh = log(relh);
	double dew = (b1 * (logrh + (a1 * ktas)/(b1 + ktas))) / (a1 - logrh - (a1 * ktas)/(b1 + ktas));
	return dew;
}


inline std::vector<double> seq(double start, double end, double increment) {
	std::vector<double> out;
	size_t s = floor((end - start) / increment);
	out.reserve(s);
	for (size_t i=0; i<=s; i++) {
		double val = start + i * increment;
		out.push_back(val);
	}
	return out;
}

inline double viscosity(const double& Tk) { 
  //viscosity of air, kg/(m s)
  // https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html?vA=290&units=K#
	double visc = (((Tk / 97.0) - 2.9) / 0.4 * (-0.034)) + 1.048;
	return 0.0000026693 * pow((28.97 * Tk), 0.5) / (pow(3.617, 2.0) * visc);
}

inline double thermal_cond(const double& viscosity) {
  //Thermal conductivity of air, W/(m K). Reference: BSL, page 257.
	return viscosity * (cp + 1.25 * r_air);
}

inline double h_evap(const double& Tk) {
	return (313.15 - Tk) / 30.0 * (-71100.0) + 2407300.0;
}

inline double h_cylinder_in_air(const double& Tk, const double& speed) {
	double density = Pair * 100. / (r_air * Tk); // 
	double okspeed = speed < min_speed ?  min_speed : speed;
  // Reynolds number,  different than the one in h_sphere_in_air
	double Re = okspeed * density * diam_wick / viscosity(Tk); 
  // Nusselt number, different than the one in h_sphere_in_air 
	double Nu = 0.281 * pow(Re, 0.6) * pow(Pr, 0.44); 
  //replaces 
	double thermal_con = thermal_cond(viscosity(Tk));
	return Nu * thermal_con / diam_wick;
}

inline double h_sphere_in_air(const double& Tk, const double& speed) { 
  //Convective heat transfer coefficient for flow around a sphere, W/(m2 K). 
  // Reference: Bird, Stewart, and Lightfoot (BSL), page 409
  
	double thermal_con = thermal_cond(viscosity(Tk));
	double density = Pair * 100. / (r_air * Tk); // 
	double okspeed = speed < min_speed ?  min_speed : speed;
  // Reynolds number for sphere
	double Re = okspeed * density * diam_globe /  viscosity(Tk);
  // Nusselt number for sphere
	double Nu = 2.0 + 0.6 * pow(Re, 0.5) * pow(Pr, 0.3333); 
	return thermal_con * Nu / diam_globe; 
}



inline double calZenith(const int& doy, const int& year_num,  double lat_deg) { 
  // return zenith in degrees
	const double DECL1 = 0.006918;
	const double DECL2 = 0.399912;
	const double DECL3 = 0.070257;
	const double DECL4 = 0.006758;
	const double DECL5 = 0.000907;
	const double DECL6 = 0.002697;
	const double DECL7 = 0.00148;
	
	const double utc_hour = 12.0;
	const double TimeOffset = 0.0;
	const double TrueSolarTime = (utc_hour * 60.0) + TimeOffset;
	double HaDeg = ((TrueSolarTime/4.0) - 180.0);
	// double HaRad = degToRad(HaDeg);
	double HaRad = toRad * HaDeg;
	
	//	double lat_rad = degToRad(lat_deg);
	double lat_rad = toRad * lat_deg;
	// get number of days in a year. Deals with leap years.
	int dpy = (year_num % 400 == 0 || year_num % 4 == 0) ? 366 : 365;
	
	//Evaluate the fractional year in radians 
	double Gamma = 2. * pi * (((doy * 1.0) - 1.0) + (utc_hour/24.0))/(dpy * 1.0); 
	
	double Decli = DECL1 - DECL2 * cos(Gamma) + DECL3 * sin(Gamma) - DECL4 * cos(2. * Gamma) 
			+ DECL5 * sin(2 * Gamma) -  DECL6 * cos(3. * Gamma) + DECL7 * sin(3. * Gamma); 
	double CosZen = (sin(lat_rad) * sin(Decli) + cos(lat_rad) * cos(Decli) * cos(HaRad));
	CosZen = CosZen > 1. ? 1. : (CosZen < -1. ? 1. : CosZen);
  
	return acos(CosZen); // returns in rads
}

// function to be minimized. Returns globe temperature in deg C.
// Tglobe_prev is the value of Tair over which the optimization occurs. The range is Tair-2, Tair+10

inline double fr_tg(const double &Tglobe_prev, const double &Tair, const double &hurs, const double &speed, 
                    const double &radiation, const double &zenith_rad, const double &viscosity_out, 
                    const double &emis_atm_out) {
	double cza = cos(zenith_rad); //# cosine of zenith angle
  // Tsfc is surface temperature; Tair is air temp.
  // Since we don't have separate values for these Liljegren, et al, set them equal. 
  // double Tsfc = Tair[i];
	double Tref_globe = 0.5 * (Tglobe_prev + Tair);
  //Convective heat transfer coefficient for flow around a sphere, W/(m2 K)
	double h_sphere = h_sphere_in_air(Tref_globe, speed);
	double Tglobe = pow((0.5 * (emis_atm_out * pow(Tair, 4.) + emis_sfc * 
                      pow(Tair, 4.)) - h_sphere / (emis_globe * stefanb) * 
                      (Tglobe_prev - Tair) + radiation / (2. * emis_globe * stefanb) * 
                      (1. - alb_globe) * (propDirect * (1. / (2. * cza) - 1.) + 1. + SurfAlbedo)), 0.25);
	return  fabs(Tglobe - Tglobe_prev); //fabs returns double abs
}

// function to be minimized for tnwb

inline double fr_tnwb(const double &Twb_prev, const double &Tair, const double &speed, 
                      double &radiation, const double &zenith_rad, const double &viscosity_out,
                      const double &emis_atm_out, const double &eair, const double &density) {
  const double irad = 1.0;
  double Tref_cylinder = 0.5 * (Twb_prev + Tair);
  double Fatm =  stefanb * emis_wick * (0.5 * (emis_atm_out * pow(Tair, 4.0) + emis_sfc *
                                        pow(Tair, 4.0)) - pow(Twb_prev, 4.0)) + (1.0 - alb_wick) * radiation * 
                                        ((1.0 - propDirect) * (1.0 + 0.25 * diam_wick/len_wick) +
                                        ((tan(zenith_rad)/3.1416) + 0.25 * diam_wick/len_wick) * propDirect + SurfAlbedo);
  double diff = diffusivity(Tref_cylinder);
  double Sc = viscosity_out / (density * diff);
  double h_cyl = h_cylinder_in_air(Twb_prev, speed);
  double ewick = esat(Twb_prev); 
  double evap = h_evap(Twb_prev);
  double Twb = Tair - evap / ratio * (ewick - eair) / (Pair - ewick) * pow((Pr / Sc), 0.56) + Fatm / h_cyl * irad;
  return fabs(Twb - Twb_prev);
}


inline double optim_Tnwb(const double &Tair, const double &hurs, const double &speed, 
	double &radiation, const double &zenith_rad, const double &viscosity, const double &emis_atm_out, 
	const double &eair, const double &density, const double &tolerance) {
  
	double tdew = Tdew(Tair - kVal, hurs) + kVal;
	std::vector<double> rng = seq(tdew-1.0, Tair+10.0, tolerance);
	size_t m = rng.size();
  
	double tst1 = fr_tnwb(rng[0], Tair, speed, radiation, zenith_rad, viscosity, emis_atm_out, eair, density);
	for (size_t i=10; i<m; i+=10) {		  
		double tst2 = fr_tnwb(rng[i], Tair, speed, radiation, zenith_rad, viscosity, emis_atm_out, eair, density);
		if (tst2 > tst1) {
			size_t off = i==10 ? 0 : i-20;
			double tstA = fr_tnwb(rng[off], Tair, speed, radiation, zenith_rad, viscosity, emis_atm_out, eair, density);
			off++;
			for (size_t j=off; j<i; j++) { 
        
				double tstB = fr_tnwb(rng[j], Tair, speed, radiation, zenith_rad, viscosity, emis_atm_out, eair, density);
				if (tstB > tstA) {
					return rng[j-1] - kVal;
				}
				tstA = tstB;
			}
			return rng[i-1] - kVal;
		}
		tst1 = tst2;
	}  
	return rng[m-1] - kVal;
}

double optim_Tg(const double &Tair, const double &hurs, const double &speed,  double &radiation, 
                        const double &zenith_rad, const double &viscosity, const double &emis_atm_out, const double &tolerance) {
  
	std::vector<double> rng = seq(Tair-2.0, Tair+10.0, tolerance);
	size_t m = rng.size(); 
	double tst1 = fr_tg(rng[0], Tair, hurs, speed, radiation, zenith_rad, viscosity, emis_atm_out); 
	for (size_t i=10; i<m; i+=10) {		  
		double tst2 = fr_tg(rng[i], Tair, hurs, speed, radiation, zenith_rad, viscosity, emis_atm_out); 
		if (tst2 > tst1) {
			size_t off = i==10 ? 0 : i-20;
			double tstA = fr_tg(rng[off], Tair, hurs, speed, radiation, zenith_rad, viscosity, emis_atm_out); 
			off++;
			for (size_t j=off; j<i; j++) { 
				double tstB = fr_tg(rng[j], Tair, hurs, speed, radiation, zenith_rad, viscosity, emis_atm_out); 
				if (tstB > tstA) {
					return rng[j-1] - kVal;
				}
				tstA = tstB;
			}
			return rng[i-1] - kVal;
		}
		tst1 = tst2;
	}  
	return rng[m-1] - kVal; //return is in C
}

// Tg output only


void fix_zrad(double &radiation, double &zenith_rad) {
    if ((radiation > 0.) & (zenith_rad > 1.57)) zenith_rad = 1.57; // 90°
    if ((radiation > 15.) & (zenith_rad > 1.54))  zenith_rad = 1.54; // 88°
    if ((radiation > 900.) & (zenith_rad > 1.52)) zenith_rad = 1.52; // 87°
	if ((radiation < 10.) & (zenith_rad == 1.57)) radiation = 0.;
}


// globe temperature Tg
// [[Rcpp::export(name = ".Tg1")]]
std::vector<double> Tg1(const Rcpp::NumericVector tas, const Rcpp::NumericVector hurs,
	const Rcpp::NumericVector wind, const Rcpp::NumericVector srad, 
	const Rcpp::NumericVector year, const Rcpp::NumericVector doy, double lat) {
	   
	size_t n = tas.size();
  
	const double& tolerance = 0.1;
  
	std::vector<double> out;
	out.reserve(n); 		
  
	for (size_t i=0; i<n; i++) {
		if (std::isnan(tas[i]) || std::isnan(hurs[i]) || std::isnan(wind[i]) || std::isnan(srad[i])) { 
			out.push_back(NAN);
			continue; 
		}
    
		double radiation = srad[i];
		double zenith_rad = calZenith(doy[i], year[i], lat); 
		//Fix up out-of bounds problems with zenith. 
		fix_zrad(radiation, zenith_rad);

		double relh = hurs[i] * 0.01;
		double Tair = tas[i] + kVal;
		double emis_atm_out = emis_atm(Tair, relh);   
		double visc = viscosity(Tair);
		
		double Tg = optim_Tg(Tair, hurs[i], wind[i], radiation, zenith_rad, visc, emis_atm_out, tolerance);
		out.push_back(Tg);
	}
	return out; // deg C
}



// [[Rcpp::export(name = ".Tg2")]]
std::vector<double> Tg2(Rcpp::NumericMatrix tas, Rcpp::NumericMatrix hurs,
         Rcpp::NumericMatrix wind, Rcpp::NumericMatrix srad, const Rcpp::NumericVector year, 
		 Rcpp::NumericVector doy, const Rcpp::NumericVector lat) {
			 
	std::vector<double> out;
	size_t n = lat.size();
	size_t rs = n * year.size();
	out.reserve(rs);
	for (size_t i=0; i<n; i++) {
		std::vector<double> x = Tg1(tas(i, Rcpp::_), hurs(i, Rcpp::_), wind(i, Rcpp::_), srad(i, Rcpp::_), year, doy, lat(i)); 
		out.insert(out.end(), x.begin(), x.end());
	}
	return out;
}


// wet bulb temperature (Tnwb) 
// [[Rcpp::export(name = ".Tnwb1")]]
std::vector<double> Tnwb1(const Rcpp::NumericVector tas, const Rcpp::NumericVector hurs,
     const Rcpp::NumericVector wind, const Rcpp::NumericVector srad, 
	 Rcpp::NumericVector year, Rcpp::NumericVector doy, double lat, bool natural=true) {
	
	const double tolerance=0.1;
	std::vector<double> out;
	size_t ds = doy.size();
	out.reserve(ds);

	for (size_t i=0; i<ds; i++) {
		if (std::isnan(tas[i]) || std::isnan(hurs[i]) || std::isnan(wind[i]) || std::isnan(srad[i])) {
			out.push_back(NAN);
			continue;
		}
		double radiation = srad[i];
		double zenith_rad = calZenith(doy[i], year[i], lat);  // zenith is in degrees
		double relh = hurs[i] * 0.01;
		double Tair = tas[i] + kVal;
		double eair = relh * esat(Tair);	
		double emis_atm_out = emis_atm(Tair, relh);
		
		//Fix out-of bounds problems with zenith
		fix_zrad(radiation, zenith_rad);
		
		double density = Pair * 100./(Tair * r_air); // 
		double visc = viscosity(Tair);
	
		double Tnwb = optim_Tnwb(Tair, hurs[i], wind[i], radiation, zenith_rad, visc, emis_atm_out, eair, density, tolerance);
		
		if (natural) {
			out.push_back(Tnwb);    
		} else {
			double Tg = optim_Tg(Tair, hurs[i], wind[i], radiation, zenith_rad, visc, emis_atm_out, tolerance);
			double wbgt = 0.7 * Tnwb + 0.2 * Tg + 0.1 * tas[i];
			out.push_back(wbgt);    
		}		
	}
	return out;
}


// [[Rcpp::export(name = ".Tnwb2")]]
std::vector<double> Tnwb2(const Rcpp::NumericMatrix tas, const Rcpp::NumericMatrix hurs,
    const Rcpp::NumericMatrix wind, const Rcpp::NumericMatrix srad, 
	const Rcpp::NumericVector year, const Rcpp::NumericVector doy, 
	const Rcpp::NumericVector lat, bool natural=true) {

	std::vector<double> out;
	size_t n = lat.size();
	size_t rs = n * year.size();
	out.reserve(rs);
	for (size_t i=0; i<n; i++) {
		std::vector<double> x = Tnwb1(tas(i, Rcpp::_), hurs(i, Rcpp::_), wind(i, Rcpp::_), srad(i, Rcpp::_),year, doy, lat(i), natural); 
		out.insert(out.end(), x.begin(), x.end());
	}
	return out;
}

