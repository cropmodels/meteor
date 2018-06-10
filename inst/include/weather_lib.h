using namespace std;
#include <cmath>
#include <vector>
#include <numeric>

#ifndef WEATHERLIB_H_
#define WEATHERLIB_H_



/*
Variable names
temp: air temperature (degrees C)
relh: relative humidity (%)
airpress: air pressure (kPa)  hPa, PA ??? (may vary--- need to standardize)
Rn: average daily net radiation [J m-2 d-1].
 */


double photoperiod(int doy, double latitude) {
// Forsythe, William C., Edward J. Rykiel Jr., Randal S. Stahl, Hsin-i Wu and Robert M. Schoolfield, 1995.
// Ecological Modeling 80: 87-95. A Model Comparison for Daylength as a Function of Latitude and Day of the Year.
	if ((latitude > 90) || (latitude < -90)) { return(-1);}
	double P = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860*(doy-186)))));
	double pi = 3.141592653589793238462643383279502884197169399375;
	double torad = pi / 180;
	latitude = latitude * torad;
	double a =  (sin(0.8333 * torad) + sin(latitude) * sin(P)) / (cos(latitude) * cos(P));
	a =  std::min(std::max(a, -1.), 1.);
	return (  24 - (24 / pi) * acos(a) );
}
 

 
//saturated vapor pressure (Pa)
// simple SVP
double ESimple(double temp) {
	return( 610.8 * exp(17.27 * temp / (temp + 237.3) ));
}

// more refined SVP (Pa)
double ES(double temp) {
/* Derived from python meteolib by Maarten J. Waterloo and J. Delsman
   For T < 0 C the saturation vapour pressure equation for ice is used accoring to Goff and Gratch (1946), whereas for T>=0 C that of Goff (1957) is used.

  References
  ----------
  - Goff, J.A.,and S. Gratch, Low-pressure properties of water from -160 to 212 F. Transactions of the American society of heating and
  ventilating engineers, p. 95-122, presented at the 52nd annual meeting of the American society of heating and ventilating engineers, New York, 1946.
  - Goff, J. A. Saturation pressure of water on the new Kelvin temperature scale, Transactions of the American society of heating and ventilating engineers, pp 347-354,
  presented at the semi-annual meeting of the American society of heating and ventilating engineers, Murray Bay, Quebec. Canada, 1957.

  Examples
  --------
  > ESc(30.0)
  4242.725994656632
 */

  double log_p;
  if (temp < 0) { // for ice
    log_p = - 9.09718 * (273.16 / (temp + 273.15) - 1.0) - 3.56654 *
         log10(273.16 / (temp + 273.15)) + 0.876793 * (1.0 - (temp + 273.15) / 273.16) + log10(6.1071);
  } else {  // for water
    log_p = 10.79574 * (1.0 - 273.16 / (temp + 273.15)) - 5.02800 *
        log10((temp + 273.15) / 273.16) + 1.50475E-4 * (1 - pow(10, (-8.2969 * ((temp + 273.15) / 273.16 - 1.0)))) +
        0.42873E-3 *  (pow(10, (4.76955 * (1.0 - 273.16 / (temp + 273.15)))) - 1) + 0.78614;

  }
  double es = pow(10, log_p);
  return es * 100;

}


// vapor pressure
double EA(double temp, double relh) {
    return ( ES(temp) *  relh / 100 );
}


double DELTA(double temp){
/* Derived from python meteolib by Maarten J. Waterloo and J. Delsman
   Compute the slope of the temperature - vapour pressure curve (Delta) from air temperatures.
   output: slope of saturated vapour curve (kPa K-1).
    References
    ----------
    Technical regulations 49, World Meteorological Organisation, 1984. Appendix A. 1-Ap-A-3.

    Examples
    --------
        > delta(30.0)
        .24334309166827094
*/
    double es = ES(temp) / 1000; // in kPa
    return ( es * 4098.0 / pow((temp + 237.3), 2) * 1000);
}




double CP(double temp, double relh, double airpress) {
/*  Derived from meteolib by Maarten J. Waterloo and J. Delsman
    Compute the specific heat of air [J kg-1 K-1].
    References:  R.G. Allen, L.S. Pereira, D. Raes and M. Smith (1998). Crop Evaporation Guidelines for computing crop water requirements,
    FAO - Food and Agriculture Organization of the United Nations.  Irrigation and drainage paper 56, Chapter 3. Rome, Italy. (http://www.fao.org/docrep/x0490e/x0490e07.htm)

    Examples
        > CP(25,60,101300)
        1014.0749457208065
*/
    // calculate vapour pressures
    double  eact = EA(temp, relh);
    return ( 0.24 * 4185.5 * (1 + 0.8 * (0.622 * eact / (airpress - eact))) ) ;
}



double LAMBDA(double temp){
/* Compute the latent heat of vapourisation from air temperature: lambda [J kg-1 K-1].
    References:   J. Bringfelt. Test of a forest evapotranspiration model. Meteorology and Climatology Reports 52, SMHI, NorrkÃ¶pping, Sweden, 1986.
    Examples
    --------
        > lambda(25)
        2440883.8804625
*/
    return (  4185.5 * (751.78 - 0.5655 * (temp + 273.15)) );
}


double GAMMA(double temp, double relh, double airpress ){
/*  Calculate the psychrometric constant gamma (Pa K-1).

    Reference: J. Bringfelt. Test of a forest evapotranspiration model. Meteorology and Climatology Reports 52, SMHI, NorrkÃ¶pping, Sweden, 1986.

    > gamma(10,50,101300)
    66.26343318657227
*/
    double spec = CP(temp, relh, airpress);
    double L = LAMBDA(temp);
    return( spec * airpress / (0.622 * L) );
 }



double pottemp(double temp, double relh, double airpress){
/* calculate the potential temperature air, theta, from air temperatures, relative humidity
    and air pressure. Reference pressure 1000 hPa.  Ouput: theta: potential air temperature [Celsius].

    > pottemp(5,45,101300)
        3.977415823848844
*/
    double spec = CP(temp, relh, airpress);
    return ((temp + 273.15) * pow((100000.0 / airpress), (287.0 / spec)) - 273.15);
}


double Rho(double temp, double relh, double airpress){
/*  calculate the density of air, rho, from air temperatures, relative humidity and air pressure.
    Returns: rho: air density data [kg m-3].

        > Rho(10,50,101300)
        1.2431927125520903
*/
    // actual vapour pressure
    double eact = EA(temp, relh);
    // Calculate density of air rho
    double rho = 1.201 * (290.0 * (airpress - 0.378 * eact)) / (1000.0 * (temp + 273.15)) / 100.0;
    return rho; // in kg/m3
}


std::vector<double> potrad_dl(int doy, double latitude, double S) {
/*
* Slightly adapted from Exercise 3.1, Seasonal course of
* radiation (L31.FST). Goudriaan and Van Laar, 1995.
*
* RJH y HY September 1995
*
*  DAYL    Photoperiod  hr
*  RDD     Daily global radiation in    kJ/m2/d
*  LAT     Latitude of the site         degree
*  ATMTR   Atmospheric transmissivity (on a clear day 0.7)}
*/

    const double pi = 3.141592653589793238463;
//	const double TC = 4.0;
//	const double P = 1.5;
//	const double ATMTR = 1;
//	const double RAD = pi/180;
//	const double SC = 1328;

	latitude = latitude * pi / 180;

	double SINLAT = sin(latitude);
    double COSLAT = cos(latitude);
    //Maximal sine of declination
    double SINDCM = sin(pi/180*23.45);
    //Sine and cosine of declination (Eqns 3.4, 3.5);
    double SINDEC = -SINDCM*cos(2*pi*(doy+10)/365);
    double COSDEC = sqrt(1-SINDEC*SINDEC);
    // The terms A and B according to Eqn 3.3;
    double A = SINLAT*SINDEC;
    double B = COSLAT*COSDEC;
    double C = A/B;
    if (C < -1) {
		C = -1;
	} else if (C > 1) {
		C = 1;
	}
 //  {Daylength according to Eqn 3.6; arcsin(c) = arctan(c/sqrt(c*c+1))}
 	double DAYL;
    if (B != 0) {
		DAYL = 12 * (1+(2/pi)* asin(C));
	} else {
		DAYL = 12;
	}
		//  {Integral of sine of solar height (Eqn 3.7)}
    double SININT = A * DAYL + (24 * B / pi) * cos((pi/2)*((DAYL/12)-1));
    //  { Daily total of radiation under a clear sky The number 3600.
    //  converts SININT from hours to seconds per day
    double RDPOT  = S * SININT * 3600;
	// RDPOT  = RDPOT  * ATMTR
    //double Ri = 0.312 * RDDPOT + 0.401 * RDDPOT * (Sunhrs / 100);
    //Ri = max(0, Ri);

	std::vector<double> out(2);
	out[0] = RDPOT;
	out[1] = DAYL;
    return out;
}



std::vector<double> sun_NR(int doy, double lat, double S) {
/*  Compute the maximum sunshine duration [h] and incoming radiation [MJ/day] at the top of the atmosphere from day of year and latitude.
    Parameters:
     - doy: day of year.
     - lat: latitude

    Returns:
    - N: (float, array) maximum sunshine hours [h].
    - Rext: (float, array) extraterrestrial radiation [J day-1].

    Notes
    -----
    Only valid for latitudes between 0 and 67 degrees (i.e. tropics and temperate zone).

    References
    ----------
    R.G. Allen, L.S. Pereira, D. Raes and M. Smith (1998). Crop Evaporation - Guidelines for computing crop water requirements,
    FAO - Food and Agriculture Organization of the United Nations. rrigation and drainage paper 56, Chapter 3. Rome, Italy.
    (http://www.fao.org/docrep/x0490e/x0490e07.htm)

    Examples
    --------
        >>> sun_NR(50,60)
        (9.1631820597268163, 9346987.824773483)
*/
    double pi = 3.14159265359;
    // Set solar constant [W/m2]
    //double S = 1367.0; // [W/m2]
    // Print warning if latitude is above 67 degrees
//   if abs(lat) > 67.{
//   print 'WARNING: Latitude outside range of application (0-67 degrees).\n)'
    // Convert latitude [degrees] to radians
    double latrad = lat * pi / 180.0;
    // calculate solar declination dt [radians]
    double dt = 0.409 * sin(2 * pi / 365 * doy - 1.39);
    // calculate sunset hour angle [radians]
    double ws = acos(-tan(latrad) * tan(dt));
    // Calculate sunshine duration N [h]
    double N = 24 / pi * ws;
    // Calculate day angle j [radians]
    double j = 2 * pi / 365.25 * doy;
    // Calculate relative distance to sun
    double dr = 1.0 + 0.03344 * cos(j - 0.048869);
    // Calculate Rext
    double Rext = S * 86400 / pi * dr * (ws * sin(latrad) * sin(dt) + sin(ws) * cos(latrad) * cos(dt));
	std::vector<double> out(2);
	out[0] = Rext;
	out[1] = N;
    return out;
}




//double windvec(double u, double D) {
/*
    Function to calculate the wind vector from time series of wind
    speed and direction.

    Parameters:
        - u: array of wind speeds [m s-1].
        - D: array of wind directions [degrees from North].

    Returns:
        - uv: Vector wind speed [m s-1].
        - Dv: Vector wind direction [degrees from North].

    Examples
    --------

        >>> u = array([[ 3.],[7.5],[2.1]])
        >>> D = array([[340],[356],[2]])
        >>> windvec(u,D)
        (4.162354202836905, array([ 353.2118882]))
        >>> uv, Dv = windvec(u,D)
        >>> uv
        4.162354202836905
        >>> Dv
        array([ 353.2118882])

*/
/*    ve = 0.0 // define east component of wind speed
    vn = 0.0 // define north component of wind speed
    D = D * math.pi / 180.0 // convert wind direction degrees to radians
    for i in range(0, len(u)){
        ve = ve + u[i] * math.sin(D[i]) // calculate sum east speed components
        vn = vn + u[i] * math.cos(D[i]) // calculate sum north speed components
    ve = - ve / len(u) // determine average east speed component
    vn = - vn / len(u) // determine average north speed component
    uv = math.sqrt(ve * ve + vn * vn) // calculate wind speed vector magnitude
    // Calculate wind speed vector direction
    vdir = arctan2(ve, vn)
    vdir = vdir * 180.0 / math.pi // Convert radians to degrees
    if vdir < 180:
        Dv = vdir + 180.0
    else:
        if vdir > 180.0:
            Dv = vdir - 180
        else:
            Dv = vdir
    return uv, Dv // uv in m/s, Dv in dgerees from North
}
*/




/*
From evaplib: Functions for calculation of potential and actual evaporation from meteorological data.
author = "Maarten J. Waterloo <maarten.waterloo@falw.vu.nl>" version 1.0 (November 2014)
*/

double Penman_E0(double temp, double relh, double airpress, double Rs, double Rext, double u, double alpha=0.25, double Z=0) {
/*  Calculate daily Penman (open) water evaporation estimates.
    Parameters:
        - temp: daily average air temperatures [Celsius].
        - relh: daily average relative humidity [%].
        - airpress: daily average air pressure data [Pa].
        - Rs: daily incoming solar radiation [J m-2 day-1].
        - Rext: daily extraterrestrial radiation [J m-2 day-1].
        - u: daily average wind speed at 2 m [m s-1].
        - alpha: albedo [-] set at 0.08 for open water by default.
        - Z: site elevation, default is 0 m a.s.l.

    Returns:
        - E0: Penman open water evaporation values [mm day-1].

    Notes
    -----
    Meteorological parameters measured at 2 m above the surface. Albedo
    alpha set by default at 0.08 for open water (Valiantzas, 2006).

    References
    ----------
    - H.L. Penman (1948). Natural evaporation from open water, bare soil and grass. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences 193: 120-145.
    - H.L. Penman (1956). Evaporation: An introductory survey. Netherlands Journal of Agricultural Science 4: 9-29.
    - J.D. Valiantzas (2006). Simplified versions for the Penman evaporation equation using routine weather data. J. Hydrology 331: 690-702.

    Examples
    --------
        // With single values and default albedo/elevation
        > E0(20.67,67.0,101300.0,22600000.,42000000.,3.2)
        6.6029208786994467
        // With albedo is 0.18 instead of default and default elevation
        > E0(20.67,67.0,101300.0,22600000.,42000000.,3.2,alpha=0.18)
        5.9664248091431968
        // With standard albedo and Z= 250.0 m
        > E0(20.67,67.0,101300.0,22600000.,42000000.,3.2,Z=250.0)
        6.6135588207586284
        // With albedo alpha = 0.18 and elevation Z = 1000 m a.s.l.
        > E0(20.67,67.0,101300.0,22600000.,42000000.,3.2,0.18,1000.)
        6.00814764682986

    */
    // Set constants
    double sigma = 4.903E-3; // Stefan Boltzmann constant J/m2/K4/d

    // Delta, gamma and lambda
    double delta = DELTA(temp); // [Pa/K]
    double gamma = GAMMA(temp, relh,airpress); // [Pa/K]
    double lambda = LAMBDA(temp); // [J/kg]

    // saturated and actual water vapour pressures
    double es = ES(temp); // [Pa]
    double ea = EA(temp,relh); // [Pa]

   // radiation components (J/m2/day)
    double Rns = (1.0 - alpha) * Rs; // Shortwave component [J/m2/d]
    double Rs0 = (0.75 + 2E-5 * Z) * Rext; // Calculate clear sky radiation Rs0
    double f = 1.35 * Rs / Rs0 - 0.35;
    double epsilom = 0.34 - 0.14 * sqrt(ea/1000);
    double Rnl = f * epsilom * sigma * pow(temp + 273.15, 4); // Longwave component [J/m2/d]
    double Rnet = Rns - Rnl; // Net radiation [J/m2/d]
    double Ea = (1 + 0.536 * u) * (es/1000 - ea/1000);
    double E0 = (delta/(delta+gamma) * Rnet/lambda + gamma/(delta+gamma) * 6430000*Ea/lambda);

    return E0;
}


double ET0pm(double temp, double relh, double airpress, double Rs, double Rext, double u, double Z=0.0){
/*  Daily Penman-Monteith reference evapotranspiration
    Parameters:
        - temp: daily average air temperatures [Celsius].
        - relh: daily average relative humidity values [%].
        - airpress: daily average air pressure data [hPa].
        - Rs: total incoming shortwave radiation [J m-2 day-1].
        - Rext: Incoming shortwave radiation at the top of the atmosphere [J m-2 day-1].
        - u: windspeed [m s-1].
        - Z: elevation [m], default is 0 m a.s.l.

    Returns:
        - ET0pm: Penman Monteith reference evaporation (short grass with optimum water supply) values [mm day-1].

    Notes
    -----
    Meteorological measuements standard at 2 m above soil surface.

    References
    ----------
    R.G. Allen, L.S. Pereira, D. Raes and M. Smith (1998). Crop evapotranspiration - Guidelines for computing crop water requirements -
    FAO Irrigation and drainage paper 56. FAO - Food and Agriculture Organization of the United Nations, Rome, 1998.
    (http://www.fao.org/docrep/x0490e/x0490e07.htm)

    Examples
    --------
    > ET0pm(20.67,67.0,101300.0,22600000.,42000000.,3.2)
      4.7235349721073039
*/
    // Set constants
    double albedo = 0.23; // short grass albedo
    double sigma = 4.903E-3; // Stefan Boltzmann constant J/m2/K4/d

    // Calculate Delta, gamma and lambda
    double delta = DELTA(temp); // [Pa/K]
    double gamma = GAMMA(temp, relh, airpress); // [Pa/K]
    double lambda = LAMBDA(temp); // [J/kg]

    // Calculate saturated and actual water vapour pressures
    double es = ES(temp); // [Pa]
    double ea = EA(temp,relh); // [Pa]

    double Rns = (1.0-albedo)*Rs; // Shortwave component [J/m2/d]
    // Calculate clear sky radiation Rs0
    double Rs0 = (0.75+2E-5*Z)*Rext; // Clear sky radiation [J/m2/d]
    double f = 1.35*Rs/Rs0-0.35;
    double epsilom = 0.34-0.14*sqrt(ea/1000);
    double Rnl = f*epsilom*sigma * pow(temp+273.15, 4); // Longwave component [J/m2/d]
    double Rnet = Rns-Rnl; // Net radiation [J/m2/d]
    double ET0pm = (delta / 1000. * Rnet / lambda+900. / (temp+273.16)*u*(es-ea)/1000 * gamma/1000)/(delta/1000.+gamma/1000*(1.+0.34*u));
    return ET0pm; // FAO reference evaporation [mm/day]
}


double Em(double temp, double relh, double airpress, double Rs){
/*  Function to calculate Makkink evaporation (in mm/day).

    The Makkink evaporation is a reference crop evaporation used in the
    Netherlands, which is combined with a crop factor to provide an
    estimate of actual crop evaporation.

    Parameters:
        - temp: daily average air temperatures [Celsius].
        - relh: daily average relative humidity values [%].
        - airpress: daily average air pressure data [Pa].
        - Rs: average daily incoming solar radiation [J m-2 day-1].

    Returns:
        - Em: Makkink evaporation values [mm day-1].

    Notes
    -----
    Meteorological measuements standard at 2 m above soil surface.

    References
    ----------
    H.A.R. de Bruin (1987). From Penman to Makkink, in Hooghart, C. (Ed.), Evaporation and Weather, 
	Proceedings and Information. Comm. Hydrological Research TNO, The Hague. pp. 5-30.

    Examples
    --------
    > Em(21.65,67.0,101300.,24200000.)
      4.503830479197991
*/
  // Calculate Delta and gamma constants
    double delta = DELTA(temp);
	// RH  -- hPa to Pa
	airpress = airpress * 100 ;
	// RH
    double gamma = GAMMA(temp,relh,airpress);
    double lambda = LAMBDA(temp);

    // calculate Em [mm/day]
    return( 0.65 * delta/(delta + gamma) * Rs / lambda ) ;
}


double Ept(double temp, double relh, double airpress, double Rn, double G) {
/* Daily Priestley - Taylor evaporation.
    Parameters:
        - Rn: average daily net radiation [J m-2 day-1].
        - G: average daily soil heat flux [J m-2 day-1].

    Returns:
        - Ept: Priestley Taylor evaporation values [mm day-1].

    References
    ----------
    Priestley, C.H.B. and R.J. Taylor, 1972. On the assessment of surface
    heat flux and evaporation using large-scale parameters. Mon. Weather Rev. 100:81-82.

    Examples
    --------
    > Ept(21.65,67.0,101300.,18200000.,600000.)
    6.349456116128078
*/
  // Calculate Delta and gamma constants
    double delta = DELTA(temp);
	// RH for consistency use hPA (mbar)
    airpress = airpress * 100; // [Pa]
	//
    double gamma = GAMMA(temp,relh,airpress);
    double lambda = LAMBDA(temp);
    // calculate Em [mm/day]
    return ( 1.26 * delta / (delta+gamma)*(Rn-G)/lambda );
}

double ra(double z, double z0, double d, double u){
/* aerodynamic resistance.
  Parameters:
        - z: measurement height [m].
        - z0: roughness length [m].
        - d: displacement length [m].
        - u: windspeed [m s-1].

    Returns:
        - ra: aerodynamic resistances [s m-1].

    References:
    A.S. Thom (1075), Momentum, mass and heat exchange of plant communities,
    In: Monteith, J.L. Vegetation and the Atmosphere, Academic Press, London.  p. 57–109.

    Examples:
        > ra(3,0.12,2.4,5.0)
        3.2378629924752942
*/
    return pow(log((z-d)/z0), 2)/(0.16*u);
}


double Epm(double temp, double relh, double airpress, double Rn, double G, double ra, double rs){
    /* Penman Monteith evaporation.
    Parameters:
        - Rn: average daily net radiation [J].
        - G: average daily soil heat flux [J].
        - ra: aerodynamic resistance [s m-1].
        - rs: surface resistance [s m-1].
    Returns:
        - Epm: Penman Monteith evaporation values [mm].
    References
    ----------
    J.L. Monteith (1965) Evaporation and environment. Symp. Soc. Exp. Biol. 19, 205-224.

    Examples
	--------
	>>> evaplib.Epm(21.67,67.0,1013.0,14100000.,500000.,104.,70.)
		3.2433411460494073
*/  // Calculate Delta, gamma and lambda
    double delta = DELTA(temp) / 100.; // [hPa/K]
    airpress = airpress * 100.; // [Pa]
    double gamma = GAMMA(temp,relh,airpress)/100.; // [hPa/K]
    double lambda = LAMBDA(temp); // [J/kg]
    double rho = Rho(temp, relh, airpress);
    double spec = CP(temp, relh, airpress);
    // Calculate saturated and actual water vapour pressures
    double es = ES(temp) / 100.; // [hPa]
    double ea = EA(temp,relh) / 100.; // [hPa]
    //Calculate Epm
	double Epm = ((delta * (Rn-G) + rho * spec *(es-ea)/ra)/(delta + gamma * (1. + rs/ra)))/lambda;
    return Epm; // actual ET in mm
}


double tvardry(double rho, double CP , double T , double sigma_t , double z, double d= 0.0){
/* 	Calculate the sensible heat flux from high frequency temperature measurements and its standard deviation.
	Parameters:
        - rho: air density values [kg m-3].
        - CP: specific heat at constant temperature values [J kg-1 K-1].
        - T: temperature data [Celsius].
        - sigma_t: standard deviation of temperature data [Celsius].
        - z: temperature measurement height above the surface [m].
        - d: displacement height due to vegetation, default is zero [m].

    Returns:
        - H: sensible heat flux [W m-2].

    Notes
    -----
    This function holds only for free convective conditions when C2*z/L >>1, where L is the Obhukov length.

    References
    ----------
    - J.E. Tillman (1972), The indirect determination of stability, heat and momentum fluxes in the atmosphere boundary layer from simple scalar
    variables during dry unstable conditions, Journal of Applied Meteorology 11: 783-792.
    - H.F. Vugts, M.J. Waterloo, F.J. Beekman, K.F.A. Frumau and L.A. Bruijnzeel. The temperature variance method: a powerful tool in the
    estimation of actual evaporation rates. In J. S. Gladwell, editor, Hydrology of Warm Humid Regions, Proc. of the Yokohama Symp., IAHS
    Publication No. 216, pages 251-260, July 1993.

    Examples
    --------
        >>> tvardry(1.25,1035.0,25.3,0.25,3.0)
        34.658669290185287
        >>> d=0.25
        >>> tvardry(1.25,1035.0,25.3,0.25,3.0,d)
        33.183149497185511
    */
    // Define constants
    double k = 0.40; // von Karman constant
    double g = 9.81; // acceleration due to gravity [m/s^2]
    double C1 =  2.9; // De Bruin et al., 1992
    double C2 = 28.4; // De Bruin et al., 1992
    // L= Obhukov-length [m]

    //Free Convection Limit
    double H = rho * CP * sqrt( pow(sigma_t/C1 , 3) * k * g * (z-d) / (T+273.15) * C2);
    //else:
    // including stability correction
    //zoverL = z/L
    //tvardry = rho * CP * sqrt((sigma_t/C1)**3 * k*g*(z-d) / (T+273.15) *
    //          (1-C2*z/L)/(-1*z/L))

    //Check if we get complex numbers (square root of negative value) and remove
    //I = find(zoL >= 0 | H.imag != 0);
    //H(I) = ones(size(I))*NaN;

    return H; // sensible heat flux
}

std::vector<double> gash79(double Pg, double ER, double S, double St, double p, double pt) {
    /* precipitation interception loss.
    Parameters:
        - Pg: daily rainfall [mm].
        - ER: evaporation percentage of total rainfall [mm h-1].
        - S: storage capacity canopy [mm].
        - St: stem storage capacity [mm].
        - p: direct throughfall [mm].
        - pt: stem precipitation [mm].

    Returns:
        - Ei: Interception [mm].
        - TF: through fall [mm].
        - SF: stemflow [mm].

    References
    ----------
    J.H.C. Gash, An analytical model of rainfall interception by forests,
    Quarterly Journal of the Royal Meteorological Society, 1979, 105, pp. 43-55.

    Examples
    --------
        >>> gash79(12.4,0.15,1.3,0.2,0.2,0.02)
        (12.4, 8.4778854123725971, 0, 3.9221145876274024)
        >>> gash79(60.0,0.15,1.3,0.2,0.2,0.02)
        (60.0, 47.033885412372598, 0, 12.966114587627404)
    */
        //PGsat calculation (for the saturation of the canopy)
        double PGsat = -(1/ER*S)* log((1-(ER/(1-p-pt))));

        //Set initial values to zero
        double Ecan= 0.;
        double Etrunk= 0.;

        // Calculate interception for different storm sizes
        if ((Pg<PGsat) & (Pg>0)) {
            Ecan=(1-p-pt) * Pg;
            if (Pg > St/pt) {
                Etrunk=St+pt*Pg;
	    }
 	} else if ((Pg > PGsat) & (Pg < St/pt)) {
            Ecan=((((1-p-pt)*PGsat)-S) + (ER*(Pg-PGsat)) + S);
            Etrunk=0.;
 	} else if ((Pg > PGsat) & (Pg > (St/pt))) {
            Ecan= ((((1-p-pt)*PGsat)-S)+ (ER*(Pg-PGsat)) + S+(St+pt*Pg));
            Etrunk = St+pt*Pg;
	}

	double Ei = Ecan + Etrunk;
    double TF = Pg-Ei;
    double SF=0;

	std::vector<double> out(3);
	out[0] = TF;
	out[1] = SF;
	out[2] = Ei;
  return out;
}



double EThornthwaiteWilmott(double temp, int doy, double latitude) {
	// Thornthwaite, C.W., 1948. An approach toward a rational classification of climate. Geogr. Rev. 38:55�94.
	// Willmott, C.J., Rowe, C.M. and Mintz, Y., 1985. Climatology of the terrestrial seasonal water cycle. J. Climatol. 5:589�606.
	temp = std::max(0., temp);
    double ET;
	if (temp > 26) {
		// Willmott for > 26C
		ET = -415.85 + 32.24*temp - 0.43 * pow(temp, 2) ;
	} else {
		double I = pow((temp/5), 1.514);
		double alpha = 6.75e-07 * pow(I, 3) - 7.71e-05 * pow(I, 2) + 1.79e-02 * I + 0.49239;
		ET = 16 * pow((10*temp / I), alpha) ;
	}
	double L = photoperiod(doy, latitude);
	ET = ET * L/(12 * 30);
	return(ET);
}


double EThornthwaiteWilmottCamargo(double tmin, double tmax, int doy, double latitude, bool Pereira=false) {
	// Camargo, A.P., Marin, F.R., Sentelhas, P.C. and Picini, A.G., 1999. Adjust of the Thornthwaite�s method to estimate the potential evapotranspiration for arid and superhumid climates, based on daily temperature amplitude. Rev. Bras. Agrometeorol. 7(2):251�257
	// Pereira, A.R. and W.O. Pruitt, 2004. Adaptation of the Thornthwaite scheme for estimating daily reference evapotranspiration. Agricultural Water Management 66: 251-257
	double temp = (tmax + tmin) / 2;
	// Camargo FFECTIVE TEMPERATURE
	double teff = 0.36 * (3 * tmax - tmin);
	// Pereira adjustment
	if (Pereira) {
		double L = photoperiod(doy, latitude);
		teff = teff * L / (24-L) ;
	}
	teff = std::max(temp, std::min(teff, tmax));
	return( EThornthwaiteWilmott(teff, doy, latitude));
}



std::vector<double> dailyToHourlyTemperature(double tmin, double tmax, int doy, double latitude) {
	double TC = 4.0;
    double P = 1.5;
	double pi = 3.141592653589793238462643383279502884197169399375;
	double daylength = photoperiod(doy, latitude);
	double nigthlength = 24 - daylength;
    double sunrise = 12 - 0.5 * daylength;
    double sunset = 12 + 0.5 * daylength;
	std::vector<double> tmp(24);
	double tsunst;

	for (int h = 0; h < 24; h++) {
		if ( h < sunrise)  {  // midnight to sunrise;
			tsunst = tmin + (tmax-tmin) * sin(pi * (daylength/(daylength + 2 * P)));
			tmp[h] = (tmin - tsunst * exp(-nigthlength/TC) + (tsunst-tmin) * exp(-(h + 24-sunset)/TC)) / (1-exp(-nigthlength/TC));
		} else if ( h < (12+P) ) { // between sunrise and time that tmax is reached
			tmp[h] = tmin + (tmax-tmin) * sin(pi * (h - sunrise)/(daylength + 2 * P));
		} else if (h < sunset) { 	// between time of tmax and sunset;
			tmp[h] = tmin + (tmax-tmin) * sin(pi * (h - sunrise)/(daylength + 2 * P));
		} else { // sunset to midnight;
			tsunst = tmin + (tmax - tmin) * sin(pi * (daylength/(daylength + 2 * P)));
			tmp[h] = (tmin - tsunst * exp(-nigthlength / TC) + (tsunst-tmin) * exp(-(h - sunset)/TC)) / (1-exp(-nigthlength/TC));
		}
	}
	return(tmp);
}


std::vector<double> dailyToHourlyRelhum(double relh, double tmin, double tmax, int doy, double latitude) {
	std::vector<double> tmp = dailyToHourlyTemperature(tmin, tmax, doy, latitude);
	double tavg = std::accumulate(tmp.begin(), tmp.end(), 0) / 24;
	double vp = EA(tavg, relh);
    std::vector<double> out(24);
	for (int i=0; i<24; i++) {
		out[i] = std::min(100., 100 * vp / ES(tmp[i]));
	}
	return(out);
}




#endif