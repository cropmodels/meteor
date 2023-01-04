#include <cmath>
#include <vector>

const double kVal = 273.15;

// thermal radiation calculation
// The formula below can be found at https://en.wikipedia.org/wiki/Mean_radiant_temperature,


inline double thermalRadiance(double tas, double wind, double Tg) {
  //Tg is the temperature of the black globe, in C. It is calculated below in fwbgTg
	double temp1 = pow(Tg + kVal, 4.0);
	double temp2 = pow(wind, 0.6) * (Tg - tas);
	double temp3 = (temp1 + 2.5e+8 * temp2);
	double temp4 = pow(temp3, 0.25);
	double tr =  temp4 - kVal;
	return  tr;
}


inline double utci(const double &ta, const double &tg, const double &wnd, const double &hurs) {
  
// Ta   : air temperature, degree Celsius
// ehPa : water vapour pressure, hPa=hecto Pascal
// Tmrt : mean radiant temperature, degree Celsius
// wnd   : wind speed 10 m (RH or 2?)  above ground level in m/s
// tg   : black globe temperature in C. 
  
 
	double Tk = ta + kVal;
	double satVapPres = 1.004 * 6.1121 * exp(17.502 * ta / (Tk - 32.18));

//	double satVapPres = esat(Tk); // in hPa
	double pa = satVapPres * hurs/1000.; // vapor pressure in kPa
	double Tmrt = thermalRadiance(ta, wnd, tg);
	double dtm = Tmrt - ta;
	double va = wnd < 0.5 ? 0.5 : wnd > 17. ? 17. : wnd;
  
	double ta2 = ta*ta;
	double ta3 = ta2*ta;
	double ta4 = ta3*ta;
	double ta5 = ta4*ta;
	double ta6 = ta5*ta;
  
	double va2 = va*va;
	double va3 = va2*va;
	double va4 = va3*va;
	double va5 = va4*va;
	double va6 = va5*va;
  
	double dtm2 = dtm*dtm;
	double dtm3 = dtm2*dtm;
	double dtm4 = dtm3*dtm;
	double dtm5 = dtm4*dtm;
	double dtm6 = dtm5*dtm;
  
	double pa2 = pa*pa;
	double pa3 = pa2*pa;
	double pa4 = pa3*pa;
	double pa5 = pa4*pa;
	double pa6 = pa5*pa;
	
	double dtm2_pa = dtm2 * pa;
	double dtm2_pa2 = dtm2 * pa2;
	double dtm3_pa = dtm3 * pa;
	double dtm3_pa2 = dtm3 * pa2;
	double dtm3_pa3 = dtm3 * pa3;
	
	double ta2_va = ta2 * va;
	double ta2_va2 = ta2 * va2;
	double ta3_va = ta3 * va;
	double va_pa = va * pa;
	double va_dtm = va * dtm;
	double utci_out = ta + 
	  0.607562052 + 
    -0.0227712343 * ta + 0.000806470249 * ta2 + -0.000154271372 * ta3 + -3.24651735e-06 * ta4 + 7.32602852e-08 * ta5 + 1.35959073e-09 * ta6 + 
    -2.2583652 * va + 0.0880326035 * ta * va + 0.00216844454 * ta2_va + -1.53347087e-05 * ta3_va + -5.72983704e-07 * ta4 * va + -2.55090145e-09 * ta5 * va + -0.751269505 * va2 + 
    -0.00408350271 * ta * va2 + -5.21670675e-05 * ta2_va2 + 1.94544667e-06 * ta3 * va2 + 1.14099531e-08 * ta4 * va2 + 
    0.158137256 * va3 + -6.57263143e-05 * ta * va3 + 
    2.22697524e-07 * ta2 * va3 + -4.16117031e-08 * ta3 * va3 + 
    -0.0127762753 * va4 + 9.66891875e-06 * ta * va4 + 2.52785852e-09 * ta2 * va4 + 
    0.000456306672 * va5 + -1.74202546e-07 * ta * va5 + -5.91491269e-06 * va6 + 
    0.398374029 * dtm + 0.000183945314 * ta * dtm + -0.00017375451 * ta2 * dtm + -7.60781159e-07 * ta3 * dtm + 3.77830287e-08 * ta4 * dtm + 5.43079673e-10 * ta5 * dtm +
    -0.0200518269 * va_dtm + 0.000892859837 * ta * va_dtm + 3.45433048e-06 * ta2 * va_dtm + -3.77925774e-07 * ta3 * va_dtm + -1.69699377e-09 * ta4 * va_dtm + 
    0.000169992415 * va2 * dtm + -4.99204314e-05 * ta * va2 * dtm + 2.47417178e-07 * ta2_va2 * dtm + 1.07596466e-08 * ta3 * va2 * dtm + 8.49242932e-05 * va3 * dtm + 
    1.35191328e-06 * ta * va3 * dtm + -6.21531254e-09 * ta2 * va3 * dtm + 
    -4.99410301e-06 * va4 * dtm + -1.89489258e-08 * ta * va4 * dtm + 
    8.15300114e-08 * va5 * dtm + 
    0.00075504309 * dtm2 + 
    -5.65095215e-05 * ta * dtm2 + 
    -4.52166564e-07 * ta2 * dtm2 + 
    2.46688878e-08 * ta3 * dtm2 + 
    2.42674348e-10 * ta4 * dtm2 + 
    0.00015454725 * va * dtm2 + 
    5.2411097e-06 * ta * va * dtm2 + 
    -8.75874982e-08 * ta2_va * dtm2 + 
    -1.50743064e-09 * ta3_va * dtm2 + 
    -1.56236307e-05 * va2 * dtm2 + 
    -1.33895614e-07 * ta * va2 * dtm2 + 
    2.49709824e-09 * ta2_va2 * dtm2 + 
    6.51711721e-07 * va3 * dtm2 + 
    1.94960053e-09 * ta * va3 * dtm2 + 
    -1.00361113e-08 * va4 * dtm2 + 
    -1.21206673e-05 * dtm3 + -2.1820366e-07 * ta * dtm3 + 7.51269482e-09 * ta2 * dtm3 + 9.79063848e-11 * ta3 * dtm3 + 
    1.25006734e-06 * va * dtm3 + -1.81584736e-09 * ta * va * dtm3 +  -3.52197671e-10 * ta2_va * dtm3 +  -3.3651463e-08 * va2 * dtm3 + 1.35908359e-10 * ta * va2 * dtm3 +   4.1703262e-10 * va3 * dtm3 + 
    -1.30369025e-09 * dtm4 + 4.13908461e-10 * ta * dtm4 + 9.22652254e-12 * ta2 * dtm4 + -5.08220384e-09 * va * dtm4 +  -2.24730961e-11 * ta * va * dtm4 +  1.17139133e-10 * va2 * dtm4 + 
    6.62154879e-10 * dtm5 +  4.0386326e-13 * ta * dtm5 + 1.95087203e-12 * va * dtm5 + 
    -4.73602469e-12 * dtm6 + 
    5.12733497 * pa + 
    -0.312788561 * ta * pa + 
    -0.0196701861 * ta2 * pa + 
    0.00099969087 * ta3 * pa + 
    9.51738512e-06 * ta4 * pa + 
    -4.66426341e-07 * ta5 * pa + 
    0.548050612 * va_pa + 
    -0.00330552823 * ta * va_pa + 
    -0.0016411944 * ta2 * va_pa + 
    -5.16670694e-06 * ta3 * va_pa + 
    9.52692432e-07 * ta4 * va_pa + 
    -0.0429223622 * va2 * pa + 
    0.00500845667 * ta * va2 * pa + 
    1.00601257e-06 * ta2_va2 * pa + 
    -1.81748644e-06 * ta3 * va2 * pa + 
    -0.00125813502 * va3 * pa + 
    -0.000179330391 * ta * va3 * pa + 
    2.34994441e-06 * ta2 * va3 * pa + 
    0.000129735808 * va4 * pa + 
    1.2906487e-06 * ta * va4 * pa + 
    -2.28558686e-06 * va5 * pa + 
    -0.0369476348 * dtm * pa + 
    0.00162325322 * ta * dtm * pa + 
    -3.1427968e-05 * ta2 * dtm * pa + 
    2.59835559e-06 * ta3 * dtm * pa + 
    -4.77136523e-08 * ta4 * dtm * pa + 
    0.0086420339 * va_dtm * pa + 
    -0.000687405181 * ta * va_dtm * pa + 
    -9.13863872e-06 * ta2 * va_dtm * pa + 
    5.15916806e-07 * ta3 * va_dtm * pa + 
    -3.59217476e-05 * va2 * dtm * pa + 
    3.28696511e-05 * ta * va2 * dtm * pa + 
    -7.10542454e-07 * ta2_va2 * dtm * pa + 
    -1.243823e-05 * va3 * dtm * pa + 
    -7.385844e-09 * ta * va3 * dtm * pa + 
    2.20609296e-07 * va4 * dtm * pa + 
    -0.00073246918 * dtm2_pa + 
    -1.87381964e-05 * ta * dtm2_pa + 
    4.80925239e-06 * ta2 * dtm2_pa + 
    -8.7549204e-08 * ta3 * dtm2_pa + 
    2.7786293e-05 * va * dtm2_pa + 
    -5.06004592e-06 * ta * va * dtm2_pa + 
    1.14325367e-07 * ta2_va * dtm2_pa + 
    2.53016723e-06 * va2 * dtm2_pa + 
    -1.72857035e-08 * ta * va2 * dtm2_pa + 
    -3.95079398e-08 * va3 * dtm2_pa + 
    -3.59413173e-07 * dtm3_pa + 
    7.04388046e-07 * ta * dtm3_pa + 
    -1.89309167e-08 * ta2 * dtm3_pa + 
    -4.79768731e-07 * va * dtm3_pa + 
    7.96079978e-09 * ta * va * dtm3_pa + 
    1.62897058e-09 * va2 * dtm3_pa + 
    3.94367674e-08 * dtm4 * pa + 
    -1.18566247e-09 * ta * dtm4 * pa + 
    3.34678041e-10 * va * dtm4 * pa + 
    -1.15606447e-10 * dtm5 * pa + 
    -2.80626406 * pa2 + 
    0.548712484 * ta * pa2 + 
    -0.0039942841 * ta2 * pa2 + 
    -0.000954009191 * ta3 * pa2 + 
    1.93090978e-05 * ta4 * pa2 + 
    -0.308806365 * va * pa2 + 
    0.0116952364 * ta * va * pa2 + 
    0.000495271903 * ta2_va * pa2 + 
    -1.90710882e-05 * ta3_va * pa2 + 
    0.00210787756 * va2 * pa2 + 
    -0.000698445738 * ta * va2 * pa2 + 
    2.30109073e-05 * ta2_va2 * pa2 + 
    0.00041785659 * va3 * pa2 + 
    -1.27043871e-05 * ta * va3 * pa2 + 
    -3.04620472e-06 * va4 * pa2 + 
    0.0514507424 * dtm * pa2 + 
    -0.00432510997 * ta * dtm * pa2 + 
    8.99281156e-05 * ta2 * dtm * pa2 + 
    -7.14663943e-07 * ta3 * dtm * pa2 + 
    -0.000266016305 * va_dtm * pa2 + 
    0.000263789586 * ta * va_dtm * pa2 + 
    -7.01199003e-06 * ta2 * va_dtm * pa2 + 
    -0.000106823306 * va2 * dtm * pa2 + 
    3.61341136e-06 * ta * va2 * dtm * pa2 + 
    2.29748967e-07 * va3 * dtm * pa2 + 
    0.000304788893 * dtm2_pa2 + 
    -6.42070836e-05 * ta * dtm2_pa2 + 
    1.16257971e-06 * ta2 * dtm2_pa2 + 
    7.68023384e-06 * va * dtm2_pa2 + 
    -5.47446896e-07 * ta * va * dtm2_pa2 + 
    -3.5993791e-08 * va2 * dtm2_pa2 + 
    -4.36497725e-06 * dtm3_pa2 + 
    1.68737969e-07 * ta * dtm3_pa2 + 
    2.67489271e-08 * va * dtm3_pa2 + 
    3.23926897e-09 * dtm4 * pa2 + 
    -0.0353874123 * pa3 + 
    -0.22120119 * ta * pa3 + 
    0.0155126038 * ta2 * pa3 + 
    -0.000263917279 * ta3 * pa3 + 
    0.0453433455 * va * pa3 + 
    -0.00432943862 * ta * va * pa3 + 
    0.000145389826 * ta2_va * pa3 + 
    0.00021750861 * va2 * pa3 + 
    -6.66724702e-05 * ta * va2 * pa3 + 
    3.3321714e-05 * va3 * pa3 + 
    -0.00226921615 * dtm * pa3 + 
    0.000380261982 * ta * dtm * pa3 + 
    -5.45314314e-09 * ta2 * dtm * pa3 + 
    -0.000796355448 * va_dtm * pa3 + 
    2.53458034e-05 * ta * va_dtm * pa3 + 
    -6.31223658e-06 * va2 * dtm * pa3 + 
    0.000302122035 * dtm2 * pa3 + 
    -4.77403547e-06 * ta * dtm2 * pa3 + 
    1.73825715e-06 * va * dtm2 * pa3 + 
    -4.09087898e-07 * dtm3_pa3 + 
    0.614155345 * pa4 + 
    -0.0616755931 * ta * pa4 + 
    0.00133374846 * ta2 * pa4 + 
    0.00355375387 * va * pa4 + 
    -0.000513027851 * ta * va * pa4 + 
    0.000102449757 * va2 * pa4 + 
    -0.00148526421 * dtm * pa4 + 
    -4.11469183e-05 * ta * dtm * pa4 + 
    -6.80434415e-06 * va_dtm * pa4 + 
    -9.77675906e-06 * dtm2 * pa4 + 
    0.0882773108 * pa5 + 
    -0.00301859306 * ta * pa5 + 
    0.00104452989 * va * pa5 + 
    0.000247090539 * dtm * pa5 + 
    0.00148348065 * pa6;

	return utci_out;
}



// [[Rcpp::export(name = ".utci")]]
std::vector<double> v_utci(const std::vector<double> &ta, const std::vector<double> &tg, const std::vector<double> &va, const std::vector<double> &hurs) {
	size_t n = ta.size();
	std::vector<double> out;
	out.reserve(n);
	for (size_t i=0; i<n; i++) {
		out[i] = utci(ta[i], tg[i], va[i], hurs[i]);
	}
	return out;
}



// [[Rcpp::export(name = ".pwc_utci")]]
std::vector<double> pwc_utci(const std::vector<double>& utci) {
	 
	const double hn = 6.;
	const double level1 = 15.8;
	const double level2 = 35.6;
	const double level3 = 42.5;
	const double level4 = 50.8;
  
	const double l4minusl3 = level4 - level3;
	const double l3minusl2 = level3 - level2;
	const double l2minusl1 = level2 - level1;

	size_t n = utci.size();
	std::vector<double> out;
	out.reserve(n);
	
	for (size_t i=0; i<n; i++) {
		double pwc = 100. / (1. + pow(45.33 / utci[i], -4.30));
		if (utci[i] >= level4) {
			pwc -=  2. * hn + 4.86;
		} else if (utci[i] >= level3) {
			pwc += ((utci[i] - level3)/l4minusl3) * (-(2. * hn + 4.86)) + (-1. * (utci[i] - level4)/l4minusl3) * (-(1.1 * hn + 0.98));
		} else if (utci[i] >= level2) {
			pwc += ((utci[i] - level2)/l3minusl2) * (-(1.1 * hn + 0.98)) + (-1. * (utci[i] - level3)/l3minusl2) * (-(0.65 * hn + 1.3));
		} else if (utci[i] > level1) {
			pwc += ((utci[i] - level1)/l2minusl1) * (-(0.65 * hn + 1.3));
		} 
		if (pwc < 0) pwc = 0.; // pwc can't be less than 0
		
		out.push_back(pwc);
	}
	return out;
}

// [[Rcpp::export(name = ".pwc_wbgt")]]
std::vector<double> pwc_wbgt(const std::vector<double>& wbgt) {

	const double level1 = 12.6;
	const double level2 = 29.4;
	const double level3 = 33.4;
	const double level4 = 36.1;
  
	const double l4minusl3 = level4 - level3;
	const double l3minusl2 = level3 - level2;
	const double l2minusl1 = level2 - level1;

	size_t n = wbgt.size();
	std::vector<double> out;
	out.reserve(n);

//	const double hn = 6.;
	
	for (size_t i=0; i<n; i++) {
		if ((std::isnan(wbgt[i])) || (wbgt[i] < 2)) {
			out.push_back(NAN);
			continue;			
		}
		double pwc = 100. / (1. + pow(33.63 / wbgt[i], -6.33));	  
		if (wbgt[i] >= level4) {
//			pwc -=  2. * hn + 4.86;
			pwc -= 16.86;
			pwc = pwc < 0 ? 0 : pwc; 
		} else if (wbgt[i] >= level3) {
//			pwc += ((wbgt[i] - level3)/l4minusl3) * (-(2.* hn + 4.86)) + (-1. * (wbgt[i] - level4)/l4minusl3) * (-(1.1 * hn + 0.98));
			pwc += ((wbgt[i] - level3)/l4minusl3) * -16.86 + ((wbgt[i] - level4)/l4minusl3) * 7.58;
		} else if (wbgt[i] >= level2) {
//			pwc += ((wbgt[i] - level2)/l3minusl2) * (-(1.1 * hn + 0.98)) + (-1. * (wbgt[i] - level3)/l3minusl2) * (-(0.65 * hn + 1.3));
			pwc += ((wbgt[i] - level2)/l3minusl2) * -7.58 + ((wbgt[i] - level3)/l3minusl2) * 5.2;
		} else if (wbgt[i] > level1) {
//			pwc += ((wbgt[i] - level1)/l2minusl1) * (-(0.65 * hn + 1.3));
			pwc += ((wbgt[i] - level1)/l2minusl1) *  -5.2;
		}
		pwc = std::round(pwc*10) / 10;
		out.push_back(pwc);
	}
	return out;
}


