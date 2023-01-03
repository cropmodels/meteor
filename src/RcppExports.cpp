// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bcppvars
std::vector<double> bcppvars(std::vector<double> prec, std::vector<double> tmin, std::vector<double> tmax);
RcppExport SEXP _meteor_bcppvars(SEXP precSEXP, SEXP tminSEXP, SEXP tmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type prec(precSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type tmin(tminSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type tmax(tmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(bcppvars(prec, tmin, tmax));
    return rcpp_result_gen;
END_RCPP
}
// v_utci
std::vector<double> v_utci(const std::vector<double>& ta, const std::vector<double>& tg, const std::vector<double>& va, const std::vector<double>& hurs);
RcppExport SEXP _meteor_v_utci(SEXP taSEXP, SEXP tgSEXP, SEXP vaSEXP, SEXP hursSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type ta(taSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type tg(tgSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type va(vaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type hurs(hursSEXP);
    rcpp_result_gen = Rcpp::wrap(v_utci(ta, tg, va, hurs));
    return rcpp_result_gen;
END_RCPP
}
// pwc_utci
std::vector<double> pwc_utci(const std::vector<double>& utci);
RcppExport SEXP _meteor_pwc_utci(SEXP utciSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type utci(utciSEXP);
    rcpp_result_gen = Rcpp::wrap(pwc_utci(utci));
    return rcpp_result_gen;
END_RCPP
}
// pwc_wbgt
std::vector<double> pwc_wbgt(const std::vector<double>& wbgt);
RcppExport SEXP _meteor_pwc_wbgt(SEXP wbgtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type wbgt(wbgtSEXP);
    rcpp_result_gen = Rcpp::wrap(pwc_wbgt(wbgt));
    return rcpp_result_gen;
END_RCPP
}
// Tg1
std::vector<double> Tg1(const Rcpp::NumericVector tas, const Rcpp::NumericVector hurs, const Rcpp::NumericVector wind, const Rcpp::NumericVector srad, const Rcpp::NumericVector year, const Rcpp::NumericVector doy, double lat);
RcppExport SEXP _meteor_Tg1(SEXP tasSEXP, SEXP hursSEXP, SEXP windSEXP, SEXP sradSEXP, SEXP yearSEXP, SEXP doySEXP, SEXP latSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type tas(tasSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type hurs(hursSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type wind(windSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type srad(sradSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type year(yearSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< double >::type lat(latSEXP);
    rcpp_result_gen = Rcpp::wrap(Tg1(tas, hurs, wind, srad, year, doy, lat));
    return rcpp_result_gen;
END_RCPP
}
// Tg2
std::vector<double> Tg2(Rcpp::NumericMatrix tas, Rcpp::NumericMatrix hurs, Rcpp::NumericMatrix wind, Rcpp::NumericMatrix srad, const Rcpp::NumericVector year, Rcpp::NumericVector doy, const Rcpp::NumericVector lat);
RcppExport SEXP _meteor_Tg2(SEXP tasSEXP, SEXP hursSEXP, SEXP windSEXP, SEXP sradSEXP, SEXP yearSEXP, SEXP doySEXP, SEXP latSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tas(tasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type hurs(hursSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type wind(windSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type srad(sradSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type year(yearSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type lat(latSEXP);
    rcpp_result_gen = Rcpp::wrap(Tg2(tas, hurs, wind, srad, year, doy, lat));
    return rcpp_result_gen;
END_RCPP
}
// Tnwb1
std::vector<double> Tnwb1(const Rcpp::NumericVector tas, const Rcpp::NumericVector hurs, const Rcpp::NumericVector wind, const Rcpp::NumericVector srad, const Rcpp::NumericVector year, const Rcpp::NumericVector doy, const double& lat, const bool& kelvin, const bool& natural);
RcppExport SEXP _meteor_Tnwb1(SEXP tasSEXP, SEXP hursSEXP, SEXP windSEXP, SEXP sradSEXP, SEXP yearSEXP, SEXP doySEXP, SEXP latSEXP, SEXP kelvinSEXP, SEXP naturalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type tas(tasSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type hurs(hursSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type wind(windSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type srad(sradSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type year(yearSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< const double& >::type lat(latSEXP);
    Rcpp::traits::input_parameter< const bool& >::type kelvin(kelvinSEXP);
    Rcpp::traits::input_parameter< const bool& >::type natural(naturalSEXP);
    rcpp_result_gen = Rcpp::wrap(Tnwb1(tas, hurs, wind, srad, year, doy, lat, kelvin, natural));
    return rcpp_result_gen;
END_RCPP
}
// Tnwb2
std::vector<double> Tnwb2(const Rcpp::NumericMatrix tas, const Rcpp::NumericMatrix hurs, const Rcpp::NumericMatrix wind, const Rcpp::NumericMatrix srad, const Rcpp::NumericVector year, const Rcpp::NumericVector doy, const Rcpp::NumericVector lat, const bool kelvin, const bool natural);
RcppExport SEXP _meteor_Tnwb2(SEXP tasSEXP, SEXP hursSEXP, SEXP windSEXP, SEXP sradSEXP, SEXP yearSEXP, SEXP doySEXP, SEXP latSEXP, SEXP kelvinSEXP, SEXP naturalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type tas(tasSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type hurs(hursSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type wind(windSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type srad(sradSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type year(yearSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type lat(latSEXP);
    Rcpp::traits::input_parameter< const bool >::type kelvin(kelvinSEXP);
    Rcpp::traits::input_parameter< const bool >::type natural(naturalSEXP);
    rcpp_result_gen = Rcpp::wrap(Tnwb2(tas, hurs, wind, srad, year, doy, lat, kelvin, natural));
    return rcpp_result_gen;
END_RCPP
}
// markov_rain
NumericMatrix markov_rain(NumericVector rain, NumericVector rainydays, int years, double markov, unsigned seed);
RcppExport SEXP _meteor_markov_rain(SEXP rainSEXP, SEXP rainydaysSEXP, SEXP yearsSEXP, SEXP markovSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rain(rainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rainydays(rainydaysSEXP);
    Rcpp::traits::input_parameter< int >::type years(yearsSEXP);
    Rcpp::traits::input_parameter< double >::type markov(markovSEXP);
    Rcpp::traits::input_parameter< unsigned >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(markov_rain(rain, rainydays, years, markov, seed));
    return rcpp_result_gen;
END_RCPP
}
// hourlyFromDailyTemp
NumericMatrix hourlyFromDailyTemp(NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector latitude);
RcppExport SEXP _meteor_hourlyFromDailyTemp(SEXP tminSEXP, SEXP tmaxSEXP, SEXP doySEXP, SEXP latitudeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tmin(tminSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    rcpp_result_gen = Rcpp::wrap(hourlyFromDailyTemp(tmin, tmax, doy, latitude));
    return rcpp_result_gen;
END_RCPP
}
// hourlyFromDailyRH
NumericMatrix hourlyFromDailyRH(NumericVector relh, NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector latitude);
RcppExport SEXP _meteor_hourlyFromDailyRH(SEXP relhSEXP, SEXP tminSEXP, SEXP tmaxSEXP, SEXP doySEXP, SEXP latitudeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmin(tminSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    rcpp_result_gen = Rcpp::wrap(hourlyFromDailyRH(relh, tmin, tmax, doy, latitude));
    return rcpp_result_gen;
END_RCPP
}
// Photoperiod
NumericVector Photoperiod(NumericVector doy, NumericVector latitude);
RcppExport SEXP _meteor_Photoperiod(SEXP doySEXP, SEXP latitudeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    rcpp_result_gen = Rcpp::wrap(Photoperiod(doy, latitude));
    return rcpp_result_gen;
END_RCPP
}
// ET0_ThornthwaiteWilmott
NumericVector ET0_ThornthwaiteWilmott(NumericVector temp, NumericVector doy, NumericVector latitude);
RcppExport SEXP _meteor_ET0_ThornthwaiteWilmott(SEXP tempSEXP, SEXP doySEXP, SEXP latitudeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    rcpp_result_gen = Rcpp::wrap(ET0_ThornthwaiteWilmott(temp, doy, latitude));
    return rcpp_result_gen;
END_RCPP
}
// ET0_ThornthwaiteWilmottCamargo
NumericVector ET0_ThornthwaiteWilmottCamargo(NumericVector tmin, NumericVector tmax, NumericVector doy, NumericVector latitude, bool Pereira);
RcppExport SEXP _meteor_ET0_ThornthwaiteWilmottCamargo(SEXP tminSEXP, SEXP tmaxSEXP, SEXP doySEXP, SEXP latitudeSEXP, SEXP PereiraSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tmin(tminSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tmax(tmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    Rcpp::traits::input_parameter< bool >::type Pereira(PereiraSEXP);
    rcpp_result_gen = Rcpp::wrap(ET0_ThornthwaiteWilmottCamargo(tmin, tmax, doy, latitude, Pereira));
    return rcpp_result_gen;
END_RCPP
}
// ET0_PenmanMonteith
NumericVector ET0_PenmanMonteith(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rn, NumericVector G, NumericVector ra, NumericVector rs);
RcppExport SEXP _meteor_ET0_PenmanMonteith(SEXP tempSEXP, SEXP relhSEXP, SEXP atmpSEXP, SEXP RnSEXP, SEXP GSEXP, SEXP raSEXP, SEXP rsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type atmp(atmpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rn(RnSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ra(raSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rs(rsSEXP);
    rcpp_result_gen = Rcpp::wrap(ET0_PenmanMonteith(temp, relh, atmp, Rn, G, ra, rs));
    return rcpp_result_gen;
END_RCPP
}
// ET0_PriestleyTaylor
NumericVector ET0_PriestleyTaylor(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rn, NumericVector G);
RcppExport SEXP _meteor_ET0_PriestleyTaylor(SEXP tempSEXP, SEXP relhSEXP, SEXP atmpSEXP, SEXP RnSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type atmp(atmpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rn(RnSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(ET0_PriestleyTaylor(temp, relh, atmp, Rn, G));
    return rcpp_result_gen;
END_RCPP
}
// ET0_Makkink
NumericVector ET0_Makkink(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rs);
RcppExport SEXP _meteor_ET0_Makkink(SEXP tempSEXP, SEXP relhSEXP, SEXP atmpSEXP, SEXP RsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type atmp(atmpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rs(RsSEXP);
    rcpp_result_gen = Rcpp::wrap(ET0_Makkink(temp, relh, atmp, Rs));
    return rcpp_result_gen;
END_RCPP
}
// E_Penman
NumericVector E_Penman(NumericVector temp, NumericVector relh, NumericVector atmp, NumericVector Rs, NumericVector Rext, NumericVector u, NumericVector alpha, NumericVector Z);
RcppExport SEXP _meteor_E_Penman(SEXP tempSEXP, SEXP relhSEXP, SEXP atmpSEXP, SEXP RsSEXP, SEXP RextSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type atmp(atmpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rs(RsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rext(RextSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(E_Penman(temp, relh, atmp, Rs, Rext, u, alpha, Z));
    return rcpp_result_gen;
END_RCPP
}
// SVP
NumericVector SVP(NumericVector temp);
RcppExport SEXP _meteor_SVP(SEXP tempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    rcpp_result_gen = Rcpp::wrap(SVP(temp));
    return rcpp_result_gen;
END_RCPP
}
// VP
NumericVector VP(NumericVector temp, NumericVector relh);
RcppExport SEXP _meteor_VP(SEXP tempSEXP, SEXP relhSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    rcpp_result_gen = Rcpp::wrap(VP(temp, relh));
    return rcpp_result_gen;
END_RCPP
}
// VPD
NumericVector VPD(NumericVector temp, NumericVector relh);
RcppExport SEXP _meteor_VPD(SEXP tempSEXP, SEXP relhSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type temp(tempSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type relh(relhSEXP);
    rcpp_result_gen = Rcpp::wrap(VPD(temp, relh));
    return rcpp_result_gen;
END_RCPP
}
// ExtraTerrestrialRadiation
NumericMatrix ExtraTerrestrialRadiation(NumericVector doy, NumericVector latitude, double sc, bool FAO);
RcppExport SEXP _meteor_ExtraTerrestrialRadiation(SEXP doySEXP, SEXP latitudeSEXP, SEXP scSEXP, SEXP FAOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type doy(doySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latitude(latitudeSEXP);
    Rcpp::traits::input_parameter< double >::type sc(scSEXP);
    Rcpp::traits::input_parameter< bool >::type FAO(FAOSEXP);
    rcpp_result_gen = Rcpp::wrap(ExtraTerrestrialRadiation(doy, latitude, sc, FAO));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_meteor_bcppvars", (DL_FUNC) &_meteor_bcppvars, 3},
    {"_meteor_v_utci", (DL_FUNC) &_meteor_v_utci, 4},
    {"_meteor_pwc_utci", (DL_FUNC) &_meteor_pwc_utci, 1},
    {"_meteor_pwc_wbgt", (DL_FUNC) &_meteor_pwc_wbgt, 1},
    {"_meteor_Tg1", (DL_FUNC) &_meteor_Tg1, 7},
    {"_meteor_Tg2", (DL_FUNC) &_meteor_Tg2, 7},
    {"_meteor_Tnwb1", (DL_FUNC) &_meteor_Tnwb1, 9},
    {"_meteor_Tnwb2", (DL_FUNC) &_meteor_Tnwb2, 9},
    {"_meteor_markov_rain", (DL_FUNC) &_meteor_markov_rain, 5},
    {"_meteor_hourlyFromDailyTemp", (DL_FUNC) &_meteor_hourlyFromDailyTemp, 4},
    {"_meteor_hourlyFromDailyRH", (DL_FUNC) &_meteor_hourlyFromDailyRH, 5},
    {"_meteor_Photoperiod", (DL_FUNC) &_meteor_Photoperiod, 2},
    {"_meteor_ET0_ThornthwaiteWilmott", (DL_FUNC) &_meteor_ET0_ThornthwaiteWilmott, 3},
    {"_meteor_ET0_ThornthwaiteWilmottCamargo", (DL_FUNC) &_meteor_ET0_ThornthwaiteWilmottCamargo, 5},
    {"_meteor_ET0_PenmanMonteith", (DL_FUNC) &_meteor_ET0_PenmanMonteith, 7},
    {"_meteor_ET0_PriestleyTaylor", (DL_FUNC) &_meteor_ET0_PriestleyTaylor, 5},
    {"_meteor_ET0_Makkink", (DL_FUNC) &_meteor_ET0_Makkink, 4},
    {"_meteor_E_Penman", (DL_FUNC) &_meteor_E_Penman, 8},
    {"_meteor_SVP", (DL_FUNC) &_meteor_SVP, 1},
    {"_meteor_VP", (DL_FUNC) &_meteor_VP, 2},
    {"_meteor_VPD", (DL_FUNC) &_meteor_VPD, 2},
    {"_meteor_ExtraTerrestrialRadiation", (DL_FUNC) &_meteor_ExtraTerrestrialRadiation, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_meteor(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
