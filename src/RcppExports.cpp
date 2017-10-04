// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dist_mtom
Rcpp::NumericMatrix dist_mtom(const Rcpp::NumericVector& xlon, const Rcpp::NumericVector& xlat, const Rcpp::NumericVector& ylon, const Rcpp::NumericVector& ylat, std::string dist_function);
RcppExport SEXP _distRcpp_dist_mtom(SEXP xlonSEXP, SEXP xlatSEXP, SEXP ylonSEXP, SEXP ylatSEXP, SEXP dist_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type xlon(xlonSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type xlat(xlatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ylon(ylonSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ylat(ylatSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_mtom(xlon, xlat, ylon, ylat, dist_function));
    return rcpp_result_gen;
END_RCPP
}
// dist_df
Rcpp::NumericVector dist_df(const Rcpp::NumericVector& xlon, const Rcpp::NumericVector& xlat, const Rcpp::NumericVector& ylon, const Rcpp::NumericVector& ylat, std::string dist_function);
RcppExport SEXP _distRcpp_dist_df(SEXP xlonSEXP, SEXP xlatSEXP, SEXP ylonSEXP, SEXP ylatSEXP, SEXP dist_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type xlon(xlonSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type xlat(xlatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ylon(ylonSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ylat(ylatSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_df(xlon, xlat, ylon, ylat, dist_function));
    return rcpp_result_gen;
END_RCPP
}
// dist_1tom
Rcpp::NumericVector dist_1tom(const double& xlon, const double& xlat, const Rcpp::NumericVector& ylon, const Rcpp::NumericVector& ylat, std::string dist_function);
RcppExport SEXP _distRcpp_dist_1tom(SEXP xlonSEXP, SEXP xlatSEXP, SEXP ylonSEXP, SEXP ylatSEXP, SEXP dist_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type xlon(xlonSEXP);
    Rcpp::traits::input_parameter< const double& >::type xlat(xlatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ylon(ylonSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ylat(ylatSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_1tom(xlon, xlat, ylon, ylat, dist_function));
    return rcpp_result_gen;
END_RCPP
}
// dist_1to1
double dist_1to1(const double& xlon, const double& xlat, const double& ylon, const double& ylat, std::string dist_function);
RcppExport SEXP _distRcpp_dist_1to1(SEXP xlonSEXP, SEXP xlatSEXP, SEXP ylonSEXP, SEXP ylatSEXP, SEXP dist_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type xlon(xlonSEXP);
    Rcpp::traits::input_parameter< const double& >::type xlat(xlatSEXP);
    Rcpp::traits::input_parameter< const double& >::type ylon(ylonSEXP);
    Rcpp::traits::input_parameter< const double& >::type ylat(ylatSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_1to1(xlon, xlat, ylon, ylat, dist_function));
    return rcpp_result_gen;
END_RCPP
}
// popdist_weighted_mean
Rcpp::DataFrame popdist_weighted_mean(Rcpp::DataFrame x_df, Rcpp::DataFrame y_df, std::string measure_col, std::string x_id, std::string x_lon_col, std::string x_lat_col, std::string y_lon_col, std::string y_lat_col, std::string pop_col, std::string dist_function, std::string dist_transform, double decay);
RcppExport SEXP _distRcpp_popdist_weighted_mean(SEXP x_dfSEXP, SEXP y_dfSEXP, SEXP measure_colSEXP, SEXP x_idSEXP, SEXP x_lon_colSEXP, SEXP x_lat_colSEXP, SEXP y_lon_colSEXP, SEXP y_lat_colSEXP, SEXP pop_colSEXP, SEXP dist_functionSEXP, SEXP dist_transformSEXP, SEXP decaySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type x_df(x_dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type y_df(y_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type measure_col(measure_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_id(x_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lon_col(x_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lat_col(x_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lon_col(y_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lat_col(y_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type pop_col(pop_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_transform(dist_transformSEXP);
    Rcpp::traits::input_parameter< double >::type decay(decaySEXP);
    rcpp_result_gen = Rcpp::wrap(popdist_weighted_mean(x_df, y_df, measure_col, x_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, pop_col, dist_function, dist_transform, decay));
    return rcpp_result_gen;
END_RCPP
}
// dist_weighted_mean
Rcpp::DataFrame dist_weighted_mean(Rcpp::DataFrame x_df, Rcpp::DataFrame y_df, std::string measure_col, std::string x_id, std::string x_lon_col, std::string x_lat_col, std::string y_lon_col, std::string y_lat_col, std::string dist_function, std::string dist_transform, double decay);
RcppExport SEXP _distRcpp_dist_weighted_mean(SEXP x_dfSEXP, SEXP y_dfSEXP, SEXP measure_colSEXP, SEXP x_idSEXP, SEXP x_lon_colSEXP, SEXP x_lat_colSEXP, SEXP y_lon_colSEXP, SEXP y_lat_colSEXP, SEXP dist_functionSEXP, SEXP dist_transformSEXP, SEXP decaySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type x_df(x_dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type y_df(y_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type measure_col(measure_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_id(x_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lon_col(x_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lat_col(x_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lon_col(y_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lat_col(y_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_transform(dist_transformSEXP);
    Rcpp::traits::input_parameter< double >::type decay(decaySEXP);
    rcpp_result_gen = Rcpp::wrap(dist_weighted_mean(x_df, y_df, measure_col, x_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function, dist_transform, decay));
    return rcpp_result_gen;
END_RCPP
}
// dist_min
Rcpp::DataFrame dist_min(Rcpp::DataFrame x_df, Rcpp::DataFrame y_df, std::string x_id, std::string y_id, std::string x_lon_col, std::string x_lat_col, std::string y_lon_col, std::string y_lat_col, std::string dist_function);
RcppExport SEXP _distRcpp_dist_min(SEXP x_dfSEXP, SEXP y_dfSEXP, SEXP x_idSEXP, SEXP y_idSEXP, SEXP x_lon_colSEXP, SEXP x_lat_colSEXP, SEXP y_lon_colSEXP, SEXP y_lat_colSEXP, SEXP dist_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type x_df(x_dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type y_df(y_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_id(x_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_id(y_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lon_col(x_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lat_col(x_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lon_col(y_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lat_col(y_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_min(x_df, y_df, x_id, y_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function));
    return rcpp_result_gen;
END_RCPP
}
// dist_max
Rcpp::DataFrame dist_max(Rcpp::DataFrame x_df, Rcpp::DataFrame y_df, std::string x_id, std::string y_id, std::string x_lon_col, std::string x_lat_col, std::string y_lon_col, std::string y_lat_col, std::string dist_function);
RcppExport SEXP _distRcpp_dist_max(SEXP x_dfSEXP, SEXP y_dfSEXP, SEXP x_idSEXP, SEXP y_idSEXP, SEXP x_lon_colSEXP, SEXP x_lat_colSEXP, SEXP y_lon_colSEXP, SEXP y_lat_colSEXP, SEXP dist_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type x_df(x_dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type y_df(y_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_id(x_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_id(y_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lon_col(x_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lat_col(x_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lon_col(y_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lat_col(y_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_max(x_df, y_df, x_id, y_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function));
    return rcpp_result_gen;
END_RCPP
}
// dist_sum_inv
Rcpp::DataFrame dist_sum_inv(Rcpp::DataFrame x_df, Rcpp::DataFrame y_df, std::string x_id, std::string y_id, std::string x_lon_col, std::string x_lat_col, std::string y_lon_col, std::string y_lat_col, std::string dist_function, std::string dist_transform, double decay, double scale_units);
RcppExport SEXP _distRcpp_dist_sum_inv(SEXP x_dfSEXP, SEXP y_dfSEXP, SEXP x_idSEXP, SEXP y_idSEXP, SEXP x_lon_colSEXP, SEXP x_lat_colSEXP, SEXP y_lon_colSEXP, SEXP y_lat_colSEXP, SEXP dist_functionSEXP, SEXP dist_transformSEXP, SEXP decaySEXP, SEXP scale_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type x_df(x_dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type y_df(y_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_id(x_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_id(y_idSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lon_col(x_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type x_lat_col(x_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lon_col(y_lon_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type y_lat_col(y_lat_colSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_function(dist_functionSEXP);
    Rcpp::traits::input_parameter< std::string >::type dist_transform(dist_transformSEXP);
    Rcpp::traits::input_parameter< double >::type decay(decaySEXP);
    Rcpp::traits::input_parameter< double >::type scale_units(scale_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_sum_inv(x_df, y_df, x_id, y_id, x_lon_col, x_lat_col, y_lon_col, y_lat_col, dist_function, dist_transform, decay, scale_units));
    return rcpp_result_gen;
END_RCPP
}
// deg_to_rad
double deg_to_rad(const double& degree);
RcppExport SEXP _distRcpp_deg_to_rad(SEXP degreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type degree(degreeSEXP);
    rcpp_result_gen = Rcpp::wrap(deg_to_rad(degree));
    return rcpp_result_gen;
END_RCPP
}
// dist_haversine
double dist_haversine(const double& xlon, const double& xlat, const double& ylon, const double& ylat);
RcppExport SEXP _distRcpp_dist_haversine(SEXP xlonSEXP, SEXP xlatSEXP, SEXP ylonSEXP, SEXP ylatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type xlon(xlonSEXP);
    Rcpp::traits::input_parameter< const double& >::type xlat(xlatSEXP);
    Rcpp::traits::input_parameter< const double& >::type ylon(ylonSEXP);
    Rcpp::traits::input_parameter< const double& >::type ylat(ylatSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_haversine(xlon, xlat, ylon, ylat));
    return rcpp_result_gen;
END_RCPP
}
// dist_vincenty
double dist_vincenty(const double& xlon, const double& xlat, const double& ylon, const double& ylat);
RcppExport SEXP _distRcpp_dist_vincenty(SEXP xlonSEXP, SEXP xlatSEXP, SEXP ylonSEXP, SEXP ylatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type xlon(xlonSEXP);
    Rcpp::traits::input_parameter< const double& >::type xlat(xlatSEXP);
    Rcpp::traits::input_parameter< const double& >::type ylon(ylonSEXP);
    Rcpp::traits::input_parameter< const double& >::type ylat(ylatSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_vincenty(xlon, xlat, ylon, ylat));
    return rcpp_result_gen;
END_RCPP
}
// inverse_value
Rcpp::NumericVector inverse_value(const Rcpp::NumericVector& d, double exp, std::string transform);
RcppExport SEXP _distRcpp_inverse_value(SEXP dSEXP, SEXP expSEXP, SEXP transformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type exp(expSEXP);
    Rcpp::traits::input_parameter< std::string >::type transform(transformSEXP);
    rcpp_result_gen = Rcpp::wrap(inverse_value(d, exp, transform));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_distRcpp_dist_mtom", (DL_FUNC) &_distRcpp_dist_mtom, 5},
    {"_distRcpp_dist_df", (DL_FUNC) &_distRcpp_dist_df, 5},
    {"_distRcpp_dist_1tom", (DL_FUNC) &_distRcpp_dist_1tom, 5},
    {"_distRcpp_dist_1to1", (DL_FUNC) &_distRcpp_dist_1to1, 5},
    {"_distRcpp_popdist_weighted_mean", (DL_FUNC) &_distRcpp_popdist_weighted_mean, 12},
    {"_distRcpp_dist_weighted_mean", (DL_FUNC) &_distRcpp_dist_weighted_mean, 11},
    {"_distRcpp_dist_min", (DL_FUNC) &_distRcpp_dist_min, 9},
    {"_distRcpp_dist_max", (DL_FUNC) &_distRcpp_dist_max, 9},
    {"_distRcpp_dist_sum_inv", (DL_FUNC) &_distRcpp_dist_sum_inv, 12},
    {"_distRcpp_deg_to_rad", (DL_FUNC) &_distRcpp_deg_to_rad, 1},
    {"_distRcpp_dist_haversine", (DL_FUNC) &_distRcpp_dist_haversine, 4},
    {"_distRcpp_dist_vincenty", (DL_FUNC) &_distRcpp_dist_vincenty, 4},
    {"_distRcpp_inverse_value", (DL_FUNC) &_distRcpp_inverse_value, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_distRcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
