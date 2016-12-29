#ifndef DISTRCPP_SHARED_H
#define DISTRCPP_SHARED_H
#include <Rcpp.h>

using namespace Rcpp;

#define a 6378137.0
#define f 1 / 298.257223563
#define b (1. - f) * a

double deg_to_rad(double degree);

double dist_haversine(double xlon, double xlat, double ylon, double ylat);

double dist_vincenty(double xlon, double xlat, double ylon, double ylat);

NumericVector inverse_value(const NumericVector& d, double exp, std::string transform);

typedef double (*funcPtr)(double xlon, double xlat, double ylon, double ylat);

XPtr<funcPtr> choose_func(std::string funcnamestr);

#endif


