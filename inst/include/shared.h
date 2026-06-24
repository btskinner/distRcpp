#ifndef DISTRCPP_SHARED_H
#define DISTRCPP_SHARED_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// length in meters of semi major axis (equatorial radius) for WGS84
#define a 6378137.0
// flattening of ellipsoid for WGS84 
#define f 1 / 298.257223563
// length in meters of semi minor axis (polar radius) for WGS84
#define b (1.0 - f) * a
// mean Earth radius for WGS84
#define r 1.0 / 3.0 * (2 * a + b)

double deg_to_rad(const double& degree);

double dist_haversine(const double& xlon,
                      const double& xlat,
                      const double& ylon,
                      const double& ylat,
                      std::string radius);

double dist_vincenty(const double& xlon,
                     const double& xlat,
                     const double& ylon,
                     const double& ylat);

Rcpp::NumericVector inverse_value(const Rcpp::NumericVector& d,
                                  double exp,
                                  std::string transform);

typedef double (*funcPtr)(const double& xlon,
                          const double& xlat,
                          const double& ylon,
                          const double& ylat);

Rcpp::XPtr<funcPtr> choose_func(std::string funcnamestr);

#endif


