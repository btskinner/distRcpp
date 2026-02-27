#ifndef DISTRCPP_SHARED_H
#define DISTRCPP_SHARED_H

#include <RcppArmadillo.h>

#define a 6378137.0
#define f 1 / 298.257223563
#define b (1. - f) * a

double deg_to_rad(const double& degree);

double dist_haversine(const double& xlon,
                      const double& xlat,
                      const double& ylon,
                      const double& ylat);

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


