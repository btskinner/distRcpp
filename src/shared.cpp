// shared.cpp
#include <shared.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Convert degrees to radians
//'
//' @param degree Degree value
//' @return Radian value (double)
//' @export
// [[Rcpp::export]]
double deg_to_rad(const double& degree) {
  return (degree * M_PI / 180.0);
}

//' Compute Haversine distance between two points
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Longitude for ending coordinate pair
//' @param ylat Latitude for ending coordinate pair
//' @return Double of distance between coordinate pairs in meters
//' @export
// [[Rcpp::export]]
double dist_haversine(const double& xlon,
                      const double& xlat,
                      const double& ylon,
                      const double& ylat) { 
    // return 0 if same point
    if (xlon == ylon && xlat == ylat) return 0;
    // formula:
    // 2r * sin^-1(sqrt(sin^2(d_lat/2) + cos(xlat)cos(ylat)sin^2(d_lon/2)))
    // where:
    //  r     := radius of sphere (using mean Earth radius)
    //  d_lat := lat_x - lat_y
    //  d_lon := lon_x - lon_y
    // identities:
    //  sin^-1   ==> arcsin
    //  sin^2(x) ==> sin(x) * sin(x)
    // rewrite:
    // 2r * arcsin(sqrt(d1^2 + cos(xlat)cos(ylat))d2^2)
    // where:
    //  d1 := sin(delta_lat/2)
    //  d2 := sin(delta_lon/2)
    // convert degrees to radians
    double xlonr = deg_to_rad(xlon);
    double xlatr = deg_to_rad(xlat);
    double ylonr = deg_to_rad(ylon);
    double ylatr = deg_to_rad(ylat);
    // compute parenthetical: sin(delta / 2)
    double d1 = sin((ylatr - xlatr) / 2.0);
    double d2 = sin((ylonr - xlonr) / 2.0);
    // compute
    double dist;
    dist = 2.0 * r * asin(sqrt(d1 * d1 + cos(xlatr) * cos(ylatr) * d2 * d2)); 
    // return
    return dist;
}

//' Compute Vincenty distance between two points
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Longitude for ending coordinate pair
//' @param ylat Latitude for ending coordinate pair
//' @return Double of distance between coordinate pairs in meters
//' @export
// [[Rcpp::export]]
double dist_vincenty(const double& xlon,
                     const double& xlat,
                     const double& ylon,
                     const double& ylat) {

    // return 0 if same point
    if (xlon == ylon && xlat == ylat) return 0;
    // convert degrees to radians
    double xlonr = deg_to_rad(xlon);
    double xlatr = deg_to_rad(xlat);
    double ylonr = deg_to_rad(ylon);
    double ylatr = deg_to_rad(ylat);
    // compute parts
    double U1 = atan((1.0 - f) * tan(xlatr));
    double U2 = atan((1.0 - f) * tan(ylatr));
    double sinU1 = sin(U1);
    double sinU2 = sin(U2);
    double cosU1 = cos(U1);
    double cosU2 = cos(U2);
    double L = ylonr - xlonr;
    double lambda = L;
    // set up iterations, tolerance, and setting for while loop
    int iters = 100;
    double tol = 1.0e-12;
    int again = 1;
    // init double variables for loop
    double sinLambda, cosLambda,
        p1, p2, sinsig, cossig, sigma, sina, cos2a, cos2sigm, C, lambdaOld,
        Usq, A, B, dsigma;
    // loop until tolerance or iteration limit is reached
    while (again) {
        // compute lambda using Vincenty formula
        // https://en.wikipedia.org/wiki/Vincenty%27s_formulae
        sinLambda = sin(lambda);
        cosLambda = cos(lambda);
        p1 = cosU2 * sinLambda;
        p2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
        sinsig = sqrt(p1 * p1 + p2 * p2);
        cossig = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma = atan2(sinsig, cossig);
        sina = cosU1 * cosU2 * sinLambda /  sinsig;
        cos2a = 1.0 - (sina * sina);
        cos2sigm = cossig - 2.0 * sinU1 * sinU2 / cos2a;
        C = f / 16.0 * cos2a * (4.0 + f * (4.0 - 3.0 * cos2a));
        // store prior lambda for comparison
        lambdaOld = lambda;
        // compute updated lambda
        lambda = L + (1.0 - C) * f * sina *
            (sigma + C * sinsig *
                (cos2sigm + C * cossig *
                    (-1.0 + 2.0 * cos2sigm * cos2sigm)
                )
            );
        // decrement iterations
        iters -= 1;
        // if difference is > tolerance AND iterations remain, do again; 
        // otherwise, move on
        again = (fabs(lambda - lambdaOld) > tol && iters > 0);
  }
    // if iterations == 0, let user know that convergence did not
    // happen; otherwise, do final calculation
    if (iters == 0) {
        throw Rcpp::exception("Failed to converge!"); 
    } else {    
        Usq = cos2a * (a * a - b * b) / (b * b);
        A = 1.0 + Usq / 16384.0 * 
            (4096.0 + Usq * (-768.0 + Usq * (320.0 - 175.0 * Usq)));
        B = Usq / 1024.0 * (256.0 + Usq * (-128.0 + Usq * (74.0 - 47.0 * Usq)));
        dsigma = B * sinsig *
            (cos2sigm + B / 4.0 *
                (cossig * (-1.0 + 2.0 * cos2sigm * cos2sigm)
                 - B / 6.0 * cos2sigm * (-3.0 + 4.0 * sinsig * sinsig)
                 * (-3.0 + 4.0 * cos2sigm * cos2sigm)
                )
            );
    double dist = b * A * (sigma - dsigma);
    return dist;
  }
}

//' Compute inverse values from vector
//'
//' @param d Vector of values (e.g., distances)
//' @param exp Rate of decay
//' @param transform == "log" if natural log transform
//' @return Vector of inverse weights
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector inverse_value(const Rcpp::NumericVector& d,
                                  double exp,
                                  std::string transform) {
    if (transform == "log")
        return 1.0 / pow(log(d), exp);
    else
        return 1.0 / pow(d, exp);
}

// function to choose distance measurement method
typedef double (*funcPtr)(const double& xlon,
                          const double& xlat,
                          const double& ylon,
                          const double& ylat);

Rcpp::XPtr<funcPtr> choose_func(std::string funcnamestr) {
    if (funcnamestr == "Haversine")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&dist_haversine)));
    else if (funcnamestr == "Vincenty")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&dist_vincenty)));
    else
        return Rcpp::XPtr<funcPtr>(R_NilValue);
}
