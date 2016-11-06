// dist.cpp
#include <Rcpp.h>

using namespace Rcpp;

// mean Earth radius in meters
#define a 6378137.0
#define f 1 / 298.257223563
#define b (1. - f) * a

// convert degrees to radians
double degToRad(double degree) {
  return (degree * M_PI / 180);
}

// compute Haversine distance
double dist_haversine(double xlon, double xlat, double ylon, double ylat) {

  // convert degrees to radians
  xlon = degToRad(xlon);
  xlat = degToRad(xlat);
  ylon = degToRad(ylon);
  ylat = degToRad(ylat);

  // compute parenthetical: sin(delta / 2)
  double d1 = sin((ylat - xlat) / 2.);
  double d2 = sin((ylon - xlon) / 2.);

  return 2.0 * a * asin(sqrt(d1 * d1 + cos(xlat) * cos(ylat) * d2 * d2));
}

// compute Vincenty distance
double dist_vincenty(double xlon, double xlat, double ylon, double ylat) {

  // convert degrees to radians
  xlon = degToRad(xlon);
  xlat = degToRad(xlat);
  ylon = degToRad(ylon);
  ylat = degToRad(ylat);

  double U1 = atan((1. - f) * tan(xlat));
  double U2 = atan((1. - f) * tan(ylat));
  double sinU1 = sin(U1);
  double sinU2 = sin(U2);
  double cosU1 = cos(U1);
  double cosU2 = cos(U2);
  double L = ylon - xlon;
  double lambda = L;

  int iters = 100;
  double tol = 1.0e-12;
  int again = 1;

  double sinLambda, cosLambda,
    p1, p2, sinsig, cossig, sigma, sina, cos2a, cos2sigm, C, lambdaOld,
    Usq, A, B, dsigma;

  while (again) {

    // https://en.wikipedia.org/wiki/Vincenty%27s_formulae

    sinLambda = sin(lambda);
    cosLambda = cos(lambda);

    p1 = cosU2 * sinLambda;
    p2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;

    sinsig = sqrt(p1 * p1 + p2 * p2);
    cossig = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

    sigma = atan2(sinsig, cossig);

    sina = cosU1 * cosU2 * sinLambda / sinsig;

    cos2a = 1. - (sina * sina);

    cos2sigm = cossig - 2. * sinU1 * sinU2 / cos2a;

    C = f / 16. * cos2a * (4. + f * (4. - 3. * cos2a));

    lambdaOld = lambda;

    lambda = L + (1 - C) * f * sina *
      (sigma + C * sinsig *
       (cos2sigm + C * cossig *
	(-1. + 2. * cos2sigm * cos2sigm)));

    iters -= 1;

    again = (fabs(lambda - lambdaOld) > tol && iters > 0);

  }

  if (iters == 0) {

    throw exception("Failed to converge!");

  }
  else {

    Usq = cos2a * (a * a - b * b) / (b * b);

    A = 1. + Usq / 16384. * (4096. + Usq * (-768. + Usq * (320. - 175. * Usq)));
    B = Usq / 1024. * (256. + Usq * (-128. + Usq * (74. - 47. * Usq)));

    dsigma = B * sinsig *
      (cos2sigm + B / 4. *
       (cossig * (-1. + 2. * cos2sigm * cos2sigm)
	- B / 6. * cos2sigm * (-3. + 4. * sinsig * sinsig)
	* (-3. + 4. * cos2sigm * cos2sigm)));

    return b * A * (sigma - dsigma);

  }
}

// function to choose distance measurement method
typedef double (*funcPtr)(double xlon, double xlat,
			  double ylon, double ylat);

XPtr<funcPtr> choose_func(std::string funcnamestr) {

  if (funcnamestr == "Haversine")
    return(XPtr<funcPtr>(new funcPtr(&dist_haversine)));
  else if (funcnamestr == "Vincenty")
    return(XPtr<funcPtr>(new funcPtr(&dist_vincenty)));
  else
    return XPtr<funcPtr>(R_NilValue);

}

//' Compute distance between each coordinate pair (many to many)
//' and return matrix.
//'
//' @param xlon Vector of longitudes for starting coordinate pairs
//' @param xlat Vector of latitudes for starting coordinate pairs
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Matrix of distances between each coordinate pair in meters
//' @export
// [[Rcpp::export]]
NumericMatrix dist_mtom(const NumericVector& xlon,
			const NumericVector& xlat,
			const NumericVector& ylon,
			const NumericVector& ylat,
			std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = choose_func(funname);
  funcPtr fun = *xpfun;

  int n = xlon.size();
  int k = ylon.size();

  NumericMatrix dist(n,k);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < k; j++) {

      // compute distance and store
      dist(i,j) = fun(xlon[i], xlat[i], ylon[j], ylat[j]);

    }
  }

  return dist;

}
//' Compute distance between corresponding coordinate pairs in data frame.
//'
//' Compute distance between corresponding coordinate pairs and return vector.
//' For use when creating a new data frame or tbl_df column.
//'
//' @param xlon Vector of longitudes for starting coordinate pairs
//' @param xlat Vector of latitudes for starting coordinate pairs
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Vector of distances between each coordinate pair in meters
//' @export
// [[Rcpp::export]]
NumericVector dist_df(const NumericVector& xlon,
		      const NumericVector& xlat,
		      const NumericVector& ylon,
		      const NumericVector& ylat,
		      std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = choose_func(funname);
  funcPtr fun = *xpfun;

  int k = ylon.size();
  NumericVector dist(k);

  for(int i = 0; i < k; i++) {

    // compute distance and store
    dist[i] = fun(xlon[i], xlat[i], ylon[i], ylat[i]);

  }

  return dist;

}

//' Compute one to many distances.
//'
//' Compute distances between single starting coordinate and vector of
//' ending coordinates (one to many) and return vector.
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Vector of distances in meters
//' @export
// [[Rcpp::export]]
NumericVector dist_1tom(const double& xlon,
			const double& xlat,
			const NumericVector& ylon,
			const NumericVector& ylat,
			std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = choose_func(funname);
  funcPtr fun = *xpfun;

  int k = ylon.size();
  NumericVector dist(k);

  for(int i = 0; i < k; i++) {

    // compute distance and store
    dist[i] = fun(xlon, xlat, ylon[i], ylat[i]);

  }

  return dist;

}

//' Compute one to one distance.
//'
//' Compute distance between two points (one to one) and return single value.
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Longitude for ending coordinate pair
//' @param ylat Latitude for ending coordinate pair
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Distance in meters
//' @export
// [[Rcpp::export]]
double dist_1to1(const double& xlon,
		 const double& xlat,
		 const double& ylon,
		 const double& ylat,
		 std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = choose_func(funname);
  funcPtr fun = *xpfun;

  // compute and return
  return fun(xlon, xlat, ylon, ylat);

}

NumericVector inverse_dist_weight(const NumericVector& d,
				  double exp,
				  std::string transform) {

  if (transform == "log")

    return 1 / pow(log(d), exp);

  else

    return 1 / pow(d, exp);

}

//' Interpolate population/inverse-distance-weighted measures.
//'
//' Interpolate population/inverse-distance-weighted measures for each \strong{x}
//' coordinate using measures taken at surrounding \strong{y} coordinates.
//' Ending measures are double weighted by population and distance so that
//' surrounding measures taken in nearby areas and those with greater
//' populations are given more weight in final average.
//'
//' @param x_df DataFrame with coordinates that need weighted measures
//' @param y_df DataFrame with coordinates at which measures were taken
//' @param measure_col String name of measure column in y_df
//' @param x_id String name of unique identifer column in x_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param pop_col String name of column in x_df with population values
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @param dist_transform String value of distance weight transform: "level" (default)
//' or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of population/distance-weighted values
//' @export
// [[Rcpp::export]]
DataFrame popdist_weighted_mean(DataFrame x_df,
				DataFrame y_df,
				std::string measure_col,
				std::string x_id = "id",
				std::string x_lon_col = "lon",
				std::string x_lat_col = "lat",
				std::string y_lon_col = "lon",
				std::string y_lat_col = "lat",
				std::string pop_col = "pop",
				std::string dist_function = "Haversine",
				std::string dist_transform = "level",
				double decay = 2) {

  // init
  CharacterVector id = x_df[x_id];
  NumericVector xlon = x_df[x_lon_col];
  NumericVector xlat = x_df[x_lat_col];
  NumericVector meas = y_df[measure_col];
  NumericVector ylon = y_df[y_lon_col];
  NumericVector ylat = y_df[y_lat_col];
  NumericVector popw = y_df[pop_col];

  int n = xlon.size();
  int k = ylon.size();
  NumericVector out(n);

  // loop to compute
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 100 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    NumericVector dist = dist_1tom(xlon[i], xlat[i], ylon, ylat, dist_function);

    // inverse distance weights
    NumericVector idw = inverse_dist_weight(dist, decay, dist_transform);

    // population adjusted idw
    NumericVector w = idw * popw;

    // weight denominator
    double w_sum = 0;
    for (int j = 0; j < k; j++) {
      w_sum += w[j];
    }

    // weight the measure_i
    NumericVector tmp = (w * meas) / (w_sum);

    // sum the wmeasure_i
    double sum = 0;
    for (int j = 0; j < k; j++) {
      sum += tmp[j];
    }

    out[i] = sum;

  }

  return DataFrame::create(_["id"] = id,
			   _["wmeasure"] = out);

}

//' Interpolate inverse-distance-weighted measures.
//'
//' Interpolate inverse-distance-weighted measures for each \strong{x}
//' coordinate using measures taken at surrounding \strong{y} coordinates.
//' Ending measures are weighted by inverse distance so that
//' surrounding measures taken in nearby areas are given more weight in final
//' average.
//'
//' @param x_df DataFrame with coordinates that need weighted measures
//' @param y_df DataFrame with coordinates at which measures were taken
//' @param measure_col String name of measure column in y_df
//' @param x_id String name of unique identifer column in x_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @param dist_transform String value of distance weight transform: "level" (default)
//' or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of distance-weighted values
//' @export
// [[Rcpp::export]]
DataFrame dist_weighted_mean(DataFrame x_df,
			     DataFrame y_df,
			     std::string measure_col,
			     std::string x_id = "id",
			     std::string x_lon_col = "lon",
			     std::string x_lat_col = "lat",
			     std::string y_lon_col = "lon",
			     std::string y_lat_col = "lat",
			     std::string dist_function = "Haversine",
			     std::string dist_transform = "level",
			     double decay = 2) {

  // init
  CharacterVector id = x_df[x_id];
  NumericVector xlon = x_df[x_lon_col];
  NumericVector xlat = x_df[x_lat_col];
  NumericVector meas = y_df[measure_col];
  NumericVector ylon = y_df[y_lon_col];
  NumericVector ylat = y_df[y_lat_col];

  int n = xlon.size();
  int k = ylon.size();
  NumericVector out(n);

  // loop to compute
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 100 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    NumericVector dist = dist_1tom(xlon[i], xlat[i], ylon, ylat, dist_function);

    // inverse distance weights
    NumericVector w = inverse_dist_weight(dist, decay, dist_transform);

    // weight denominator
    double w_sum = 0;
    for (int j = 0; j < k; j++) {
      w_sum += w[j];
    }

    // weight the measure_i
    NumericVector tmp = (w * meas) / (w_sum);

    // sum the wmeasure_i
    double sum = 0;
    for (int j = 0; j < k; j++) {
      sum += tmp[j];
    }

    out[i] = sum;

  }

  return DataFrame::create(_["id"] = id,
			   _["wmeasure"] = out);

}

//' Find minimum distance.
//'
//' Find minimum distance between each starting point in \strong{x} and
//' possible end points, \strong{y}.
//'
//' @param x_df DataFrame with coordinates that need weighted measures
//' @param y_df DataFrame with coordinates at which measures were taken
//' @param x_id String name of unique identifer column in x_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @return DataFrame with minimum distance in meters
//' @export
// [[Rcpp::export]]
DataFrame dist_min(DataFrame x_df,
		   DataFrame y_df,
		   std::string x_id = "id",
		   std::string x_lon_col = "lon",
		   std::string x_lat_col = "lat",
		   std::string y_lon_col = "lon",
		   std::string y_lat_col = "lat",
		   std::string dist_function = "Haversine") {

   // init
  CharacterVector id = x_df[x_id];
  NumericVector xlon = x_df[x_lon_col];
  NumericVector xlat = x_df[x_lat_col];
  NumericVector ylon = y_df[y_lon_col];
  NumericVector ylat = y_df[y_lat_col];

  int n = xlon.size();
  NumericVector out(n);

  // loop
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 100 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    NumericVector dist = dist_1tom(xlon[i], xlat[i], ylon, ylat, dist_function);

    // add min to out
    out[i] = min(dist);

  }

  return DataFrame::create(_["id"] = id,
			   _["mindist"] = out);

}

