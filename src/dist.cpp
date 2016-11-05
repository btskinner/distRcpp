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
double distHaversine(double xlon, double xlat, double ylon, double ylat) {

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
double distVincenty(double xlon, double xlat, double ylon, double ylat) {

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

XPtr<funcPtr> chooseFunc(std::string funcnamestr) {

  if (funcnamestr == "Haversine")
    return(XPtr<funcPtr>(new funcPtr(&distHaversine)));
  else if (funcnamestr == "Vincenty")
    return(XPtr<funcPtr>(new funcPtr(&distVincenty)));
  else
    return XPtr<funcPtr>(R_NilValue);

}

//' Compute distance between each coordinate pair and return matrix
//'
//' @param xlon Vector of longitudes for starting coordinate pairs
//' @param xlat Vector of latitudes for starting coordinate pairs
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Matrix of distances between each coordinate pair in meters
// [[Rcpp::export]]
NumericMatrix getDistMat(const NumericVector& xlon,
			 const NumericVector& xlat,
			 const NumericVector& ylon,
			 const NumericVector& ylat,
			 std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = chooseFunc(funname);
  funcPtr fun = *xpfun;

  int n = xlon.size();
  int k = ylon.size();

  NumericMatrix dist(n,k);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < k; j++) {
      dist(i,j) = fun(xlon[i], xlat[i], ylon[j], ylat[j]);
    }
  }

  return dist;

}
//' Compute distance between corresponding coordinate pairs and return vector.
//' For use in a data frame.
//'
//' @param xlon Vector of longitudes for starting coordinate pairs
//' @param xlat Vector of latitudes for starting coordinate pairs
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Vector of distances between each coordinate pair in meters
// [[Rcpp::export]]
NumericVector getDistDF(const NumericVector& xlon,
			const NumericVector& xlat,
			const NumericVector& ylon,
			const NumericVector& ylat,
			std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = chooseFunc(funname);
  funcPtr fun = *xpfun;

  int k = ylon.size();
  NumericVector dist(k);

  for(int i = 0; i < k; i++) {
    dist[i] = fun(xlon[i], xlat[i], ylon[i], ylat[i]);
  }

  return dist;

}

//' Compute distances between single starting coordinate and vector of
//' ending coordinates.
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Vector of distances in meters
// [[Rcpp::export]]
NumericVector getDistVec(const double& xlon,
			 const double& xlat,
			 const NumericVector& ylon,
			 const NumericVector& ylat,
			 std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = chooseFunc(funname);
  funcPtr fun = *xpfun;

  int k = ylon.size();
  NumericVector dist(k);

  for(int i = 0; i < k; i++) {
    dist[i] = fun(xlon, xlat, ylon[i], ylat[i]);
  }

  return dist;

}

//' Compute distance between two points.
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Longitude for ending coordinate pair
//' @param ylat Latitude for ending coordinate pair
//' @param funname String name of distance function: Haversine, Vincenty
//' @return Distance in meters
// [[Rcpp::export]]
double getDist(const double& xlon,
	       const double& xlat,
	       const double& ylon,
	       const double& ylat,
	       std::string funname) {

  // select function
  XPtr<funcPtr> xpfun = chooseFunc(funname);
  funcPtr fun = *xpfun;

  // compute and return
  return fun(xlon, xlat, ylon, ylat);

}

NumericVector invDistWeight(const NumericVector& d,
			    double exp,
			    std::string transform) {

  if (transform == "Log")

    return 1 / pow(log(d), exp);

  else

    return 1 / pow(d, exp);

}

//' Compute weighted measures for each coordinate using measures taken at
//' surrounding coordinates. Ending measures are double weighted by population
//' and distance so that surrounding measures taken in nearby areas and those
//' with greater populations are given more weight.
//'
//' @param fromDF DataFrame with coordinates that need weighted measures
//' @param toDF DataFrame with coordinates at which measures were taken
//' @param measureName String name of measure column in toDF
//' @param fromID String name of unique identifer column in fromDF
//' @param fromYear String name of year column in fromDF (assumes multiple years)
//' @param fromLonName String name of column in fromDF with longitude values
//' @param fromLatName String name of column in fromDF with latitude values
//' @param toLonName String name of column in toDF with longitude values
//' @param toLatName String name of column in toDF with latitude values
//' @param popName String name of column in fromDF with population values
//' @param distFuncName String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @param distTransform String value of distance weight transform: " " (default)
//' or "Log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of population/distance-weighted values
// [[Rcpp::export]]
DataFrame popDistWMean(DataFrame fromDF,
		       DataFrame toDF,
		       std::string measureName,
		       std::string fromID = "unitid",
		       std::string fromYear = "year",
		       std::string fromLonName = "lon",
		       std::string fromLatName = "lat",
		       std::string toLonName = "lon",
		       std::string toLatName = "lat",
		       std::string popName = "pop",
		       std::string distFuncName = "Haversine",
		       std::string distTransform = " ",
		       double decay = 2) {

  // init
  CharacterVector id = fromDF[fromID];
  NumericVector year = fromDF[fromYear];
  NumericVector xlon = fromDF[fromLonName];
  NumericVector xlat = fromDF[fromLatName];
  NumericVector meas = toDF[measureName];
  NumericVector ylon = toDF[toLonName];
  NumericVector ylat = toDF[toLatName];
  NumericVector popw = toDF[popName];

  int n = xlon.size();
  int k = ylon.size();
  NumericVector out(n);

  // loop to compute
  for (int i = 0; i < n; i++) {

    // distance vector
    NumericVector dist = getDistVec(xlon[i], xlat[i], ylon, ylat, distFuncName);

    // inverse distance weights
    NumericVector idw = invDistWeight(dist, decay, distTransform);

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
			   _["year"] = year,
			   _["wmeasure"] = out);

}

//' Compute distance-weighted measures for each coordinate using measures taken at
//' surrounding coordinates. Ending measures are inverse distance-weighted so that
//' surrounding measures taken in nearby areas are given more weight.
//'
//' @param fromDF DataFrame with coordinates that need weighted measures
//' @param toDF DataFrame with coordinates at which measures were taken
//' @param measureName String name of measure column in toDF
//' @param fromID String name of unique identifer column in fromDF
//' @param fromYear String name of year column in fromDF (assumes multiple years)
//' @param fromLonName String name of column in fromDF with longitude values
//' @param fromLatName String name of column in fromDF with latitude values
//' @param toLonName String name of column in toDF with longitude values
//' @param toLatName String name of column in toDF with latitude values
//' @param distFuncName String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @param distTransform String value of distance weight transform: " " (default)
//' or "Log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of distance-weighted values
// [[Rcpp::export]]
DataFrame distWMean(DataFrame fromDF,
		    DataFrame toDF,
		    std::string measureName,
		    std::string fromID = "unitid",
		    std::string fromYear = "year",
		    std::string fromLonName = "lon",
		    std::string fromLatName = "lat",
		    std::string toLonName = "lon",
		    std::string toLatName = "lat",
		    std::string distFuncName = "Haversine",
		    std::string distTransform = " ",
		    double decay = 2) {

  // init
  CharacterVector id = fromDF[fromID];
  NumericVector year = fromDF[fromYear];
  NumericVector xlon = fromDF[fromLonName];
  NumericVector xlat = fromDF[fromLatName];
  NumericVector meas = toDF[measureName];
  NumericVector ylon = toDF[toLonName];
  NumericVector ylat = toDF[toLatName];

  int n = xlon.size();
  int k = ylon.size();
  NumericVector out(n);

  // loop to compute
  for (int i = 0; i < n; i++) {

    // distance vector
    NumericVector dist = getDistVec(xlon[i], xlat[i], ylon, ylat, distFuncName);

    // inverse distance weights
    NumericVector w = invDistWeight(dist, decay, distTransform);

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
			   _["year"] = year,
			   _["wmeasure"] = out);

}

//' Find minimum distance between starting point and possible end points.
//'
//' @param fromDF DataFrame with coordinates that need weighted measures
//' @param toDF DataFrame with coordinates at which measures were taken
//' @param fromID String name of unique identifer column in fromDF
//' @param fromYear String name of year column in fromDF (assumes multiple years)
//' @param fromLonName String name of column in fromDF with longitude values
//' @param fromLatName String name of column in fromDF with latitude values
//' @param toLonName String name of column in toDF with longitude values
//' @param toLatName String name of column in toDF with latitude values
//' @param distFuncName String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @return DataFrame with minimum distance in meters
// [[Rcpp::export]]
DataFrame minDist(DataFrame fromDF,
		  DataFrame toDF,
		  std::string fromID = "unitid",
		  std::string fromYear = "year",
		  std::string fromLonName = "lon",
		  std::string fromLatName = "lat",
		  std::string toLonName = "lon",
		  std::string toLatName = "lat",
		  std::string distFuncName = "Haversine") {

   // init
  CharacterVector id = fromDF[fromID];
  NumericVector year = fromDF[fromYear];
  NumericVector xlon = fromDF[fromLonName];
  NumericVector xlat = fromDF[fromLatName];
  NumericVector ylon = toDF[toLonName];
  NumericVector ylat = toDF[toLatName];

  int n = xlon.size();
  NumericVector out(n);

  // loop
  for (int i = 0; i < n; i++) {

    // distance vector
    NumericVector dist = getDistVec(xlon[i], xlat[i], ylon, ylat, distFuncName);

    // add min to out
    out[i] = min(dist);

  }

  return DataFrame::create(_["id"] = id,
			   _["year"] = year,
			   _["mindist"] = out);

}

