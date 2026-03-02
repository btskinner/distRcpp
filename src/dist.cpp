// dist.cpp
#include <shared.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute one to one distance.
//'
//' Compute distance between two points (one to one) and return single value.
//'
//' @param xlon Longitude for starting coordinate pair
//' @param xlat Latitude for starting coordinate pair
//' @param ylon Longitude for ending coordinate pair
//' @param ylat Latitude for ending coordinate pair
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Distance in meters
//' @export
// [[Rcpp::export]]
double dist_1to1(const double& xlon,
                 const double& xlat,
                 const double& ylon,
                 const double& ylat,
                 const std::string dist_function="Haversine") {
  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;
  // compute and return
  return fun(xlon, xlat, ylon, ylat);
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
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Vector of distances in meters
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_1tom(const double& xlon,
                              const double& xlat,
                              const Rcpp::NumericVector& ylon,
                              const Rcpp::NumericVector& ylat,
                              const std::string dist_function="Haversine") {
  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;
  // size of end points
  const int k = ylon.size();
  // init vector that is ylon.size() since there's just one start point
  Rcpp::NumericVector dist(k);
  // loop through each end point (same starting point each time)
  for(int i = 0; i < k; i++) {
    // compute distance and store
    dist[i] = fun(xlon, xlat, ylon[i], ylat[i]);
  }
  // return vector
  return dist;
}

//' Compute distance between each coordinate pair (many to many)
//' and return matrix.
//'
//' @param xlon Vector of longitudes for starting coordinate pairs
//' @param xlat Vector of latitudes for starting coordinate pairs
//' @param ylon Vector of longitudes for ending coordinate pairs
//' @param ylat Vector of latitudes for ending coordinate pairs
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Matrix of distances between each coordinate pair in meters
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_mtom(const Rcpp::NumericVector& xlon,
                              const Rcpp::NumericVector& xlat,
                              const Rcpp::NumericVector& ylon,
                              const Rcpp::NumericVector& ylat,
                              const std::string dist_function="Haversine") {
  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;
  // get vector sizes for setting matrix size
  const int n = xlon.size();
  const int k = ylon.size();
  // init output matrix that is x.size() by y.size()
  Rcpp::NumericMatrix dist(n,k);
  // double loop: for each start, for each end...
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < k; j++) {
      // compute distance and store in matrix cell
      dist(i,j) = fun(xlon[i], xlat[i], ylon[j], ylat[j]);
    }
  }
  // return matrix
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
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Vector of distances between each coordinate pair in meters
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_df(const Rcpp::NumericVector& xlon,
                            const Rcpp::NumericVector& xlat,
                            const Rcpp::NumericVector& ylon,
                            const Rcpp::NumericVector& ylat,
                            const std::string dist_function="Haversine") {
  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;
  // get the size of end vector (assumed since in a dataframe)
  const int k = ylon.size();
  // init distance vector of ylon.size()
  Rcpp::NumericVector dist(k);
  // loop through each row (just one by one measures)
  for(int i = 0; i < k; i++) {
    // compute distance and store
    dist[i] = fun(xlon[i], xlat[i], ylon[i], ylat[i]);
  }
  // return vector
  return dist;
}

//' Find maximum distance.
//'
//' Find maximum distance between each starting point in \strong{x} and
//' possible end points, \strong{y}.
//'
//' @param x_df DataFrame with starting coordinates
//' @param y_df DataFrame with ending coordinates
//' @param x_id String name of unique identifer column in x_df
//' @param y_id String name of unique identifer column in y_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param dist_function String name of distance function: "Haversine" (default)
//'     or "Vincenty"
//' @return DataFrame with id of farthest point and distance in meters
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_max(Rcpp::DataFrame x_df,
                         Rcpp::DataFrame y_df,
                         const std::string& x_id="id",
                         const std::string& y_id="id",
                         const std::string& x_lon_col="lon",
                         const std::string& x_lat_col="lat",
                         const std::string& y_lon_col="lon",
                         const std::string& y_lat_col="lat",
                         const std::string& dist_function="Haversine") {
  // init vectors
  Rcpp::CharacterVector idx = x_df[x_id];     // id: start
  Rcpp::CharacterVector idy = y_df[y_id];     // id: end
  Rcpp::NumericVector xlon = x_df[x_lon_col]; // x: lon
  Rcpp::NumericVector xlat = x_df[x_lat_col]; // x: lat
  Rcpp::NumericVector ylon = y_df[y_lon_col]; // y: lon
  Rcpp::NumericVector ylat = y_df[y_lat_col]; // y: lat
  // get size for loop
  const int n = xlon.size();
  Rcpp::NumericVector dist(n);  // vector of distances
  Rcpp::CharacterVector end(n); // vector ids for those distances
  // loop through start points
  for (int i = 0; i < n; i++) {
    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // compute distance vector
    Rcpp::NumericVector distvec = dist_1tom(xlon[i], xlat[i],
                                            ylon, ylat, dist_function);
    // add maximum to distance output
    dist[i] = max(distvec);
    // get id of maximum
    end[i] = idy[which_max(distvec)];
  }
  // return data frame
  return Rcpp::DataFrame::create(Rcpp::Named("id_start") = idx,
                                 Rcpp::Named("id_end") = end,
                                 Rcpp::Named("meters") = dist,
                                 Rcpp::Named("stringsAsFactors") = false);
}

//' Find minimum distance.
//'
//' Find minimum distance between each starting point in \strong{x} and
//' possible end points, \strong{y}.
//'
//' @param x_df DataFrame with starting coordinates
//' @param y_df DataFrame with ending coordinates
//' @param x_id String name of unique identifer column in x_df
//' @param y_id String name of unique identifer column in y_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param dist_function String name of distance function: "Haversine" (default)
//'     or "Vincenty"
//' @return DataFrame with id of closest point and distance in meters
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_min(Rcpp::DataFrame x_df,
                         Rcpp::DataFrame y_df,
                         const std::string& x_id="id",
                         const std::string& y_id="id",
                         const std::string& x_lon_col="lon",
                         const std::string& x_lat_col="lat",
                         const std::string& y_lon_col="lon",
                         const std::string& y_lat_col="lat",
                         const std::string& dist_function="Haversine") {
  // init vectors
  Rcpp::CharacterVector idx = x_df[x_id];     // id: start
  Rcpp::CharacterVector idy = y_df[y_id];     // id: end
  Rcpp::NumericVector xlon = x_df[x_lon_col]; // x: lon
  Rcpp::NumericVector xlat = x_df[x_lat_col]; // x: lat
  Rcpp::NumericVector ylon = y_df[y_lon_col]; // y: lon
  Rcpp::NumericVector ylat = y_df[y_lat_col]; // y: lat
  // get size for loop
  const int n = xlon.size();
  Rcpp::NumericVector dist(n);  // vector of distances
  Rcpp::CharacterVector end(n); // vector ids for those distances
  // loop through start points
  for (int i = 0; i < n; i++) {
    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // compute distance vector
    Rcpp::NumericVector distvec = dist_1tom(xlon[i], xlat[i],
                                            ylon, ylat, dist_function);
    // add minimum to distance output
    dist[i] = min(distvec);
    // get id of minimum
    end[i] = idy[which_min(distvec)];
  }
  // return data frame
  return Rcpp::DataFrame::create(Rcpp::Named("id_start") = idx,
                                 Rcpp::Named("id_end") = end,
                                 Rcpp::Named("meters") = dist,
                                 Rcpp::Named("stringsAsFactors") = false);
}

//' Find nearest X number of points.
//'
//' Find nearest X values between each starting point in \strong{x} and
//' possible end points, \strong{y}.
//'
//' @param x_df DataFrame with starting coordinates
//' @param y_df DataFrame with ending coordinates
//' @param num_nearest The number of closest points to return
//' @param x_id String name of unique identifer column in x_df
//' @param y_id String name of unique identifer column in y_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param dist_function String name of distance function: "Haversine" (default)
//'     or "Vincenty"
//' @return DataFrame with id of X closest points and distance in meters
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_nearest_x(Rcpp::DataFrame& x_df,
                               Rcpp::DataFrame& y_df,
                               const int num_nearest=10,
                               const std::string x_id="id",
                               const std::string y_id="id",
                               const std::string x_lon_col="lon",
                               const std::string x_lat_col="lat",
                               const std::string y_lon_col="lon",
                               const std::string y_lat_col="lat",
                               const std::string dist_function="Haversine") {
  // init
  Rcpp::CharacterVector idx = x_df[x_id];
  Rcpp::CharacterVector idy = y_df[y_id];
  Rcpp::NumericVector xlon = x_df[x_lon_col];
  Rcpp::NumericVector xlat = x_df[x_lat_col];
  Rcpp::NumericVector ylon = y_df[y_lon_col];
  Rcpp::NumericVector ylat = y_df[y_lat_col];
  // size of starting points
  const auto n = idx.size();
  // adjusted number nearest: num nearest <= length of idy
  const auto adjnn = (num_nearest > idy.size()) ? idy.size() : num_nearest;
  // output vectors need to be size of idx * adjusted number nearest
  Rcpp::CharacterVector out_idx(n * adjnn);
  Rcpp::CharacterVector out_idy(n * adjnn);
  Rcpp::NumericVector out_met(n * adjnn);
  // loop through all starting points
  for (int i = 0; i < n; i++) {
    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // distance vector
    arma::vec distvec = dist_1tom(xlon[i], xlat[i], ylon, ylat, dist_function);
    // sort distance output
    arma::vec sort_distvec = arma::sort(distvec);
    arma::uvec sort_distids = arma::sort_index(distvec);
    // loop through distance vectors to grab adjnn closest
    for (int j = 0; j < adjnn; j++) {
      // i = 0, adjnn = 3: (0) + j --> k = 0, 1, 2
      // i = 1, adjnn = 3: (3) + j --> k = 3, 4, 5
      int k = (i * adjnn) + j;
      out_idx[k] = idx[i];                // repeat starting idx
      out_idy[k] = idy[sort_distids[j]];  // get idy of closest value
      out_met[k] = sort_distvec[j];       // get measure of closest value
    }
  }
  // return tidy data frame
  return Rcpp::DataFrame::create(Rcpp::Named("id_x") = out_idx,
                                 Rcpp::Named("id_y") = out_idy,
                                 Rcpp::Named("meters") = out_met,
                                 Rcpp::Named("stringsAsFactors") = false);

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
//' @param dist_function String name of distance function: "Haversine" (default)
//'     or "Vincenty"
//' @param dist_transform String value of distance weight transform: "level"
//'     (default) or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of distance-weighted values
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_weighted_mean(Rcpp::DataFrame x_df,
                                   Rcpp::DataFrame y_df,
                                   const std::string& measure_col,
                                   const std::string& x_id="id",
                                   const std::string& x_lon_col="lon",
                                   const std::string& x_lat_col="lat",
                                   const std::string& y_lon_col="lon",
                                   const std::string& y_lat_col="lat",
                                   const std::string& dist_function="Haversine",
                                   const std::string& dist_transform="level",
                                   const double decay=2) {
  // init vectors
  Rcpp::CharacterVector id = x_df[x_id];        // id of start
  Rcpp::NumericVector xlon = x_df[x_lon_col];   // x: lon vector
  Rcpp::NumericVector xlat = x_df[x_lat_col];   // x: lat vector
  Rcpp::NumericVector meas = y_df[measure_col]; // what is being averaged
  Rcpp::NumericVector ylon = y_df[y_lon_col];   // y: lon vector
  Rcpp::NumericVector ylat = y_df[y_lat_col];   // y: lat vector
  // get sizes for output vectors and loop
  const int n = xlon.size();
  const int k = ylon.size();
  // init output vector for weighted measure
  Rcpp::NumericVector out(n);
  // loop to compute
  for (int i = 0; i < n; i++) {
    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // compute distance vector
    Rcpp::NumericVector dist = dist_1tom(xlon[i], xlat[i],
                                         ylon, ylat, dist_function);
    // compute inverse distance weights with decay
    Rcpp::NumericVector w = inverse_value(dist, decay, dist_transform);
    // compute weight denominator
    double w_sum = 0.0;
    for (int j = 0; j < k; j++) {
      w_sum += w[j];
    }
    // compute weighted measure_i
    Rcpp::NumericVector tmp = (w * meas) / (w_sum);
    // sum the wmeasure_i
    double sum = 0.0;
    for (int j = 0; j < k; j++) {
      sum += tmp[j];
    }
    // store weighted measure in output vector
    out[i] = sum;
  }
  // return data frame
  return Rcpp::DataFrame::create(Rcpp::Named("id") = id,
                                 Rcpp::Named("wmeasure") = out,
                                 Rcpp::Named("stringsAsFactors") = false);
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
//' @param dist_function String name of distance function: "Haversine" (default)
//'     or "Vincenty"
//' @param dist_transform String value of distance weight transform: "level"
//'     (default) or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of population/distance-weighted values
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame popdist_weighted_mean(Rcpp::DataFrame x_df,
                                      Rcpp::DataFrame y_df,
                                      const std::string& measure_col,
                                      const std::string& x_id="id",
                                      const std::string& x_lon_col="lon",
                                      const std::string& x_lat_col="lat",
                                      const std::string& y_lon_col="lon",
                                      const std::string& y_lat_col="lat",
                                      const std::string& pop_col="pop",
                                      const std::string& dist_function="Haversine",
                                      const std::string& dist_transform="level",
                                      const double decay=2) {

  // init vectors
  Rcpp::CharacterVector id = x_df[x_id];        // id of start
  Rcpp::NumericVector xlon = x_df[x_lon_col];   // x: lon vector
  Rcpp::NumericVector xlat = x_df[x_lat_col];   // x: lat vector
  Rcpp::NumericVector meas = y_df[measure_col]; // what is being averaged
  Rcpp::NumericVector ylon = y_df[y_lon_col];   // y: lon vector
  Rcpp::NumericVector ylat = y_df[y_lat_col];   // y: lat vector
  Rcpp::NumericVector popw = y_df[pop_col];     // population number for weight
  // compute sizes for loops and out vector
  const int n = xlon.size();
  const int k = ylon.size();
  // init out vector for the weighted measure
  Rcpp::NumericVector out(n);
  // loop through each starting point
  for (int i = 0; i < n; i++) {
    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // distance vector from starting point to many end points
    Rcpp::NumericVector dist = dist_1tom(xlon[i], xlat[i],
                                         ylon, ylat, dist_function);
    // compute inverse distance weights with decay
    Rcpp::NumericVector idw = inverse_value(dist, decay, dist_transform);
    // population adjusted inverse distance weights
    Rcpp::NumericVector w = idw * popw;
    // weight denominator
    double w_sum = 0.0;
    for (int j = 0; j < k; j++) {
      w_sum += w[j];
    }
    // apply weight to the measure_i
    Rcpp::NumericVector tmp = (w * meas) / (w_sum);
    // sum the wmeasure_i
    double sum = 0.0;
    for (int j = 0; j < k; j++) {
      sum += tmp[j];
    }
    // store weighted value in vector
    out[i] = sum;
  }
  // return data frame
  return Rcpp::DataFrame::create(Rcpp::Named("id") = id,
                                 Rcpp::Named("wmeasure") = out,
                                 Rcpp::Named("stringsAsFactors") = false);
}

//' Sum inverse distances.
//'
//' Find sum of inverse distances between each starting point in \strong{x}
//' and possible end points, \strong{y}.
//'
//' @param x_df DataFrame with starting coordinates
//' @param y_df DataFrame with ending coordinates
//' @param x_id String name of unique identifer column in x_df
//' @param y_id String name of unique identifer column in y_df
//' @param x_lon_col String name of column in x_df with longitude values
//' @param x_lat_col String name of column in x_df with latitude values
//' @param y_lon_col String name of column in y_df with longitude values
//' @param y_lat_col String name of column in y_df with latitude values
//' @param dist_function String name of distance function: "Haversine" (default)
//'     or "Vincenty"
//' @param dist_transform String value of distance transform: "level" (default)
//'     or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @param scale_units Double value to divide return value by (e.g., 1000 == km)
//' @return DataFrame with sum of distances
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_sum_inv(Rcpp::DataFrame x_df,
                             Rcpp::DataFrame y_df,
                             const std::string& x_id="id",
                             const std::string& y_id="id",
                             const std::string& x_lon_col="lon",
                             const std::string& x_lat_col="lat",
                             const std::string& y_lon_col="lon",
                             const std::string& y_lat_col="lat",
                             const std::string& dist_function="Haversine",
                             const std::string& dist_transform="level",
                             const double decay=2,
                             const double scale_units=1) {

  // init vectors
  Rcpp::CharacterVector idx = x_df[x_id];     // id of start
  Rcpp::CharacterVector idy = y_df[y_id];     // id of end
  Rcpp::NumericVector xlon = x_df[x_lon_col]; // x: lon
  Rcpp::NumericVector xlat = x_df[x_lat_col]; // x: lat
  Rcpp::NumericVector ylon = y_df[y_lon_col]; // y: lon
  Rcpp::NumericVector ylat = y_df[y_lat_col]; // y: lat
  // get sizes for output vector and loop
  int n = xlon.size();
  // init output vector for inverse distances
  Rcpp::NumericVector idist(n);
  // loop through each starting point
  for (int i = 0; i < n; i++) {
    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // compute distance vector
    Rcpp::NumericVector distvec = dist_1tom(xlon[i], xlat[i],
                                            ylon, ylat, dist_function);
    // scale units: e.g., 1000 to convert meters to kilometers
    distvec = distvec / scale_units;
    // compute inverse distance weights
    Rcpp::NumericVector inv_distvec = inverse_value(distvec, decay,
                                                    dist_transform);
    // replace infinite values with 0 since this means dist is very small
    inv_distvec[is_infinite(inv_distvec)] = 0;
    // add summed vector to inverse distance output
    idist[i] = sum(inv_distvec);
  }
  // return data frame
  return Rcpp::DataFrame::create(Rcpp::Named("id") = idx,
                                 Rcpp::Named("inv_distance") = idist,
                                 Rcpp::Named("stringsAsFactors") = false);
}
