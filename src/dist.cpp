// dist.cpp
#include <shared.h>
#include <Rcpp.h>

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
			      std::string dist_function="Haversine") {

  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;

  int n = xlon.size();
  int k = ylon.size();

  Rcpp::NumericMatrix dist(n,k);

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
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Vector of distances between each coordinate pair in meters
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_df(const Rcpp::NumericVector& xlon,
			    const Rcpp::NumericVector& xlat,
			    const Rcpp::NumericVector& ylon,
			    const Rcpp::NumericVector& ylat,
			    std::string dist_function="Haversine") {

  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;

  int k = ylon.size();
  Rcpp::NumericVector dist(k);

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
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Vector of distances in meters
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dist_1tom(const double& xlon,
			      const double& xlat,
			      const Rcpp::NumericVector& ylon,
			      const Rcpp::NumericVector& ylat,
			      std::string dist_function="Haversine") {

  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;

  int k = ylon.size();
  Rcpp::NumericVector dist(k);

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
//' @param dist_function String name of distance function: Haversine, Vincenty
//' @return Distance in meters
//' @export
// [[Rcpp::export]]
double dist_1to1(const double& xlon,
		 const double& xlat,
		 const double& ylon,
		 const double& ylat,
		 std::string dist_function="Haversine") {

  // select function
  Rcpp::XPtr<funcPtr> xpfun = choose_func(dist_function);
  funcPtr fun = *xpfun;

  // compute and return
  return fun(xlon, xlat, ylon, ylat);

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
Rcpp::DataFrame popdist_weighted_mean(Rcpp::DataFrame x_df,
				      Rcpp::DataFrame y_df,
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
  Rcpp::CharacterVector id = x_df[x_id];
  Rcpp::NumericVector xlon = x_df[x_lon_col];
  Rcpp::NumericVector xlat = x_df[x_lat_col];
  Rcpp::NumericVector meas = y_df[measure_col];
  Rcpp::NumericVector ylon = y_df[y_lon_col];
  Rcpp::NumericVector ylat = y_df[y_lat_col];
  Rcpp::NumericVector popw = y_df[pop_col];

  int n = xlon.size();
  int k = ylon.size();
  Rcpp::NumericVector out(n);

  // loop to compute
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    Rcpp::NumericVector dist = dist_1tom(xlon[i], xlat[i],
					 ylon, ylat, dist_function);

    // inverse distance weights
    Rcpp::NumericVector idw = inverse_value(dist, decay, dist_transform);

    // population adjusted idw
    Rcpp::NumericVector w = idw * popw;

    // weight denominator
    double w_sum = 0;
    for (int j = 0; j < k; j++) {
      w_sum += w[j];
    }

    // weight the measure_i
    Rcpp::NumericVector tmp = (w * meas) / (w_sum);

    // sum the wmeasure_i
    double sum = 0;
    for (int j = 0; j < k; j++) {
      sum += tmp[j];
    }

    out[i] = sum;

  }

  return Rcpp::DataFrame::create(Rcpp::Named("id") = id,
				 Rcpp::Named("wmeasure") = out,
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
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @param dist_transform String value of distance weight transform: "level" (default)
//' or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @return Dataframe of distance-weighted values
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_weighted_mean(Rcpp::DataFrame x_df,
				   Rcpp::DataFrame y_df,
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
  Rcpp::CharacterVector id = x_df[x_id];
  Rcpp::NumericVector xlon = x_df[x_lon_col];
  Rcpp::NumericVector xlat = x_df[x_lat_col];
  Rcpp::NumericVector meas = y_df[measure_col];
  Rcpp::NumericVector ylon = y_df[y_lon_col];
  Rcpp::NumericVector ylat = y_df[y_lat_col];

  int n = xlon.size();
  int k = ylon.size();
  Rcpp::NumericVector out(n);

  // loop to compute
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    Rcpp::NumericVector dist = dist_1tom(xlon[i], xlat[i],
					 ylon, ylat, dist_function);

    // inverse distance weights
    Rcpp::NumericVector w = inverse_value(dist, decay, dist_transform);

    // weight denominator
    double w_sum = 0;
    for (int j = 0; j < k; j++) {
      w_sum += w[j];
    }

    // weight the measure_i
    Rcpp::NumericVector tmp = (w * meas) / (w_sum);

    // sum the wmeasure_i
    double sum = 0;
    for (int j = 0; j < k; j++) {
      sum += tmp[j];
    }

    out[i] = sum;

  }

  return Rcpp::DataFrame::create(Rcpp::Named("id") = id,
				 Rcpp::Named("wmeasure") = out,
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
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @return DataFrame with id of closest point and distance in meters
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_min(Rcpp::DataFrame x_df,
			 Rcpp::DataFrame y_df,
			 std::string x_id = "id",
			 std::string y_id = "id",
			 std::string x_lon_col = "lon",
			 std::string x_lat_col = "lat",
			 std::string y_lon_col = "lon",
			 std::string y_lat_col = "lat",
			 std::string dist_function = "Haversine") {

  // init
  Rcpp::CharacterVector idx = x_df[x_id];
  Rcpp::CharacterVector idy = y_df[y_id];
  Rcpp::NumericVector xlon = x_df[x_lon_col];
  Rcpp::NumericVector xlat = x_df[x_lat_col];
  Rcpp::NumericVector ylon = y_df[y_lon_col];
  Rcpp::NumericVector ylat = y_df[y_lat_col];

  int n = xlon.size();
  Rcpp::NumericVector dist(n);
  Rcpp::CharacterVector end(n);

  // loop
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    Rcpp::NumericVector distvec = dist_1tom(xlon[i], xlat[i],
					    ylon, ylat, dist_function);

    // add minimum to distance output
    dist[i] = min(distvec);

    // get ID of minimum
    int j;
    j = which_min(distvec);
    end[i] = idy[j];

  }

  return Rcpp::DataFrame::create(Rcpp::Named("id_start") = idx,
				 Rcpp::Named("id_end") = end,
				 Rcpp::Named("meters") = dist,
				 Rcpp::Named("stringsAsFactors") = false);

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
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @return DataFrame with id of farthest point and distance in meters
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_max(Rcpp::DataFrame x_df,
			 Rcpp::DataFrame y_df,
			 std::string x_id = "id",
			 std::string y_id = "id",
			 std::string x_lon_col = "lon",
			 std::string x_lat_col = "lat",
			 std::string y_lon_col = "lon",
			 std::string y_lat_col = "lat",
			 std::string dist_function = "Haversine") {

  // init
  Rcpp::CharacterVector idx = x_df[x_id];
  Rcpp::CharacterVector idy = y_df[y_id];
  Rcpp::NumericVector xlon = x_df[x_lon_col];
  Rcpp::NumericVector xlat = x_df[x_lat_col];
  Rcpp::NumericVector ylon = y_df[y_lon_col];
  Rcpp::NumericVector ylat = y_df[y_lat_col];

  int n = xlon.size();
  Rcpp::NumericVector dist(n);
  Rcpp::CharacterVector end(n);

  // loop
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    Rcpp::NumericVector distvec = dist_1tom(xlon[i], xlat[i],
					    ylon, ylat, dist_function);

    // add maximum to distance output
    dist[i] = max(distvec);

    // get ID of maximum
    int j;
    j = which_max(distvec);
    end[i] = idy[j];

  }

  return Rcpp::DataFrame::create(Rcpp::Named("id_start") = idx,
				 Rcpp::Named("id_end") = end,
				 Rcpp::Named("meters") = dist,
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
//' @param dist_function String name of distance function: "Haversine" (default) or
//' "Vincenty"
//' @param dist_transform String value of distance transform: "level" (default)
//' or "log"
//' @param decay Numeric value of distance weight decay: 2 (default)
//' @param scale_units Double value to divide return value by (e.g., 1000 == km)
//' @return DataFrame with sum of distances
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dist_sum_inv(Rcpp::DataFrame x_df,
			     Rcpp::DataFrame y_df,
			     std::string x_id = "id",
			     std::string y_id = "id",
			     std::string x_lon_col = "lon",
			     std::string x_lat_col = "lat",
			     std::string y_lon_col = "lon",
			     std::string y_lat_col = "lat",
			     std::string dist_function = "Haversine",
			     std::string dist_transform = "level",
			     double decay = 2, 
			     double scale_units = 1) {

  // init
  Rcpp::CharacterVector idx = x_df[x_id];
  Rcpp::CharacterVector idy = y_df[y_id];
  Rcpp::NumericVector xlon = x_df[x_lon_col];
  Rcpp::NumericVector xlat = x_df[x_lat_col];
  Rcpp::NumericVector ylon = y_df[y_lon_col];
  Rcpp::NumericVector ylat = y_df[y_lat_col];

  int n = xlon.size();
  Rcpp::NumericVector dist(n);

  // loop
  for (int i = 0; i < n; i++) {

    // check for interrupt
    if(i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    // distance vector
    Rcpp::NumericVector distvec = dist_1tom(xlon[i], xlat[i],
					    ylon, ylat, dist_function);

    // scale units
    distvec = distvec / scale_units;

    // inverse distance weights
    Rcpp::NumericVector inv_distvec = inverse_value(distvec, decay,
						    dist_transform);

    // replace infinite values with 0
    inv_distvec[is_infinite(inv_distvec)] = 0;

    // add maximum to distance output
    dist[i] = sum(inv_distvec);

  }

  return Rcpp::DataFrame::create(Rcpp::Named("id") = idx,
				 Rcpp::Named("inv_distance") = dist,
				 Rcpp::Named("stringsAsFactors") = false);

}
