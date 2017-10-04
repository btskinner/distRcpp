#ifndef DISTRCPP_DIST_H
#define DISTRCPP_DIST_H

Rcpp::NumericMatrix dist_mtom(const Rcpp::NumericVector& xlon,
			      const Rcpp::NumericVector& xlat,
			      const Rcpp::NumericVector& ylon,
			      const Rcpp::NumericVector& ylat,
			      std::string dist_function);

Rcpp::NumericVector dist_df(const Rcpp::NumericVector& xlon,
			    const Rcpp::NumericVector& xlat,
			    const Rcpp::NumericVector& ylon,
			    const Rcpp::NumericVector& ylat,
			    std::string dist_function);

Rcpp::NumericVector dist_1tom(const double& xlon,
			      const double& xlat,
			      const Rcpp::NumericVector& ylon,
			      const Rcpp::NumericVector& ylat,
			      std::string dist_function);

double dist_1to1(const double& xlon,
		 const double& xlat,
		 const double& ylon,
		 const double& ylat,
		 std::string dist_function);

Rcpp::NumericVector inverse_dist_weight(const Rcpp::NumericVector& d,
					double exp,
					std::string transform);

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
				      double decay = 2);

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
				   double decay = 2);

Rcpp::DataFrame dist_min(Rcpp::DataFrame x_df,
			 Rcpp::DataFrame y_df,
			 std::string x_id = "id",
			 std::string y_id = "id",
			 std::string x_lon_col = "lon",
			 std::string x_lat_col = "lat",
			 std::string y_lon_col = "lon",
			 std::string y_lat_col = "lat",
			 std::string dist_function = "Haversine");

Rcpp::DataFrame dist_max(Rcpp::DataFrame x_df,
			 Rcpp::DataFrame y_df,
			 std::string x_id = "id",
			 std::string y_id = "id",
			 std::string x_lon_col = "lon",
			 std::string x_lat_col = "lat",
			 std::string y_lon_col = "lon",
			 std::string y_lat_col = "lat",
			 std::string dist_function = "Haversine");

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
			     double scale_units = 1);

#endif

