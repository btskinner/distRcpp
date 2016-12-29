#ifndef DISTRCPP_DIST_H
#define DISTRCPP_DIST_H

#include <inst/shared.h>
#include <Rcpp.h>

using namespace Rcpp;

NumericMatrix dist_mtom(const NumericVector& xlon,
			const NumericVector& xlat,
			const NumericVector& ylon,
			const NumericVector& ylat,
			std::string funname);

NumericVector dist_df(const NumericVector& xlon,
		      const NumericVector& xlat,
		      const NumericVector& ylon,
		      const NumericVector& ylat,
		      std::string funname);

NumericVector dist_1tom(const double& xlon,
			const double& xlat,
			const NumericVector& ylon,
			const NumericVector& ylat,
			std::string funname);

double dist_1to1(const double& xlon,
		 const double& xlat,
		 const double& ylon,
		 const double& ylat,
		 std::string funname);

NumericVector inverse_dist_weight(const NumericVector& d,
				  double exp,
				  std::string transform);

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
				double decay = 2);

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
			     double decay = 2);

DataFrame dist_min(DataFrame x_df,
		   DataFrame y_df,
		   std::string x_id = "id",
		   std::string y_id = "id",
		   std::string x_lon_col = "lon",
		   std::string x_lat_col = "lat",
		   std::string y_lon_col = "lon",
		   std::string y_lat_col = "lat",
		   std::string dist_function = "Haversine");

DataFrame dist_max(DataFrame x_df,
		   DataFrame y_df,
		   std::string x_id = "id",
		   std::string y_id = "id",
		   std::string x_lon_col = "lon",
		   std::string x_lat_col = "lat",
		   std::string y_lon_col = "lon",
		   std::string y_lat_col = "lat",
		   std::string dist_function = "Haversine");

#endif

