#ifndef _HermiteInterpolation
#define _HermiteInterpolation 1

#include <Rcpp.h>

double pchip(Rcpp::NumericVector knot_x, Rcpp::NumericVector knot_y, double x);
Rcpp::NumericVector pchip(Rcpp::NumericVector knot_x, Rcpp::NumericVector knot_y, Rcpp::NumericVector x);

#endif

