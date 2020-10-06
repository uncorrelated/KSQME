#ifndef _HermiteInterpolation
#define _HermiteInterpolation 1

#include <Rcpp.h>
#include <RcppParallel.h>

double pchip(Rcpp::NumericVector knot_x, Rcpp::NumericVector knot_y, double x);
double pchip(const double *knot_x, const double *knot_y, const int kn, const double x);
Rcpp::NumericVector pchip(Rcpp::NumericVector knot_x, Rcpp::NumericVector knot_y, Rcpp::NumericVector x);

#endif

