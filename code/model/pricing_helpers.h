#pragma once
#include <Rcpp.h>
Rcpp::List y_fitting_jfec_cpp(const Rcpp::NumericMatrix& X0,
                              const Rcpp::NumericVector& A_X_for,
                              const Rcpp::NumericMatrix& B_X_for,
                              const Rcpp::NumericVector& A_X_exp,
                              const double r_lb,
                              const Rcpp::NumericVector& s_n,
                              const Rcpp::NumericVector& A_X_for_pi,
                              const Rcpp::NumericMatrix& B_X_pi,
                              const Rcpp::NumericMatrix& Sigma2_X,
                              const Rcpp::NumericMatrix& B_X_cum,
                              const Rcpp::NumericMatrix& B_X_cum_pi);
