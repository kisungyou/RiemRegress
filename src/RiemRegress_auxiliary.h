#ifndef RiemRegress_AUXILIARY_H
#define RiemRegress_AUXILIARY_H

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"

arma::mat cpp_equivariant(arma::cube X, std::string mfdname);

#endif