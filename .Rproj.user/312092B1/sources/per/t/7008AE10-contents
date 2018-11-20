// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"
#include "RiemRegress_auxiliary.h"


/*
* (01) cpp_equivariant : convert manifold-valued data into embedded ones
*/

////////////////////////////////////////////////////////////////////////
// (01) cpp_equivariant : convert manifold-valued data into embedded ones.
// [[Rcpp::export]]
arma::mat cpp_equivariant(arma::cube X, std::string mfdname){
  // taking one exempler data
  arma::mat testdata = X.slice(0);
  // get some relevant integer parameters
  int N = X.n_slices;
  int dnrow = testdata.n_rows;
  int dncol = testdata.n_cols;
  // apply equivariant embedding for the testdata to get the size.
  arma::vec testequiv = riemfunc_equiv(testdata,dnrow,dncol,mfdname);
  int P = testequiv.n_elem;
  
  // now apply for each data, stacked as rows
  arma::mat output(N,P,fill::zeros);
  for (int i=0;i<N;i++){
    testdata  = X.slice(i);
    testequiv = riemfunc_equiv(testdata,dnrow,dncol,mfdname);
    output.row(i) = testequiv.t();
  }
  return(output);
}
