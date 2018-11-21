// Extrinsic Local Regression 
// (1) elocal_fit_cv  : measures MSE only 
// (2) elocal_fit     : fit/smoothed list
// (3) elocal_predict : prediction, very similar to (1)

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

#include "RcppArmadillo.h"
#include "riemfactory.hpp"
#include "RiemRegress_auxiliary.h"

//////////////////////////////////////////////////////////
// (1) elocal_fit_cv only measures MSE for a given h value.
//      - cube fitted : fitted points invequiv'ed 
//      - vec  fscore : MSE
// [[Rcpp::export]]
double elocal_fit_cv(double h, arma::cube ytrain, arma::cube ytest, arma::mat xtrain, arma::mat xtest, std::string mfdname){
  // get size parameters
  int ntrain = ytrain.n_slices;
  int ntest  = ytest.n_slices;
  int dnrow  = ytrain.n_rows;
  int dncol  = ytrain.n_cols;
  double h2  = (2.0*h*h); 

  // KH computation
  arma::mat KH(ntrain,ntest,fill::zeros);
  for (int i=0;i<ntrain;i++){
    for (int j=0;j<ntest;j++){
      KH(i,j) = static_cast<double>(std::exp(static_cast<float>(-std::pow(arma::norm(xtrain.row(i)-xtest.row(j), 2), 2.0)/h2)));
    }
  }

  // conversion of 'ytrain' cube
  arma::mat ytrainmat = cpp_equivariant(ytrain, mfdname);
  
  // find the fitted matrices
  arma::cube ytestfit(dnrow,dncol,ntest,fill::zeros);
  arma::rowvec yvectmp(ytrainmat.n_cols,fill::zeros);
  for (int j=0;j<ntest;j++){ // iterate over test samples
    // initialize
    yvectmp.fill(0.0);
    for (int i=0;i<ntrain;i++){
      yvectmp += KH(i,j)*ytrainmat.row(i);
    }
    yvectmp /= arma::sum(KH.col(j));
    ytestfit.slice(j) = riemfunc_invequiv(yvectmp.t(),dnrow,dncol,mfdname);
  }

  // compute the squared intrinsic distances and sum it E_k(h)
  double cvscore = 0.0;
  double tmpdist = 0.0;
  arma::mat mat1(dnrow,dncol,fill::zeros);
  arma::mat mat2(dnrow,dncol,fill::zeros);
  for (int j=0;j<ntest;j++){
    mat1 = ytestfit.slice(j);
    mat2 = ytest.slice(j);
    tmpdist = riemfunc_dist(mat1,mat2,mfdname);
    cvscore += (tmpdist*tmpdist);
  }
  
  // return
  return(cvscore);
}
  
// (2) elocal_fit    : fit/smoothed list
// [[Rcpp::export]]
Rcpp::List elocal_fit(double h, arma::cube y, arma::mat x, std::string mfdname){
  // get parameters
  int n = x.n_rows;
  int dnrow = y.n_rows;
  int dncol = y.n_cols;
  double h2 = (2.0*h*h);
  
  // KH comutation : weight matrix
  arma::mat KH(n,n,fill::zeros);
  for (int i=0;i<n;i++){
    KH(i,i) = 1.0;
  }
  double tmpval;
  for (int i=0;i<(n-1);i++){
    for (int j=(i+1);j<n;j++){
      tmpval  = static_cast<double>(std::exp(static_cast<float>(-std::pow(arma::norm(x.row(i)-x.row(j), 2), 2.0)/h2)));
      KH(i,j) = tmpval;
      KH(j,i) = tmpval;
    }
  }
  
  // conversion : matrix stacked as rows
  arma::mat ymat = cpp_equivariant(y, mfdname);
  int ymatdim = ymat.n_cols;
  arma::rowvec yfit(ymatdim, fill::zeros);
  arma::cube yfitted(dnrow,dncol,n,fill::zeros);
  for (int i=0;i<n;i++){
    yfit.fill(0.0);
    for (int j=0;j<n;j++){
      yfit += KH(i,j)*ymat.row(j);
    }
    yfit /= arma::sum(KH.col(i));
    yfitted.slice(i) = riemfunc_invequiv(yfit.t(),dnrow,dncol,mfdname);
  }
  
  // compute MSE
  double mse = 0.0;
  for (int i=0;i<n;i++){
    mse += static_cast<double>(std::pow(static_cast<float>(riemfunc_dist(y.slice(i),yfitted.slice(i),mfdname)),2.0));
  }
  mse /= static_cast<double>(n);
  
  // return output
  Rcpp::List output;
  output["fitted"] = yfitted;
  output["mse"] = mse;
  return(output);
}
  

// (3) elocal_predict : prediction, very similar to (1)
// [[Rcpp::export]]
arma::cube elocal_predict(double h, arma::cube y, arma::mat x, arma::mat xpredict, std::string mfdname){
  // get size parameters
  int dnrow = y.n_rows;
  int dncol = y.n_cols;
  int nold  = x.n_rows;
  int nnew  = xpredict.n_rows;
  double h2 = (2.0*h*h);
  
  // KH : weight matrix computation
  arma::mat KH(nold,nnew,fill::zeros);
  for (int i=0;i<nold;i++){
    for (int j=0;j<nnew;j++){
      KH(i,j) = static_cast<double>(std::exp(static_cast<float>(-std::pow(arma::norm(x.row(i)-xpredict.row(j), 2), 2.0)/h2)));
    }
  }
  
  // conversion : matrix stacked as rows
  arma::mat ymat = cpp_equivariant(y, mfdname);
  int ydim = ymat.n_cols;
  
  // prediction
  arma::cube yfitted(dnrow,dncol,nnew,fill::zeros);
  arma::rowvec ytmp(ydim,fill::zeros);
  
  for (int j=0;j<nnew;j++){
    ytmp.fill(0.0);
    for (int i=0;i<nold;i++){
      ytmp += KH(i,j)*ymat.row(i);
    }
    ytmp /= arma::sum(KH.col(j));
    yfitted.slice(j) = riemfunc_invequiv(ytmp.t(),dnrow,dncol,mfdname);
  }
  return(yfitted);
}