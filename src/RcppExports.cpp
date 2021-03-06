// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cpp_equivariant
arma::mat cpp_equivariant(arma::cube X, std::string mfdname);
RcppExport SEXP _RiemRegress_cpp_equivariant(SEXP XSEXP, SEXP mfdnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type X(XSEXP);
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_equivariant(X, mfdname));
    return rcpp_result_gen;
END_RCPP
}
// elocal_fit_cv
double elocal_fit_cv(double h, arma::cube ytrain, arma::cube ytest, arma::mat xtrain, arma::mat xtest, std::string mfdname);
RcppExport SEXP _RiemRegress_elocal_fit_cv(SEXP hSEXP, SEXP ytrainSEXP, SEXP ytestSEXP, SEXP xtrainSEXP, SEXP xtestSEXP, SEXP mfdnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ytrain(ytrainSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ytest(ytestSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xtrain(xtrainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xtest(xtestSEXP);
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    rcpp_result_gen = Rcpp::wrap(elocal_fit_cv(h, ytrain, ytest, xtrain, xtest, mfdname));
    return rcpp_result_gen;
END_RCPP
}
// elocal_fit
Rcpp::List elocal_fit(double h, arma::cube y, arma::mat x, std::string mfdname);
RcppExport SEXP _RiemRegress_elocal_fit(SEXP hSEXP, SEXP ySEXP, SEXP xSEXP, SEXP mfdnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    rcpp_result_gen = Rcpp::wrap(elocal_fit(h, y, x, mfdname));
    return rcpp_result_gen;
END_RCPP
}
// elocal_predict
arma::cube elocal_predict(double h, arma::cube y, arma::mat x, arma::mat xpredict, std::string mfdname);
RcppExport SEXP _RiemRegress_elocal_predict(SEXP hSEXP, SEXP ySEXP, SEXP xSEXP, SEXP xpredictSEXP, SEXP mfdnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xpredict(xpredictSEXP);
    Rcpp::traits::input_parameter< std::string >::type mfdname(mfdnameSEXP);
    rcpp_result_gen = Rcpp::wrap(elocal_predict(h, y, x, xpredict, mfdname));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _RiemRegress_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _RiemRegress_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _RiemRegress_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _RiemRegress_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RiemRegress_cpp_equivariant", (DL_FUNC) &_RiemRegress_cpp_equivariant, 2},
    {"_RiemRegress_elocal_fit_cv", (DL_FUNC) &_RiemRegress_elocal_fit_cv, 6},
    {"_RiemRegress_elocal_fit", (DL_FUNC) &_RiemRegress_elocal_fit, 4},
    {"_RiemRegress_elocal_predict", (DL_FUNC) &_RiemRegress_elocal_predict, 5},
    {"_RiemRegress_rcpparma_hello_world", (DL_FUNC) &_RiemRegress_rcpparma_hello_world, 0},
    {"_RiemRegress_rcpparma_outerproduct", (DL_FUNC) &_RiemRegress_rcpparma_outerproduct, 1},
    {"_RiemRegress_rcpparma_innerproduct", (DL_FUNC) &_RiemRegress_rcpparma_innerproduct, 1},
    {"_RiemRegress_rcpparma_bothproducts", (DL_FUNC) &_RiemRegress_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RiemRegress(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
