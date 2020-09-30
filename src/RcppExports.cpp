// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ADMM_ENrun
List ADMM_ENrun(arma::vec tildelogY, arma::mat X, arma::mat D, arma::mat tildedelta, double rho, double eta, double tau, double lambda, double alpha, arma::vec w, arma::vec Gamma, arma::vec Beta, arma::vec Theta, unsigned int max_iter, double tol_abs, double tol_rel, double gamma, double euc_tildelogY, arma::vec Xbeta, arma::vec tXB, unsigned int n, unsigned int l, unsigned int p);
RcppExport SEXP _penAFT_ADMM_ENrun(SEXP tildelogYSEXP, SEXP XSEXP, SEXP DSEXP, SEXP tildedeltaSEXP, SEXP rhoSEXP, SEXP etaSEXP, SEXP tauSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP wSEXP, SEXP GammaSEXP, SEXP BetaSEXP, SEXP ThetaSEXP, SEXP max_iterSEXP, SEXP tol_absSEXP, SEXP tol_relSEXP, SEXP gammaSEXP, SEXP euc_tildelogYSEXP, SEXP XbetaSEXP, SEXP tXBSEXP, SEXP nSEXP, SEXP lSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type tildelogY(tildelogYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tildedelta(tildedeltaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol_abs(tol_absSEXP);
    Rcpp::traits::input_parameter< double >::type tol_rel(tol_relSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type euc_tildelogY(euc_tildelogYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xbeta(XbetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tXB(tXBSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type l(lSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(ADMM_ENrun(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, Gamma, Beta, Theta, max_iter, tol_abs, tol_rel, gamma, euc_tildelogY, Xbeta, tXB, n, l, p));
    return rcpp_result_gen;
END_RCPP
}
// ADMM_SGrun
List ADMM_SGrun(arma::vec tildelogY, arma::mat X, arma::mat D, arma::mat tildedelta, double rho, double eta, double tau, double lambda, double alpha, arma::vec w, arma::vec v, arma::vec borderIndexes, arma::vec Gamma, arma::vec Beta, arma::vec Theta, unsigned int max_iter, double tol_abs, double tol_rel, double gamma, double euc_tildelogY, arma::vec Xbeta, arma::vec tXB, unsigned int n, unsigned int l, unsigned int p, int G);
RcppExport SEXP _penAFT_ADMM_SGrun(SEXP tildelogYSEXP, SEXP XSEXP, SEXP DSEXP, SEXP tildedeltaSEXP, SEXP rhoSEXP, SEXP etaSEXP, SEXP tauSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP wSEXP, SEXP vSEXP, SEXP borderIndexesSEXP, SEXP GammaSEXP, SEXP BetaSEXP, SEXP ThetaSEXP, SEXP max_iterSEXP, SEXP tol_absSEXP, SEXP tol_relSEXP, SEXP gammaSEXP, SEXP euc_tildelogYSEXP, SEXP XbetaSEXP, SEXP tXBSEXP, SEXP nSEXP, SEXP lSEXP, SEXP pSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type tildelogY(tildelogYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tildedelta(tildedeltaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type borderIndexes(borderIndexesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol_abs(tol_absSEXP);
    Rcpp::traits::input_parameter< double >::type tol_rel(tol_relSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type euc_tildelogY(euc_tildelogYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xbeta(XbetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tXB(tXBSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type l(lSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(ADMM_SGrun(tildelogY, X, D, tildedelta, rho, eta, tau, lambda, alpha, w, v, borderIndexes, Gamma, Beta, Theta, max_iter, tol_abs, tol_rel, gamma, euc_tildelogY, Xbeta, tXB, n, l, p, G));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_penAFT_ADMM_ENrun", (DL_FUNC) &_penAFT_ADMM_ENrun, 23},
    {"_penAFT_ADMM_SGrun", (DL_FUNC) &_penAFT_ADMM_SGrun, 26},
    {NULL, NULL, 0}
};

RcppExport void R_init_penAFT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
