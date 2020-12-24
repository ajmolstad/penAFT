#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
double softThreshold(double z,double gam){
	double sparse=0;
	if(z>0 && gam<fabs(z)) return(z-gam);
	if(z<0 && gam<fabs(z)) return(z+gam);
	if(gam>=fabs(z)) return(sparse);
	else return(0);
}


//[[Rcpp::export]]
arma::sp_mat SpMatMult(const arma::sp_mat& a, const arma::sp_mat& b) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a * b);
    return result;

}

//[[Rcpp::export]]
arma::sp_mat SpMatMultDZ(const arma::sp_mat& D, const arma::mat& Z) {
    // sparse x sparse -> sparse
    arma::sp_mat result(D.t() * Z);
    return result;

}

//[[Rcpp::export]]
arma::sp_mat SpMatMultDS(const arma::sp_mat& D, const arma::sp_mat& S) {
    // sparse x sparse -> sparse
    arma::sp_mat result(D.t() * S);
    return result;

}


//[[Rcpp::export]]
arma::sp_mat SpMatMultXB(const arma::sp_mat& B, const arma::mat& X) {
    // sparse x dense -> sparse
    arma::sp_mat result(X * B);
    return result;
}



//[[Rcpp::export]]
arma::vec ThetaUpdate(const arma::mat& temp, const arma::mat& tildedelta_nrho) {
   	
   	arma::vec Theta(tildedelta_nrho.n_rows); Theta.zeros(); 

   	for (unsigned int m = 0; m < tildedelta_nrho.n_rows; m++) {
   		if (temp(m) > tildedelta_nrho(m,1)) {
        	Theta(m) = temp(m) - tildedelta_nrho(m,1);
   		} else {
   			if (temp(m) < -tildedelta_nrho(m,0)) {
   				Theta(m) = temp(m) + tildedelta_nrho(m,0);
   			}
   		}
   	}

   	return Theta;

}


// List BetaUpdate(const arma::sp_mat& D, arma::mat X, arma::mat Z, arma::sp_mat Beta, double eta)






