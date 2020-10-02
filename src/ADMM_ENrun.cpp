#include "RcppArmadillo.h"

#include <cmath>
using namespace Rcpp;

double updateThetaEntry(double temp, double tildedelta_1, double tildedelta_2, double nrho) {
    double theta_entry = 0;
    if (temp > tildedelta_2 / nrho)
    {
        theta_entry = temp - (tildedelta_2 / nrho);
    }
    else if (temp < - (tildedelta_1 / nrho))
    {
        theta_entry = temp + (tildedelta_1 / nrho);
    }
    return theta_entry;
}


int signum(double num) {
  if (num < 0) {
    return -1;
  }
  else if (num > 0) {
    return 1;
  }
  return 0;
}

double maximum(double num1, double num2) {
  if (num1 > num2) {
    return num1;
  }
  else {
    return num2;
  }
}

double maximum(double num1, double num2, double num3) {
    if (num1 >= num2) {
        return(maximum(num1, num3));
    }
    else {
        return(maximum(num2, num3));
    }
}

// [[Rcpp::export]]
List ADMM_ENrun(arma::vec tildelogY, arma::mat X, arma::mat D, arma::mat tildedelta, double rho, double eta, double tau, double lambda,
    double alpha, arma::vec w, arma::vec Gamma, arma::vec Beta, arma::vec Theta, unsigned int max_iter, double tol_abs, double tol_rel, double gamma,
    double euc_tildelogY, arma::vec Xbeta, arma::vec tXB, unsigned int n, unsigned int l, unsigned int p)
{

    int updateStep = 1;
    unsigned int lll_counter = 0;
    arma::sp_mat BetaSp(Beta);
    arma::sp_mat DSp(D);

    arma::sp_mat BetaPrev(BetaSp);

    arma::sp_mat DSp_t = DSp.t();
    arma::mat X_t = X.t();

    arma::vec tTheta(l);
    arma::mat w_fact_alpha(w);
    arma::mat w_fact_1_alpha(w);
    arma::sp_mat s0(p, 1);
    arma::mat signMatrix(w);
    arma::vec temp(l);

    for (unsigned int lll = 1; lll <= max_iter; lll++)
    {
        lll_counter++;

        // ---------------------------------
        // Theta update
        // ---------------------------------
        tTheta = Theta;
        double nrho = pow(n, 2 - gamma) * rho;
        temp = (tildelogY - tXB - ((1/rho) * Gamma));

        for (unsigned int m = 0; m < l; m++)
        {
            Theta(m) = updateThetaEntry(temp(m), tildedelta(m,0), tildedelta(m,1), nrho);
        }

        // -------------------------------------
        // Beta update
        // -------------------------------------

        do
        {
            BetaPrev = BetaSp;
            double fact_alpha = pow(n, gamma) * lambda * alpha / (rho * eta);
            double fact_1_alpha = pow(n, gamma) * lambda * (1 - alpha) / (rho * eta);

            w_fact_alpha = fact_alpha * w;
            w_fact_1_alpha = (fact_1_alpha * w) + 1;
            s0 = (X_t * (DSp_t * (tildelogY - Theta - ((1/rho) * Gamma) - tXB)));
            BetaSp = ((1/eta) * s0) + BetaSp;

            BetaSp = abs(BetaSp) - w_fact_alpha;
            signMatrix = w;
            signMatrix.for_each([](arma::mat::elem_type& val) { val = signum(val); } );

            BetaSp.for_each( [](arma::sp_mat::elem_type& val) { val = maximum(val, 0.0); } );

            BetaSp = BetaSp % signMatrix;
            BetaSp = BetaSp / w_fact_1_alpha;

            tXB = DSp * (X * BetaSp);

        } while(norm(BetaPrev - BetaSp, 2) > 1e-2);

        // ----------------------------------
        // Gamma update
        // ----------------------------------
        Gamma = Gamma + tau*rho*(Theta - tildelogY + tXB);


        //-----------------------------------------------------------
        // Step size update and convergence conditions check
        //-----------------------------------------------------------
        if ((lll <= 2000) && (lll % (int)(1.6*updateStep) == 0))
        {
            double s = rho * norm(X_t * (DSp_t * (Theta - tTheta)), 2);
            double r = norm(Theta - tildelogY + tXB, 2);

            double eprim = sqrt(l) * tol_abs + tol_rel * maximum(norm(tXB, 2), norm(Theta, 2), euc_tildelogY);
            double edual = sqrt(p) * tol_abs + tol_rel * norm(X_t * (DSp_t * Gamma), 2);

            if (r/eprim > 10*s/edual){
                rho = rho*2;
            }

            if (s/edual > 10*r/eprim){
                rho = rho/2;
            }

            if(r < eprim && s < edual){
                break;
            }
            updateStep += 1;
        }
    }

    arma::vec BetaOut(BetaSp);
    arma::vec ThetaOut(l);
    arma::vec GammaOut(l);

    for (unsigned int i = 0; i < l; i++)
    {
        ThetaOut(i) = Theta(i);
    }

    for (unsigned int i = 0; i < l; i++)
    {
        GammaOut(i) = Gamma(i);
    }

    return List::create(Named("Beta") = wrap(BetaOut),
                      Named("Theta") = wrap(ThetaOut),
                      Named("Gamma") = wrap(GammaOut),
                      Named("rho") = wrap(rho),
                      Named("iter.counter") = wrap(lll_counter));
}
