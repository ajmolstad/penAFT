#include "RcppArmadillo.h"

#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

int signumSG(double num) {
  if (num < 0) {
    return -1;
  }
  else if (num > 0) {
    return 1;
  }
  return 0;
}

double maximumSG(double num1, double num2) {
 if (num1 > num2) {
   return num1;
  }
  else {
    return num2;
  }
}

double maximumSG(double num1, double num2, double num3) {
    if (num1 >= num2) {
        return(maximumSG(num1, num3));
    }
    else {
        return(maximumSG(num2, num3));
    }
}

// [[Rcpp::export]]
List ADMM_SGrun(arma::vec tildelogY, arma::mat X, arma::sp_mat D, arma::mat tildedelta, double rho, double eta, double tau, double lambda,
    double alpha, arma::vec w, arma::vec v, arma::vec borderIndexes, arma::vec Gamma, arma::vec Beta, arma::vec Theta, unsigned int max_iter, double tol_abs, double tol_rel, double gamma,
    double euc_tildelogY, unsigned int n, unsigned int l, unsigned int p, int G)
{

    double updateStep = 1.0;

    unsigned int lll_counter = 0;
    arma::sp_mat BetaSp(Beta);
    arma::sp_mat DSp(D);

    arma::sp_mat DSp_t = DSp.t();
    arma::mat X_t = X.t();

    arma::sp_mat BetaPrev(BetaSp);

    arma::vec tTheta(l);
    arma::vec tGamma(l);

    arma::mat signMatrix(w);

    arma::sp_mat tXB(DSp * (X * BetaSp));

    double lam = (pow(n, gamma)*lambda*alpha) / eta;

    arma::vec outTheta(l);

    for (unsigned int lll = 1; lll <= max_iter; lll++)
    {
        lll_counter++;

        // ---------------------------------
        // Theta update
        // ---------------------------------
        tTheta = Theta;
        double nrho = pow(n, 2 - gamma) * rho;
        arma::mat t0(tildelogY - tXB - ((1/rho) * Gamma));

        arma::mat tildedelta_nrho = (tildedelta / pow(n, 2 - gamma)) / rho;

        for (unsigned int m = 0; m < tildedelta_nrho.n_rows; m++) 
        {
            if (t0(m) > tildedelta_nrho(m,1)) 
            {
                Theta(m) = t0(m) - tildedelta_nrho(m,1);
            } 
            else 
            {
                if (t0(m) < -tildedelta_nrho(m,0)) 
                {
                    Theta(m) = t0(m) + tildedelta_nrho(m,0);
                }
                else
                {
                    Theta(m) = 0.0;
                }
            }
        }

        outTheta = Theta;

        BetaPrev = BetaSp;

        arma::sp_mat W((X_t * (DSp_t * (t0 - Theta)))/eta + BetaSp);

        int i = 0;
        int j = 0;
        for (int g = 0; g < G; g++)
        {
            i = borderIndexes(g) - 1;
            j = borderIndexes(g+1) - 2;

            arma::mat h(arma::abs(W.submat(i, 0, j, 0)) - (lam / rho) * w.submat(i, 0, j, 0));
            h.for_each( [](arma::sp_mat::elem_type& val) { val = maximumSG(val, 0.0); } );
            arma::mat signMatrix(W.submat(i, 0, j, 0));
            signMatrix.for_each([](arma::mat::elem_type& val) { val = signumSG(val); } );

            arma::mat h0(h % signMatrix);
            double h1 = norm(h0, 2);

            if (h1 > 0)
            {
                BetaSp.submat(i, 0, j, 0) = maximumSG(1 - v(g) * lambda * (1 - alpha) / (eta * rho * h1), 0) * h0;
            }

        }

        tXB = DSp * (X * BetaSp);

        // ----------------------------------
        // Gamma update
        // ----------------------------------
        Gamma = Gamma + tau*rho*(Theta - tildelogY + tXB);

        tGamma = Gamma;

        //-----------------------------------------------------------
        // Step size update and convergence conditions check
        //-----------------------------------------------------------
        if (lll % (int)(updateStep) == 0)
        {
            double s = rho * norm(X_t * (DSp_t * (Theta - tTheta)), 2);
            double r = norm(Theta - tildelogY + tXB, 2);

            double eprim = sqrt(l) * tol_abs + tol_rel * maximumSG(norm(tXB, 2), norm(Theta, 2), euc_tildelogY);
            double edual = sqrt(p) * tol_abs + tol_rel * norm(X_t * (DSp_t * Gamma), 2);

            if (r/eprim > 10*s/edual)
            {
                rho = rho*2;
            }

            if (s/edual > 10*r/eprim)
            {
                rho = rho/2;
            }

            if (lll > 10)
            {
                if (r < eprim && s < edual)
                {
                    break;
                }
            }
            updateStep = (updateStep + 1) * 1.1;
        }

    }

    arma::vec BetaOut(BetaSp);
    arma::vec ThetaOut(l);
    arma::vec GammaOut(l);

    for (unsigned int i = 0; i < l; i++)
    {
        ThetaOut(i) = outTheta(i);
    }

    for (unsigned int i = 0; i < l; i++)
    {
        GammaOut(i) = tGamma(i);
    }

    return List::create(Named("Beta") = wrap(BetaOut),
                      Named("Theta") = wrap(ThetaOut),
                      Named("Gamma") = wrap(GammaOut),
                      Named("rho") = wrap(rho),
                      Named("iter.counter") = wrap(lll_counter));
}