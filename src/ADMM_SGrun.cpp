#include "RcppArmadillo.h"

#include <cmath>
// #include <math.h>
using namespace Rcpp;

double updateThetaEntrySG(double temp, double tildedelta_1, double tildedelta_2, double nrho) {
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


double absoluteSG(double num) {
   if (num >= 0) {
       return num;
   }
   else {
       return -num;
   }
}


// [[Rcpp::export]]
List ADMM_SGrun(arma::vec tildelogY, arma::mat X, arma::mat D, arma::mat tildedelta, double rho, double eta, double tau, double lambda,
    double alpha, arma::vec w, arma::vec v, arma::vec borderIndexes, arma::vec Gamma, arma::vec Beta, arma::vec Theta, int max_iter, double tol_abs, double tol_rel, double gamma,
    double euc_tildelogY, arma::vec Xbeta, arma::vec tXB, int n, int l, int p, int G)
{

    int updateStep = 1;
    //double obj_min = 1e300; // or std::numeric_limits<double>::max();

    unsigned int lll_counter = 0;
    arma::sp_mat BetaSp(Beta);
    arma::sp_mat DSp(D);

    for (unsigned int lll = 1; lll <= max_iter; lll++)
    {
        lll_counter++;
        // ---------------------------------
        // Theta update
        // ---------------------------------
        arma::vec tTheta = Theta;
        double nrho = pow(n, 2 - gamma) * rho;
        arma::vec temp = tildelogY - tXB - ((1/rho) * Gamma);

        for (unsigned int m = 0; m < l; m++)
        {
            Theta(m) = updateThetaEntrySG(temp(m), tildedelta(m,0), tildedelta(m,1), nrho);
        }

        // -------------------------------------
        // Beta update
        // -------------------------------------

        double fact_alpha = pow(n, gamma) * lambda * alpha / (rho * eta);
        double fact_1_alpha = pow(n, gamma) * lambda * (1 - alpha) / (rho * eta);

        arma::mat w_fact_alpha = fact_alpha * w;


        arma::sp_mat s0(X.t() * (DSp.t() * (tildelogY - Theta - ((1/rho) * Gamma) - tXB)));
        arma::sp_mat A = ((1/eta) * s0) + BetaSp;

        arma::sp_mat softVec(arma::abs(A) - w_fact_alpha);
        softVec.for_each( [](arma::sp_mat::elem_type& val) { val = maximumSG(val, 0.0); } );

        arma::mat signMatrix(A);
        signMatrix.for_each([](arma::mat::elem_type& val) { val = signumSG(val); } );

        softVec = softVec % signMatrix; //ELEMENTWISE MULTIPLICATION

        int i = 0;
        int j = 0;
        double softDenom = 0.0;
        double r0 = 0.0;
        double r1 = 0.0;

        for (int g = 0; g < G - 1; g++)
        {
            i = borderIndexes(g) - 1;
            j = borderIndexes(g+1) - 2;

            softDenom = norm(softVec.submat(i, 0, j, 0), 2);

            if (abs(softDenom) < 1e-30)
            {
                r1 = 0.0;
            }
            else
            {
                r0 = (v(g) * fact_1_alpha) / softDenom;
                r1 = maximumSG(1 - r0, 0);
            }

            BetaSp.submat(i, 0, j, 0) = r1 * softVec.submat(i, 0, j, 0);
        }


        tXB = DSp * (X * BetaSp);
        // ----------------------------------
        // Gamma update
        // ----------------------------------
        Gamma = Gamma + tau*rho*(Theta - tildelogY + tXB);

        if ((lll <= 1000) && (lll % (int)(2*updateStep) == 0))
        {
            double s = rho * norm(X.t() * (DSp.t() * (Theta - tTheta)), 2);
            double r = norm(Theta - tildelogY + tXB, 2);

            double eprim = sqrt(l) * tol_abs + tol_rel * maximumSG(norm(tXB, 2), norm(Theta, 2), euc_tildelogY);
            double edual = sqrt(p) * tol_abs + tol_rel * norm(X.t() * (DSp.t() * Gamma), 2);

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