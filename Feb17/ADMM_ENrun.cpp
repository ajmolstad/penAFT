#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
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

double absolute(double num) {
	if (num >= 0) {
		return num;
	}
	else {
		return -num;
	}
}

// [[Rcpp::export]]
List ADMM_ENrun(const arma::vec& tildelogY, const arma::mat& X, 
    const arma::mat& D_pos, const arma::mat& D_vert_1, const arma::mat& D_vert_neg1, const arma::mat& tildedelta, 
    double rho, double eta, double tau, double lambda,
    double alpha, const arma::vec& w, arma::vec Gamma, const arma::vec& Beta, arma::vec Theta, 
    unsigned int max_iter, double tol_abs, double tol_rel, double gamma,
    double euc_tildelogY)
{

    int p = size(X)(1);
    int n = size(X)(0);
    int l = size(tildelogY)(0);

    double updateStep = 1.0;

    arma::mat lam((pow(n, gamma)*lambda*alpha*w)/(eta));
    arma::mat lam2((pow(n, gamma)*lambda*(1-alpha)*w)/(eta));


    unsigned int lll_counter = 0;

    arma::sp_mat BetaSp(Beta);
    //arma::sp_mat DSp(D);

    //arma::sp_mat DSp_t = D.t();
    arma::mat X_t = X.t();

    arma::vec tXB(l);
    arma::vec XBeta(X * BetaSp);

    arma::vec t0_Theta(l);
    arma::vec D_t0_Theta(n);

    arma::vec Theta_tTheta(l);
    arma::vec D_Theta_tTheta(n);

    arma::vec D_Gamma(n);

    //arma::mat tXB(DSp * (X * BetaSp));
    int ii = 0;
    int jj = 0;

    for (unsigned int k = 0; k < l; k++) 
    {
        ii = D_pos(k, 0) - 1;
        jj = D_pos(k, 1) - 1;
        tXB(k) = XBeta(ii) - XBeta(jj);

    }
    arma::vec tTheta(l);


    for (unsigned int lll = 1; lll <= max_iter; lll++)
    {
        lll_counter++;

        // ---------------------------------
        // Theta update
        // ---------------------------------
        tTheta = Theta;
        double nrho = pow(n, 2 - gamma) * rho;
        arma::vec t0(tildelogY - tXB - ((1/rho) * Gamma));

        for (unsigned int m = 0; m < l; m++)
        {
            Theta(m) = updateThetaEntry(t0(m), tildedelta(m,0), tildedelta(m,1), nrho);
        }

        // -------------------------------------
        // Beta update
        // -------------------------------------


        //BetaSp = ((1 / eta) * X_t * (DSp_t * (t0 - Theta))) + BetaSp;

        t0_Theta = t0 - Theta;

        for (unsigned int ee = 0; ee < n; ee++) 
    	{

        	double sumRow = 0.0;
        	for (unsigned int ff = 0; ff < n; ff++) 
        	{
        		ii = D_vert_1(ee,ff) - 1;
        		jj = D_vert_neg1(ee,ff) - 1;

        		if (ii >= 0) 
        		{
        			sumRow += t0_Theta(ii);
        		}
        		if (jj >= 0)
        		{
        			sumRow -= t0_Theta(jj);
        		}
        	}

        	D_t0_Theta(ee) = sumRow;

    	}

    	BetaSp = ((1 / eta) * X_t * D_t0_Theta) + BetaSp;

        arma::sp_mat signMatrix(BetaSp);

        signMatrix.for_each([](arma::mat::elem_type& val) { val = signum(val); } );

        BetaSp = abs(BetaSp);

        BetaSp = BetaSp - (lam / rho);

        BetaSp.for_each( [](arma::sp_mat::elem_type& val) { val = maximum(val, 0.0); } );

        BetaSp = BetaSp % signMatrix;

        BetaSp = BetaSp / (1 + (lam2 / rho));

        //tXB = DSp * (X * BetaSp);
        XBeta = X * BetaSp;
        for (unsigned int k = 0; k < l; k++) 
        {
            ii = D_pos(k, 0) - 1;
            jj = D_pos(k, 1) - 1;
            tXB(k) = XBeta(ii) - XBeta(jj);
        }

        // ----------------------------------
        // Gamma update
        // ----------------------------------
        Gamma = Gamma + tau*rho*(Theta - tildelogY + tXB);

        //-----------------------------------------------------------
        // Step size update and convergence conditions check
        //-----------------------------------------------------------
        
        if (lll % (int)updateStep == 0)
        {

        	Theta_tTheta = Theta - tTheta;


	        for (unsigned int ee = 0; ee < n; ee++) 
	    	{

	        	double sumRow = 0.0;
	        	for (unsigned int ff = 0; ff < n; ff++) 
	        	{
	        		ii = D_vert_1(ee,ff) - 1;
	        		jj = D_vert_neg1(ee,ff) - 1;

	        		if (ii >= 0) 
	        		{
	        			sumRow += Theta_tTheta(ii);
	        		}
	        		if (jj >= 0)
	        		{
	        			sumRow -= Theta_tTheta(jj);
	        		}
	        	}

	        	D_Theta_tTheta(ee) = sumRow;

	    	}


            //double s = rho * norm(X_t * (DSp_t * (Theta - tTheta)), 2);
            double s = rho * norm(X_t * D_Theta_tTheta, 2);
            double r = norm(Theta - tildelogY + tXB, 2);

            double eprim = sqrt(l) * tol_abs + tol_rel * maximum(norm(tXB, 2), norm(Theta, 2), euc_tildelogY);

            for (unsigned int ee = 0; ee < n; ee++) 
	    	{

	        	double sumRow = 0.0;
	        	for (unsigned int ff = 0; ff < n; ff++) 
	        	{
	        		ii = D_vert_1(ee,ff) - 1;
	        		jj = D_vert_neg1(ee,ff) - 1;

	        		if (ii >= 0) 
	        		{
	        			sumRow += Gamma(ii);
	        		}
	        		if (jj >= 0)
	        		{
	        			sumRow -= Gamma(jj);
	        		}
	        	}

	        	D_Gamma(ee) = sumRow;

	    	}

            double edual = sqrt(p) * tol_abs + tol_rel * norm(X_t * D_Gamma, 2);

            if (r/eprim > 10*s/edual){
                rho = rho*2;
            }

            if (s/edual > 10*r/eprim){
                rho = rho/2;
            }
            if (lll > 10) {
                if (r < eprim && s < edual){
                    break;
                }
            }
            updateStep = (updateStep + 1)*1.1;
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
