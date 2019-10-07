
//[[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;

//----------- Internal functions -----------//

// rgig is imported from GIGrvg
double rgig(double lambda, double chi, double psi) {
 PutRNGstate();
 SEXP (*fun)(SEXP, SEXP, SEXP, SEXP) = NULL;
 if (!fun) {
  fun = (SEXP(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("GIGrvg", "rgig");
 }
 GetRNGstate();
 return as<double>(fun(wrap(1), wrap(lambda), wrap(chi), wrap(psi)));
}

//----------- External functions -----------//

// Bayesian ridge regression
// Future extension for bridge():
// - consider empirical Bayes estimation of \tau^2.
// ( -> add argument prior = c("bp", "ig", "eb")?)

// [[Rcpp::export(.bridge)]]
Rcpp::List bridge(arma::colvec y, arma::mat X, const int prior, const double a = 0.5, const double b = 0.5, const int mcmc = 1000, const int  burnin = 1000, const int thin = 10, bool verbose = true){

	// Dimension data
	const int n = X.n_rows;
	const int p = X.n_cols;

	// Initialisation
	const double aStar = a + 0.5*p;
	const double cStar = 0.5*(n+p);
	const double rStar = b + 0.5*p;
	const double ustar = a - 0.5*p;
	double gammaminus2 = 1;
	double deltaminus2 = 1;
	double tauminus2 = 1;
	double sigmaminus2 = 1;
	const int nruns = mcmc + burnin;

	// SVD decomposition
	arma::mat u;
	arma::colvec d;
	arma::mat v;
	if(n >= p){
	  arma::mat XTX = X.t()*X;
		svd_econ(u, d, v, XTX, "right");
		d = sqrt(d);
		u = X*v*diagmat(1/d);
	}
	else{
	  arma::mat XXT = X*X.t();
		svd_econ(u, d, v, XXT, "left");
		d = sqrt(d);
		v = X.t()*u*diagmat(1/d);
	}

	// Useful quantities
	arma::colvec duy = d % (u.cols(0, d.n_elem-1).t()*y);
	arma::mat ud = u * diagmat(d);
	u.clear();
	
	// MCMC samples
	arma::mat thetasamp = zeros(mcmc, d.n_elem);
	arma::colvec tau2samp = zeros(mcmc);
	arma::colvec sigma2samp = zeros(mcmc);

	arma::colvec thetavar;
	arma::colvec thetamean;
	arma::colvec theta(d.n_elem);
	arma::mat Xbeta;
	double btb;
	double dStar;
	int k = 0;

	//  Gibbs algorithm
	for(int i = 0; i < nruns; i++){

		for(int j = 0; j < thin; j++){

			// Sample from P(\theta | ...)
			thetavar = 1/(square(d) + tauminus2);
			thetamean = thetavar % duy;
			thetavar /= sigmaminus2;
			for(uword l = 0; l < thetamean.n_elem; l++){
				theta(l) = R::rnorm(thetamean(l), sqrt(thetavar(l)));
			}

			// Sample from P(\tau^2 | ...)
			btb = sum(square(theta));
			switch(prior){
			case 1:
  			gammaminus2 = 1/rgig(ustar, deltaminus2*sigmaminus2*btb, 2);
  			deltaminus2 = R::rgamma(rStar, 1/(0.5*gammaminus2*sigmaminus2*btb + 1));
  			tauminus2 = gammaminus2*deltaminus2;
			case 2:
			  tauminus2 = R::rgamma(aStar, 1/(b + 0.5*sigmaminus2*btb));
			}
			
			// Sample from P(\sigma^2| ...)
			Xbeta = ud * theta;
			dStar = 0.5*( sum(square(y-Xbeta)) + tauminus2*btb );
			sigmaminus2 = R::rgamma(cStar, 1/dStar);

		}

		// Save samples
		if(i >= burnin){

			thetasamp.row(k) += theta.t();
			tau2samp(k) += 1/(tauminus2);
			sigma2samp(k) += 1/sigmaminus2;
			k++;

			if(verbose && (k%1000==0)){
				Rcpp::Rcout << k << " samples generated" << std::endl;
			}
		}

	}

	if(verbose){
		Rcpp::Rcout << "Thanks, that's enough samples." << std::endl;
	}

	// Get sample from P(\beta | ...)
	arma::mat betasamp = trans(v * thetasamp.t());
	v.clear();

	// Output
	Rcpp::List out = Rcpp::List::create(
			Rcpp::Named("betas") = betasamp,
			Rcpp::Named("tau2s") = tau2samp,
			Rcpp::Named("sigma2s") = sigma2samp
		);

	return out;
}



