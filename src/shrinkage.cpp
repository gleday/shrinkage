//[[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;

////////////////////////////////////////////////////
//-------------- Internal functions --------------//
////////////////////////////////////////////////////

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

colvec rmnorm(const int d) {
  colvec out = zeros(d);
  for(int j = 0; j < d; j++){
    out(j) += R::rnorm(0, 1);
  }
  return out;
}


////////////////////////////////////////////////////
//-------------- External functions --------------//
////////////////////////////////////////////////////

// Future extension for Bayesian ridge:
// - implement empirical Bayes ML-based estimation of \tau^2.
// ( -> add argument prior = c("bp", "ig", "ml")?)

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



// [[Rcpp::export(.bgridge)]]
Rcpp::List bgridge(arma::colvec y, arma::mat X, arma::colvec g, const double a = 0.5, const double b = 0.5, const double c = 1, const int mcmc = 1000, const int  burnin = 1000, const int thin = 10, bool verbose = true, bool light = true){
  
  // Dimension data
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  // Group sizes
  arma::colvec gr = unique(g);
  const int K = gr.n_elem;
  arma::colvec pk = zeros(gr.n_elem);
  for(int u = 0; u < K; u++){
    pk(u) += sum(g == gr(u));
  }

  // Initialisation
  const double cStar = 0.5*(n+p);
  const double ustar = a - 0.5*p;
  double tauminus2 = 1;
  double sigmaminus2 = 1;
  const int nruns = mcmc + burnin;
  arma::colvec wk = ones(K)*1/K;
  arma::colvec d = zeros(p);
  for(int u = 0; u < K; u++){
    d(find(g == gr(u))).fill(pk(u)/wk(u));
  }
  arma::mat D = diagmat(d);

  // Useful quantities
  arma::mat XTX = X.t()*X;
  arma::mat XTy = X.t()*y;
  
  // MCMC samples
  arma::mat betasamp = zeros(mcmc, p);
  arma::mat wksamp = zeros(mcmc, K);
  arma::colvec tau2samp = zeros(mcmc);
  arma::colvec sigma2samp = zeros(mcmc);
  arma::mat betavar, cholvar;
  arma::colvec betamean, resid;
  arma::colvec beta(p);
  double btDb, bktDbk, dStar;
  uvec indk;
  int k = 0;
  arma::mat V;
  
  //  Gibbs algorithm
  for(int i = 0; i < nruns; i++){
    for(int j = 0; j < thin; j++){

      // Sample from P(\beta | ...)
      V = sigmaminus2*(XTX + tauminus2*D);
      betavar = arma::inv_sympd(V);
      betamean = sigmaminus2 * betavar * XTy;
      beta = rmnorm(p);
      cholvar = arma::chol(betavar);
      beta = trans(beta.t() * cholvar);
      beta += betamean;
      
      // Sample from P(\tau^2 | ...)
      btDb = sum(d % square(beta));
      tauminus2 = 1/rgig(ustar, sigmaminus2*btDb, 2*b);
      
      // Sample from P(\omega_k | ...)
      for(int u = 0; u < K; u++){
        indk = find(g == gr(u));
        bktDbk = sum(square(beta(indk)));
        wk(u) = rgig(c - 0.5 * pk(u), sigmaminus2*pk(u)*bktDbk, 2*b);
      }
      wk /= sum(wk);

      // Sample from P(\sigma^2| ...)
      for(int u = 0; u < K; u++){
        d(find(g == gr(u))).fill(pk(u)/wk(u));
      }
      if(d.max() > std::numeric_limits<float>::max()){
        d(find(d > std::numeric_limits<float>::max())).fill(std::numeric_limits<float>::max());
      }
      D = diagmat(d);
      btDb = sum(d % square(beta));
      resid = y-X*beta;
      dStar = 0.5*( sum(square(resid)) + tauminus2*btDb );
      sigmaminus2 = R::rgamma(cStar, 1/dStar);
    }
    
    // Save samples
    if(i >= burnin){
      betasamp.row(k) += beta.t();
      tau2samp(k) += 1/(tauminus2);
      wksamp.row(k) += wk.t();
      sigma2samp(k) += 1/sigmaminus2;
      k++;
      
      if(verbose && (k%100==0)){
        Rcpp::Rcout << k << " samples generated" << std::endl;
      }
    }
    
  }
  
  if(verbose){
    Rcpp::Rcout << "Thanks, that's enough." << std::endl;
  }
  
  // Output
  Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("betas") = betasamp,
      Rcpp::Named("tau2s") = tau2samp,
      Rcpp::Named("ws") = wksamp,
      Rcpp::Named("sigma2s") = sigma2samp
    );

  return out;
}
