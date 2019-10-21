//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(BH)]]
#include <R.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/tools/precision.hpp>

using namespace arma;
using namespace Rcpp;
using boost::math::tools::brent_find_minima;
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::eps_tolerance;


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

double logML(const double tauminus2, const int p, const int n, const double yTy, colvec eigs, colvec thet){
  double out = -0.5*n*log(datum::pi);
  out += 0.5*p*log(tauminus2);
  out += lgamma(0.5*n);
  colvec egs = zeros(p);
  egs.subvec(0, eigs.n_elem-1) += eigs;
  out -= 0.5*sum(log(egs + tauminus2));
  out -= 0.5*n*log(0.5*yTy - 0.5*sum((square(thet)%square(eigs))/(eigs + tauminus2)));
  return(out);
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
	const double ustar = a - 0.5*p;
	const double vstar = -0.5 - 0.5*p;
	const double wstar = b/(a*a);
  const double xstar = 2*b;
	double gamma2 = 1;
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

		  if((i == 0) && (j == 0) && (burnin>0)){
		    Rcpp::Rcout << "burnin..." << std::endl;
		  }
		  
			// Sample from P(\theta | ...)
			thetavar = 1/(square(d) + tauminus2);
			thetamean = thetavar % duy;
			thetavar /= sigmaminus2;
			for(uword l = 0; l < thetamean.n_elem; l++){
				theta(l) = R::rnorm(thetamean(l), sqrt(thetavar(l)));
			}

			// Sample from P(\tau^2 | ...)
			btb = sum(square(theta));
			if(prior == 1){ // if \tau^2 ~ invGamma(a, b)
			  tauminus2 = R::rgamma(aStar, 1/(b + 0.5*sigmaminus2*btb));
			}
			if(prior == 2){ // if \tau^2 ~ BetaPrime(a, b)
			  tauminus2 = 1/rgig(ustar, sigmaminus2*btb, 2*gamma2);
			  gamma2 = R::rgamma(a+b, 1/(1/tauminus2 + 1));
			}
			if(prior == 3){ // if \tau^2 ~ invGaussian(a, b)
			  tauminus2 = 1/rgig(vstar, b + sigmaminus2*btb, wstar);
			}
			if(prior == 4){ // if \tau^2 ~ Gamma(a, b)
			  tauminus2 = 1/rgig(ustar, sigmaminus2*btb, xstar);
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

// [[Rcpp::export(.bridge_fixed)]]
Rcpp::List bridge_fixed(arma::colvec y, arma::mat X){ 
  
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  // SVD decomposition
  arma::mat u;
  arma::colvec d, d2;
  arma::mat v;
  if(n >= p){
    arma::mat XTX = X.t()*X;
    svd_econ(u, d2, v, XTX, "right");
    d = sqrt(d2);
    u = X*v*diagmat(1/d);
  }
  else{
    arma::mat XXT = X*X.t();
    svd_econ(u, d2, v, XXT, "left");
    d = sqrt(d2);
    v = X.t()*u*diagmat(1/d);
  }
  
  // Useful quantities
  arma::colvec thetahat = (1/d) % (u.cols(0, d.n_elem-1).t()*y);
  u.clear();
  const double yTy = sum(square(y));
  
  // Values of log-ML for a grid of tauminus2
  mat gridlambda(100, 2, fill::zeros);
  gridlambda.col(0) += logspace(-5, 20, 100);
  for(int j=0; j<100; j++){
    gridlambda(j, 1) = logML(gridlambda(j, 0), p, n, yTy, d2, thetahat);
  }
  
  // Optimal shrinkage
  uword idx = gridlambda.col(1).index_max();
  double lowerVal;
  double upperVal;
  if(idx==0){
    idx = 1;
  }
  if(idx==gridlambda.n_rows){
    idx = idx - 1;
  }
  lowerVal = gridlambda(idx-1, 0);
  upperVal = gridlambda(idx+1, 0);
  const auto obj = [p, n, yTy, d2, thetahat](double x) { return -logML(x, p, n, yTy, d2, thetahat); };
  boost::uintmax_t it = 1000;
  const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  auto lambdaOpt = 0.0, valOpt = 0.0;
  std::tie(lambdaOpt, valOpt) = result;
  
  // Compute marginal posterior means and variances
  colvec vartheta = 1/(d2+lambdaOpt);
  colvec thetabar = (d2 % vartheta) % thetahat;
  colvec betabar = v*thetabar;
  mat vv = v;
  v.each_row() %= sqrt(vartheta.t());
  colvec varbeta = arma::sum(square(v), 1);
  const double sigma2scale = yTy - sum((square(thetahat)%square(d2))/(d2 + lambdaOpt));
  varbeta *= sigma2scale / n;
  
  // Output object
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("grid") = gridlambda,
                                      Rcpp::Named("tauminus2") = lambdaOpt,
                                      Rcpp::Named("logML") = valOpt,
                                      Rcpp::Named("thetabar") = thetabar,
                                      Rcpp::Named("vartheta") = vartheta,
                                      Rcpp::Named("betabar") = betabar,
                                      Rcpp::Named("varbeta") = varbeta,
                                      Rcpp::Named("sigma2scale") = 0.5*sigma2scale,
                                      Rcpp::Named("v") = vv,
                                      Rcpp::Named("n") = n);
  
  return out;
}


// [[Rcpp::export(.bgridge)]]
Rcpp::List bgridge(arma::colvec y, arma::mat X, arma::colvec g, const double a = 0.00001, const double b = 0.00001, double c = 1, const int mcmc = 1000, const int  burnin = 1000, const int thin = 10, bool verbose = true, const int step = 1000){
  
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
      //Rcpp::Rcout << "#### i = " << i << std::endl;
      //Rcpp::Rcout << "## j = " << j << std::endl;
      // Sample from P(\beta | ...)
      //Rcpp::Rcout << "tauminus2 = " << tauminus2 << std::endl;
      //Rcpp::Rcout << "d = " << d << std::endl;
      //Rcpp::Rcout << "wk = " << wk << std::endl;
      V = sigmaminus2*(XTX + tauminus2*D);
      betavar = arma::inv_sympd(V);
      betamean = sigmaminus2 * betavar * XTy;
      beta = rmnorm(p);
      cholvar = arma::chol(betavar);
      beta = trans(beta.t() * cholvar);
      beta += betamean;
      btDb = sum(d % square(beta));

      // Sample from P(\tau^2 | ...)
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
      
      // empirical Bayes
      if(k%step==0){
        double e = sum(mean(log(wksamp.rows(0, k-1)), 0));
        const auto obj = [K, e](double c) { return K*R::digamma(c) - R::digamma(K*c) - e; };
        boost::uintmax_t it = 1000;
        eps_tolerance<double> tol(20);
        const auto result = bracket_and_solve_root(obj, c, 2.0, true, tol, it);
        auto val1 = 0.0, val2 = 0.0;
        std::tie(val1, val2) = result;
        c = val1;
      }
      
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
      Rcpp::Named("sigma2s") = sigma2samp,
      Rcpp::Named("c") = c
    );

  return out;
}


