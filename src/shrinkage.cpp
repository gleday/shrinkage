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

// [[Rcpp::export(.ridge)]]
Rcpp::List ridge(arma::colvec y, arma::mat X, const int prior, const double a = 0.5, const double b = 0.5, const int mcmc = 1000, const int  burnin = 1000, const int thin = 10, bool verbose = true){

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
	int th = 0;

	//  Gibbs algorithm
	for(int i = 0; i < nruns; i++){
		for(int j = 0; j < th; j++){

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

		  if(k == 0){
		    th = thin;
		  }
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

// [[Rcpp::export(.ridge_fixed)]]
Rcpp::List ridge_fixed(arma::colvec y, arma::mat X){ 
  
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


// [[Rcpp::export(.gridge)]]
Rcpp::List gridge(arma::colvec y, arma::mat X, arma::colvec g, const int prior = 1, double a = 0.00001, double b = 0.00001, const int mcmc = 1000, const int  burnin = 1000, const int thin = 10, bool verbose = true, const int step = 1000){
  
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
  arma::colvec aStar = ones(K);
  arma::colvec ustar = ones(K);
  arma::colvec vstar = ones(K);
  const double wstar = b/(a*a);
  const double xstar = 2*b;
  const double cStar = 0.5*(n+p);
  arma::colvec tauminus2 = ones(K);
  double sigmaminus2 = 1;
  arma::colvec gamma2 = ones(K);
  const int nruns = mcmc + burnin;
  arma::colvec wk = ones(K)*1/K;
  arma::colvec d = zeros(p);
  for(int u = 0; u < K; u++){
    d(find(g == gr(u))).fill(tauminus2(u));
    aStar(u) = a + 0.5*pk(u);
    ustar(u) = a - 0.5*pk(u);
    vstar(u) = -0.5 - 0.5*pk(u);
  }
  arma::mat D = diagmat(d);

  // Useful quantities
  arma::mat XTX = X.t()*X;
  arma::mat XTy = X.t()*y;
  
  // MCMC samples
  arma::mat betasamp = zeros(mcmc, p);
  arma::mat tau2samp = zeros(mcmc, K);
  arma::colvec sigma2samp = zeros(mcmc);
  arma::mat betavar, cholvar;
  arma::colvec betamean, resid;
  arma::colvec beta(p);
  double btDb, bktbk, dStar;
  uvec indk;
  int k = 0;
  int th = 1;
  arma::mat V;
  
  //  Gibbs algorithm
  for(int i = 0; i < nruns; i++){
    for(int j = 0; j < th; j++){
      //Rcpp::Rcout << "#### i = " << i << std::endl;
      //Rcpp::Rcout << "## j = " << j << std::endl;
      // Sample from P(\beta | ...)
      //Rcpp::Rcout << "k = " << k << std::endl;
      //Rcpp::Rcout << "c = " << c << std::endl;
      //Rcpp::Rcout << "dStar = " << dStar << std::endl;
      //Rcpp::Rcout << "sigmaminus2 = " << sigmaminus2 << std::endl;
      //Rcpp::Rcout << "tauminus2 = " << tauminus2 << std::endl;
      //Rcpp::Rcout << "d.min() = " << d.min() << std::endl;
      //Rcpp::Rcout << "d.max() = " << d.max() << std::endl;
      //Rcpp::Rcout << "wk = " << wk << std::endl;
      
      V = sigmaminus2*(XTX + D);
      betavar = arma::inv_sympd(V);
      betamean = sigmaminus2 * betavar * XTy;
      beta = rmnorm(p);
      cholvar = arma::chol(betavar);
      beta = trans(beta.t() * cholvar);
      beta += betamean;

      // Sample from P(\tau_k^2 | ...)
      if(prior == 1){ // if \tau_k^2 ~ invGamma(a, b)
        for(int u = 0; u < K; u++){
          indk = find(g == gr(u));
          bktbk = sum(square(beta(indk)));
          tauminus2(u) = R::rgamma(aStar(u), 1/(b + 0.5*sigmaminus2*bktbk));
        }
      }
      if(prior == 2){ // if \tau_k^2 ~ BetaPrime(a, b)
        for(int u = 0; u < K; u++){
          indk = find(g == gr(u));
          bktbk = sum(square(beta(indk)));
          tauminus2(u) = 1/rgig(ustar(u), sigmaminus2*bktbk, 2*gamma2(u));
          gamma2(u) = R::rgamma(a+b, 1/(1/tauminus2(u) + 1));
        }
      }
      if(prior == 3){ // if \tau_k^2 ~ invGaussian(a, b)
        for(int u = 0; u < K; u++){
          indk = find(g == gr(u));
          bktbk = sum(square(beta(indk)));
          tauminus2(u) = 1/rgig(vstar(u), b + sigmaminus2*bktbk, wstar);
        }
      }
      if(prior == 4){ // if \tau_k^2 ~ Gamma(a, b)
        for(int u = 0; u < K; u++){
          indk = find(g == gr(u));
          bktbk = sum(square(beta(indk)));
          tauminus2(u) = 1/rgig(ustar(u), sigmaminus2*bktbk, xstar);
        }
      }

      // Sample from P(\sigma^2| ...)
      for(int u = 0; u < K; u++){
        d(find(g == gr(u))).fill(tauminus2(u));
      }
      if(d.max() > std::numeric_limits<float>::max()){
        d(find(d > std::numeric_limits<float>::max())).fill(std::numeric_limits<float>::max());
      }
      D = diagmat(d);
      btDb = sum(d % square(beta));
      resid = y-X*beta;
      dStar = 0.5*( sum(square(resid)) + btDb );
      if(dStar > std::numeric_limits<float>::max()){
        dStar = std::numeric_limits<float>::max();
      }
      sigmaminus2 = R::rgamma(cStar, 1/dStar);
      
    }
    
    // Save samples
    if(i >= burnin){
      if(k == 0){
        th = thin;
      }
      betasamp.row(k) += beta.t();
      tau2samp.row(k) += 1/(tauminus2.t());
      sigma2samp(k) += 1/sigmaminus2;
      k++;
      
      // empirical Bayes
      if(k%step==0){
        // estimate c
        //double elogwk = sum(mean(log(wksamp.rows(0, k-1)), 0));
        //double elogtau2 = mean(log(tau2samp.subvec(0, k-1)));
        //double etau2 = mean(tau2samp.subvec(0, k-1));// K*R::digamma(c) - R::digamma(K*c) - e
        //const auto obj = [K, elogwk](double c) { return K*R::digamma(c) - R::digamma(K*c) - elogwk; };
        //const auto obj = [K, elogwk, elogtau2, etau2](double c) { return K*R::digamma(c) - elogwk + 0.00001*0.5*(1/sqrt(c))*etau2 - K*elogtau2 - 0.5*(K/c) - 0.5*K*log(c) + 5*K*log(10); };
        //boost::uintmax_t it = 1000;
        //eps_tolerance<double> tol(20);
        //const auto result = bracket_and_solve_root(obj, c, 2.0, true, tol, it);
        //auto val1 = 0.0, val2 = 0.0;
        //std::tie(val1, val2) = result;
        //c = val1;
        
        //Rcpp::Rcout << "val1 = " << val1 << std::endl;
        //Rcpp::Rcout << "val2 = " << val2 << std::endl;
        
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
      Rcpp::Named("sigma2s") = sigma2samp,
      Rcpp::Named("a") = a,
      Rcpp::Named("b") = b
    );

  return out;
}

// [[Rcpp::export(.gd)]]
Rcpp::List gd(arma::colvec y, arma::mat X, arma::colvec g, const int prior = 1, double a = 0.00001, double b = 0.00001, double c = 1, const int mcmc = 1000, const int  burnin = 1000, const int thin = 10, bool verbose = true, const int step = 1000){
  
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
  double ustar = a - 0.5*p;
  double tauminus2 = 1;
  double sigmaminus2 = 1;
  double gamma2 = 1;
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
  double btDb, bktbk, dStar;
  uvec indk;
  int k = 0;
  int th = 1;
  arma::mat V;
  
  //  Gibbs algorithm
  for(int i = 0; i < nruns; i++){
    for(int j = 0; j < th; j++){
      //Rcpp::Rcout << "#### i = " << i << std::endl;
      //Rcpp::Rcout << "## j = " << j << std::endl;
      // Sample from P(\beta | ...)
      //Rcpp::Rcout << "k = " << k << std::endl;
      //Rcpp::Rcout << "c = " << c << std::endl;
      //Rcpp::Rcout << "dStar = " << dStar << std::endl;
      //Rcpp::Rcout << "sigmaminus2 = " << sigmaminus2 << std::endl;
      //Rcpp::Rcout << "tauminus2 = " << tauminus2 << std::endl;
      //Rcpp::Rcout << "d.min() = " << d.min() << std::endl;
      //Rcpp::Rcout << "d.max() = " << d.max() << std::endl;
      //Rcpp::Rcout << "wk = " << wk << std::endl;
      
      V = sigmaminus2*(XTX + tauminus2*D);
      betavar = arma::inv_sympd(V);
      betamean = sigmaminus2 * betavar * XTy;
      beta = rmnorm(p);
      cholvar = arma::chol(betavar);
      beta = trans(beta.t() * cholvar);
      beta += betamean;
      btDb = sum(d % square(beta));
      
      if(prior == 1){
        
        // Sample from P(\tau^2 | ...)
        tauminus2 = 1/rgig(ustar, sigmaminus2*btDb, 2*b);
        
        // Sample from P(\omega_k | ...)
        for(int u = 0; u < K; u++){
          indk = find(g == gr(u));
          bktbk = sum(square(beta(indk)));
          wk(u) = rgig(c - 0.5 * pk(u), sigmaminus2*pk(u)*bktbk, 2*b);
        }
        wk /= sum(wk);
        
      }else{
        
        // Sample from P(\tau^2 | ...)
        tauminus2 = 1/rgig(ustar, sigmaminus2*btDb, 2*gamma2);
        
        // Sample from P(\gamma^2 | ...)
        gamma2 = R::rgamma(a+b, 1/(1/tauminus2 + 1));
        
        // Sample from P(\omega_k | ...)
        for(int u = 0; u < K; u++){
          indk = find(g == gr(u));
          bktbk = sum(square(beta(indk)));
          wk(u) = rgig(c - 0.5 * pk(u), sigmaminus2*pk(u)*bktbk, 2*gamma2);
        }
        wk /= sum(wk);
        
      }
      
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
      if(dStar > std::numeric_limits<float>::max()){
        dStar = std::numeric_limits<float>::max();
      }
      sigmaminus2 = R::rgamma(cStar, 1/dStar);
    }
    
    // Save samples
    if(i >= burnin){
      if(k == 0){
        th = thin;
      }
      betasamp.row(k) += beta.t();
      tau2samp(k) += 1/(tauminus2);
      wksamp.row(k) += wk.t();
      sigma2samp(k) += 1/sigmaminus2;
      k++;
      
      // empirical Bayes
      if(k%step==0){
        // estimate c
        double elogwk = sum(mean(log(wksamp.rows(0, k-1)), 0));
        //double elogtau2 = mean(log(tau2samp.subvec(0, k-1)));
        //double etau2 = mean(tau2samp.subvec(0, k-1));// K*R::digamma(c) - R::digamma(K*c) - e
        const auto obj = [K, elogwk](double c) { return K*R::digamma(c) - R::digamma(K*c) - elogwk; };
        //const auto obj = [K, elogwk, elogtau2, etau2](double c) { return K*R::digamma(c) - elogwk + 0.00001*0.5*(1/sqrt(c))*etau2 - K*elogtau2 - 0.5*(K/c) - 0.5*K*log(c) + 5*K*log(10); };
        boost::uintmax_t it = 1000;
        eps_tolerance<double> tol(20);
        const auto result = bracket_and_solve_root(obj, c, 2.0, true, tol, it);
        auto val1 = 0.0, val2 = 0.0;
        std::tie(val1, val2) = result;
        c = val1;
        
        Rcpp::Rcout << "val1 = " << val1 << std::endl;
        Rcpp::Rcout << "val2 = " << val2 << std::endl;
        
        // update a and b
        a = K*c;
        if(prior == 1){
          b = 0.00001*sqrt(c);
        }else{
          b = a;
        }
        ustar = a - 0.5*p;
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

