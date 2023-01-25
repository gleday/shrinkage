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

// Sample from:
//   N(\mu, \Sigma),
// where
// \Sigma = (\Phi^T \Phi + D^{-1})^{-1}
// \mu = \Sigma \Phi^T \alpha
// using Rue (2001) algorithm
colvec r_beta_rue(mat phiTphi, colvec Ddiag, colvec phiTalpha){
  
  arma::mat Q = phiTphi;
  Q.diag() += 1/Ddiag;
  arma::mat R;
  bool success;
  success = chol(R, Q);
  //Rcpp::Rcout << "success 1 = " << success << std::endl;
  if(success == false){
    Q.diag() += 1e-3;
    success = chol(R, Q);
  }
  //Rcpp::Rcout << "success 2 = " << success << std::endl;
  arma::sp_mat L = sp_mat(R);
  arma::colvec v = spsolve(L.t(), phiTalpha, "lapack");
  arma::colvec m = spsolve(L, v, "lapack");
  arma::colvec s = rmnorm(phiTphi.n_cols);
  arma::colvec w = spsolve(L, s, "lapack");
  arma::colvec b = m + w;
  
  return b;
}

// Sample from:
//   N(\mu, \Sigma),
// where
// \Sigma = (\Phi^T \Phi + D^{-1})^{-1}
// \mu = \Sigma \Phi^T \alpha
// using Bhattacharya et al. (2016) algorithm
// to be generalized for logistic regression, see Cong et al. (2017, Algprithm 4)
arma::colvec r_beta_bha(mat phi, colvec Ddiag, colvec alpha){
  
  arma::colvec u = rmnorm(phi.n_cols)%sqrt(Ddiag);
  arma::colvec delta = rmnorm(phi.n_rows);
  arma::colvec v = phi*u + delta;
  arma::mat DphiT = phi.t();
  DphiT.each_col() %= Ddiag;
  arma::mat M = phi*DphiT + eye(phi.n_rows, phi.n_rows);
  arma::sp_mat L = sp_mat(chol(M));
  arma::colvec r = spsolve(L.t(), alpha-v, "lapack");
  arma::colvec w = spsolve(L, r, "lapack");
  arma::colvec b = u + DphiT*w;
  
  return b;
}

double logML(const double tauminus2, const int p, const int n, const double yTy, colvec eigvals, colvec thetahat){
  double out = -0.5*n*log(datum::pi);
  out += 0.5*p*log(tauminus2);
  out += lgamma(0.5*n);
  colvec egs = zeros(p);
  egs.subvec(0, eigvals.n_elem-1) += eigvals;
  out -= 0.5*sum(log(egs + tauminus2));
  out -= 0.5*n*log(0.5*yTy - 0.5*sum((square(thetahat)%square(eigvals))/(eigvals + tauminus2)));
  return(out);
}

double logML(arma::colvec tauminus2, const int p, const int n, const double yTy, colvec eigvals, colvec thetahat){
  double out = -0.5*n*log(datum::pi);
  out += 0.5*p*sum(log(tauminus2));
  out += lgamma(0.5*n);
  colvec eigs = zeros(p);
  eigs.subvec(0, eigvals.n_elem-1) += eigvals;
  out -= 0.5*sum(log(eigs + tauminus2));
  out -= 0.5*n*log(0.5*yTy - 0.5*sum((square(thetahat)%square(eigvals))/(eigvals + tauminus2)));
  return(out);
}

void fast_svd(arma::mat* X, arma::mat* u, arma::colvec* d, arma::mat* v){
  
  // SVD decomposition
  if(X->n_rows >= X->n_cols){
    arma::mat XTX = X->t()* *X;
    svd_econ(*u, *d, *v, XTX, "right");
    *d = sqrt(*d);
    *u = *X * *v * diagmat(1/ *d);
  }
  else{
    arma::mat XXT = *X * X->t();
    svd_econ(*u, *d, *v, XXT, "left");
    *d = sqrt(*d);
    *v = X->t() * *u *diagmat(1/ *d);
  }
}

// [[Rcpp::export(.fast_svd)]]
Rcpp::List fast_svd_list(arma::mat X){
  
  // initialization
  arma::mat u, v;
  arma::colvec d;

  fast_svd(&X, &u, &d, &v);
  
  // Output object
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("u") = u,
    Rcpp::Named("d") = d,
    Rcpp::Named("v") = v
  );
  
  return out;
}

void center(arma::mat* x){
  arma::rowvec cm = mean(*x);
  x->each_row() -= cm;
}

void scale(arma::mat* x){
  rowvec cs = 1/sqrt(sum(square(*x)) / x->n_rows);
  x->each_row() %= cs;
}

double get_brg_opt_tauminus2(const double* yTy, arma::mat* u, arma::colvec* d, arma::mat* v, arma::colvec* thetahat){
  
  const int n = u->n_rows;
  const int p = v->n_rows;
  boost::uintmax_t it = 1000;
  
  // Useful quantities
  arma::colvec d2 = square(*d);
  
  // Values of log-ML for a grid of tauminus2
  mat gridlambda(100, 2, fill::zeros);
  gridlambda.col(0) += logspace(-5, 20, 100);
  for(int j=0; j<100; j++){
    gridlambda(j, 1) = logML(gridlambda(j, 0), p, n, *yTy, d2, *thetahat);
  }
  
  // Optimal shrinkage
  uword idx = gridlambda.col(1).index_max();
  double lowerVal;
  double upperVal;
  if(idx==0){
    idx = 1;
  }
  if(idx == gridlambda.n_rows){
    idx = idx - 1;
  }
  lowerVal = gridlambda(idx-1, 0);
  upperVal = gridlambda(idx+1, 0);
  const auto obj = [p, n, yTy, d2, thetahat](double x) { return -logML(x, p, n, *yTy, d2, *thetahat); };
  const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
  auto lambdaOpt = 0.0, valOpt = 0.0;
  std::tie(lambdaOpt, valOpt) = result;
  
  return(lambdaOpt);
}

// [[Rcpp::export(.get_brg_opt_tauminus2)]]
double brg_eb_tauminus2(arma::colvec y, arma::mat X){
  
  // SVD decomposition
  arma::mat u, v;
  arma::colvec d;
  fast_svd(&X, &u, &d, &v);
  arma::colvec d2 = square(d);
  
  // quantities needed
  const double yTy = sum(square(y));
  arma::colvec thetahat = (1/ d) % (u.cols(0, d.n_elem - 1).t() * y);
  
  // maximizer of log-ML
  double lambdaOpt = get_brg_opt_tauminus2(&yTy, &u, &d, &v, &thetahat);
  
  return lambdaOpt;
}


// // [[Rcpp::export(.get_brg_opt_tauminus2)]]
// Rcpp::List get_brg_opt_tauminus2(arma::colvec y, arma::mat X){
// 
//   // svd decomposition
//   Rcpp::List mysvd = fast_svd(X);
//   arma::colvec d2 = square(mysvd[1]);
// 
//   // Useful quantities
//   arma::colvec thetahat = (1/mysvd[1]) % (mysvd[0].cols(0, mysvd[1].n_elem-1).t()*y);
//   const double yTy = sum(square(y));
//   boost::uintmax_t it = 1000;
// 
//   // Values of log-ML for a grid of tauminus2
//   mat gridlambda(100, 2, fill::zeros);
//   gridlambda.col(0) += logspace(-5, 20, 100);
//   for(int j=0; j<100; j++){
//     gridlambda(j, 1) = logML(gridlambda(j, 0), X.n_cols, X.n_rows, yTy, d2, thetahat);
//   }
// 
//   // // Optimal shrinkage
//   // uword idx = gridlambda.col(1).index_max();
//   // double lowerVal;
//   // double upperVal;
//   // if(idx==0){
//   //   idx = 1;
//   // }
//   // if(idx==gridlambda.n_rows){
//   //   idx = idx - 1;
//   // }
//   // lowerVal = gridlambda(idx-1, 0);
//   // upperVal = gridlambda(idx+1, 0);
//   // const auto obj = [p, n, yTy, d2, thetahat](double x) { return -logML(x, p, n, yTy, d2, thetahat); };
//   // const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
//   // lambdaOpt = 0.0, valOpt = 0.0;
//   // std::tie(lambdaOpt, valOpt) = result;
// 
//   // Output object
//   Rcpp::List out = Rcpp::List::create(
//     Rcpp::Named("gridlambda") = gridlambda
//   );
//   
//   return out;
// }

////////////////////////////////////////////////////
//-------------- External functions --------------//
////////////////////////////////////////////////////

// [[Rcpp::export(.brg_gibbs)]]
Rcpp::List brg_gibbs(arma::colvec y, arma::mat X, const int prior,
                     const double a = 0.5, const double b = 0.5,
                     const int mcmc = 1000, const int  burnin = 1000,
                     const int thin = 10, bool verbose = true,
                     const int bp = 1){

	// Dimension data
	const int n = X.n_rows;
	const int p = X.n_cols;

	// Constants
	const double aStar = a + 0.5*p;
	const double bStar = b + 0.5*p;
	const double cStar = 0.5*(n+p);
	const double ustar = a - 0.5*p;
	const double vstar = -0.5 - 0.5*p;
	const double wstar = b/(a*a);
  const double xstar = 2*b;
  const int nruns = mcmc + burnin;
  
  // Starting values of parameters
	double gamma2 = 1;
	double tauminus2 = 1;
	double sigmaminus2 = 1;

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
	arma::colvec thetahat = (1/d) % (u.cols(0, d.n_elem-1).t()*y);
	arma::mat ud = u * diagmat(d);
	u.clear();
  
	// MML as starting value
	boost::uintmax_t it = 1000;
	const double yTy = sum(square(y));
	mat gridlambda(100, 2, fill::zeros);
	gridlambda.col(0) += logspace(-5, 20, 100);
	for(int j=0; j<100; j++){
	  gridlambda(j, 1) = logML(gridlambda(j, 0), p, n, yTy, square(d), thetahat);
	}
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
	const auto obj = [p, n, yTy, d, thetahat](double x) { return -logML(x, p, n, yTy, square(d), thetahat); };
	const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
	auto lambdaOpt = 0.0, valOpt = 0.0;
	std::tie(lambdaOpt, valOpt) = result;
	tauminus2 = lambdaOpt;
	
	// MCMC samples
	arma::mat thetasamp = zeros(d.n_elem, mcmc);
	arma::colvec tau2samp = zeros(mcmc);
	arma::colvec sigma2samp = zeros(mcmc);
	
	// Initialization
	arma::colvec thetavar;
	arma::colvec thetamean;
	arma::colvec theta(d.n_elem);
	arma::mat Xbeta;
	double btb;
	double dStar;
	int k = 0;
	int th = 1;
	
	//  Gibbs algorithm
	for(int i = 0; i < nruns; i++){
		for(int j = 0; j < th; j++){
		  
		  if(verbose && (i == 0) && (j == 0) && (burnin>0)){
		    Rcpp::Rcout << "burnin..." << std::endl;
		  }

		  // Sample from P(\theta | ...)
			thetavar = 1/(square(d) + tauminus2);
			thetamean = thetavar % duy;
			thetavar /= sigmaminus2;
			//thetavar %= sigma2;
			for(uword l = 0; l < thetamean.n_elem; l++){
				theta(l) = R::rnorm(thetamean(l), sqrt(thetavar(l)));
			}

			// Sample from P(\tau^2 | ...)
			btb = sum(square(theta));
			//Rcpp::Rcout << "btb = " << btb << std::endl;
			if(prior == 1){ // if \tau^2 ~ invGamma(a, b)
			  tauminus2 = R::rgamma(aStar, 1/(b + 0.5*sigmaminus2*btb));
			}
			if(prior == 2){ // if \tau^2 ~ BetaPrime(a, b)
			  if(bp == 1){// using gamma-gamma representation
			    tauminus2 = 1/rgig(ustar, sigmaminus2*btb, 2*gamma2);
			    gamma2 = R::rgamma(a+b, 1/(1/tauminus2 + 1));
			    //tau2 = rgig(ustar, (1/sigma2)*btb, 2*gamma2);
			    //gamma2 = R::rgamma(a+b, 1/(tau2 + 1));
			  }
			  if(bp == 2){// using inverse gamma-inverse gamma representation
			    tauminus2 = R::rgamma(bStar, 1/(0.5*sigmaminus2*btb + 1/gamma2));
			    gamma2 = 1/R::rgamma(a+b, 1/(tauminus2 + 1));
			    //tauminus2 = R::rgamma(bStar, 1/(0.5*sigmaminus2*btb + gammaminus2));
			  }
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
			thetasamp.col(k) = theta;
			tau2samp(k) = 1/tauminus2;
			sigma2samp(k) = 1/sigmaminus2;
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
	arma::mat betasamp = v * thetasamp;
	v.clear();
	thetasamp.clear();

	// Output
	Rcpp::List out = Rcpp::List::create(
			Rcpp::Named("betas") = betasamp,
			Rcpp::Named("tau2s") = tau2samp,
			Rcpp::Named("sigma2s") = sigma2samp
		);

	return out;
}

// [[Rcpp::export(.brg_closedform)]]
Rcpp::List brg_closedform(arma::colvec y, arma::mat X, bool fixed = false,
                          const double tau2 = 1){
  
  // SVD decomposition
  arma::mat u, v;
  arma::colvec d;
  fast_svd(&X, &u, &d, &v);
  arma::colvec d2 = square(d);
  
  // quantities needed
  const double yTy = sum(square(y));
  arma::colvec thetahat = (1/ d) % (u.cols(0, d.n_elem - 1).t() * y);
  
  // maximizer of log-ML
  double lambdaOpt;
  if(fixed){
    lambdaOpt = 1/tau2;
  }
  else{
    lambdaOpt = get_brg_opt_tauminus2(&yTy, &u, &d, &v, &thetahat);
  }
  
  // posterior mean and variance of \theta
  colvec thetavar = 1/(d2+lambdaOpt);
  colvec thetabar = (d2 % thetavar) % thetahat;
  colvec betabar = v*thetabar;
  mat vv = v;
  v.each_row() %= sqrt(thetavar.t());
  colvec betavar = arma::sum(square(v), 1);
  const double sigma2scale = yTy - sum((square(thetahat)%square(d2))/(d2 + lambdaOpt));
  betavar *= sigma2scale / X.n_rows;
  thetavar *= sigma2scale / X.n_rows;
  
  // Output object
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("beta") = Rcpp::List::create(
      Rcpp::Named("bar") = betabar,
      Rcpp::Named("var") = betavar
    ),
    Rcpp::Named("tau2s") = 1/lambdaOpt,
    Rcpp::Named("svd") = Rcpp::List::create(
      Rcpp::Named("d") = d,
      Rcpp::Named("v") = vv
    ),
    Rcpp::Named("theta") = Rcpp::List::create(
      Rcpp::Named("hat") = thetahat,
      Rcpp::Named("bar") = thetabar,
      Rcpp::Named("var") = thetavar
    ),
    Rcpp::Named("sigma2scale") = 0.5*sigma2scale,
    Rcpp::Named("n") = X.n_rows
    );
  
  return out;
}

// [[Rcpp::export(.brl_gibbs)]]
Rcpp::List brl_gibbs(arma::colvec y, arma::mat X, arma::colvec g,
                     const int prior = 1, double a = 0.00001,
                     double b = 0.00001, const int mcmc = 1000,
                     const int  burnin = 1000, const int thin = 10,
                     bool verbose = true, const int bp = 2){
  
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
  arma::colvec bStar = ones(K);
  arma::colvec ustar = ones(K);
  arma::colvec vstar = ones(K);
  const double wstar = b/(a*a);
  const double xstar = 2*b;
  const double sigshape = 0.5*(n+p);
  arma::colvec tauminus2 = ones(K);
  double sigmaminus2 = 1;
  arma::colvec gamma2 = ones(K);
  arma::colvec wk = ones(K)*1/K;
  arma::colvec Ddiag = zeros(p);
  for(int u = 0; u < K; u++){
    Ddiag(find(g == gr(u))).fill(tauminus2(u));
    aStar(u) = a + 0.5*pk(u);
    bStar(u) = b + 0.5*pk(u);
    ustar(u) = a - 0.5*pk(u);
    vstar(u) = -0.5 - 0.5*pk(u);
  }
  
  // Useful quantities
  arma::mat XTX;
  arma::colvec XTy;
  if(n > p){
    XTX = X.t()*X;
    XTy = X.t()*y;
  }
  
  // MCMC samples
  const int nruns = mcmc + burnin;
  arma::mat betasamp = zeros(p, mcmc);
  arma::mat tau2samp = zeros(K, mcmc);
  arma::colvec sigma2samp = zeros(mcmc);
  arma::colvec resid;
  arma::colvec beta(p);
  double btDb, bktbk, sigscale;
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
      if(n > p){
        beta = r_beta_rue(XTX*sigmaminus2, 1/(Ddiag*sigmaminus2), XTy*sigmaminus2);
      }else{
        beta = r_beta_bha(X*sqrt(sigmaminus2), 1/(Ddiag*sigmaminus2), y*sqrt(sigmaminus2));
      }
      
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
          if(bp == 1){// using gamma-gamma representation
            //tau2(u) = rgig(ustar(u), (1/sigma2)*bktbk, 2*gamma2(u));
            //gamma2(u) = R::rgamma(a+b, 1/(tau2(u) + 1));
            tauminus2(u) = 1/rgig(ustar(u), sigmaminus2*bktbk, 2*gamma2(u));
            gamma2(u) = R::rgamma(a+b, 1/((1/tauminus2(u)) + 1));
          }
          if(bp == 2){// using inverse gamma-inverse gamma representation
            tauminus2(u) = R::rgamma(bStar(u), 1/(0.5*sigmaminus2*bktbk + (gamma2(u))));
            gamma2(u) = 1/R::rgamma(a+b, 1/(tauminus2(u) + 1));
          }
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
        Ddiag(find(g == gr(u))).fill(tauminus2(u));
      }
      if(Ddiag.max() > std::numeric_limits<float>::max()){
        Ddiag(find(Ddiag > std::numeric_limits<float>::max())).fill(std::numeric_limits<float>::max());
      }
      btDb = sum(Ddiag % square(beta));
      resid = y-X*beta;
      sigscale = 0.5*( sum(square(resid)) + btDb );
      if(sigscale > std::numeric_limits<float>::max()){
        sigscale = std::numeric_limits<float>::max();
      }
      sigmaminus2 = R::rgamma(sigshape, 1/sigscale);
      
    }
    
    // Save samples
    if(i >= burnin){
      if(k == 0){
        th = thin;
      }
      betasamp.col(k) += beta;
      tau2samp.col(k) += 1/tauminus2;
      sigma2samp(k) += 1/sigmaminus2;
      k++;
      
      // empirical Bayes
      //if(k%step==0){
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
        
      //}
      
      if(verbose && (k%1000==0)){
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
      Rcpp::Named("sigma2s") = sigma2samp
    );

  return out;
}

/*
// [[Rcpp::export(.brl_closedform)]]
Rcpp::List brl_closedform(arma::colvec y, arma::mat X, arma::colvec tauminus2,
 bool fixed = false, const double tau2 = 1){
  
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  // scale X with tau2
  X.each_row() %= 1/sqrt(tauminus2);

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
  Rcpp::List out = Rcpp::List::create(
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
*/



