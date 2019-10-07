\name{bridge}
\alias{bridge}

\title{
Bayesian linear ridge regression
}

\description{
This function implements fast Gibbs algorithms for the Bayesian ridge regression model.
}

\usage{
bridge(y, X, prior = "bp", a = 0.5, b = 0.5,
mcmc = 5000L, burnin = 1000L, thin = 10L,
verbose = TRUE, light = FALSE)
}

\arguments{
  \item{y}{response vector of length n}
  \item{X}{n by p matrix of predictors}
  \item{prior}{character. Type of prior on the variance of regression parameters.}
  \item{a}{hyperparameter}
  \item{b}{hyperparameter}
  \item{mcmc}{integer. Number of desired samples}
  \item{burnin}{integer. Number of burn-in samples}
  \item{thin}{integer. Number of burn-in samples}
  \item{verbose}{logical. Whether information on progress should be be printed.}
  \item{light}{logical. Whether to return a summary of samples only.}
}

\examples{
\dontrun{
# Generate data
set.seed(2019)
n <- 20
p <- 20
X <- matrix(rnorm(n*p), n, p)
trueBeta <- rep(0, p)
trueBeta[1:10] <- 1
y <- X%*%trueBeta + rnorm(n, sd=0.5)

# Run gibbs sampler
res <- bridge(y, X, prior="bp", a=0.5, b=0.5, mcmc=1000, burnin=1000, thin=10, verbose=TRUE)
res$summary

}
}
